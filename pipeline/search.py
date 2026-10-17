"""Period searches: astropy BoxLeastSquares and transitleastsquares.

BLS grid: frequency spacing df = q_ref / (oversample * baseline * P), i.e. the phase drift of a transit of
reference duration q_ref over the full baseline is at most 1/oversample of the duration. For a 1426-d
baseline, q_ref = 2 h and oversample = 3 that is 211,097 trial periods between 1 and 60 d. Durations
1-6 h. power = objective 'snr' (depth / depth_err of the best box at each period, the quantity usually
called "BLS SNR"); SDE = (power - median) / std over the periodogram.

TLS: Hippke & Heller 2019 (A&A 623, A39); period grid from transitleastsquares.period_grid for the given
stellar mass/radius (Ofir 2014 optimal sampling), restricted to [period_min, period_max] (default 1-30 d,
~181k periods for this baseline at oversampling_factor=3; measured 1.6 ms/period with use_threads=3 on
46k points => ~5 min). Quadratic limb darkening u=(0.2125, 0.3984) = TIC/TLS catalog_info for TIC
120686465 (Claret 2017 TESS-band tables; adequate for a 3633 K dwarf in the Kepler band for a search).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from astropy.timeseries import BoxLeastSquares

from .fetch import BKJD_OFFSET

STAR = {"R_star": 0.507, "M_star": 0.505, "R_star_min": 0.45, "R_star_max": 0.56,
        "M_star_min": 0.45, "M_star_max": 0.56, "u": [0.2125, 0.3984]}


def bls_period_grid(baseline: float, pmin: float = 1.0, pmax: float = 60.0, q_ref_days: float = 2.0 / 24,
                    oversample: float = 3.0) -> np.ndarray:
    fmax, fmin = 1.0 / pmin, 1.0 / pmax
    freqs = [fmax]
    while freqs[-1] > fmin:
        P = 1.0 / freqs[-1]
        freqs.append(freqs[-1] - q_ref_days / (oversample * baseline * P))
    return np.sort(1.0 / np.array(freqs))


@dataclass
class BLSResult:
    period: float; t0: float; duration: float; depth: float; depth_err: float; power: float; sde: float
    depth_snr: float; n_periods: int; runtime_s: float
    periods: np.ndarray = field(repr=False); powers: np.ndarray = field(repr=False)
    stats: dict = field(default_factory=dict, repr=False)
    top_peaks: list = field(default_factory=list)

    def summary(self) -> dict:
        return {"period_d": self.period, "t0_bkjd": self.t0, "t0_bjd_tdb": self.t0 + BKJD_OFFSET,
                "duration_h": self.duration * 24, "depth_ppm": self.depth * 1e6, "depth_err_ppm": self.depth_err * 1e6,
                "power_snr_objective": self.power, "sde": self.sde, "depth_snr": self.depth_snr,
                "n_periods": self.n_periods, "runtime_s": self.runtime_s, "top_peaks": self.top_peaks,
                "stats": self.stats}


def _top_peaks(periods, power, n=6, tol=0.01):
    o = np.argsort(power)[::-1]; shown = []
    sde = (power - np.median(power)) / np.std(power)
    for j in o:
        p = periods[j]
        if all(abs(p / q - 1) > tol and abs(p / q - 2) > tol and abs(p / q - 0.5) > tol
               and abs(p / q - 3) > tol and abs(p / q - 1 / 3) > tol and abs(p / q - 4) > tol for q in shown):
            shown.append(p)
            yield {"period_d": float(p), "power": float(power[j]), "sde": float(sde[j])}
        if len(shown) >= n:
            break


def bls_search(t, f, e, pmin=1.0, pmax=60.0, durations_h=(1, 1.5, 2, 2.5, 3, 4, 5, 6), oversample=3.0,
               q_ref_h=2.0, refine=True, verbose=True) -> BLSResult:
    periods = bls_period_grid(np.ptp(t), pmin, pmax, q_ref_h / 24, oversample)
    durs = np.asarray(durations_h, float) / 24
    bls = BoxLeastSquares(t, f, dy=e)
    t0 = time.time()
    r = bls.power(periods, durs, objective="snr")
    rt = time.time() - t0
    pw = np.asarray(r.power); sde = (pw - np.median(pw)) / np.std(pw)
    i = int(np.argmax(pw))
    P, T0, D = float(r.period[i]), float(r.transit_time[i]), float(r.duration[i])
    if verbose:
        print(f"[bls] {len(periods)} periods x {len(durs)} durations in {rt:.1f}s; best P={P:.6f} d power={pw[i]:.2f} SDE={sde[i]:.2f}")
    if refine:
        fine = np.linspace(P * (1 - 0.003), P * (1 + 0.003), 1500)
        fdurs = np.arange(1.0, 4.01, 0.25) / 24
        t1 = time.time(); rr = bls.power(fine, fdurs, objective="snr"); rt += time.time() - t1
        j = int(np.argmax(rr.power))
        P, T0, D = float(rr.period[j]), float(rr.transit_time[j]), float(rr.duration[j])
        if verbose:
            print(f"[bls] refined: P={P:.6f} d T0={T0:.5f} BKJD dur={D*24:.2f} h depth={rr.depth[j]*1e6:.0f}+/-{rr.depth_err[j]*1e6:.0f} ppm (+{time.time()-t1:.1f}s)")
    st = bls.compute_stats(P, D, T0)
    stats = {k: (float(v) if np.ndim(v) == 0 else np.asarray(v).tolist()) for k, v in st.items()
             if k not in ("transit_times", "per_transit_count", "per_transit_log_likelihood")}
    stats["n_transits_with_data"] = int(np.sum(np.asarray(st["per_transit_count"]) > 0))
    stats["transit_times_bkjd"] = np.asarray(st["transit_times"]).tolist()
    stats["per_transit_count"] = np.asarray(st["per_transit_count"]).tolist()
    depth, depth_err = float(st["depth"][0]), float(st["depth"][1])
    return BLSResult(P, T0, D, depth, depth_err, float(pw[i]), float(sde[i]), depth / depth_err, len(periods), rt,
                     np.asarray(r.period), pw, stats, list(_top_peaks(np.asarray(r.period), pw)))


def bls_refine(t, f, e, period, t0_guess=None, half_width_frac=0.003, n=1500, durations_h=np.arange(1.0, 4.01, 0.25),
               verbose=True) -> dict:
    """Fine BLS around a known period (used on the masked-detrended light curve): refined P, T0, duration,
    depth +/- err and astropy compute_stats (odd/even, phase-0.5 depth, harmonic test, per-transit counts)."""
    bls = BoxLeastSquares(t, f, dy=e)
    fine = np.linspace(period * (1 - half_width_frac), period * (1 + half_width_frac), n)
    t0 = time.time(); r = bls.power(fine, np.asarray(durations_h) / 24, objective="snr"); rt = time.time() - t0
    j = int(np.argmax(r.power))
    P, T0, D = float(r.period[j]), float(r.transit_time[j]), float(r.duration[j])
    st = bls.compute_stats(P, D, T0)
    out = {"period_d": P, "t0_bkjd": T0, "t0_bjd_tdb": T0 + BKJD_OFFSET, "duration_h": D * 24,
           "depth_ppm": float(st["depth"][0]) * 1e6, "depth_err_ppm": float(st["depth"][1]) * 1e6,
           "depth_snr": float(st["depth"][0] / st["depth"][1]), "power_snr_objective": float(r.power[j]),
           "depth_odd_ppm": [float(st["depth_odd"][0]) * 1e6, float(st["depth_odd"][1]) * 1e6],
           "depth_even_ppm": [float(st["depth_even"][0]) * 1e6, float(st["depth_even"][1]) * 1e6],
           "depth_phase_0p5_ppm": [float(st["depth_phased"][0]) * 1e6, float(st["depth_phased"][1]) * 1e6],
           "depth_half_period_ppm": [float(st["depth_half"][0]) * 1e6, float(st["depth_half"][1]) * 1e6],
           "harmonic_delta_log_likelihood": float(st["harmonic_delta_log_likelihood"]),
           "n_transits_with_data": int(np.sum(np.asarray(st["per_transit_count"]) > 0)), "runtime_s": rt}
    if verbose:
        print(f"[bls-refine] P={P:.6f} d T0={T0:.5f} BKJD dur={D*24:.2f} h depth={out['depth_ppm']:.0f}+/-{out['depth_err_ppm']:.0f} ppm "
              f"(SNR {out['depth_snr']:.1f}); odd {out['depth_odd_ppm'][0]:.0f} even {out['depth_even_ppm'][0]:.0f}; "
              f"phase0.5 {out['depth_phase_0p5_ppm'][0]:.0f}+/-{out['depth_phase_0p5_ppm'][1]:.0f} ppm; {out['n_transits_with_data']} transits with data")
    return out


def null_sde_excluding(periods, power, period, band=0.02, harmonics=(1, 2, 3, 4, 0.5, 1 / 3, 0.25)):
    """Best SDE of a (null) periodogram outside +/-band around the candidate period and its harmonics."""
    sde = (power - np.median(power)) / np.std(power)
    keep = np.ones(len(periods), bool)
    for h in harmonics:
        keep &= np.abs(periods / (period * h) - 1) > band
    k = int(np.argmax(np.where(keep, sde, -np.inf)))
    return {"best_sde_excluding_candidate_band": float(sde[k]), "period_d": float(periods[k]), "power": float(power[k]),
            "best_sde_overall": float(np.max(sde)), "period_overall_d": float(periods[int(np.argmax(sde))])}


@dataclass
class TLSResult:
    period: float; period_unc: float; t0: float; duration: float; depth: float; sde: float; snr: float
    rp_rs: float; transit_count: int; distinct_transit_count: int; odd_even_mismatch: float
    depth_mean_odd: tuple; depth_mean_even: tuple; n_periods: int; runtime_s: float; fap: float
    periods: np.ndarray = field(repr=False); powers: np.ndarray = field(repr=False)
    folded_phase: np.ndarray = field(repr=False); folded_y: np.ndarray = field(repr=False)
    model_phase: np.ndarray = field(repr=False); model_y: np.ndarray = field(repr=False)
    transit_times: list = field(default_factory=list, repr=False); transit_depths: list = field(default_factory=list, repr=False)
    per_transit_count: list = field(default_factory=list, repr=False); snr_per_transit: list = field(default_factory=list, repr=False)
    settings: dict = field(default_factory=dict)

    def summary(self) -> dict:
        return {"period_d": self.period, "period_unc_d": self.period_unc, "t0_bkjd": self.t0,
                "t0_bjd_tdb": self.t0 + BKJD_OFFSET, "duration_h": self.duration * 24, "depth_ppm": (1 - self.depth) * 1e6,
                "sde": self.sde, "snr": self.snr, "fap": self.fap, "rp_rs": self.rp_rs,
                "transit_count": self.transit_count, "distinct_transit_count": self.distinct_transit_count,
                "odd_even_mismatch_sigma": self.odd_even_mismatch,
                "depth_mean_odd_ppm": [(1 - self.depth_mean_odd[0]) * 1e6, self.depth_mean_odd[1] * 1e6],
                "depth_mean_even_ppm": [(1 - self.depth_mean_even[0]) * 1e6, self.depth_mean_even[1] * 1e6],
                "n_periods": self.n_periods, "runtime_s": self.runtime_s, "settings": self.settings,
                "transit_times_bkjd": self.transit_times, "transit_depths": self.transit_depths,
                "per_transit_count": self.per_transit_count, "snr_per_transit": self.snr_per_transit}


def tls_search(t, f, e, pmin=1.0, pmax=30.0, oversampling_factor=3, duration_grid_step=1.1, use_threads=3,
               star=STAR, verbose=True) -> TLSResult:
    from transitleastsquares import transitleastsquares
    model = transitleastsquares(t, f, e)
    kw = dict(period_min=pmin, period_max=pmax, oversampling_factor=oversampling_factor,
              duration_grid_step=duration_grid_step, use_threads=use_threads, show_progress_bar=False, **star)
    t0 = time.time()
    r = model.power(**kw)
    rt = time.time() - t0
    if verbose:
        print(f"[tls] {len(r.periods)} periods ({pmin}-{pmax} d, osf={oversampling_factor}) in {rt:.1f}s with {use_threads} threads; "
              f"P={r.period:.6f} d SDE={r.SDE:.2f} SNR={r.snr:.2f} depth={(1-r.depth)*1e6:.0f} ppm dur={r.duration*24:.2f} h "
              f"T0={r.T0:.5f} BKJD transits={r.transit_count} distinct={r.distinct_transit_count}")

    def _f(x):
        try:
            return float(x)
        except Exception:
            return float("nan")
    return TLSResult(float(r.period), _f(r.period_uncertainty), float(r.T0), float(r.duration), float(r.depth),
                     float(r.SDE), _f(r.snr), _f(r.rp_rs), int(r.transit_count), int(r.distinct_transit_count),
                     _f(r.odd_even_mismatch), (_f(r.depth_mean_odd[0]), _f(r.depth_mean_odd[1])),
                     (_f(r.depth_mean_even[0]), _f(r.depth_mean_even[1])), len(r.periods), rt, _f(r.FAP),
                     np.asarray(r.periods), np.asarray(r.power), np.asarray(r.folded_phase), np.asarray(r.folded_y),
                     np.asarray(r.model_folded_phase), np.asarray(r.model_folded_model),
                     np.asarray(r.transit_times).tolist(), np.asarray(r.transit_depths, float).tolist(),
                     np.asarray(r.per_transit_count).tolist(), np.asarray(r.snr_per_transit, float).tolist(),
                     {k: (v if not isinstance(v, np.ndarray) else v.tolist()) for k, v in kw.items()})
