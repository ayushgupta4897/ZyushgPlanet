"""Stitch, clean and detrend a KBONUS-BKG light curve.

Steps (in this order):
1. per-quarter quality mask: SAP_QUALITY == 0 by default ("hard"; this is the Kepler TPF QUALITY
   bitmask copied into the HLSP) or lightkurve's KeplerQualityFlags.DEFAULT_BITMASK ("default"),
   plus finite time/flux/err and err > 0;
2. per-quarter normalisation by the quarter median (PSF flux levels differ quarter to quarter
   because the PSF model and pixel footprint change with season);
3. wotan biweight detrending with a window ~3x the longest plausible transit duration
   (see `default_window_days`), segments broken at gaps > `break_tolerance` days;
4. asymmetric sigma clipping on the *detrended* residuals: only upward outliers (flares, cosmic rays,
   bad-cadence spikes) above +n_sigma robust sigma are removed; nothing below the trend is ever
   clipped, so transits of any depth survive.

Why the detrending window is ~0.5 d
-----------------------------------
For a planet on a P = 13.8 d circular orbit around a 0.507 Rsun / 0.505 Msun M dwarf:
    a = (M P^2)^(1/3) AU = 0.0897 AU,  a/R* = 38.0,
    T14(b=0, k=0.05) = P/pi * asin(sqrt((1+k)^2 - b^2)/(a/R*)) = 2.91 h.
Allowing eccentricity up to e = 0.3 at apastron lengthens that by sqrt((1+e)/(1-e)) = 1.36 to ~3.97 h,
which is the longest plausible transit for this period/star. The biweight window should be about three
times the transit duration so that the robust location estimate is dominated by out-of-transit points
and does not eat the dip (Hippke et al. 2019, AJ 158, 143, Sect. 4.2 recommends ~3x T14): 3 x 3.97 h
= 11.9 h = 0.50 d. For a blind search extending to 60 d the longest T14 would be ~4.8 h (6.5 h with e = 0.3) and one
would use ~0.8 d; the default here is tuned to the 13.8-d candidate as the task requires.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

AU_IN_RSUN = 215.032  # 1 au / 1 Rsun


def t14_central_hours(period_days: float, r_star: float, m_star: float, k: float = 0.05, b: float = 0.0) -> float:
    """Full transit duration (hours) for a circular orbit, radius ratio k and impact parameter b."""
    a_au = (m_star * (period_days / 365.25) ** 2) ** (1.0 / 3.0)
    a_over_r = a_au * AU_IN_RSUN / r_star
    return period_days / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / a_over_r) * 24.0


def default_window_days(period_ref: float = 13.8, r_star: float = 0.507, m_star: float = 0.505,
                        e_max: float = 0.3, factor: float = 3.0) -> float:
    t14 = t14_central_hours(period_ref, r_star, m_star) * np.sqrt((1 + e_max) / (1 - e_max))
    return factor * t14 / 24.0


@dataclass
class LightCurve:
    time: np.ndarray        # BKJD
    flux: np.ndarray        # normalised (per quarter), before detrending
    flux_err: np.ndarray    # normalised errors
    quarter: np.ndarray
    sap_flux: np.ndarray    # normalised SAP flux (same mask), for dilution checks
    flat: np.ndarray | None = None   # detrended flux (≈1)
    trend: np.ndarray | None = None

    def __len__(self):
        return len(self.time)


def quality_mask(quality: np.ndarray, mode: str = "hard") -> np.ndarray:
    if mode == "hard":
        return quality == 0
    if mode == "default":
        from lightkurve.utils import KeplerQualityFlags
        return (quality & KeplerQualityFlags.DEFAULT_BITMASK) == 0
    if mode == "none":
        return np.ones(len(quality), bool)
    raise ValueError(mode)


def stitch(quarters: Iterable, quality_mode: str = "hard", min_points: int = 300, verbose: bool = True) -> LightCurve:
    """Mask, normalise per quarter and concatenate. Returns a LightCurve sorted in time."""
    T, F, E, Q, S = [], [], [], [], []
    for q in quarters:
        good = quality_mask(q.quality, quality_mode) & np.isfinite(q.time) & np.isfinite(q.flux) \
            & np.isfinite(q.flux_err) & (q.flux_err > 0)
        n = int(good.sum())
        if n < min_points:
            if verbose:
                print(f"[prep] Q{q.quarter}: only {n} good points, skipped")
            continue
        med = np.median(q.flux[good])
        smed = np.nanmedian(q.sap_flux[good])
        T.append(q.time[good]); F.append(q.flux[good] / med); E.append(q.flux_err[good] / med)
        Q.append(np.full(n, q.quarter)); S.append(q.sap_flux[good] / smed if smed > 0 else np.full(n, np.nan))
        if verbose:
            print(f"[prep] Q{q.quarter:2d}: {n:5d}/{len(q.time):5d} cadences kept (quality={quality_mode}), "
                  f"median PSF flux {med:8.1f} e-/s, PSFFRAC {q.header.get('PSFFRAC', float('nan')):.3f}, "
                  f"CROWDSAP {q.header.get('CROWDSAP', float('nan')):.3f}")
    t = np.concatenate(T); o = np.argsort(t)
    return LightCurve(t[o], np.concatenate(F)[o], np.concatenate(E)[o], np.concatenate(Q)[o], np.concatenate(S)[o])


def transit_mask(time, period, t0, duration, factor=1.5, min_width_h=3.0):
    """Boolean mask, True inside +/- max(factor*duration, min_width)/2 of every predicted transit."""
    from wotan import transit_mask as _tm
    width = max(factor * duration, min_width_h / 24.0)
    return _tm(time=time, period=period, duration=width, T0=t0)


def detrend(lc: LightCurve, window_length: float | None = None, break_tolerance: float = 0.5,
            method: str = "biweight", mask: np.ndarray | None = None, verbose: bool = True) -> LightCurve:
    """wotan detrending; writes lc.flat and lc.trend (NaN-free: NaN trend points are dropped).

    mask: optional boolean array (True = in-transit) excluded from the trend estimate (second, "masked"
    pass). Without it a biweight window that contains a 1700-ppm, 2-h transit is pulled down by
    ~(n_in/n_window) x depth ~ 250 ppm, which biases the depth low and creates positive "shoulders"
    around the transit; with the mask the trend is interpolated across the transit from out-of-transit
    points only.
    """
    from wotan import flatten
    if window_length is None:
        window_length = default_window_days()
    flat, trend = flatten(lc.time, lc.flux, method=method, window_length=window_length,
                          break_tolerance=break_tolerance, edge_cutoff=0.0, return_trend=True, mask=mask)
    ok = np.isfinite(flat)
    if verbose:
        print(f"[prep] wotan {method} window={window_length:.3f} d ({window_length*24:.1f} h), "
              f"break_tolerance={break_tolerance} d, masked points={0 if mask is None else int(mask.sum())}; "
              f"{(~ok).sum()} NaN points dropped")
    lc = _subset(lc, ok)
    lc.flat = flat[ok]; lc.trend = trend[ok]
    return lc


def clip_upward(lc: LightCurve, n_sigma: float = 4.0, verbose: bool = True) -> LightCurve:
    """Remove only points ABOVE the trend by more than n_sigma robust sigma. Never clips dips."""
    resid = lc.flat - 1.0
    sig = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    bad = resid > n_sigma * sig
    if verbose:
        print(f"[prep] upward-only clip at +{n_sigma} x {sig*1e6:.0f} ppm: removed {bad.sum()} points "
              f"(no downward clipping)")
    return _subset(lc, ~bad)


def _subset(lc: LightCurve, m: np.ndarray) -> LightCurve:
    return LightCurve(lc.time[m], lc.flux[m], lc.flux_err[m], lc.quarter[m], lc.sap_flux[m],
                      None if lc.flat is None else lc.flat[m], None if lc.trend is None else lc.trend[m])


def noise_metrics(lc: LightCurve, cadence_min: float = 29.4244, verbose: bool = True) -> dict:
    """Point-to-point scatter and a 6-hour CDPP-like metric on the detrended flux (ppm)."""
    f = lc.flat
    resid = f - 1.0
    std = float(np.std(resid))
    mad_sigma = float(1.4826 * np.median(np.abs(resid - np.median(resid))))
    p2p = float(np.std(np.diff(f)) / np.sqrt(2))
    # 6-h CDPP-like: rms of the mean of the flux in sliding 6-h windows (12 x 29.4-min cadences), robust
    n_cad = int(round(6.0 * 60 / cadence_min))
    bins = np.floor((lc.time - lc.time[0]) / (n_cad * cadence_min / 1440.0)).astype(int)
    uniq, inv, cnt = np.unique(bins, return_inverse=True, return_counts=True)
    sums = np.bincount(inv, weights=resid)
    full = cnt >= max(3, n_cad // 2)
    bmeans = sums[full] / cnt[full]
    cdpp6_binned = float(1.4826 * np.median(np.abs(bmeans - np.median(bmeans))))
    out = {"n_points": int(len(f)), "n_quarters": int(len(np.unique(lc.quarter))),
           "quarters": [int(q) for q in np.unique(lc.quarter)],
           "baseline_days": float(np.ptp(lc.time)), "std_ppm": std * 1e6, "mad_sigma_ppm": mad_sigma * 1e6,
           "point_to_point_ppm": p2p * 1e6, "median_flux_err_ppm": float(np.median(lc.flux_err)) * 1e6,
           "cdpp6h_binned_robust_ppm": cdpp6_binned * 1e6, "cdpp6h_n_bins": int(full.sum()),
           "white_noise_6h_expectation_ppm": p2p / np.sqrt(n_cad) * 1e6}
    try:
        import lightkurve as lk
        klc = lk.LightCurve(time=lc.time, flux=f, flux_err=lc.flux_err)
        out["cdpp6h_lightkurve_ppm"] = float(klc.estimate_cdpp(transit_duration=n_cad).value)
    except Exception as ex:  # pragma: no cover
        out["cdpp6h_lightkurve_ppm"] = None
        out["cdpp_error"] = repr(ex)
    if verbose:
        print(f"[prep] N={out['n_points']} points in {out['n_quarters']} quarters, baseline {out['baseline_days']:.1f} d; "
              f"std {out['std_ppm']:.0f} ppm, robust sigma {out['mad_sigma_ppm']:.0f} ppm, point-to-point {out['point_to_point_ppm']:.0f} ppm; "
              f"6-h CDPP-like: binned {out['cdpp6h_binned_robust_ppm']:.0f} ppm, lightkurve {out['cdpp6h_lightkurve_ppm']} ppm, "
              f"white-noise expectation {out['white_noise_6h_expectation_ppm']:.0f} ppm")
    return out


def prepare(quarters, quality_mode="hard", window_length=None, break_tolerance=0.5, n_sigma_up=4.0,
            mask_ephemeris=None, verbose=True):
    """stitch -> detrend -> upward clip. mask_ephemeris=(period, t0, duration) enables the masked pass."""
    lc = stitch(quarters, quality_mode=quality_mode, verbose=verbose)
    mask = None
    if mask_ephemeris is not None:
        P, T0, D = mask_ephemeris
        mask = transit_mask(lc.time, P, T0, D)
    lc = detrend(lc, window_length=window_length, break_tolerance=break_tolerance, mask=mask, verbose=verbose)
    lc = clip_upward(lc, n_sigma=n_sigma_up, verbose=verbose)
    return lc


def invert(lc: LightCurve) -> LightCurve:
    """Null-test copy with flat -> 2 - flat (dips become bumps, noise/systematics statistics preserved)."""
    out = _subset(lc, np.ones(len(lc), bool))
    out.flat = 2.0 - lc.flat
    out.flux = 2.0 - lc.flux
    return out
