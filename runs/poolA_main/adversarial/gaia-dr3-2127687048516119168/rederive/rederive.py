#!/usr/bin/env python
"""Independent re-derivation for gaia-dr3-2127687048516119168.

Written from scratch (numpy/scipy/astropy only; no pipeline/ code, no wotan).
Steps: load KBONUS-BKG stitched LC -> quality==0 -> per-quarter normalise ->
running-median detrend (1 d window) -> blind BLS 1-50 d -> box statistics at best
ephemeris (depth w/ random-phase empirical error, odd/even, secondary, per-year,
n_transits) -> same depth on SAP_FLUX -> inverted-flux null BLS.
"""
import json, sys, time
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FNAME = "gaia-dr3-2127687048516119168"
FITS = f"/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_{FNAME}_kepler_v1.0_lc.fits"
OUT = f"/home/user/ZyushgPlanet/runs/poolA_main/adversarial/{FNAME}/rederive"
RNG = np.random.default_rng(20261017)
WIN = 1.0           # running-median window (days)
PMIN, PMAX = 1.0, 50.0
DURS = np.array([0.04, 0.06, 0.08, 0.10, 0.13, 0.17, 0.22])  # days
PIPE_P, PIPE_DEPTH, PIPE_DUR = 32.86870, 351.0, 1.8


def running_median(t, f, win):
    """Time-windowed running median (own implementation; O(N*w))."""
    half = win / 2.0
    lo = np.searchsorted(t, t - half)
    hi = np.searchsorted(t, t + half, side="right")
    out = np.empty_like(f)
    for i in range(len(t)):
        out[i] = np.median(f[lo[i]:hi[i]])
    return out


def load(col="FLUX"):
    with fits.open(FITS) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d[col], float)
        q = np.asarray(d["QUARTER"], int)
        qual = np.asarray(d["SAP_QUALITY"], int)
    m = np.isfinite(t) & np.isfinite(f) & (qual == 0) & (f > 0)
    t, f, q = t[m], f[m], q[m]
    o = np.argsort(t)
    t, f, q = t[o], f[o], q[o]
    # per-quarter normalise
    fn = np.empty_like(f)
    for qq in np.unique(q):
        s = q == qq
        fn[s] = f[s] / np.nanmedian(f[s])
    # detrend per quarter with running median so quarter edges don't leak
    trend = np.empty_like(fn)
    for qq in np.unique(q):
        s = q == qq
        trend[s] = running_median(t[s], fn[s], WIN)
    fd = fn / trend
    # sigma clip high outliers only (keep dips), 4 sigma via MAD
    mad = 1.4826 * np.median(np.abs(fd - np.median(fd)))
    keep = fd < 1 + 4 * mad
    return t[keep], fd[keep], q[keep], mad


def period_grid(t, pmin, pmax, dur_min, os_factor=3.0):
    T = t.max() - t.min()
    periods = [pmin]
    while periods[-1] < pmax:
        p = periods[-1]
        periods.append(p + p * p * dur_min / (os_factor * T))
    return np.array(periods)


def sde(power):
    return (power.max() - np.median(power)) / np.std(power)


def run_bls(t, f, periods, durs, label):
    t0 = time.time()
    bls = BoxLeastSquares(t, f)
    res = bls.power(periods, durs, objective="snr")
    pw = np.asarray(res.power)
    # detrend power spectrum trend vs log P with a running median to suppress red noise in the periodogram
    k = max(51, len(pw) // 200)
    k += (k % 2 == 0)
    from scipy.ndimage import median_filter
    pw_d = pw - median_filter(pw, size=k, mode="nearest")
    # top peaks (distinct, separated by >2% in period)
    order = np.argsort(pw_d)[::-1]
    peaks = []
    for i in order:
        p = periods[i]
        if all(abs(p - pp["period"]) / pp["period"] > 0.02 for pp in peaks):
            peaks.append(dict(period=float(p), power=float(pw[i]), power_detrended=float(pw_d[i]),
                              t0=float(res.transit_time[i]), duration=float(res.duration[i]),
                              depth=float(res.depth[i]), depth_snr=float(res.depth_snr[i])))
        if len(peaks) >= 5:
            break
    out = dict(label=label, n_points=int(len(t)), n_periods=int(len(periods)),
               sde_raw=float(sde(pw)), sde_detrended=float((pw_d.max() - np.median(pw_d)) / np.std(pw_d)),
               top_peaks=peaks, runtime_s=time.time() - t0)
    print(f"[{label}] {len(periods)} periods, {time.time()-t0:.0f}s, SDE raw={out['sde_raw']:.2f} det={out['sde_detrended']:.2f}")
    for pk in peaks:
        print(f"   P={pk['period']:.5f} pow={pk['power']:.2f} powd={pk['power_detrended']:.2f} dur={pk['duration']*24:.2f}h depth={pk['depth']*1e6:.0f}ppm t0={pk['t0']:.4f}")
    return out, pw, pw_d


def refine(t, f, p_guess, label):
    """Fine grid around a peak: period, t0, duration via BLS."""
    bls = BoxLeastSquares(t, f)
    T = t.max() - t.min()
    dp = p_guess ** 2 * 0.02 / (10 * T)
    periods = np.arange(p_guess * 0.995, p_guess * 1.005, dp)
    durs = np.linspace(0.03, 0.25, 23)
    res = bls.power(periods, durs, objective="snr")
    i = np.argmax(res.power)
    P, t0, dur = float(res.period[i]), float(res.transit_time[i]), float(res.duration[i])
    # stats
    st = bls.compute_stats(P, dur, t0)
    print(f"[{label} refine] P={P:.6f} t0={t0:.5f} dur={dur*24:.2f}h depth={res.depth[i]*1e6:.0f}+-{res.depth_err[i]*1e6:.0f} ppm snr={res.depth_snr[i]:.2f}")
    return P, t0, dur, st


def box_depth(t, f, P, t0, dur, phase_offset=0.0):
    ph = ((t - t0) / P + 0.5 - phase_offset) % 1.0 - 0.5
    inn = np.abs(ph) < 0.5 * dur / P
    # local out-of-transit: within 3 durations but not within 1.5 durations
    oot = (np.abs(ph) < 3 * dur / P) & (np.abs(ph) > 0.75 * dur / P)
    if inn.sum() < 3:
        return np.nan, np.nan, int(inn.sum())
    d = np.median(f[oot]) - np.mean(f[inn])
    err = np.std(f[oot]) / np.sqrt(inn.sum())
    return float(d), float(err), int(inn.sum())


def empirical_depth(t, f, P, t0, dur, nrand=400):
    d, e_formal, n = box_depth(t, f, P, t0, dur)
    rand = []
    for _ in range(nrand):
        off = RNG.uniform(0.1, 0.9)  # avoid primary within +-0.1 phase
        dd, _, nn = box_depth(t, f, P, t0, dur, off)
        if np.isfinite(dd) and nn >= 0.5 * n:
            rand.append(dd)
    rand = np.array(rand)
    return d, float(np.std(rand)), e_formal, n, float(np.mean(rand))


def per_event(t, f, P, t0, dur):
    n0 = int(np.floor((t.min() - t0) / P)) - 1
    n1 = int(np.ceil((t.max() - t0) / P)) + 1
    evs = []
    for n in range(n0, n1 + 1):
        tc = t0 + n * P
        inn = np.abs(t - tc) < 0.5 * dur
        oot = (np.abs(t - tc) < 3 * dur) & (np.abs(t - tc) > 0.75 * dur)
        if inn.sum() >= 2 and oot.sum() >= 10:
            d = np.median(f[oot]) - np.mean(f[inn])
            e = np.std(f[oot]) / np.sqrt(inn.sum())
            evs.append(dict(n=n, tc=float(tc), n_in=int(inn.sum()), n_oot=int(oot.sum()), depth_ppm=float(d * 1e6), err_ppm=float(e * 1e6)))
    return evs


def main():
    results = dict(fname=FNAME, date="2026-10-17", method=dict(
        quality="SAP_QUALITY==0 & finite & flux>0", normalise="per-quarter median",
        detrend=f"per-quarter time-windowed running median, window={WIN} d (own code)",
        clip="upper 4*MAD only", bls="astropy BoxLeastSquares objective=snr, own period grid",
        grid=f"P {PMIN}-{PMAX} d, dP = P^2*dur_min/(3*T), durations {list(DURS)} d",
        depth="mean in-box vs median of local (0.75-3 dur) out-of-transit; empirical error = std of 400 random-phase boxes"))

    t, f, q, mad = load("FLUX")
    print(f"PSF: {len(t)} pts, baseline {t.min():.1f}-{t.max():.1f}, MAD scatter {mad*1e6:.0f} ppm")
    results["psf_n_points"] = int(len(t)); results["psf_point_scatter_ppm"] = float(mad * 1e6)

    periods = period_grid(t, PMIN, PMAX, DURS.min())
    blind, pw, pw_d = run_bls(t, f, periods, DURS, "PSF blind")
    results["blind_bls"] = blind
    top = blind["top_peaks"][0]["period"]
    # is the pipeline period or an alias in the top 3?
    def is_alias(p, ref):
        for k in [1, 2, 3, 0.5, 1 / 3]:
            if abs(p / (ref * k) - 1) < 0.005:
                return k
        return None
    results["pipeline_period_in_top3"] = [dict(period=pk["period"], alias_factor=is_alias(pk["period"], PIPE_P)) for pk in blind["top_peaks"][:3]]
    recovered = any(is_alias(pk["period"], PIPE_P) == 1 for pk in blind["top_peaks"][:1])
    results["recovered_blind_top1"] = bool(recovered)
    results["recovered_blind_top3_any_alias"] = bool(any(d["alias_factor"] for d in results["pipeline_period_in_top3"]))

    # inverted null
    inv, _, _ = run_bls(t, 2 - f, periods, DURS, "PSF inverted null")
    results["inverted_null_bls"] = inv

    # pick ephemeris: best peak if it matches pipeline, else also evaluate pipeline period
    P, t0, dur, st = refine(t, f, top, "PSF top peak")
    results["best"] = dict(period_d=P, t0_bkjd=t0, t0_bjd=t0 + 2454833.0, duration_h=dur * 24)
    d, e_emp, e_form, n_in, rand_mean = empirical_depth(t, f, P, t0, dur)
    results["best"].update(depth_ppm=d * 1e6, depth_err_empirical_ppm=e_emp * 1e6, depth_err_formal_ppm=e_form * 1e6,
                           n_in_transit_points=n_in, random_box_mean_ppm=rand_mean * 1e6,
                           empirical_snr=d / e_emp if e_emp > 0 else np.nan)
    # odd/even
    ph_n = np.round((t - t0) / P).astype(int)
    odd = ph_n % 2 == 1
    do, eo, no_ = box_depth(t[odd], f[odd], P, t0, dur)
    de, ee, ne = box_depth(t[~odd], f[~odd], P, t0, dur)
    # empirical errors for odd/even subsets
    _, eo_emp, _, _, _ = empirical_depth(t[odd], f[odd], P, t0, dur, 200)
    _, ee_emp, _, _, _ = empirical_depth(t[~odd], f[~odd], P, t0, dur, 200)
    results["best"]["odd_even"] = dict(odd_ppm=do * 1e6, odd_err_ppm=eo_emp * 1e6, even_ppm=de * 1e6, even_err_ppm=ee_emp * 1e6,
                                       sigma=float(abs(do - de) / np.hypot(eo_emp, ee_emp)))
    # secondary at phase 0.5
    ds, _, ns = box_depth(t, f, P, t0, dur, 0.5)
    results["best"]["secondary"] = dict(depth_ppm=ds * 1e6, err_ppm=e_emp * 1e6, sigma=float(ds / e_emp), n_points=ns)
    # scan for strongest secondary at any phase in 0.2-0.8 (same box)
    phs = np.linspace(0.15, 0.85, 141)
    sec = np.array([box_depth(t, f, P, t0, dur, p)[0] for p in phs])
    i = np.nanargmax(sec)
    results["best"]["max_other_phase_dip"] = dict(phase=float(phs[i]), depth_ppm=float(sec[i] * 1e6), sigma=float(sec[i] / e_emp))
    # per year
    years = {"Y1(Q1-4)": (0, 540), "Y2(Q5-8)": (540, 906), "Y3(Q9-12)": (906, 1274), "Y4(Q13-17)": (1274, 1700)}
    py = {}
    for k, (a, b) in years.items():
        s = (t >= a) & (t < b)
        dd, ee_, nn = box_depth(t[s], f[s], P, t0, dur)
        _, ee_emp2, _, _, _ = empirical_depth(t[s], f[s], P, t0, dur, 150)
        py[k] = dict(depth_ppm=dd * 1e6, err_ppm=ee_emp2 * 1e6, n_in=nn)
    results["best"]["per_year"] = py
    evs = per_event(t, f, P, t0, dur)
    results["best"]["events"] = evs
    results["best"]["n_transits_with_data"] = len(evs)
    results["best"]["n_transits_expected_in_baseline"] = int(np.floor((t.max() - t.min()) / P)) + 1
    # per-event: how many are individually >1 sigma deep; single-event dominance
    dep = np.array([e["depth_ppm"] for e in evs]); err = np.array([e["err_ppm"] for e in evs])
    w = 1 / err ** 2
    results["best"]["weighted_mean_event_depth_ppm"] = float(np.sum(w * dep) / np.sum(w))
    results["best"]["max_single_event_fraction_of_chi2"] = float(np.max(w * dep ** 2) / np.sum(w * dep ** 2))
    results["best"]["events_positive_fraction"] = float(np.mean(dep > 0))

    # V-shape proxy: compare depth with half-duration central box vs full box
    dc, _, _ = box_depth(t, f, P, t0, dur * 0.5)
    results["best"]["central_half_box_depth_ppm"] = dc * 1e6

    # Pipeline ephemeris check (if different)
    if abs(P / PIPE_P - 1) > 0.002:
        Pp, t0p, durp, _ = refine(t, f, PIPE_P, "PSF pipeline period")
        dp_, ep_, _, np_, _ = empirical_depth(t, f, Pp, t0p, durp)
        results["at_pipeline_period"] = dict(period_d=Pp, t0_bkjd=t0p, duration_h=durp * 24, depth_ppm=dp_ * 1e6, depth_err_ppm=ep_ * 1e6, snr=dp_ / ep_)

    # SAP flux
    ts, fs, qs, mads = load("SAP_FLUX")
    print(f"SAP: {len(ts)} pts, MAD scatter {mads*1e6:.0f} ppm")
    dsap, esap, _, nsap, _ = empirical_depth(ts, fs, P, t0, dur)
    # also allow SAP its own small t0 refinement? keep same ephemeris (fair comparison)
    results["sap"] = dict(n_points=int(len(ts)), point_scatter_ppm=mads * 1e6, depth_ppm=dsap * 1e6, depth_err_ppm=esap * 1e6,
                          snr=dsap / esap, n_in=nsap, ratio_sap_over_psf=dsap / d if d else np.nan)
    # SAP blind BLS at narrow window around P for comparison
    bls_s = BoxLeastSquares(ts, fs)
    rs = bls_s.power(np.linspace(P * 0.99, P * 1.01, 400), np.linspace(0.04, 0.2, 9), objective="snr")
    j = np.argmax(rs.power)
    results["sap"]["bls_local_peak"] = dict(period=float(rs.period[j]), depth_ppm=float(rs.depth[j] * 1e6), snr=float(rs.depth_snr[j]), duration_h=float(rs.duration[j] * 24))

    # ---- figure ----
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    ax = ax.ravel()
    ax[0].plot(periods, pw_d, lw=0.5, color="0.3")
    ax[0].axvline(PIPE_P, color="C3", ls="--", lw=1, label=f"pipeline P={PIPE_P}")
    ax[0].set_xscale("log"); ax[0].set_xlabel("Period [d]"); ax[0].set_ylabel("BLS SNR power (median-detrended)")
    ax[0].set_title(f"Blind BLS 1-50 d (own grid), SDE={blind['sde_detrended']:.1f}"); ax[0].legend()
    ph = ((t - t0) / P + 0.5) % 1.0 - 0.5
    x = ph * P * 24
    s = np.abs(x) < 12
    ax[1].plot(x[s], (f[s] - 1) * 1e6, ".", ms=2, color="0.7", alpha=0.6)
    bins = np.arange(-12, 12.01, 0.5)
    idx = np.digitize(x[s], bins)
    bx = [0.5 * (bins[i - 1] + bins[i]) for i in range(1, len(bins))]
    by = [np.mean((f[s] - 1)[idx == i]) * 1e6 if np.any(idx == i) else np.nan for i in range(1, len(bins))]
    ax[1].plot(bx, by, "o-", color="C0", ms=4)
    ax[1].axvspan(-dur * 12, dur * 12, color="C3", alpha=0.1)
    ax[1].set_ylim(-3000, 2000); ax[1].set_xlabel("hours from mid-transit"); ax[1].set_ylabel("ppm")
    ax[1].set_title(f"PSF fold P={P:.5f} d, depth={d*1e6:.0f}+-{e_emp*1e6:.0f} ppm, dur={dur*24:.2f} h")
    # odd/even
    for sel, c, lab in [(odd, "C1", "odd"), (~odd, "C2", "even")]:
        xx = x[sel & s]; yy = (f[sel & s] - 1) * 1e6
        idx = np.digitize(xx, bins)
        by = [np.mean(yy[idx == i]) if np.any(idx == i) else np.nan for i in range(1, len(bins))]
        ax[2].plot(bx, by, "o-", color=c, ms=4, label=lab)
    ax[2].axvspan(-dur * 12, dur * 12, color="C3", alpha=0.1); ax[2].legend(); ax[2].set_ylim(-1500, 1000)
    ax[2].set_title(f"odd {do*1e6:.0f}+-{eo_emp*1e6:.0f} / even {de*1e6:.0f}+-{ee_emp*1e6:.0f} ppm ({results['best']['odd_even']['sigma']:.1f} sigma)")
    ax[2].set_xlabel("hours from mid-transit")
    # secondary
    ph2 = ((t - t0) / P) % 1.0
    xb = np.linspace(0, 1, 201)
    idx = np.digitize(ph2, xb)
    by = [np.mean((f - 1)[idx == i]) * 1e6 if np.any(idx == i) else np.nan for i in range(1, len(xb))]
    ax[3].plot(0.5 * (xb[1:] + xb[:-1]), by, "-", color="0.3")
    ax[3].axvline(0.5, color="C3", ls="--"); ax[3].set_ylim(-600, 400)
    ax[3].set_title(f"full phase (200 bins); phase-0.5 box depth {ds*1e6:.0f}+-{e_emp*1e6:.0f} ppm"); ax[3].set_xlabel("phase")
    # per event
    ax[4].errorbar([e["tc"] for e in evs], dep, err, fmt="o", ms=4, color="C0")
    ax[4].axhline(d * 1e6, color="C3", ls="--"); ax[4].axhline(0, color="k", lw=0.5)
    ax[4].set_xlabel("event mid-time [BKJD]"); ax[4].set_ylabel("event depth [ppm]"); ax[4].set_title(f"{len(evs)} events with data")
    # SAP fold
    phs_ = ((ts - t0) / P + 0.5) % 1.0 - 0.5; xs = phs_ * P * 24; ss = np.abs(xs) < 12
    idx = np.digitize(xs[ss], bins)
    bys = [np.mean((fs[ss] - 1)[idx == i]) * 1e6 if np.any(idx == i) else np.nan for i in range(1, len(bins))]
    ax[5].plot(xs[ss], (fs[ss] - 1) * 1e6, ".", ms=2, color="0.7", alpha=0.6)
    ax[5].plot(bx, bys, "o-", color="C4", ms=4); ax[5].set_ylim(-3000, 2000)
    ax[5].axvspan(-dur * 12, dur * 12, color="C3", alpha=0.1)
    ax[5].set_title(f"SAP_FLUX fold, depth={dsap*1e6:.0f}+-{esap*1e6:.0f} ppm (ratio {dsap/d:.2f})"); ax[5].set_xlabel("hours from mid-transit")
    fig.suptitle(f"{FNAME} independent re-derivation (2026-10-17)")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fold_rederive.png", dpi=110)

    def clean(o):
        if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [clean(v) for v in o]
        if isinstance(o, (np.floating, float)): return None if not np.isfinite(o) else float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.ndarray): return clean(o.tolist())
        return o
    with open(f"{OUT}/rederive_results.json", "w") as fh:
        json.dump(clean(results), fh, indent=2)
    print(json.dumps(clean(results["best"]), indent=1)[:4000])
    print("SAP:", json.dumps(clean(results["sap"]), indent=1))
    if "at_pipeline_period" in results: print("AT PIPE:", json.dumps(clean(results["at_pipeline_period"]), indent=1))


if __name__ == "__main__":
    main()
