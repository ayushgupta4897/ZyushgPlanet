#!/usr/bin/env python
"""
Independent re-derivation for kic-006515744 (KBONUS-BKG PSF light curve).

Written from scratch with numpy / scipy / astropy only. No code under pipeline/ is
imported or read. Nothing here is a planet; this is a check of a first-screen signal.

Steps
 1. Load LIGHTCURVE_STITCHED, keep SAP_QUALITY == 0 and finite cadences, normalise per
    quarter, detrend with a 1-d time-windowed running median (own implementation),
    clip only *positive* outliers.
 2. Blind BLS 1-50 d (astropy BoxLeastSquares, own frequency grid), weighted by a
    per-quarter robust scatter. Same search on the inverted light curve (2 - f) as a null.
 3. At the best ephemeris: refined P/t0/duration, box depth with empirical
    (random-phase) error, odd/even, secondary at phase 0.5, per-year depths, per-event
    depths, number of transits with data.
 4. Same depth measurement on SAP_FLUX (dilution check).
 5. Writes results JSON + a fold figure.  Comparison with pipeline/vet.json is done
    afterwards by hand (not in this script).
"""
import json
import os
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import BoxLeastSquares

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-006515744_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006515744/rederive/"
FNAME = "kic-006515744"
PIPE_P = 36.20750          # pipeline claim, used ONLY to flag aliases in the blind spectrum
RNG = np.random.default_rng(20261017)

os.makedirs(OUT, exist_ok=True)


# --------------------------------------------------------------------------------------
# 1. loader + detrend
# --------------------------------------------------------------------------------------
def load(col="FLUX"):
    with fits.open(FITS) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d[col], float)
        fe = np.asarray(d[col + "_ERR"], float)
        q = np.asarray(d["QUARTER"], int)
        qual = np.asarray(d["SAP_QUALITY"], int)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (f > 0)
    t, f, fe, q = t[m], f[m], fe[m], q[m]
    o = np.argsort(t)
    t, f, fe, q = t[o], f[o], fe[o], q[o]
    fn = np.empty_like(f)
    fen = np.empty_like(fe)
    for qq in np.unique(q):
        s = q == qq
        med = np.median(f[s])
        fn[s] = f[s] / med
        fen[s] = fe[s] / med
    return t, fn, fen, q


def running_median_time(t, y, window=1.0):
    """Median of y within +-window/2 in time around each point (own implementation).
    Called per quarter so the window never bridges inter-quarter gaps."""
    hw = window / 2.0
    lo = np.searchsorted(t, t - hw, side="left")
    hi = np.searchsorted(t, t + hw, side="right")
    out = np.empty_like(y)
    for i in range(len(t)):
        out[i] = np.median(y[lo[i]:hi[i]])
    return out


def detrend(t, y, q, window=1.0, n_iter=2, clip_hi=4.0):
    """Per-quarter running-median detrend. Positive outliers (> clip_hi sigma above the
    trend) are masked from the trend estimate and from the output; negative excursions
    are never clipped (they could be transits). Returns detrended flux, mask of kept
    points, and a per-quarter robust scatter used as the BLS weight."""
    trend = np.ones_like(y)
    keep = np.ones(len(y), bool)
    for _ in range(n_iter):
        for qq in np.unique(q):
            s = np.where((q == qq) & keep)[0]
            tr = running_median_time(t[s], y[s], window)
            # interpolate trend onto clipped points too, so they can be re-evaluated
            s_all = np.where(q == qq)[0]
            trend[s_all] = np.interp(t[s_all], t[s], tr)
        resid = y / trend - 1.0
        sig = np.zeros_like(y)
        for qq in np.unique(q):
            s = q == qq
            mad = np.median(np.abs(resid[s & keep] - np.median(resid[s & keep])))
            sig[s] = 1.4826 * mad
        keep = resid < clip_hi * sig            # clip only upward outliers
    yd = y / trend
    sig_q = np.zeros_like(y)
    for qq in np.unique(q):
        s = q == qq
        r = yd[s & keep] - 1.0
        sig_q[s] = 1.4826 * np.median(np.abs(r - np.median(r)))
    return yd, keep, sig_q


# --------------------------------------------------------------------------------------
# 2. blind BLS
# --------------------------------------------------------------------------------------
def period_grid(pmin, pmax, baseline, dur_ref=0.045, oversample=3.0):
    """Uniform-in-phase-drift grid: df = dur_ref * f / (oversample * baseline), i.e. the
    accumulated transit-time error over the full baseline is <= dur_ref/oversample."""
    fmin, fmax = 1.0 / pmax, 1.0 / pmin
    freqs = [fmax]
    f = fmax
    while f > fmin:
        f -= dur_ref * f / (oversample * baseline)
        freqs.append(f)
    return 1.0 / np.array(freqs[::-1])


def sde(power):
    p = np.asarray(power, float)
    return float((np.nanmax(p) - np.nanmean(p)) / np.nanstd(p))


def flatten_spectrum(periods, power, nwin=2001):
    """Remove the slow period-dependent floor of the BLS statistic with a running
    median over the (log-spaced) period grid, so SDE is not dominated by the long-
    period rise. Returns flattened power (power - running median)."""
    n = len(power)
    half = nwin // 2
    # pad with edge values, then compute medians over strided windows
    padded = np.concatenate([np.full(half, power[0]), power, np.full(half, power[-1])])
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(padded, nwin), axis=-1)
    return power - med


def top_peaks(periods, power, n=3, excl_frac=0.03):
    """Distinct local maxima, excluding neighbours within excl_frac in period of a
    stronger peak."""
    order = np.argsort(power)[::-1]
    peaks = []
    for i in order:
        p = periods[i]
        if all(abs(p / pk["period"] - 1.0) > excl_frac for pk in peaks):
            peaks.append({"period": float(p), "power": float(power[i]), "index": int(i)})
        if len(peaks) >= n:
            break
    return peaks


def alias_tag(p, ref=PIPE_P, tol=0.01):
    for num in range(1, 7):
        for den in range(1, 7):
            if abs(p / ref - num / den) < tol * (num / den):
                return f"{num}/{den} x {ref}"
    return None


def blind_bls(t, y, dy, pmin=1.0, pmax=50.0, label=""):
    baseline = t.max() - t.min()
    periods = period_grid(pmin, pmax, baseline)
    durs = np.array([0.04, 0.055, 0.075, 0.10, 0.13, 0.17, 0.22, 0.30])
    bls = BoxLeastSquares(t, y, dy=dy)
    t0 = time.time()
    res = bls.power(periods, durs, objective="likelihood")
    dt = time.time() - t0
    pw = np.asarray(res.power, float)
    pw_flat = flatten_spectrum(periods, pw)
    snr = np.asarray(res.depth_snr, float)
    out = {
        "label": label,
        "n_periods": int(len(periods)),
        "durations_d": durs.tolist(),
        "runtime_s": round(dt, 1),
        "sde_loglike_raw": sde(pw),
        "sde_loglike_flattened": sde(pw_flat),
        "sde_depth_snr": sde(snr),
        "top3_loglike_raw": top_peaks(periods, pw),
        "top3_loglike_flattened": top_peaks(periods, pw_flat),
        "top3_depth_snr": top_peaks(periods, snr),
    }
    for key in ["top3_loglike_raw", "top3_loglike_flattened", "top3_depth_snr"]:
        for pk in out[key]:
            i = pk["index"]
            pk["duration_d"] = float(res.duration[i])
            pk["t0_bkjd"] = float(res.transit_time[i])
            pk["depth_ppm"] = float(res.depth[i] * 1e6)
            pk["depth_snr"] = float(res.depth_snr[i])
            pk["alias_of_pipeline_P"] = alias_tag(pk["period"])
    return out, periods, pw, pw_flat, snr, res


# --------------------------------------------------------------------------------------
# 3. ephemeris statistics
# --------------------------------------------------------------------------------------
def phase_days(t, P, t0):
    return (t - t0 + 0.5 * P) % P - 0.5 * P


def box_depth(t, y, w, P, t0, dur, phase_offset=0.0, sel=None, local_win=3.0):
    """Weighted mean deficit of points inside a box of width dur centred at phase
    offset `phase_offset` (days), relative to the weighted mean of out-of-box points
    within +- local_win*dur of the box centre. Returns (depth, n_in, n_out)."""
    if sel is None:
        sel = np.ones(len(t), bool)
    ph = phase_days(t - phase_offset, P, t0)
    inb = sel & (np.abs(ph) < dur / 2)
    outb = sel & (np.abs(ph) >= dur / 2) & (np.abs(ph) < (local_win + 0.5) * dur)
    if inb.sum() == 0 or outb.sum() < 5:
        return np.nan, int(inb.sum()), int(outb.sum())
    m_in = np.sum(w[inb] * y[inb]) / np.sum(w[inb])
    m_out = np.sum(w[outb] * y[outb]) / np.sum(w[outb])
    return float(m_out - m_in), int(inb.sum()), int(outb.sum())


def empirical_depth(t, y, w, P, t0, dur, sel=None, n_rand=3000, avoid=None):
    """Depth at phase 0 plus an empirical error: the scatter (1.4826*MAD and std) of the
    same statistic evaluated at random phases that avoid the primary (and optional
    secondary) windows."""
    d0, n_in, n_out = box_depth(t, y, w, P, t0, dur, 0.0, sel)
    offs = RNG.uniform(-P / 2, P / 2, size=n_rand * 2)
    bad = np.abs(offs) < 1.5 * dur
    if avoid is not None:
        for a in avoid:
            bad |= np.abs(((offs - a + P / 2) % P) - P / 2) < 1.5 * dur
    offs = offs[~bad][:n_rand]
    vals = np.array([box_depth(t, y, w, P, t0, dur, o, sel)[0] for o in offs])
    vals = vals[np.isfinite(vals)]
    mad = 1.4826 * np.median(np.abs(vals - np.median(vals)))
    sd = np.std(vals)
    # p-value-like: fraction of random boxes at least as deep
    frac = float(np.mean(vals >= d0)) if np.isfinite(d0) else np.nan
    return {
        "depth_ppm": d0 * 1e6 if np.isfinite(d0) else None,
        "err_mad_ppm": mad * 1e6,
        "err_std_ppm": sd * 1e6,
        "snr_mad": d0 / mad if np.isfinite(d0) and mad > 0 else None,
        "snr_std": d0 / sd if np.isfinite(d0) and sd > 0 else None,
        "n_in": n_in,
        "n_out_local": n_out,
        "n_random": int(len(vals)),
        "frac_random_deeper": frac,
    }


def per_event(t, y, w, P, t0, dur):
    """Per-epoch local depth: in-box weighted mean vs weighted mean of +-3 dur either
    side of the box within the same epoch."""
    ep = np.round((t - t0) / P).astype(int)
    rows = []
    for e in np.unique(ep):
        s = ep == e
        tc = t0 + e * P
        dtc = t[s] - tc
        inb = np.abs(dtc) < dur / 2
        outb = (np.abs(dtc) >= dur / 2) & (np.abs(dtc) < 3.5 * dur)
        if inb.sum() == 0:
            continue
        yy, ww = y[s], w[s]
        m_in = np.sum(ww[inb] * yy[inb]) / np.sum(ww[inb])
        if outb.sum() >= 4:
            m_out = np.sum(ww[outb] * yy[outb]) / np.sum(ww[outb])
            loc_sig = np.std(yy[outb])
        else:
            m_out, loc_sig = np.nan, np.nan
        # distance to nearest data gap (> 0.5 d without data) from the event centre
        tt = t[s]
        gaps = np.where(np.diff(t) > 0.5)[0]
        edges = np.concatenate([[t.min()], t[gaps], t[gaps + 1], [t.max()]])
        gap_dist = float(np.min(np.abs(edges - tc)))
        rows.append({
            "epoch": int(e),
            "tc_bkjd": float(tc),
            "n_in": int(inb.sum()),
            "n_out_local": int(outb.sum()),
            "depth_ppm": float((m_out - m_in) * 1e6) if np.isfinite(m_out) else None,
            "depth_err_formal_ppm": float(loc_sig / np.sqrt(inb.sum()) * 1e6) if np.isfinite(loc_sig) else None,
            "min_in_flux_ppm": float((1 - yy[inb].min()) * 1e6),
            "gap_dist_d": gap_dist,
            "quarter": int(np.round(np.median(QUARTER_LOOKUP[s]))),
        })
    return rows


def refine(t, y, dy, P0, span=0.004, n=4001):
    """Fine BLS around P0 with a dense duration set; returns best P, t0, dur, depth."""
    periods = np.linspace(P0 * (1 - span), P0 * (1 + span), n)
    durs = np.arange(0.02, 0.42, 0.01)
    bls = BoxLeastSquares(t, y, dy=dy)
    res = bls.power(periods, durs, objective="likelihood", oversample=20)
    i = int(np.argmax(res.power))
    return (float(res.period[i]), float(res.transit_time[i]), float(res.duration[i]),
            float(res.depth[i]), float(res.depth_err[i]), float(res.depth_snr[i]))


def duration_scan(t, y, dy, P, t0):
    """Depth SNR vs duration at fixed P (t0 free within +-0.1 d)."""
    durs = np.arange(0.02, 0.5, 0.01)
    bls = BoxLeastSquares(t, y, dy=dy)
    res = bls.power([P], durs, objective="snr", oversample=20)
    # astropy returns best duration only; scan manually with compute_stats
    out = []
    for d in durs:
        st = bls.compute_stats(P, d, t0)
        out.append((float(d), float(st["depth"][0] * 1e6), float(st["depth"][0] / st["depth"][1])))
    return out


def year_of(bkjd):
    return Time(bkjd + 2454833.0, format="jd").to_datetime().year


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
if __name__ == "__main__":
    results = {"fname": FNAME, "date": "2026-10-17",
               "method": "independent numpy/astropy re-derivation; no pipeline/ code used",
               "fits": FITS}

    # ---- PSF flux ----
    t, y_raw, ye, q = load("FLUX")
    y, keep, sig_q = detrend(t, y_raw, q, window=1.0)
    QUARTER_LOOKUP = q
    results["n_cadences_quality0"] = int(len(t))
    results["n_after_positive_clip"] = int(keep.sum())
    results["baseline_d"] = float(t.max() - t.min())
    results["quarters"] = [int(x) for x in np.unique(q)]
    results["per_quarter_scatter_ppm"] = {int(qq): float(1e6 * sig_q[q == qq][0]) for qq in np.unique(q)}
    tk, yk, dyk, qk = t[keep], y[keep], sig_q[keep], q[keep]
    wk = 1.0 / dyk ** 2

    # ---- blind BLS, real ----
    blind, periods, pw, pw_flat, snr, res = blind_bls(tk, yk, dyk, 1.0, 50.0, "PSF flux, real")
    results["blind_bls_real"] = blind
    # ---- blind BLS, inverted null ----
    blind_inv, _, pw_inv, pw_inv_flat, snr_inv, _ = blind_bls(tk, 2.0 - yk, dyk, 1.0, 50.0, "PSF flux, inverted (2-f)")
    results["blind_bls_inverted"] = blind_inv

    # which peak to follow: top raw-loglike peak
    P_top = blind["top3_loglike_raw"][0]["period"]
    results["top_peak_is_pipeline_period_or_alias"] = (
        abs(P_top / PIPE_P - 1) < 0.01 or alias_tag(P_top) is not None)
    results["top_peak_within_1pct_of_pipeline_P"] = bool(abs(P_top / PIPE_P - 1) < 0.01)

    # ---- refine around the top peak (and, separately, around the pipeline P if different) ----
    P, t0, dur, dep, dep_err, dep_snr = refine(tk, yk, dyk, P_top)
    results["refined"] = {"period_d": P, "t0_bkjd": t0, "t0_bjd": t0 + 2454833.0,
                          "duration_d": dur, "duration_h": dur * 24,
                          "bls_depth_ppm": dep * 1e6, "bls_depth_err_ppm": dep_err * 1e6,
                          "bls_depth_snr_formal": dep_snr}
    if abs(P_top / PIPE_P - 1) > 0.01:
        Pp, t0p, durp, depp, depp_err, snrp = refine(tk, yk, dyk, PIPE_P)
        results["refined_at_pipeline_period"] = {
            "period_d": Pp, "t0_bkjd": t0p, "duration_d": durp, "bls_depth_ppm": depp * 1e6,
            "bls_depth_snr_formal": snrp}

    # ---- duration scan ----
    results["duration_scan_P_fixed"] = [{"dur_d": a, "depth_ppm": b, "snr_formal": c}
                                        for a, b, c in duration_scan(tk, yk, dyk, P, t0)]

    # ---- box depth with empirical error ----
    prim = empirical_depth(tk, yk, wk, P, t0, dur, avoid=[P / 2])
    results["primary_box"] = prim

    # ---- secondary at phase 0.5 ----
    sec0, _, _ = box_depth(tk, yk, wk, P, t0, dur, P / 2)
    sec = empirical_depth(tk, yk, wk, P, t0 + P / 2, dur, avoid=[P / 2])
    results["secondary_phase05"] = sec
    # also scan the full phase for the deepest "other" box (secondary at any phase)
    offs = np.linspace(-P / 2, P / 2, 2000)
    sc = np.array([box_depth(tk, yk, wk, P, t0, dur, o)[0] for o in offs])
    ok = np.abs(offs) > 1.5 * dur
    j = np.nanargmax(np.where(ok, sc, -np.inf))
    results["deepest_other_phase_box"] = {"phase_offset_d": float(offs[j]), "phase": float(offs[j] / P),
                                          "depth_ppm": float(sc[j] * 1e6),
                                          "sigma_vs_random": float(sc[j] / (prim["err_mad_ppm"] * 1e-6))}

    # ---- odd / even ----
    ep = np.round((tk - t0) / P).astype(int)
    odd = empirical_depth(tk, yk, wk, P, t0, dur, sel=(ep % 2 == 1), n_rand=1500, avoid=[P / 2])
    even = empirical_depth(tk, yk, wk, P, t0, dur, sel=(ep % 2 == 0), n_rand=1500, avoid=[P / 2])
    results["odd"] = odd
    results["even"] = even
    if odd["depth_ppm"] is not None and even["depth_ppm"] is not None:
        diff = odd["depth_ppm"] - even["depth_ppm"]
        err = np.hypot(odd["err_mad_ppm"], even["err_mad_ppm"])
        results["odd_even_sigma"] = float(abs(diff) / err)
        results["odd_even_diff_ppm"] = float(diff)

    # ---- per-year ----
    yrs = np.array([year_of(x) for x in tk])
    per_year = {}
    for yr in np.unique(yrs):
        s = yrs == yr
        r = empirical_depth(tk, yk, wk, P, t0, dur, sel=s, n_rand=1000, avoid=[P / 2])
        r["n_events_with_data"] = int(len(np.unique(ep[s & (np.abs(phase_days(tk, P, t0)) < dur / 2)])))
        per_year[int(yr)] = r
    results["per_year"] = per_year

    # ---- per event ----
    ev = per_event(tk, yk, wk, P, t0, dur)
    results["events"] = ev
    results["n_transits_with_data"] = int(len(ev))
    results["n_transits_ge2pts"] = int(sum(1 for e in ev if e["n_in"] >= 2))
    results["n_epochs_in_baseline"] = int(np.floor((tk.max() - t0) / P) - np.ceil((tk.min() - t0) / P) + 1)
    # single-event dominance: share of the summed weighted deficit coming from the deepest event
    deps = np.array([e["depth_ppm"] if e["depth_ppm"] is not None else 0.0 for e in ev])
    nin = np.array([e["n_in"] for e in ev])
    contrib = np.clip(deps, 0, None) * nin
    results["single_event_dominance"] = float(contrib.max() / contrib.sum()) if contrib.sum() > 0 else None
    results["n_events_depth_gt_half_mean"] = int(np.sum(deps > 0.5 * np.nanmean(deps)))

    # ---- robustness: different detrend window ----
    rob = {}
    for win in [0.5, 2.0]:
        y2, keep2, sig2 = detrend(t, y_raw, q, window=win)
        w2 = 1 / sig2[keep2] ** 2
        r = empirical_depth(t[keep2], y2[keep2], w2, P, t0, dur, n_rand=1000, avoid=[P / 2])
        rob[f"window_{win}d"] = {"depth_ppm": r["depth_ppm"], "err_mad_ppm": r["err_mad_ppm"], "snr_mad": r["snr_mad"]}
    # and raw (per-quarter normalised, no detrend)
    r = empirical_depth(t, y_raw, 1 / sig_q ** 2, P, t0, dur, n_rand=1000, avoid=[P / 2])
    rob["no_detrend"] = {"depth_ppm": r["depth_ppm"], "err_mad_ppm": r["err_mad_ppm"], "snr_mad": r["snr_mad"]}
    results["depth_robustness"] = rob

    # ---- SAP flux ----
    ts, ys_raw, yse, qs = load("SAP_FLUX")
    ys, keeps, sigs = detrend(ts, ys_raw, qs, window=1.0)
    ws = 1 / sigs[keeps] ** 2
    sap = empirical_depth(ts[keeps], ys[keeps], ws, P, t0, dur, avoid=[P / 2])
    sap_ev = per_event(ts[keeps], ys[keeps], ws, P, t0, dur) if False else None
    results["sap_box"] = sap
    results["sap_per_quarter_scatter_ppm"] = {int(qq): float(1e6 * sigs[qs == qq][0]) for qq in np.unique(qs)}
    if sap["depth_ppm"] and prim["depth_ppm"]:
        results["sap_over_psf_depth_ratio"] = float(sap["depth_ppm"] / prim["depth_ppm"])
    # SAP per-event (for the figure/table)
    QUARTER_LOOKUP = qs[keeps]
    results["sap_events"] = per_event(ts[keeps], ys[keeps], ws, P, t0, dur)
    QUARTER_LOOKUP = qk

    # ---- also: depth on PSF_FLAT_FLUX column (the HLSP's own flattened PSF flux) ----
    try:
        tf, yf_raw, yfe, qf = load("PSF_FLAT_FLUX")
        yf, keepf, sigf = detrend(tf, yf_raw, qf, window=1.0)
        wf = 1 / sigf[keepf] ** 2
        results["psf_flat_box"] = empirical_depth(tf[keepf], yf[keepf], wf, P, t0, dur, n_rand=1000, avoid=[P / 2])
    except Exception as e:  # noqa
        results["psf_flat_box"] = {"error": str(e)}

    # ---- figure ----
    fig, ax = plt.subplots(3, 2, figsize=(14, 12))
    ax = ax.ravel()
    ax[0].plot(periods, pw, lw=0.5, color="0.3")
    ax[0].plot(periods, pw_inv, lw=0.5, color="tab:red", alpha=0.5, label="inverted (2-f)")
    ax[0].axvline(PIPE_P, color="tab:blue", ls="--", lw=0.8, label=f"pipeline P={PIPE_P}")
    ax[0].set_xscale("log"); ax[0].set_xlabel("period [d]"); ax[0].set_ylabel("BLS log-likelihood")
    ax[0].set_title(f"blind BLS 1-50 d (own grid, {blind['n_periods']} periods); top={P_top:.4f} d")
    ax[0].legend(fontsize=8)
    ax[1].plot(periods, snr, lw=0.5, color="0.3")
    ax[1].plot(periods, snr_inv, lw=0.5, color="tab:red", alpha=0.5)
    ax[1].axvline(PIPE_P, color="tab:blue", ls="--", lw=0.8)
    ax[1].set_xscale("log"); ax[1].set_xlabel("period [d]"); ax[1].set_ylabel("BLS depth SNR")
    ax[1].set_title(f"SDE(loglike raw)={blind['sde_loglike_raw']:.1f}, flattened={blind['sde_loglike_flattened']:.1f}; "
                    f"inverted: {blind_inv['sde_loglike_raw']:.1f}/{blind_inv['sde_loglike_flattened']:.1f}", fontsize=9)
    ph = phase_days(tk, P, t0) * 24
    s = np.abs(ph) < 12
    ax[2].scatter(ph[s], (yk[s] - 1) * 1e6, s=4, c=qk[s], cmap="tab10", alpha=0.7)
    # binned
    bins = np.arange(-12, 12.01, 0.5)
    idx = np.digitize(ph[s], bins)
    bx = 0.5 * (bins[1:] + bins[:-1])
    by = np.array([np.average((yk[s][idx == i] - 1) * 1e6, weights=wk[s][idx == i]) if np.any(idx == i) else np.nan
                   for i in range(1, len(bins))])
    ax[2].plot(bx, by, "k-", lw=1.5)
    ax[2].axvspan(-dur * 12, dur * 12, color="tab:blue", alpha=0.15)
    ax[2].set_xlabel("hours from mid-transit"); ax[2].set_ylabel("detrended PSF flux - 1 [ppm]")
    ax[2].set_title(f"fold P={P:.5f} d, t0={t0:.4f} BKJD, dur={dur*24:.2f} h, "
                    f"box depth={prim['depth_ppm']:.0f}+-{prim['err_mad_ppm']:.0f} ppm", fontsize=9)
    ax[2].set_ylim(-4 * max(prim["depth_ppm"] or 1000, 3000), 3 * max(prim["err_mad_ppm"], 1000))
    # phase 0.5
    ph5 = phase_days(tk - P / 2, P, t0) * 24
    s5 = np.abs(ph5) < 12
    ax[3].scatter(ph5[s5], (yk[s5] - 1) * 1e6, s=4, color="0.4", alpha=0.6)
    ax[3].axvspan(-dur * 12, dur * 12, color="tab:orange", alpha=0.15)
    ax[3].set_xlabel("hours from phase 0.5"); ax[3].set_ylabel("ppm")
    ax[3].set_title(f"secondary box: {sec['depth_ppm']:.0f}+-{sec['err_mad_ppm']:.0f} ppm", fontsize=9)
    ax[3].set_ylim(ax[2].get_ylim())
    # per-event depths
    ec = [e["tc_bkjd"] for e in ev]
    ed = [e["depth_ppm"] if e["depth_ppm"] is not None else np.nan for e in ev]
    ee = [e["depth_err_formal_ppm"] if e["depth_err_formal_ppm"] is not None else np.nan for e in ev]
    ax[4].errorbar(ec, ed, yerr=ee, fmt="o", color="k", ms=4)
    ax[4].axhline(prim["depth_ppm"], color="tab:blue", ls="--")
    ax[4].axhline(0, color="0.6", lw=0.5)
    ax[4].set_xlabel("event centre [BKJD]"); ax[4].set_ylabel("local event depth [ppm]")
    ax[4].set_title(f"{len(ev)} epochs with in-window data; odd-even {results.get('odd_even_sigma', np.nan):.1f} sigma", fontsize=9)
    # odd / even stacked
    for par, col, lab in [(1, "tab:red", "odd"), (0, "tab:green", "even")]:
        ss = s & (ep % 2 == par)
        idx2 = np.digitize(ph[ss], bins)
        by2 = np.array([np.average((yk[ss][idx2 == i] - 1) * 1e6, weights=wk[ss][idx2 == i]) if np.any(idx2 == i) else np.nan
                        for i in range(1, len(bins))])
        ax[5].plot(bx, by2, "-", color=col, label=f"{lab}: {results[lab]['depth_ppm']:.0f}+-{results[lab]['err_mad_ppm']:.0f} ppm")
    ax[5].axvspan(-dur * 12, dur * 12, color="tab:blue", alpha=0.15)
    ax[5].legend(fontsize=8); ax[5].set_xlabel("hours from mid-transit"); ax[5].set_ylabel("binned ppm")
    ax[5].set_title("odd vs even (0.5 h bins)", fontsize=9)
    fig.suptitle(f"{FNAME} independent re-derivation (PSF flux, Q{results['quarters']})", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT + "fold_rederive.png", dpi=130)

    with open(OUT + "rederive_results.json", "w") as fh:
        json.dump(results, fh, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
    print(json.dumps({k: v for k, v in results.items() if k not in ("events", "sap_events", "duration_scan_P_fixed")},
                     indent=1, default=str))
    print("events:")
    for e in ev:
        print(e)
