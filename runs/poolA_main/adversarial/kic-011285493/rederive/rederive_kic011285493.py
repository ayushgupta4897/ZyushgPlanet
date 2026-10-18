#!/usr/bin/env python
"""
Independent re-derivation of the kic-011285493 first-screen signal.

Written from scratch for the adversarial step (CLAUDE.md "Before anything goes
public", item 4).  It deliberately imports NOTHING from pipeline/ and uses no
wotan / TLS: only numpy, scipy, astropy, matplotlib.

Steps
  1. load LIGHTCURVE_STITCHED, SAP_QUALITY == 0, per-quarter median normalisation,
     running-median detrend (time-based window, segment-aware), robust clipping
  2. blind BLS 1-50 d with astropy BoxLeastSquares on an own geometric grid
  3. at the best ephemeris: box depth with random-phase empirical error,
     odd/even, secondary at phase 0.5, per-year depths, per-transit depths,
     number of transits, trapezoid duration
  4. same measurement on SAP_FLUX for a dilution check
  5. inverted-flux (2 - f) null with the identical search
Outputs (all under this directory): results.json, fold.png, bls_periodogram.csv
"""
import json
import os
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.optimize import least_squares

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-011285493_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-011285493/rederive"
FNAME = "kic-011285493"
PIPE_P = 10.79100          # pipeline value we are testing against (from the task text)
RNG = np.random.default_rng(20261017)

# ----------------------------------------------------------------------------
# 1. loading and detrending
# ----------------------------------------------------------------------------

def load(col="FLUX"):
    """Return t (BKJD), normalised flux, normalised error, quarter for good cadences."""
    with fits.open(FITS) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        t = np.array(d["TIME"], float)
        q = np.array(d["QUARTER"], int)
        f = np.array(d[col], float)
        ecol = col + "_ERR"
        fe = np.array(d[ecol], float)
        sq = np.array(d["SAP_QUALITY"], int)
    good = (sq == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (f > 0)
    t, q, f, fe = t[good], q[good], f[good], fe[good]
    o = np.argsort(t)
    t, q, f, fe = t[o], q[o], f[o], fe[o]
    fn = np.empty_like(f)
    en = np.empty_like(f)
    for Q in np.unique(q):
        m = q == Q
        med = np.median(f[m])
        fn[m] = f[m] / med
        en[m] = fe[m] / med
    return t, fn, en, q


def segments(t, gap=0.5):
    """Indices splitting the time series at gaps longer than `gap` days."""
    br = np.where(np.diff(t) > gap)[0] + 1
    return np.split(np.arange(len(t)), br)


def running_median(t, f, window=1.0, mask=None, gap=0.5):
    """Time-based running median trend. Points with mask==True are excluded from
    the median (used to protect in-transit points) but still receive a trend."""
    trend = np.full_like(f, np.nan)
    use = np.ones_like(f, bool) if mask is None else ~mask
    for idx in segments(t, gap):
        ts, fs, us = t[idx], f[idx], use[idx]
        tu, fu = ts[us], fs[us]
        if len(tu) < 5:
            trend[idx] = np.median(fs)
            continue
        lo = np.searchsorted(tu, ts - window / 2)
        hi = np.searchsorted(tu, ts + window / 2)
        out = np.empty(len(ts))
        for k in range(len(ts)):
            a, b = lo[k], hi[k]
            if b - a < 3:          # widen at segment edges / sparse regions
                a, b = max(0, a - 3), min(len(tu), b + 3)
            out[k] = np.median(fu[a:b])
        trend[idx] = out
    return trend


def detrend(t, f, window=1.0, mask=None, clip=5.0):
    tr = running_median(t, f, window, mask)
    r = f / tr
    resid = r - 1
    mad = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    keep = np.abs(resid) < clip * mad
    if mask is not None:          # never clip protected (in-transit) points downward
        keep |= (mask & (resid < 0))
    return r, keep, mad


# ----------------------------------------------------------------------------
# 2. blind BLS
# ----------------------------------------------------------------------------

def period_grid(T, pmin=1.0, pmax=50.0, dur_min=0.05, oversample=2.5):
    """Geometric frequency grid: successive frequencies differ by dur_min/(T*oversample)
    so that the phase drift across the baseline is < dur_min/oversample at every P."""
    ratio = 1.0 + dur_min / (T * oversample)
    n = int(np.ceil(np.log(pmax / pmin) / np.log(ratio))) + 1
    f = (1.0 / pmax) * ratio ** np.arange(n)
    return 1.0 / f[f <= 1.0 / pmin][::-1]


def sde_of(periods, power, nwin=2001):
    """Remove a slow trend (running median in grid index) and return the SDE
    spectrum plus the robust std used."""
    from scipy.ndimage import median_filter
    base = median_filter(power, size=nwin, mode="nearest")
    det = power - base
    mad = 1.4826 * np.median(np.abs(det - np.median(det)))
    return det / mad, det, mad


def top_peaks(periods, sde, n=3, sep=0.02):
    """Local maxima sorted by SDE, distinct in period by a fraction `sep`."""
    loc = np.where((sde[1:-1] > sde[:-2]) & (sde[1:-1] >= sde[2:]))[0] + 1
    loc = loc[np.argsort(sde[loc])[::-1]]
    chosen = []
    for i in loc:
        if all(abs(periods[i] / periods[j] - 1) > sep for j in chosen):
            chosen.append(i)
        if len(chosen) == n:
            break
    return chosen


def alias_of(P, Pref, tol=0.003):
    for k, name in [(1, "1x"), (2, "2x"), (3, "3x"), (0.5, "1/2x"), (1 / 3, "1/3x"),
                    (1.5, "3/2x"), (2 / 3, "2/3x")]:
        if abs(P / (k * Pref) - 1) < tol:
            return name
    return None


def run_bls(t, f, e, periods, durations):
    bls = BoxLeastSquares(t, f, dy=e)
    res = bls.power(periods, durations, objective="likelihood")
    return bls, res


# ----------------------------------------------------------------------------
# 3. measurements at an ephemeris
# ----------------------------------------------------------------------------

def phase_days(t, P, t0):
    """Time from nearest transit centre, in days, in [-P/2, P/2)."""
    return ((t - t0 + P / 2) % P) - P / 2


def box_depth(t, f, P, t0, dur, exclude=None, centre=0.0):
    """Mean(out) - mean(in) for a box of width dur centred at phase offset `centre`
    (days).  `exclude` = list of (centre, halfwidth) in days removed from 'out'."""
    dt = phase_days(t, P, t0 - centre)
    inn = np.abs(dt) < dur / 2
    out = np.abs(dt) > dur          # 1 duration either side of the box
    if exclude:
        dt0 = phase_days(t, P, t0)
        for c, hw in exclude:
            out &= np.abs(((dt0 - c + P / 2) % P) - P / 2) > hw
    if inn.sum() < 3:
        return np.nan, inn.sum()
    return np.mean(f[out]) - np.mean(f[inn]), int(inn.sum())


def empirical_depth_error(t, f, P, t0, dur, ntrial=3000, avoid=None):
    """Std of box depths at random phases, avoiding the (centre, halfwidth) regions
    in `avoid` (days).  Returns the robust std and the null sample."""
    vals = []
    tries = 0
    while len(vals) < ntrial and tries < 20 * ntrial:
        tries += 1
        c = RNG.uniform(-P / 2, P / 2)
        if avoid and any(abs(((c - a + P / 2) % P) - P / 2) < hw for a, hw in avoid):
            continue
        d, n = box_depth(t, f, P, t0, dur, centre=c)
        if np.isfinite(d):
            vals.append(d)
    vals = np.array(vals)
    sig = 1.4826 * np.median(np.abs(vals - np.median(vals)))
    return sig, vals


def trapezoid(p, dt):
    """depth, T14, T23/T14, dt0 -> model flux (1 - dip)."""
    depth, T14, frac, dt0 = p
    frac = np.clip(frac, 0.0, 1.0)
    T23 = T14 * frac
    x = np.abs(dt - dt0)
    m = np.ones_like(dt)
    full = x <= T23 / 2
    ing = (x > T23 / 2) & (x < T14 / 2)
    m[full] = 1 - depth
    if T14 > T23:
        m[ing] = 1 - depth * (T14 / 2 - x[ing]) / ((T14 - T23) / 2)
    return m


def fit_trapezoid(dt, f, depth0, dur0):
    sel = np.abs(dt) < 3 * dur0
    x, y = dt[sel], f[sel]

    def resid(p):
        return trapezoid(p, x) - y
    best = None
    for frac0 in (0.3, 0.6, 0.9):
        r = least_squares(resid, [depth0, dur0, frac0, 0.0],
                          bounds=([0, 0.3 * dur0, 0, -dur0 / 2], [0.05, 3 * dur0, 1, dur0 / 2]))
        if best is None or r.cost < best.cost:
            best = r
    return best.x, best


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    os.makedirs(OUT, exist_ok=True)
    log = {}
    t0wall = time.time()

    # ---------- PSF flux ----------
    t, f_raw, e_raw, q = load("FLUX")
    T = t.max() - t.min()
    f1, keep1, mad1 = detrend(t, f_raw, window=1.0)
    t_s, f_s, e_s = t[keep1], f1[keep1], e_raw[keep1]
    log["n_good_quality0"] = int(len(t))
    log["n_after_clip"] = int(len(t_s))
    log["baseline_d"] = float(T)
    log["robust_scatter_ppm_30min"] = float(mad1 * 1e6)
    print(f"loaded {len(t)} q=0 cadences, baseline {T:.1f} d, kept {len(t_s)} after 5-MAD clip, "
          f"robust scatter {mad1*1e6:.0f} ppm", flush=True)

    # ---------- blind BLS 1-50 d ----------
    durations = np.array([0.05, 0.0667, 0.0833, 0.125, 0.1667, 0.2083, 0.25, 0.3333])
    periods = period_grid(T, 1.0, 50.0, dur_min=0.05, oversample=2.5)
    print(f"BLS grid: {len(periods)} periods, {len(durations)} durations", flush=True)
    tt = time.time()
    bls, res = run_bls(t_s, f_s, e_s, periods, durations)
    print(f"BLS done in {time.time()-tt:.0f} s", flush=True)
    sde, det, madp = sde_of(res.period, res.power)
    peaks = top_peaks(res.period, sde, n=3)
    ib = peaks[0]
    P_bls = float(res.period[ib])
    t0_bls = float(res.transit_time[ib])
    dur_bls = float(res.duration[ib])
    dep_bls = float(res.depth[ib])
    sde_best = float(sde[ib])
    sde_raw = float((res.power[ib] - np.median(res.power)) / np.std(res.power))
    print(f"top peak: P={P_bls:.5f} d  t0={t0_bls:.4f}  dur={dur_bls*24:.2f} h  depth={dep_bls*1e6:.0f} ppm  "
          f"SDE={sde_best:.1f} (raw {sde_raw:.1f})", flush=True)
    top3 = []
    for i in peaks:
        top3.append({"period_d": float(res.period[i]), "sde": float(sde[i]),
                     "power": float(res.power[i]), "depth_ppm": float(res.depth[i] * 1e6),
                     "duration_h": float(res.duration[i] * 24), "t0_bkjd": float(res.transit_time[i]),
                     "alias_of_pipeline_P": alias_of(res.period[i], PIPE_P)})
        print("  peak", top3[-1], flush=True)
    recovered_blind = alias_of(P_bls, PIPE_P) == "1x"
    recovered_alias = alias_of(P_bls, PIPE_P)
    np.savetxt(os.path.join(OUT, "bls_periodogram.csv"),
               np.c_[res.period, res.power, sde], delimiter=",", header="period_d,power,sde", comments="")

    # ---------- refine around the top peak ----------
    fine_p = np.linspace(P_bls * (1 - 3e-3), P_bls * (1 + 3e-3), 3001)
    fine_d = np.linspace(0.04, 0.4, 37)
    bls_f, resf = run_bls(t_s, f_s, e_s, fine_p, fine_d)
    j = int(np.argmax(resf.power))
    P = float(resf.period[j]); t0 = float(resf.transit_time[j]); dur = float(resf.duration[j])
    # move t0 to the middle of the data span to minimise P/t0 covariance in reporting
    n_mid = np.round((np.median(t_s) - t0) / P)
    t0 = t0 + n_mid * P
    print(f"refined: P={P:.6f} d  t0={t0:.5f} BKJD  dur_bls={dur*24:.2f} h", flush=True)

    # ---------- re-detrend with the transit masked ----------
    dt_all = phase_days(t, P, t0)
    tmask = np.abs(dt_all) < 1.0 * dur          # protect +-1 duration
    f2, keep2, mad2 = detrend(t, f_raw, window=1.0, mask=tmask)
    tm, fm, em, qm = t[keep2], f2[keep2], e_raw[keep2], q[keep2]
    dtm = phase_days(tm, P, t0)

    # ---------- trapezoid duration ----------
    dep0, _ = box_depth(tm, fm, P, t0, dur)
    (tdepth, T14, frac, dt0), fitres = fit_trapezoid(dtm, fm, max(dep0, 1e-4), dur)
    T23 = T14 * frac
    print(f"trapezoid: depth={tdepth*1e6:.0f} ppm T14={T14*24:.2f} h T23={T23*24:.2f} h dt0={dt0*24*60:.1f} min",
          flush=True)
    dur_use = T14 if 0.3 * dur < T14 < 3 * dur else dur

    # ---------- box depth + empirical error ----------
    depth, n_in = box_depth(tm, fm, P, t0, dur_use, exclude=[(P / 2, dur_use)])
    err, null = empirical_depth_error(tm, fm, P, t0, dur_use, ntrial=3000,
                                      avoid=[(0.0, 1.5 * dur_use)])
    snr = depth / err
    print(f"box depth {depth*1e6:.0f} +- {err*1e6:.0f} ppm (empirical, {len(null)} random phases)  SNR={snr:.1f}  "
          f"n_in={n_in}", flush=True)
    # formal error for reference
    inn = np.abs(dtm) < dur_use / 2
    formal = np.std(fm[~inn]) / np.sqrt(inn.sum())
    print(f"formal error {formal*1e6:.0f} ppm (white-noise)", flush=True)

    # ---------- secondary at phase 0.5 ----------
    sec, n_sec = box_depth(tm, fm, P, t0, dur_use, exclude=[(0.0, dur_use)], centre=P / 2)
    sec_sig = sec / err
    # also scan the deepest secondary at any phase away from primary (worst case)
    sc_phases = np.linspace(-P / 2, P / 2, 400)
    sc = np.array([box_depth(tm, fm, P, t0, dur_use, centre=c)[0] if abs(c) > 1.5 * dur_use else np.nan
                   for c in sc_phases])
    k = int(np.nanargmax(sc))
    print(f"secondary @0.5: {sec*1e6:.0f} ppm ({sec_sig:.1f} sigma, n={n_sec}); deepest off-transit box "
          f"{sc[k]*1e6:.0f} ppm at phase {sc_phases[k]/P:+.3f}", flush=True)

    # ---------- odd / even ----------
    ep = np.round((tm - t0) / P).astype(int)
    odd = ep % 2 == 1
    d_odd, n_odd = box_depth(tm[odd], fm[odd], P, t0, dur_use)
    d_even, n_even = box_depth(tm[~odd], fm[~odd], P, t0, dur_use)
    e_odd = err * np.sqrt(n_in / max(n_odd, 1))
    e_even = err * np.sqrt(n_in / max(n_even, 1))
    oe_sig = abs(d_odd - d_even) / np.hypot(e_odd, e_even)
    print(f"odd {d_odd*1e6:.0f}+-{e_odd*1e6:.0f} (n={n_odd})  even {d_even*1e6:.0f}+-{e_even*1e6:.0f} (n={n_even})  "
          f"-> {oe_sig:.2f} sigma", flush=True)

    # ---------- per-year ----------
    year_edges = [(131, 540, "Y1 Q1-Q4"), (540, 910, "Y2 Q5-Q8"), (910, 1275, "Y3 Q9-Q12"), (1275, 1600, "Y4 Q13-Q17")]
    per_year = []
    for a, b, name in year_edges:
        m = (tm >= a) & (tm < b)
        d_y, n_y = box_depth(tm[m], fm[m], P, t0, dur_use)
        e_y = err * np.sqrt(n_in / max(n_y, 1))
        per_year.append({"year": name, "depth_ppm": float(d_y * 1e6) if np.isfinite(d_y) else None,
                         "err_ppm": float(e_y * 1e6), "n_in": int(n_y)})
        print(f"  {name}: {d_y*1e6:.0f} +- {e_y*1e6:.0f} ppm (n_in={n_y})", flush=True)
    yd = np.array([p["depth_ppm"] for p in per_year if p["depth_ppm"] is not None])
    ye = np.array([p["err_ppm"] for p in per_year if p["depth_ppm"] is not None])
    chi2_year = float(np.sum(((yd - depth * 1e6) / ye) ** 2))

    # ---------- per transit ----------
    expected_in = dur_use / (29.4 / 60 / 24)
    per_transit = []
    n_transits_any, n_transits_half = 0, 0
    for n in np.unique(ep):
        m = ep == n
        inn_n = m & (np.abs(dtm) < dur_use / 2)
        loc = m & (np.abs(dtm) > dur_use) & (np.abs(dtm) < 4 * dur_use)
        if inn_n.sum() == 0:
            continue
        n_transits_any += 1
        if inn_n.sum() >= 0.5 * expected_in:
            n_transits_half += 1
        d_n = (np.mean(fm[loc]) if loc.sum() > 5 else 1.0) - np.mean(fm[inn_n])
        e_n = mad2 / np.sqrt(inn_n.sum())
        per_transit.append({"epoch": int(n), "tc_bkjd": float(t0 + n * P), "depth_ppm": float(d_n * 1e6),
                            "err_ppm": float(e_n * 1e6), "n_in": int(inn_n.sum()), "n_local_oot": int(loc.sum())})
    pt_d = np.array([p["depth_ppm"] for p in per_transit])
    pt_e = np.array([p["err_ppm"] for p in per_transit])
    full = np.array([p["n_in"] >= 0.5 * expected_in for p in per_transit])
    # single-event dominance: fraction of total (depth/err^2) weight carried by the strongest event
    w = pt_d / pt_e ** 2
    dom = float(np.max(w[full]) / np.sum(w[full])) if full.sum() and np.sum(w[full]) > 0 else None
    # leave-one-out: depth when the strongest full event is removed
    strongest = int(np.argmax((pt_d / pt_e) * full))
    m_loo = ep != per_transit[strongest]["epoch"]
    d_loo, n_loo = box_depth(tm[m_loo], fm[m_loo], P, t0, dur_use)
    print(f"transits with data: {n_transits_any}, with >=50% coverage: {n_transits_half} "
          f"(expected {int(np.floor(T/P))+1} in baseline); strongest event epoch {per_transit[strongest]['epoch']} "
          f"depth {pt_d[strongest]:.0f}+-{pt_e[strongest]:.0f}; depth without it {d_loo*1e6:.0f} ppm; "
          f"dominance {dom}", flush=True)

    # ---------- SAP flux ----------
    ts, s_raw, se_raw, qs = load("SAP_FLUX")
    tmask_s = np.abs(phase_days(ts, P, t0)) < 1.0 * dur
    s2, keeps, mads = detrend(ts, s_raw, window=1.0, mask=tmask_s)
    tsm, ssm = ts[keeps], s2[keeps]
    sdepth, n_in_s = box_depth(tsm, ssm, P, t0, dur_use, exclude=[(P / 2, dur_use)])
    serr, _ = empirical_depth_error(tsm, ssm, P, t0, dur_use, ntrial=2000, avoid=[(0.0, 1.5 * dur_use)])
    with fits.open(FITS) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        g = d["SAP_QUALITY"] == 0
        med_psf = float(np.nanmedian(d["FLUX"][g])); med_sap = float(np.nanmedian(d["SAP_FLUX"][g]))
    ratio = sdepth / depth
    print(f"SAP depth {sdepth*1e6:.0f} +- {serr*1e6:.0f} ppm; PSF depth {depth*1e6:.0f}; ratio SAP/PSF = {ratio:.2f}; "
          f"median PSF/SAP flux = {med_psf/med_sap:.3f}", flush=True)

    # ---------- inverted-flux null (same search) ----------
    tt = time.time()
    f_inv = 2.0 - f_s
    bls_i, res_i = run_bls(t_s, f_inv, e_s, periods, durations)
    sde_i, _, _ = sde_of(res_i.period, res_i.power)
    pk_i = top_peaks(res_i.period, sde_i, n=3)
    inv_top = [{"period_d": float(res_i.period[i]), "sde": float(sde_i[i]),
                "depth_ppm": float(res_i.depth[i] * 1e6)} for i in pk_i]
    print(f"inverted null: top SDE {inv_top[0]['sde']:.1f} at P={inv_top[0]['period_d']:.4f} "
          f"({time.time()-tt:.0f} s)", flush=True)

    # ---------- is the pipeline period in the real periodogram? ----------
    kp = int(np.argmin(np.abs(res.period - PIPE_P)))
    w = slice(max(0, kp - 400), kp + 400)
    k_near = w.start + int(np.argmax(sde[w]))
    sde_at_pipeP = float(np.max(sde[w]))

    # ---------- figure ----------
    make_figure(res.period, sde, peaks, P, t0, dur_use, tm, fm, odd, tdepth, T14, T23, dt0, depth, sec, OUT)

    # ---------- results ----------
    results = {
        "fname": FNAME, "date": "2026-10-17", "method": "independent re-derivation; no pipeline/ code",
        "data": {"fits": FITS, "column": "FLUX (PSF)", "quality_cut": "SAP_QUALITY==0",
                 "normalisation": "per-quarter median", "detrend": "running median, 1.0 d time window, "
                 "segment-aware (gap>0.5 d), 5-MAD clip; second pass with +-1 duration masked around transits",
                 **log},
        "blind_bls": {"period_range_d": [1, 50], "n_periods": int(len(periods)),
                      "durations_h": (durations * 24).round(2).tolist(), "objective": "likelihood",
                      "sde_definition": "(power - running median over 2001 grid points) / 1.4826*MAD",
                      "top3": top3, "top_is_pipeline_period": bool(recovered_blind),
                      "top_alias_of_pipeline_period": recovered_alias,
                      "sde_within_0.4pct_of_pipeline_P": sde_at_pipeP,
                      "sde_raw_mean_std_definition": sde_raw},
        "ephemeris": {"period_d": P, "t0_bkjd": t0, "t0_bjd": t0 + 2454833.0,
                      "duration_bls_h": dur * 24, "duration_trapezoid_T14_h": T14 * 24,
                      "duration_trapezoid_T23_h": T23 * 24, "T23_over_T14": float(frac),
                      "trapezoid_depth_ppm": float(tdepth * 1e6), "duration_used_for_boxes_h": dur_use * 24},
        "depth": {"box_depth_ppm": float(depth * 1e6), "empirical_err_ppm": float(err * 1e6),
                  "formal_err_ppm": float(formal * 1e6), "empirical_snr": float(snr), "n_in_transit": int(n_in),
                  "n_random_phase_trials": int(len(null))},
        "odd_even": {"odd_ppm": float(d_odd * 1e6), "odd_err_ppm": float(e_odd * 1e6), "n_odd": int(n_odd),
                     "even_ppm": float(d_even * 1e6), "even_err_ppm": float(e_even * 1e6), "n_even": int(n_even),
                     "sigma": float(oe_sig)},
        "secondary": {"phase": 0.5, "depth_ppm": float(sec * 1e6), "sigma": float(sec_sig), "n_in": int(n_sec),
                      "deepest_offtransit_box_ppm": float(sc[k] * 1e6), "deepest_offtransit_phase": float(sc_phases[k] / P),
                      "deepest_offtransit_sigma": float(sc[k] / err)},
        "per_year": per_year, "per_year_chi2_vs_mean": chi2_year, "per_year_dof": int(len(yd) - 1),
        "transits": {"n_with_any_data": int(n_transits_any), "n_with_ge50pct_coverage": int(n_transits_half),
                     "n_expected_in_baseline": int(np.floor(T / P)) + 1,
                     "single_event_dominance": dom, "strongest_epoch": per_transit[strongest],
                     "depth_without_strongest_ppm": float(d_loo * 1e6), "per_transit": per_transit},
        "sap_check": {"sap_depth_ppm": float(sdepth * 1e6), "sap_err_ppm": float(serr * 1e6),
                      "sap_over_psf_depth_ratio": float(ratio), "median_psf_flux_e_s": med_psf,
                      "median_sap_flux_e_s": med_sap, "psf_over_sap_flux": med_psf / med_sap,
                      "n_in_sap": int(n_in_s)},
        "inverted_null": {"top3": inv_top, "real_top_sde": sde_best},
        "pipeline_reference": {"period_d": PIPE_P, "depth_ppm": 796, "duration_h": 4.5},
        "wall_s": time.time() - t0wall,
    }
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(results, fh, indent=1, default=float)
    print(f"wrote {OUT}/results.json  ({time.time()-t0wall:.0f} s)", flush=True)


def make_figure(periods, sde, peaks, P, t0, dur, tm, fm, odd, tdepth, T14, T23, dt0, depth, sec, out):
    # palette: text-neutral ink, one data hue (blue), one model hue (orange); recessive hairline grid
    INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e2e1dc", "#fcfcfb"
    RAW, BIN, MODEL, ALT = "#c3c2b7", "#2a78d6", "#eb6834", "#1baf7a"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "text.color": INK, "axes.facecolor": SURF, "figure.facecolor": SURF,
                         "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.2, 1])

    ax = fig.add_subplot(gs[0, :])
    ax.plot(periods, sde, color=BIN, lw=0.6)
    for i in peaks:
        ax.annotate(f"{periods[i]:.3f} d", (periods[i], sde[i]), textcoords="offset points", xytext=(4, 2), color=INK2)
    ax.set_xscale("log"); ax.set_xlabel("period [d]"); ax.set_ylabel("BLS SDE")
    ax.set_title(f"{FNAME}: own blind BLS 1-50 d on PSF flux (running-median detrend)", loc="left", color=INK)

    dt = phase_days(tm, P, t0)
    ax = fig.add_subplot(gs[1, :])
    sel = np.abs(dt) < 4 * dur
    ax.plot(dt[sel] * 24, (fm[sel] - 1) * 1e6, ".", ms=2, color=RAW, rasterized=True, label="30-min cadences")
    bw = 20 / 60 / 24
    edges = np.arange(-4 * dur, 4 * dur + bw, bw)
    idx = np.digitize(dt[sel], edges)
    bx, by, be = [], [], []
    for b in range(1, len(edges)):
        m = idx == b
        if m.sum() > 3:
            bx.append(0.5 * (edges[b - 1] + edges[b])); by.append(np.mean(fm[sel][m]) - 1)
            be.append(np.std(fm[sel][m]) / np.sqrt(m.sum()))
    bx, by, be = np.array(bx) * 24, np.array(by) * 1e6, np.array(be) * 1e6
    ax.errorbar(bx, by, be, fmt="o", ms=4, color=BIN, ecolor=BIN, elinewidth=1, label="20-min bins")
    xx = np.linspace(-4 * dur, 4 * dur, 800)
    ax.plot(xx * 24, (trapezoid([tdepth, T14, T23 / T14 if T14 > 0 else 0.5, dt0], xx) - 1) * 1e6, color=MODEL, lw=2,
            label=f"trapezoid: {tdepth*1e6:.0f} ppm, T14={T14*24:.1f} h")
    ax.axhline(-depth * 1e6, color=MODEL, lw=1, ls=":", label=f"box depth {depth*1e6:.0f} ppm")
    ax.set_ylim(-4000, 3000); ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("rel. flux [ppm]")
    ax.set_title(f"fold at P={P:.5f} d, t0={t0:.4f} BKJD", loc="left", color=INK)
    ax.legend(loc="lower right", frameon=False, ncol=2)

    def binned(dtx, fx, lab, col, axx):
        selx = np.abs(dtx) < 4 * dur
        idx = np.digitize(dtx[selx], edges)
        xs, ys, es = [], [], []
        for b in range(1, len(edges)):
            m = idx == b
            if m.sum() > 3:
                xs.append(0.5 * (edges[b - 1] + edges[b]) * 24); ys.append((np.mean(fx[selx][m]) - 1) * 1e6)
                es.append(np.std(fx[selx][m]) / np.sqrt(m.sum()) * 1e6)
        axx.errorbar(xs, ys, es, fmt="o", ms=4, color=col, elinewidth=1, label=lab)

    ax = fig.add_subplot(gs[2, 0])
    binned(dt[odd], fm[odd], "odd", BIN, ax); binned(dt[~odd], fm[~odd], "even", MODEL, ax)
    ax.set_ylim(-2500, 1500); ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("rel. flux [ppm]")
    ax.set_title("odd vs even (20-min bins)", loc="left", color=INK); ax.legend(frameon=False)

    ax = fig.add_subplot(gs[2, 1])
    dts = phase_days(tm, P, t0 - P / 2)
    binned(dts, fm, "phase 0.5", ALT, ax)
    ax.axhline(-sec * 1e6, color=INK2, lw=1, ls=":", label=f"box at 0.5: {sec*1e6:.0f} ppm")
    ax.set_ylim(-2500, 1500); ax.set_xlabel("hours from phase 0.5"); ax.set_ylabel("rel. flux [ppm]")
    ax.set_title("secondary window", loc="left", color=INK); ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(os.path.join(out, "fold.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
