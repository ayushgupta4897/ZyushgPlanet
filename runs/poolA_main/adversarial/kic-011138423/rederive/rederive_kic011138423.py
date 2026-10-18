#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Independent re-derivation of the kic-011138423 first-screen signal
(ZyushgPlanet, run poolA_main, adversarial stage, 2026-10-17).

Written from scratch with numpy / scipy / astropy only.  Nothing under pipeline/ is
imported, read or copied; the only pipeline numbers used are the ones quoted in the
task statement (P = 4.89450 d) so that the blind result can be compared against them.

Steps
  1. KBONUS-BKG stitched light curve -> SAP_QUALITY == 0 -> per-quarter median
     normalisation -> gap-aware running median (1 d window) -> one-sided sigma clip
     (positive outliers > 4 sigma, negatives only beyond 10 sigma).
  2. Blind BLS 1-50 d on PSF flux with an own log-spaced period grid (phase drift
     < 1/3 of a 2 h box over the baseline), 7 durations 1-6 h; SDE raw and with the
     red-noise trend of the spectrum removed; top-3 distinct peaks; alias check
     against 4.89450 d.  Identical search on the inverted flux (2 - f) as a null.
  3. Narrow refinement around the best peak; transit-masked re-detrend; box depth with
     empirical error from identical boxes at random phases; odd/even; secondary at
     phase 0.5 plus a full phase scan; per-year depths; per-event depths with gap
     proximity; number of transits; trapezoid fit for duration/shape.
  4. Same depth on SAP_FLUX (same ephemeris, same detrend) -> dilution check.
Outputs (all in OUT): results.json, fold_rederive.png, periodogram.npz, events.csv
"""
import json
import os
import time

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter
from scipy.optimize import minimize
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-011138423_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-011138423/rederive"
FNAME = "kic-011138423"
P_PIPE = 4.89450          # pipeline period quoted in the task statement (comparison only)
BKJD0 = 2454833.0
CADENCE = 0.0204337       # d, Kepler long cadence
WINDOW = 1.0              # d, running-median window
GAP = 0.5                 # d, break the light curve into segments at gaps longer than this
PMIN, PMAX = 1.0, 50.0
DUR_GRID = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]) / 24.0
OS = 3                    # phase oversampling of the log-period grid (reference box = 2 h)
N_NULL = 4000             # random-phase boxes per empirical error
SEED = 20261017
rng = np.random.default_rng(SEED)
YEARS = {"Y1": [0, 1, 2, 3, 4], "Y2": [5, 6, 7, 8], "Y3": [9, 10, 11, 12],
         "Y4": [13, 14, 15, 16, 17]}

# figure palette (dataviz reference palette, light surface)
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
MUTED, GRID, INK, INK2, SURF = "#898781", "#e1e0d9", "#0b0b0b", "#52514e", "#fcfcfb"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ----------------------------------------------------------------------------- loading
def load():
    d = fits.open(FITS)["LIGHTCURVE_STITCHED"].data
    cols = dict(t="TIME", f="FLUX", fe="FLUX_ERR", s="SAP_FLUX", se="SAP_FLUX_ERR",
                q="SAP_QUALITY", Q="QUARTER")
    a = {}
    for k, v in cols.items():
        arr = np.asarray(d[v])
        a[k] = arr.astype(int) if k in ("q", "Q") else arr.astype(float)
    m = a["q"] == 0
    for k in ("t", "f", "fe", "s", "se"):
        m &= np.isfinite(a[k])
    m &= (a["f"] > 0) & (a["s"] > 0)
    o = np.argsort(a["t"][m])
    return {k: v[m][o] for k, v in a.items()}, int(len(a["t"])), int(m.sum())


def normalise(f, fe, Q):
    fn = np.empty_like(f)
    fen = np.empty_like(fe)
    for qq in np.unique(Q):
        k = Q == qq
        med = np.median(f[k])
        fn[k] = f[k] / med
        fen[k] = fe[k] / med
    return fn, fen


# --------------------------------------------------------------------------- detrending
def segments(t, gap=GAP):
    br = np.where(np.diff(t) > gap)[0] + 1
    return np.split(np.arange(len(t)), br)


def running_median(t, f, use=None, window=WINDOW):
    """Gap-aware running median over a +-window/2 time window.  Points with use==False
    are not used in the medians (transit-masked re-detrend) but still get a trend value."""
    if use is None:
        use = np.ones(len(t), bool)
    trend = np.full(len(t), np.nan)
    hw = window / 2.0
    for seg in segments(t):
        ts, fs, us = t[seg], f[seg], use[seg]
        tu, fu = ts[us], fs[us]
        if len(tu) < 5:
            trend[seg] = np.median(fs)
            continue
        lo = np.searchsorted(tu, ts - hw)
        hi = np.searchsorted(tu, ts + hw, side="right")
        for i in range(len(seg)):
            if hi[i] - lo[i] >= 3:
                trend[seg[i]] = np.median(fu[lo[i]:hi[i]])
            else:
                j = int(np.clip(np.searchsorted(tu, ts[i]), 1, len(tu) - 1))
                trend[seg[i]] = np.median(fu[max(0, j - 12):j + 12])
    return trend


def robust_std(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def clip_mask(r, hi=4.0, lo=10.0):
    med = np.median(r)
    s = robust_std(r)
    return (r - med < hi * s) & (r - med > -lo * s)


# --------------------------------------------------------------------------------- BLS
def period_grid(pmin, pmax, baseline, dref=2.0 / 24.0, os_=OS):
    dlnp = dref / (os_ * baseline)
    n = int(np.ceil(np.log(pmax / pmin) / dlnp))
    return pmin * np.exp(dlnp * np.arange(n + 1))


def bls_search(t, y, dy, periods, durs):
    bls = BoxLeastSquares(t, y, dy=dy)
    res = bls.power(periods, durs, objective="likelihood", method="fast")
    pw = np.asarray(res.power, float)
    size = min(2001, (len(pw) // 10) * 2 + 1)
    trend = median_filter(pw, size=size, mode="nearest")
    pdt = pw - trend
    stats = dict(sde_raw=float((pw.max() - pw.mean()) / pw.std()),
                 sde_detrended=float((pdt.max() - pdt.mean()) / pdt.std()),
                 sde_robust=float(pdt.max() / robust_std(pdt)),
                 i_raw=int(np.argmax(pw)), i_det=int(np.argmax(pdt)))
    return res, pw, pdt, stats


def top_peaks(periods, pdt, n=3, excl=0.01):
    m = np.ones(len(pdt), bool)
    out = []
    for _ in range(n):
        i = int(np.argmax(np.where(m, pdt, -np.inf)))
        out.append(i)
        m &= np.abs(periods / periods[i] - 1) > excl
    return out


def alias(p, pref, tol=0.005):
    for a, b in [(1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (3, 2), (2, 3), (4, 1), (1, 4),
                 (5, 1), (1, 5), (5, 2), (2, 5)]:
        if abs(p / (pref * a / b) - 1) < tol:
            return f"{a}/{b}"
    return None


def refine(t, y, dy, p0, baseline):
    dln = (2.0 / 24.0) / (OS * baseline)
    periods = p0 * np.exp(np.linspace(-6 * dln, 6 * dln, 2401))
    durs = np.arange(0.75, 6.01, 0.125) / 24.0
    bls = BoxLeastSquares(t, y, dy=dy)
    res = bls.power(periods, durs, objective="likelihood", method="fast", oversample=20)
    i = int(np.argmax(res.power))
    return dict(period=float(res.period[i]), t0=float(res.transit_time[i]),
                duration=float(res.duration[i]), depth=float(res.depth[i]),
                depth_err=float(res.depth_err[i]), depth_snr=float(res.depth_snr[i]),
                power=float(res.power[i]), periods=np.asarray(res.period),
                power_arr=np.asarray(res.power))


# ------------------------------------------------------------------- depth statistics
def phase_days(t, P, t0, offset=0.0):
    ph = ((t - t0) / P - offset + 0.5) % 1.0 - 0.5
    return ph * P


def box_stats(t, y, P, t0, dur, sel=None, offset=0.0, n_null=N_NULL, excl=1.5):
    """Box depth at phase `offset` (0 primary, 0.5 secondary) on subset `sel`.
    Empirical error = scatter of identical boxes placed at random phases that avoid the
    primary transit (+-excl*dur), i.e. it includes red noise and gaps."""
    if sel is None:
        sel = np.ones(len(t), bool)
    tt, yy = t[sel], y[sel]
    dp0 = phase_days(tt, P, t0, 0.0)
    oot = np.abs(dp0) > excl * dur
    base = float(np.median(yy[oot]))
    sig_pt = robust_std(yy[oot])
    inb = np.abs(phase_days(tt, P, t0, offset)) < dur / 2
    n_in = int(inb.sum())
    depth = base - float(np.mean(yy[inb])) if n_in > 0 else np.nan
    err_formal = sig_pt / np.sqrt(max(n_in, 1))
    u = rng.uniform(0, 1, size=n_null * 3)
    dist = np.minimum(u, 1 - u) * P
    u = u[dist > excl * dur + dur / 2][:n_null]
    null = np.full(len(u), np.nan)
    for k, uk in enumerate(u):
        ib = np.abs(phase_days(tt, P, t0, uk)) < dur / 2
        if ib.sum() > 0:
            null[k] = base - np.mean(yy[ib])
    null = null[np.isfinite(null)]
    err_emp = float(np.std(null))
    return dict(depth_ppm=1e6 * depth, err_emp_ppm=1e6 * err_emp,
                err_emp_robust_ppm=1e6 * float(robust_std(null)),
                err_formal_ppm=1e6 * err_formal,
                snr_emp=float(depth / err_emp) if err_emp > 0 else np.nan,
                snr_formal=float(depth / err_formal), n_in=n_in, n_points=int(len(tt)),
                n_null=int(len(null)), null_mean_ppm=1e6 * float(np.mean(null)),
                null_p_ge_depth=float(np.mean(null >= depth)) if n_in > 0 else np.nan,
                baseline=base, sigma_point_ppm=1e6 * sig_pt)


def phase_scan(t, y, P, t0, dur, err_emp, step=0.005, excl=1.5):
    dp0 = phase_days(t, P, t0, 0.0)
    oot = np.abs(dp0) > excl * dur
    base = np.median(y[oot])
    rows = []
    for o in np.arange(0.0, 1.0, step):
        if min(o, 1 - o) * P < excl * dur + dur / 2:
            continue
        ib = np.abs(phase_days(t, P, t0, o)) < dur / 2
        if ib.sum() < 3:
            continue
        rows.append((o, 1e6 * (base - np.mean(y[ib])), int(ib.sum())))
    arr = np.array(rows)
    k = int(np.argmax(np.abs(arr[:, 1])))
    return dict(max_abs_phase=float(arr[k, 0]), max_abs_depth_ppm=float(arr[k, 1]),
                max_abs_sigma=float(arr[k, 1] / (1e6 * err_emp)), n_boxes=int(len(arr)),
                scan=arr)


def gap_edges(t, gap=0.1):
    i = np.where(np.diff(t) > gap)[0]
    return np.concatenate([[t[0]], t[i], t[i + 1], [t[-1]]])


def per_event(t, y, P, t0, dur):
    n = np.round((t - t0) / P).astype(int)
    edges = gap_edges(t)
    exp_n = dur / CADENCE
    rows = []
    for e in np.unique(n):
        tc = t0 + e * P
        d = t - tc
        inb = np.abs(d) < dur / 2
        loc = (np.abs(d) < 1.0) & (np.abs(d) > 1.5 * dur)
        n_in = int(inb.sum())
        if n_in == 0 or loc.sum() < 10:
            continue
        base = np.median(y[loc])
        sig = robust_std(y[loc])
        rows.append(dict(epoch=int(e), tc_bkjd=float(tc), n_in=n_in,
                         coverage=float(n_in / exp_n),
                         depth_ppm=1e6 * float(base - np.mean(y[inb])),
                         err_ppm=1e6 * float(sig / np.sqrt(n_in)), n_local=int(loc.sum()),
                         gap_dist_d=float(np.min(np.abs(edges - tc)))))
    return rows


def trapezoid(x, depth, T14, T12, dt):
    xx = np.abs(x - dt)
    half = T14 / 2.0
    T12 = min(max(T12, 1e-4), half)
    s = np.zeros_like(xx)
    s[xx <= half - T12] = 1.0
    ramp = (xx > half - T12) & (xx < half)
    s[ramp] = (half - xx[ramp]) / T12
    return 1.0 - depth * s


def fit_trapezoid(t, y, P, t0, dur, depth0):
    dp = phase_days(t, P, t0)
    w = np.abs(dp) < max(4 * dur, 0.3)
    x, yy = dp[w], y[w]
    sig = robust_std(yy)

    def chi2(p):
        depth, T14, T12, dt = p
        if depth < 0 or T14 <= 0.01 or T14 > 0.5 or T12 < 0 or T12 > T14 / 2 or abs(dt) > dur:
            return 1e30
        return float(np.sum(((yy - trapezoid(x, depth, T14, T12, dt)) / sig) ** 2))

    best = None
    for f14 in (0.7, 1.0, 1.3):
        for f12 in (0.05, 0.25, 0.5):
            p0 = [max(depth0, 1e-4), f14 * dur, f12 * f14 * dur / 2, 0.0]
            r = minimize(chi2, p0, method="Nelder-Mead",
                         options=dict(xatol=1e-6, fatol=1e-4, maxiter=6000))
            if best is None or r.fun < best.fun:
                best = r
    depth, T14, T12, dt = best.x
    chi2_box = chi2([depth, T14, 1e-4, dt])
    chi2_flat = float(np.sum(((yy - 1.0) / sig) ** 2))
    return dict(depth_ppm=1e6 * depth, T14_h=24 * T14, T12_h=24 * T12,
                T23_h=24 * (T14 - 2 * T12), shape_T23_over_T14=float((T14 - 2 * T12) / T14),
                dt_min=1440 * dt, chi2=float(best.fun), chi2_box_same_T14=chi2_box,
                chi2_flat=chi2_flat, n_points=int(w.sum()),
                delta_chi2_vs_flat=chi2_flat - float(best.fun))


def binned(x, y, width):
    edges = np.arange(x.min(), x.max() + width, width)
    idx = np.digitize(x, edges) - 1
    xb, yb, eb = [], [], []
    for i in range(len(edges) - 1):
        k = idx == i
        if k.sum() >= 3:
            xb.append(x[k].mean())
            yb.append(y[k].mean())
            eb.append(np.std(y[k]) / np.sqrt(k.sum()))
    return np.array(xb), np.array(yb), np.array(eb)


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


# -------------------------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()
    d, n_all, n_good = load()
    t, Q = d["t"], d["Q"]
    baseline = float(t.max() - t.min())
    log(f"loaded {n_all} cadences, {n_good} with SAP_QUALITY==0 and finite flux; "
        f"baseline {baseline:.1f} d, quarters {sorted(set(Q.tolist()))}")

    fn, fen = normalise(d["f"], d["fe"], Q)
    sn, sen = normalise(d["s"], d["se"], Q)

    # --- pass 1 detrend + one-sided clip, pass 2 detrend on kept points
    tr1 = running_median(t, fn)
    keep = clip_mask(fn / tr1, 4.0, 10.0)
    tr2 = running_median(t, fn, use=keep)
    y_blind = fn / tr2
    t, Q, fn, fen, sn, sen, y_blind = [a[keep] for a in (t, Q, fn, fen, sn, sen, y_blind)]
    n_used = int(len(t))
    dy = np.empty_like(y_blind)
    for qq in np.unique(Q):
        k = Q == qq
        dy[k] = robust_std(y_blind[k])
    log(f"clipped {int((~keep).sum())} cadences; {n_used} used; "
        f"median point scatter {1e6*np.median(dy):.0f} ppm")

    # --- blind BLS 1-50 d
    periods = period_grid(PMIN, PMAX, baseline)
    log(f"blind BLS: {len(periods)} periods x {len(DUR_GRID)} durations")
    res, pw, pdt, st = bls_search(t, y_blind, dy, periods, DUR_GRID)
    log(f"blind BLS done in {time.time()-t_start:.0f}s; SDE raw {st['sde_raw']:.2f}, "
        f"detrended {st['sde_detrended']:.2f}, robust {st['sde_robust']:.2f}")
    pk = top_peaks(periods, pdt, n=3)
    peaks = []
    for rank, i in enumerate(pk, 1):
        peaks.append(dict(rank=rank, period_d=float(periods[i]), power=float(pw[i]),
                          power_detrended=float(pdt[i]),
                          sde_detrended=float(pdt[i] / pdt.std()),
                          sde_robust=float(pdt[i] / robust_std(pdt)),
                          duration_h=24 * float(res.duration[i]),
                          depth_ppm=1e6 * float(res.depth[i]),
                          depth_snr_bls=float(res.depth_snr[i]),
                          t0_bkjd=float(res.transit_time[i]),
                          alias_of_pipeline_P=alias(float(periods[i]), P_PIPE)))
        log(f"  peak {rank}: P={periods[i]:.5f} d dur={24*res.duration[i]:.2f} h "
            f"depth={1e6*res.depth[i]:.0f} ppm pdt={pdt[i]:.1f} alias={peaks[-1]['alias_of_pipeline_P']}")
    i_raw = st["i_raw"]
    raw_top = dict(period_d=float(periods[i_raw]), power=float(pw[i_raw]),
                   alias_of_pipeline_P=alias(float(periods[i_raw]), P_PIPE))
    top_P = peaks[0]["period_d"]
    recovered_exact = abs(top_P / P_PIPE - 1) < 1e-3
    recovered_alias = peaks[0]["alias_of_pipeline_P"] is not None
    # power at the pipeline period itself (nearest grid point)
    ip = int(np.argmin(np.abs(periods - P_PIPE)))
    at_pipe = dict(period_d=float(periods[ip]), power_detrended=float(pdt[ip]),
                   sde_robust=float(pdt[ip] / robust_std(pdt)),
                   rank_among_detrended_peaks=int((pdt > pdt[ip]).sum()) + 1)

    # --- inverted-flux null with the identical search
    log("inverted-flux null BLS ...")
    res_i, pw_i, pdt_i, st_i = bls_search(t, 2.0 - y_blind, dy, periods, DUR_GRID)
    pk_i = top_peaks(periods, pdt_i, n=3)
    inv = dict(sde_raw=st_i["sde_raw"], sde_detrended=st_i["sde_detrended"],
               sde_robust=st_i["sde_robust"],
               top_peaks=[dict(period_d=float(periods[i]), power_detrended=float(pdt_i[i]),
                               sde_robust=float(pdt_i[i] / robust_std(pdt_i)),
                               depth_ppm=1e6 * float(res_i.depth[i]),
                               duration_h=24 * float(res_i.duration[i])) for i in pk_i])
    log(f"inverted null: SDE raw {st_i['sde_raw']:.2f} det {st_i['sde_detrended']:.2f} "
        f"rob {st_i['sde_robust']:.2f} top P {periods[pk_i[0]]:.4f}")
    np.savez_compressed(os.path.join(OUT, "periodogram.npz"), period=periods, power=pw,
                        power_detrended=pdt, power_inverted=pw_i,
                        power_inverted_detrended=pdt_i, duration=np.asarray(res.duration),
                        depth=np.asarray(res.depth))

    # --- refinement around the top peak
    rf = refine(t, y_blind, dy, top_P, baseline)
    P, t0_raw, dur = rf["period"], rf["t0"], rf["duration"]
    n_mid = np.round((0.5 * (t.min() + t.max()) - t0_raw) / P)
    t0 = float(t0_raw + n_mid * P)
    log(f"refined: P={P:.6f} d, dur={24*dur:.2f} h, t0={t0:.4f} BKJD, box depth "
        f"{1e6*rf['depth']:.0f} ppm (BLS err {1e6*rf['depth_err']:.0f})")

    # --- transit-masked re-detrend (PSF and SAP) at the refined ephemeris
    intr = np.abs(phase_days(t, P, t0)) < 1.0 * dur
    y_psf = fn / running_median(t, fn, use=~intr)
    y_sap = sn / running_median(t, sn, use=~intr)
    # trapezoid fit for duration / shape and a t0 correction
    trap = fit_trapezoid(t, y_psf, P, t0, dur, rf["depth"])
    t0 += trap["dt_min"] / 1440.0
    log(f"trapezoid: T14={trap['T14_h']:.2f} h, T23/T14={trap['shape_T23_over_T14']:.2f}, "
        f"depth={trap['depth_ppm']:.0f} ppm, dt={trap['dt_min']:.1f} min")

    # --- depth statistics on PSF flux
    prim = box_stats(t, y_psf, P, t0, dur)
    prim_blind = box_stats(t, y_blind, P, t0, dur, n_null=1000)
    sec = box_stats(t, y_psf, P, t0, dur, offset=0.5)
    scan = phase_scan(t, y_psf, P, t0, dur, prim["err_emp_ppm"] / 1e6)
    epoch = np.round((t - t0) / P).astype(int)
    odd = box_stats(t, y_psf, P, t0, dur, sel=(epoch % 2 == 1))
    even = box_stats(t, y_psf, P, t0, dur, sel=(epoch % 2 == 0))
    oe_sigma = abs(odd["depth_ppm"] - even["depth_ppm"]) / np.hypot(odd["err_emp_ppm"],
                                                                   even["err_emp_ppm"])
    log(f"PSF depth {prim['depth_ppm']:.0f} +- {prim['err_emp_ppm']:.0f} ppm (emp), "
        f"SNR {prim['snr_emp']:.1f}; odd {odd['depth_ppm']:.0f}+-{odd['err_emp_ppm']:.0f}, "
        f"even {even['depth_ppm']:.0f}+-{even['err_emp_ppm']:.0f} -> {oe_sigma:.2f} sigma; "
        f"secondary {sec['depth_ppm']:.0f}+-{sec['err_emp_ppm']:.0f}")

    # --- per-year
    ev = per_event(t, y_psf, P, t0, dur)
    years = {}
    for yname, qs in YEARS.items():
        sel = np.isin(Q, qs)
        if sel.sum() < 200:
            continue
        s = box_stats(t, y_psf, P, t0, dur, sel=sel, n_null=2000)
        qsel = np.isin(Q, qs)
        n_tr = sum(1 for e in ev if np.any(qsel & (np.abs(t - e["tc_bkjd"]) < dur / 2))
                   and e["coverage"] >= 0.5)
        s["n_transits_cov50"] = int(n_tr)
        s["quarters"] = [int(q) for q in np.unique(Q[sel])]
        years[yname] = s
        log(f"  {yname}: depth {s['depth_ppm']:.0f} +- {s['err_emp_ppm']:.0f} ppm, "
            f"{n_tr} transits")

    # --- events / single-event dominance
    n_any = len(ev)
    n_cov50 = sum(1 for e in ev if e["coverage"] >= 0.5)
    snr2 = np.array([(e["depth_ppm"] / e["err_ppm"]) ** 2 for e in ev])
    pos = np.array([e["depth_ppm"] > 0 for e in ev])
    dom = float(snr2[pos].max() / snr2[pos].sum()) if pos.any() else np.nan
    strongest = ev[int(np.argmax(np.where(pos, snr2, -1)))] if pos.any() else None
    n_pred = int(np.floor((t.max() - t0) / P) - np.ceil((t.min() - t0) / P) + 1)
    with open(os.path.join(OUT, "events.csv"), "w") as fh:
        fh.write("epoch,tc_bkjd,tc_bjd,n_in,coverage,depth_ppm,err_ppm,n_local,gap_dist_d\n")
        for e in ev:
            fh.write(f"{e['epoch']},{e['tc_bkjd']:.5f},{e['tc_bkjd']+BKJD0:.5f},{e['n_in']},"
                     f"{e['coverage']:.2f},{e['depth_ppm']:.1f},{e['err_ppm']:.1f},"
                     f"{e['n_local']},{e['gap_dist_d']:.3f}\n")
    log(f"transits: {n_pred} predicted in baseline, {n_any} with data, {n_cov50} with >=50% "
        f"coverage; strongest event carries {100*dom:.0f}% of sum SNR^2")

    # --- SAP flux depth at the same ephemeris
    sap = box_stats(t, y_sap, P, t0, dur)
    sap_sec = box_stats(t, y_sap, P, t0, dur, offset=0.5, n_null=1000)
    ratio = sap["depth_ppm"] / prim["depth_ppm"]
    ratio_err = abs(ratio) * np.hypot(sap["err_emp_ppm"] / sap["depth_ppm"],
                                      prim["err_emp_ppm"] / prim["depth_ppm"])
    flux_ratio = float(np.median(d["s"]) / np.median(d["f"]))
    log(f"SAP depth {sap['depth_ppm']:.0f} +- {sap['err_emp_ppm']:.0f} ppm; "
        f"SAP/PSF depth ratio {ratio:.2f} +- {ratio_err:.2f}; median SAP/PSF flux {flux_ratio:.2f}")

    # --- odd/even and per-year on SAP too (cheap, useful cross-check)
    sap_odd = box_stats(t, y_sap, P, t0, dur, sel=(epoch % 2 == 1), n_null=1000)
    sap_even = box_stats(t, y_sap, P, t0, dur, sel=(epoch % 2 == 0), n_null=1000)

    results = dict(
        fname=FNAME, date="2026-10-17", script=os.path.abspath(__file__),
        fits=FITS, pipeline_period_tested=P_PIPE,
        data=dict(n_cadences_total=n_all, n_quality0=n_good, n_used_after_clip=n_used,
                  baseline_d=baseline, quarters=sorted(set(Q.tolist())),
                  median_point_scatter_ppm=1e6 * float(np.median(dy)),
                  median_psf_flux_e_s=float(np.median(d["f"])),
                  median_sap_flux_e_s=float(np.median(d["s"])), sap_over_psf_flux=flux_ratio),
        method=dict(quality="SAP_QUALITY==0", normalisation="per-quarter median",
                    detrend=f"gap-aware running median, {WINDOW} d window, segments split at "
                            f"gaps > {GAP} d; clip > +4 sigma and < -10 sigma (robust); "
                            "transit-masked re-detrend (+-1 duration) for depth statistics",
                    bls="astropy BoxLeastSquares, objective=likelihood, method=fast, "
                        f"log-period grid {PMIN}-{PMAX} d with dlnP = 2h/({OS} x baseline) "
                        f"({len(periods)} periods), durations {sorted((24*DUR_GRID).tolist())} h",
                    sde="raw: (max-mean)/std of power; detrended: same after subtracting a "
                        "2001-point running median of the spectrum; robust: max/(1.4826 MAD)",
                    empirical_error=f"std of the depth of {N_NULL} identical boxes at random "
                                    "phases avoiding +-1.5 durations of the primary",
                    seed=SEED),
        blind_search=dict(n_periods=int(len(periods)), sde_raw=st["sde_raw"],
                          sde_detrended=st["sde_detrended"], sde_robust=st["sde_robust"],
                          top3_peaks=peaks, raw_power_top=raw_top,
                          at_pipeline_period=at_pipe,
                          recovered_exact_top_peak=bool(recovered_exact),
                          recovered_alias_top_peak=bool(recovered_alias)),
        inverted_null=inv,
        ephemeris=dict(period_d=P, t0_bkjd=t0, t0_bjd=t0 + BKJD0, duration_h=24 * dur,
                       refinement_window_d=[float(rf["periods"].min()), float(rf["periods"].max())],
                       bls_box_depth_ppm=1e6 * rf["depth"], bls_depth_err_ppm=1e6 * rf["depth_err"],
                       bls_depth_snr=rf["depth_snr"],
                       period_err_est_d=float(dur / (3 * baseline) * P)),
        trapezoid=trap,
        primary_psf=prim, primary_psf_blind_detrend=prim_blind,
        secondary_psf=sec,
        secondary_scan=dict(max_abs_phase=scan["max_abs_phase"],
                            max_abs_depth_ppm=scan["max_abs_depth_ppm"],
                            max_abs_sigma=scan["max_abs_sigma"], n_boxes=scan["n_boxes"]),
        odd_psf=odd, even_psf=even, odd_even_sigma=float(oe_sigma),
        per_year=years,
        transits=dict(n_predicted_in_baseline=n_pred, n_with_data=n_any,
                      n_coverage_ge_50pct=n_cov50, strongest_event=strongest,
                      strongest_event_snr2_fraction=dom, n_events_positive=int(pos.sum())),
        sap=dict(primary=sap, secondary=sap_sec, odd=sap_odd, even=sap_even,
                 depth_ratio_sap_over_psf=float(ratio), depth_ratio_err=float(ratio_err)),
        runtime_s=float(time.time() - t_start),
    )
    results = jsonable(results)
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    log("results.json written")

    # ------------------------------------------------------------------------ figure
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                         "axes.titlecolor": INK, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.facecolor": SURF, "figure.facecolor": SURF})
    fig, axs = plt.subplots(3, 2, figsize=(12.5, 12.5))
    fig.suptitle(f"{FNAME}: independent re-derivation (own loader/detrend/BLS, no pipeline code)"
                 f"   P = {P:.5f} d, t0 = {t0+BKJD0:.4f} BJD, box duration {24*dur:.2f} h",
                 fontsize=10.5, color=INK)

    ax = axs[0, 0]
    ax.plot(periods, pdt_i, color=MUTED, lw=0.5, alpha=0.7, label="inverted flux (2 - f), same search")
    ax.plot(periods, pdt, color=BLUE, lw=0.6, label="PSF flux, trend-removed BLS power")
    for p_ in peaks:
        ax.plot(p_["period_d"], p_["power_detrended"], "v", color=ORANGE, ms=7, mec=SURF, mew=0.8)
        ax.annotate(f"{p_['rank']}: {p_['period_d']:.4f} d", (p_["period_d"], p_["power_detrended"]),
                    xytext=(4, 4), textcoords="offset points", fontsize=8, color=INK2)
    ax.axvline(P_PIPE, color=INK2, ls="--", lw=0.8, label=f"pipeline P = {P_PIPE} d")
    ax.set_xscale("log")
    ax.set_xlabel("period (d)")
    ax.set_ylabel("BLS power, red-noise trend removed")
    ax.set_title(f"blind BLS 1-50 d: SDE raw {st['sde_raw']:.1f} / robust {st['sde_robust']:.1f}; "
                 f"inverted null robust {st_i['sde_robust']:.1f}", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(True, which="both", axis="x")

    def fold_panel(ax, y, title, color, label, show_models=True, win_h=8.0):
        dp = phase_days(t, P, t0) * 24.0
        w = np.abs(dp) < win_h
        ax.plot(dp[w], 1e6 * (y[w] - 1), ".", color=MUTED, ms=2, alpha=0.25, rasterized=True)
        xb, yb, eb = binned(dp[w], y[w], 1.0 / 3.0)
        ax.errorbar(xb, 1e6 * (yb - 1), 1e6 * eb, fmt="o", color=color, ms=4, lw=1,
                    mec=SURF, mew=0.6, label=label)
        if show_models:
            xx = np.linspace(-win_h, win_h, 800)
            ax.plot(xx, 1e6 * (trapezoid(xx / 24, trap["depth_ppm"] / 1e6, trap["T14_h"] / 24,
                                         trap["T12_h"] / 24, 0.0) - 1), color=ORANGE, lw=1.4,
                    label=f"trapezoid T14 {trap['T14_h']:.2f} h")
            box = np.where(np.abs(xx) < 12 * dur, -prim["depth_ppm"], 0.0)
            ax.plot(xx, box, color=INK2, lw=1.0, ls="--", label=f"box {24*dur:.2f} h")
        ax.set_xlim(-win_h, win_h)
        ax.set_ylim(-6000, 4000)
        ax.set_xlabel("hours from mid-transit")
        ax.set_ylabel("relative flux (ppm)")
        ax.set_title(title, fontsize=9)
        ax.grid(True, axis="y")
        ax.legend(frameon=False, fontsize=8, loc="lower right")

    fold_panel(axs[0, 1], y_psf,
               f"PSF flux fold (20-min bins): depth {prim['depth_ppm']:.0f} +- "
               f"{prim['err_emp_ppm']:.0f} ppm, empirical SNR {prim['snr_emp']:.1f}",
               BLUE, "PSF flux, binned")

    ax = axs[1, 0]
    dp = phase_days(t, P, t0) * 24.0
    w = np.abs(dp) < 8
    for sel, col, lab in ((epoch % 2 == 1, AQUA, f"odd  {odd['depth_ppm']:.0f} +- {odd['err_emp_ppm']:.0f} ppm"),
                          (epoch % 2 == 0, VIOLET, f"even {even['depth_ppm']:.0f} +- {even['err_emp_ppm']:.0f} ppm")):
        xb, yb, eb = binned(dp[w & sel], y_psf[w & sel], 0.5)
        ax.errorbar(xb, 1e6 * (yb - 1), 1e6 * eb, fmt="o", color=col, ms=4, lw=1, mec=SURF,
                    mew=0.6, label=lab)
    ax.axhline(0, color=GRID, lw=0.8)
    ax.set_xlim(-8, 8)
    ax.set_ylim(-4000, 2500)
    ax.set_xlabel("hours from mid-transit")
    ax.set_ylabel("relative flux (ppm)")
    ax.set_title(f"odd vs even transits (30-min bins): difference {oe_sigma:.2f} sigma", fontsize=9)
    ax.grid(True, axis="y")
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax = axs[1, 1]
    dps = phase_days(t, P, t0, 0.5) * 24.0
    w = np.abs(dps) < 8
    ax.plot(dps[w], 1e6 * (y_psf[w] - 1), ".", color=MUTED, ms=2, alpha=0.25, rasterized=True)
    xb, yb, eb = binned(dps[w], y_psf[w], 1.0 / 3.0)
    ax.errorbar(xb, 1e6 * (yb - 1), 1e6 * eb, fmt="o", color=BLUE, ms=4, lw=1, mec=SURF, mew=0.6,
                label="PSF flux at phase 0.5")
    ax.axhline(0, color=GRID, lw=0.8)
    ax.set_xlim(-8, 8)
    ax.set_ylim(-6000, 4000)
    ax.set_xlabel("hours from phase 0.5")
    ax.set_ylabel("relative flux (ppm)")
    ax.set_title(f"secondary at phase 0.5: {sec['depth_ppm']:.0f} +- {sec['err_emp_ppm']:.0f} ppm "
                 f"({sec['snr_emp']:.1f} sigma); max |box| elsewhere {scan['max_abs_sigma']:.1f} sigma "
                 f"at phase {scan['max_abs_phase']:.3f}", fontsize=9)
    ax.grid(True, axis="y")
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax = axs[2, 0]
    w = np.abs(dp) < 8
    for y_, col, lab in ((y_psf, BLUE, f"PSF  {prim['depth_ppm']:.0f} +- {prim['err_emp_ppm']:.0f} ppm"),
                         (y_sap, ORANGE, f"SAP  {sap['depth_ppm']:.0f} +- {sap['err_emp_ppm']:.0f} ppm")):
        xb, yb, eb = binned(dp[w], y_[w], 0.5)
        ax.errorbar(xb, 1e6 * (yb - 1), 1e6 * eb, fmt="o", color=col, ms=4, lw=1, mec=SURF,
                    mew=0.6, label=lab)
    ax.axhline(0, color=GRID, lw=0.8)
    ax.set_xlim(-8, 8)
    ax.set_ylim(-3000, 1500)
    ax.set_xlabel("hours from mid-transit")
    ax.set_ylabel("relative flux (ppm)")
    ax.set_title(f"PSF vs SAP flux (30-min bins): depth ratio SAP/PSF = {ratio:.2f} +- {ratio_err:.2f}; "
                 f"median SAP/PSF flux = {flux_ratio:.2f}", fontsize=9)
    ax.grid(True, axis="y")
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax = axs[2, 1]
    names = list(years.keys())
    xs = np.arange(len(names))
    ax.errorbar(xs, [years[n]["depth_ppm"] for n in names], [years[n]["err_emp_ppm"] for n in names],
                fmt="o", color=BLUE, ms=6, lw=1.2, mec=SURF, mew=0.6, label="per-year box depth")
    ax.axhline(prim["depth_ppm"], color=INK2, ls="--", lw=0.9, label=f"all data {prim['depth_ppm']:.0f} ppm")
    ax.axhline(0, color=GRID, lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{n}\nQ{years[n]['quarters'][0]}-Q{years[n]['quarters'][-1]}\n"
                        f"{years[n]['n_transits_cov50']} transits" for n in names])
    ax.set_ylabel("box depth (ppm)")
    ax.set_title(f"per-year depths; {n_cov50} transits with >=50% coverage of {n_pred} predicted; "
                 f"strongest event = {100*dom:.0f}% of sum SNR^2", fontsize=9)
    ax.grid(True, axis="y")
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    for a in axs.ravel():
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(OUT, "fold_rederive.png"), dpi=150)
    log(f"figure written; total runtime {time.time()-t_start:.0f}s")


if __name__ == "__main__":
    main()
