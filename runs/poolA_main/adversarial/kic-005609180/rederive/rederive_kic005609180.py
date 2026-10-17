#!/usr/bin/env python
"""
Independent re-derivation for kic-005609180 (KBONUS-BKG PSF light curve).

Written from scratch with numpy/scipy/astropy only. Deliberately does NOT import
or read anything under pipeline/ and does NOT use wotan or TLS. The purpose is an
adversarial, shared-code-free check of the period, depth, duration, odd/even,
secondary and dilution numbers reported by the main pipeline for this star.

Steps
  1. load LIGHTCURVE_STITCHED, keep SAP_QUALITY==0 and finite cadences,
     normalise per quarter, detrend with a time-windowed running median (1.0 d)
  2. blind BLS 1-50 d on my own period grid (astropy BoxLeastSquares)
  3. refine ephemeris, re-detrend with the transit masked, measure box depth with
     an empirical (random-phase) error, odd/even, secondary at phase 0.5,
     per-year depths, N transits, duration
  4. same depth measurement on SAP_FLUX for a dilution check
  5. write results JSON + one fold figure
"""
import json
import os
import sys
import time as _time

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FNAME = "kic-005609180"
FITS_PATH = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
             "hlsp_kbonus-bkg_kepler_kepler_kic-005609180_kepler_v1.0_lc.fits")
OUT_DIR = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005609180/rederive"
PIPELINE_P = 21.10630      # quoted by the task; only used for *comparison*, never as a prior
PIPELINE_DEPTH_PPM = 656.0
PIPELINE_DUR_H = 3.4

RNG = np.random.default_rng(20261017)


# --------------------------------------------------------------------------- #
# 1. loading / normalisation / detrending
# --------------------------------------------------------------------------- #
def load(col="FLUX"):
    with fits.open(FITS_PATH) as h:
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
    # per-quarter normalisation
    fn = np.empty_like(f)
    fen = np.empty_like(fe)
    for qq in np.unique(q):
        s = q == qq
        med = np.nanmedian(f[s])
        fn[s] = f[s] / med
        fen[s] = fe[s] / med
    return t, fn, fen, q


def segments(t, gap=0.5):
    """indices of contiguous segments, broken at gaps > `gap` days"""
    br = np.where(np.diff(t) > gap)[0] + 1
    edges = np.concatenate([[0], br, [len(t)]])
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def running_median(t, f, window=1.0, mask=None):
    """time-windowed running median trend; points with mask==True are excluded
    from the median (used to protect in-transit points on the second pass).
    Trend for masked points is the median of unmasked neighbours in the window."""
    trend = np.full_like(f, np.nan)
    if mask is None:
        mask = np.zeros_like(f, bool)
    half = window / 2.0
    for (a, b) in segments(t):
        ts, fs, ms = t[a:b], f[a:b], mask[a:b]
        tg, fg = ts[~ms], fs[~ms]
        if len(tg) < 5:
            trend[a:b] = np.nanmedian(fs) if len(fs) else 1.0
            continue
        lo = np.searchsorted(tg, ts - half, side="left")
        hi = np.searchsorted(tg, ts + half, side="right")
        for i in range(len(ts)):
            if hi[i] - lo[i] >= 3:
                trend[a + i] = np.median(fg[lo[i]:hi[i]])
            else:
                # tiny window at an edge: widen until >= 3 points
                l2, h2 = lo[i], hi[i]
                while h2 - l2 < 3 and (l2 > 0 or h2 < len(fg)):
                    l2 = max(l2 - 1, 0)
                    h2 = min(h2 + 1, len(fg))
                trend[a + i] = np.median(fg[l2:h2])
    return trend


def detrend(t, f, fe, window=1.0, mask=None, clip_hi=4.0):
    tr = running_median(t, f, window, mask)
    fd = f / tr
    fed = fe / tr
    # clip positive outliers only (never clip downward dips)
    r = fd - 1.0
    mad = 1.4826 * np.median(np.abs(r - np.median(r)))
    keep = r < clip_hi * mad
    return fd, fed, keep, mad


# --------------------------------------------------------------------------- #
# 2. BLS on my own grid
# --------------------------------------------------------------------------- #
def period_grid(t, pmin=1.0, pmax=50.0, oversample=3.0, dur_scale=0.5):
    """Frequency grid whose step at each period is set by the shortest transit
    duration plausibly seen there (0.5 x central-transit duration of a Sun-like
    star), divided by the baseline and an oversampling factor."""
    T = t.max() - t.min()
    fmin, fmax = 1.0 / pmax, 1.0 / pmin
    freqs = [fmax]
    f = fmax
    while f > fmin:
        P = 1.0 / f
        dmin = dur_scale * 0.076 * P ** (1.0 / 3.0)   # days
        df = dmin / (T * oversample)
        f = f - df
        freqs.append(f)
    freqs = np.array(freqs)
    freqs = freqs[freqs >= fmin]
    return 1.0 / freqs[::-1]   # ascending periods


def sde_of(power):
    p = power[np.isfinite(power)]
    return (p.max() - np.median(p)) / np.std(p)


def top_peaks(periods, power, n=3, excl_frac=0.02):
    """greedy top-N peaks, excluding +-excl_frac in period around each pick"""
    pw = power.copy()
    out = []
    for _ in range(n):
        i = int(np.nanargmax(pw))
        P = periods[i]
        out.append((P, float(power[i]), i))
        pw[np.abs(periods - P) / P < excl_frac] = -np.inf
    return out


def run_bls(t, f, fe, pmin=1.0, pmax=50.0, durations_h=(1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 9.0, 12.0)):
    periods = period_grid(t, pmin, pmax)
    durs = np.array(durations_h) / 24.0
    bls = BoxLeastSquares(t, f, dy=fe)
    t0 = _time.time()
    res = bls.power(periods, durs, objective="likelihood", method="fast")
    dt = _time.time() - t0
    return bls, periods, res, dt


# --------------------------------------------------------------------------- #
# 3. measurements at an ephemeris
# --------------------------------------------------------------------------- #
def phase_fold(t, P, t0):
    return ((t - t0) / P + 0.5) % 1.0 - 0.5   # phase in [-0.5, 0.5), transit at 0


def box_depth(t, f, fe, P, t0, dur, phase_center=0.0, oot_window=3.0):
    """depth = (local out-of-transit mean) - (in-box mean), in flux units.
    oot = points within +-oot_window*dur of centre but outside +-0.75 dur."""
    ph = phase_fold(t, P, t0) - phase_center
    ph = (ph + 0.5) % 1.0 - 0.5
    dtp = ph * P
    inb = np.abs(dtp) < dur / 2.0
    oot = (np.abs(dtp) > 0.75 * dur) & (np.abs(dtp) < oot_window * dur)
    if inb.sum() < 3 or oot.sum() < 10:
        return np.nan, np.nan, int(inb.sum())
    w_in = 1.0 / fe[inb] ** 2
    w_out = 1.0 / fe[oot] ** 2
    m_in = np.sum(w_in * f[inb]) / np.sum(w_in)
    m_out = np.sum(w_out * f[oot]) / np.sum(w_out)
    # formal error from scatter
    e_in = np.std(f[inb], ddof=1) / np.sqrt(inb.sum())
    e_out = np.std(f[oot], ddof=1) / np.sqrt(oot.sum())
    return m_out - m_in, np.hypot(e_in, e_out), int(inb.sum())


def empirical_depth_error(t, f, fe, P, t0, dur, n=400, avoid=(0.0, 0.5), avoid_w=None):
    """std of box depths at random phases, avoiding the transit and secondary."""
    if avoid_w is None:
        avoid_w = 2.0 * dur / P
    deps = []
    tries = 0
    while len(deps) < n and tries < 20 * n:
        tries += 1
        pc = RNG.uniform(-0.5, 0.5)
        if any(np.abs(((pc - a) + 0.5) % 1.0 - 0.5) < avoid_w for a in avoid):
            continue
        d, _, nin = box_depth(t, f, fe, P, t0, dur, phase_center=pc)
        if np.isfinite(d):
            deps.append(d)
    deps = np.array(deps)
    return np.std(deps), np.mean(deps), len(deps)


def transit_epochs(t, P, t0, dur, min_frac=0.5):
    """list of (epoch_n, t_mid, n_in_box, n_expected) for epochs with data."""
    n_all = np.arange(np.floor((t.min() - t0) / P) - 1, np.ceil((t.max() - t0) / P) + 2).astype(int)
    cad = 29.4244 / 60.0 / 24.0
    n_exp = dur / cad
    out = []
    for n in n_all:
        tm = t0 + n * P
        if tm < t.min() - dur or tm > t.max() + dur:
            continue
        nin = np.sum(np.abs(t - tm) < dur / 2.0)
        nloc = np.sum((np.abs(t - tm) > 0.75 * dur) & (np.abs(t - tm) < 3 * dur))
        out.append((int(n), float(tm), int(nin), float(n_exp), int(nloc)))
    return out


def per_event_depth(t, f, fe, tm, dur):
    inb = np.abs(t - tm) < dur / 2.0
    oot = (np.abs(t - tm) > 0.75 * dur) & (np.abs(t - tm) < 3 * dur)
    if inb.sum() < 2 or oot.sum() < 8:
        return np.nan, np.nan, int(inb.sum())
    d = np.mean(f[oot]) - np.mean(f[inb])
    e = np.hypot(np.std(f[oot], ddof=1) / np.sqrt(oot.sum()),
                 np.std(f[oot], ddof=1) / np.sqrt(inb.sum()))
    return d, e, int(inb.sum())


def refine_ephemeris(t, f, fe, P0, t00, dur0, dP_frac=2e-4, nP=201, nT=201, durs_h=None):
    """fine grid over P, T0 and duration maximising box SNR (depth / formal err)"""
    if durs_h is None:
        durs_h = np.arange(1.0, 8.01, 0.25)
    best = (-np.inf, P0, t00, dur0)
    Ps = P0 * (1.0 + np.linspace(-dP_frac, dP_frac, nP))
    T0s = t00 + np.linspace(-0.5, 0.5, nT) * dur0
    # first P/T0 at fixed duration
    for P in Ps:
        for T0 in T0s:
            d, e, nin = box_depth(t, f, fe, P, T0, dur0)
            if np.isfinite(d) and e > 0 and d / e > best[0]:
                best = (d / e, P, T0, dur0)
    _, P, T0, _ = best
    # then duration at fixed P/T0, re-centering T0 per duration
    bestd = (-np.inf, dur0, T0)
    for dh in durs_h:
        dur = dh / 24.0
        for T0b in T0 + np.linspace(-0.25, 0.25, 41) * dur:
            d, e, nin = box_depth(t, f, fe, P, T0b, dur)
            if np.isfinite(d) and e > 0 and d / e > bestd[0]:
                bestd = (d / e, dur, T0b)
    return P, bestd[2], bestd[1], bestd[0]


def duration_profile(t, f, fe, P, t0, durs_h):
    """box depth and snr versus assumed duration; also trapezoid-ish shape check"""
    rows = []
    for dh in durs_h:
        d, e, nin = box_depth(t, f, fe, P, t0, dh / 24.0)
        rows.append((dh, d, e, (d / e) if (np.isfinite(d) and e > 0) else np.nan, nin))
    return rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = {"fname": FNAME, "date": "2026-10-17", "fits": FITS_PATH,
               "method": "independent numpy/astropy re-derivation (no pipeline code, no wotan, no TLS)"}

    # ------------------------------------------------------------------ PSF
    t, f, fe, q = load("FLUX")
    results["n_cadences_quality0"] = int(len(t))
    results["baseline_d"] = float(t.max() - t.min())
    results["quarters"] = [int(x) for x in np.unique(q)]

    fd, fed, keep, mad = detrend(t, f, fe, window=1.0)
    results["detrend"] = {"method": "time-windowed running median, 1.0 d window, per contiguous segment (gap>0.5 d)",
                          "positive_clip_sigma": 4.0, "n_clipped": int((~keep).sum()),
                          "point_to_point_mad_ppm": float(mad * 1e6)}
    t1, f1, fe1, q1 = t[keep], fd[keep], fed[keep], q[keep]
    results["n_cadences_used"] = int(len(t1))

    # ------------------------------------------------------------------ BLS
    print("running blind BLS 1-50 d ...", flush=True)
    bls, periods, res, dt = run_bls(t1, f1, fe1)
    print(f"  {len(periods)} periods in {dt:.1f} s", flush=True)
    power = np.asarray(res.power)
    sde = sde_of(power)
    peaks = top_peaks(periods, power, n=3)
    ibest = peaks[0][2]
    P_bls = float(res.period[ibest]); t0_bls = float(res.transit_time[ibest])
    dur_bls = float(res.duration[ibest]); dep_bls = float(res.depth[ibest])
    results["bls"] = {
        "period_range_d": [1.0, 50.0], "n_periods": int(len(periods)),
        "durations_h": [1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 9.0, 12.0],
        "objective": "likelihood", "runtime_s": dt,
        "sde_top": float(sde),
        "top3": [{"period_d": float(res.period[i]), "power": float(pw),
                  "sde": float((pw - np.median(power)) / np.std(power)),
                  "depth_ppm": float(res.depth[i] * 1e6), "duration_h": float(res.duration[i] * 24),
                  "t0_bkjd": float(res.transit_time[i])} for (P, pw, i) in peaks],
        "best": {"period_d": P_bls, "t0_bkjd": t0_bls, "duration_h": dur_bls * 24, "depth_ppm": dep_bls * 1e6},
    }
    # alias check against the pipeline period
    ratio = P_bls / PIPELINE_P
    alias = None
    for num, den in [(1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (3, 2), (2, 3)]:
        if abs(ratio - num / den) < 2e-3:
            alias = f"{num}/{den}"
    results["bls"]["top_peak_vs_pipeline_period"] = {"pipeline_P": PIPELINE_P, "ratio": float(ratio),
                                                     "alias": alias, "is_same_period": alias == "1/1"}
    results["recovered_blind"] = bool(alias == "1/1")
    # also: power of the BLS periodogram at the pipeline period, for the record
    j = int(np.argmin(np.abs(periods - PIPELINE_P)))
    results["bls"]["power_near_pipeline_P"] = {"period_d": float(periods[j]), "power": float(power[j]),
                                               "sde": float((power[j] - np.median(power)) / np.std(power))}
    for k, pk in enumerate(results["bls"]["top3"]):
        print(f"  peak {k+1}: P={pk['period_d']:.5f} d  SDE={pk['sde']:.2f}  depth={pk['depth_ppm']:.0f} ppm  dur={pk['duration_h']:.2f} h")

    # ------------------------------------------------------------------ refine + masked re-detrend
    P1, T01, dur1, snr1 = refine_ephemeris(t1, f1, fe1, P_bls, t0_bls, dur_bls)
    # mask in-transit (1.5x duration) and re-detrend the raw normalised flux
    ph = phase_fold(t, P1, T01)
    inmask = np.abs(ph * P1) < 0.75 * dur1
    fd2, fed2, keep2, mad2 = detrend(t, f, fe, window=1.0, mask=inmask)
    keep2 = keep2 | inmask     # never clip in-transit points
    t2, f2, fe2, q2 = t[keep2], fd2[keep2], fed2[keep2], q[keep2]
    P2, T02, dur2, snr2 = refine_ephemeris(t2, f2, fe2, P1, T01, dur1)
    results["ephemeris"] = {"period_d": float(P2), "t0_bkjd": float(T02), "t0_bjd": float(T02 + 2454833.0),
                            "duration_h": float(dur2 * 24), "note": "box-SNR maximised on masked re-detrended PSF flux"}
    print(f"refined: P={P2:.5f}  T0={T02:.4f} BKJD  dur={dur2*24:.2f} h", flush=True)

    # duration profile
    prof = duration_profile(t2, f2, fe2, P2, T02, np.arange(1.0, 10.01, 0.5))
    results["duration_profile"] = [{"dur_h": a, "depth_ppm": b * 1e6, "err_ppm": c * 1e6, "snr": d, "n_in": e}
                                   for (a, b, c, d, e) in prof]

    # box depth + empirical error (PSF)
    dep, dep_formal, nin = box_depth(t2, f2, fe2, P2, T02, dur2)
    emp_sig, emp_mean, nemp = empirical_depth_error(t2, f2, fe2, P2, T02, dur2)
    results["depth_psf"] = {"depth_ppm": float(dep * 1e6), "formal_err_ppm": float(dep_formal * 1e6),
                            "empirical_err_ppm": float(emp_sig * 1e6), "random_phase_mean_ppm": float(emp_mean * 1e6),
                            "n_random_phases": int(nemp), "n_in_box": int(nin),
                            "empirical_snr": float(dep / emp_sig), "formal_snr": float(dep / dep_formal)}
    print(f"PSF depth = {dep*1e6:.0f} +- {emp_sig*1e6:.0f} ppm (empirical), SNR_emp={dep/emp_sig:.2f}", flush=True)

    # odd / even
    ph2 = phase_fold(t2, P2, T02)
    epoch = np.round((t2 - T02) / P2).astype(int)
    odd = (epoch % 2) == 1
    d_odd, e_odd, n_odd = box_depth(t2[odd], f2[odd], fe2[odd], P2, T02, dur2)
    d_even, e_even, n_even = box_depth(t2[~odd], f2[~odd], fe2[~odd], P2, T02, dur2)
    # empirical errors for each half
    s_odd, _, _ = empirical_depth_error(t2[odd], f2[odd], fe2[odd], P2, T02, dur2, n=200)
    s_even, _, _ = empirical_depth_error(t2[~odd], f2[~odd], fe2[~odd], P2, T02, dur2, n=200)
    oe_sig = abs(d_odd - d_even) / np.hypot(s_odd, s_even)
    results["odd_even"] = {"odd_depth_ppm": float(d_odd * 1e6), "odd_err_ppm": float(s_odd * 1e6), "n_odd_in": int(n_odd),
                           "even_depth_ppm": float(d_even * 1e6), "even_err_ppm": float(s_even * 1e6), "n_even_in": int(n_even),
                           "sigma_diff": float(oe_sig)}
    print(f"odd {d_odd*1e6:.0f}+-{s_odd*1e6:.0f}  even {d_even*1e6:.0f}+-{s_even*1e6:.0f}  -> {oe_sig:.2f} sigma", flush=True)

    # secondary at phase 0.5 (same duration), with error from random phases
    d_sec, e_sec, n_sec = box_depth(t2, f2, fe2, P2, T02, dur2, phase_center=0.5)
    results["secondary"] = {"depth_ppm": float(d_sec * 1e6), "empirical_err_ppm": float(emp_sig * 1e6),
                            "sigma": float(d_sec / emp_sig), "n_in_box": int(n_sec),
                            "max_allowed_3sig_ppm": float(3 * emp_sig * 1e6)}
    # also scan the whole phase range for the deepest non-primary box
    scan_ph = np.linspace(-0.5, 0.5, 400, endpoint=False)
    scan = []
    for pc in scan_ph:
        if abs(pc) < 2.0 * dur2 / P2:
            continue
        d, e, n = box_depth(t2, f2, fe2, P2, T02, dur2, phase_center=pc)
        scan.append((pc, d))
    scan = np.array(scan)
    k = int(np.nanargmax(scan[:, 1]))
    results["secondary"]["deepest_non_primary_box"] = {"phase": float(scan[k, 0]), "depth_ppm": float(scan[k, 1] * 1e6),
                                                       "sigma": float(scan[k, 1] / emp_sig)}
    print(f"secondary @0.5: {d_sec*1e6:.0f} ppm ({d_sec/emp_sig:.2f} sigma); deepest off-primary box {scan[k,1]*1e6:.0f} ppm at phase {scan[k,0]:.3f}", flush=True)

    # per-year depths. Kepler "years": Q1-Q4 (yr1), Q5-Q8 (yr2), Q9-Q12 (yr3), Q13-Q17 (yr4)
    yr_of_q = {1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 2, 7: 2, 8: 2, 9: 3, 10: 3, 11: 3, 12: 3, 13: 4, 14: 4, 15: 4, 16: 4, 17: 4}
    yr = np.array([yr_of_q[int(x)] for x in q2])
    per_year = []
    for y in (1, 2, 3, 4):
        s = yr == y
        d, e, n = box_depth(t2[s], f2[s], fe2[s], P2, T02, dur2)
        sy, _, _ = empirical_depth_error(t2[s], f2[s], fe2[s], P2, T02, dur2, n=150)
        per_year.append({"year": y, "quarters": sorted(set(int(x) for x in q2[s])), "depth_ppm": float(d * 1e6),
                         "empirical_err_ppm": float(sy * 1e6), "n_in_box": int(n),
                         "sigma": float(d / sy) if sy > 0 else None})
        print(f"  year {y}: {d*1e6:.0f} +- {sy*1e6:.0f} ppm (n_in={n})", flush=True)
    results["per_year"] = per_year
    yd = np.array([p["depth_ppm"] for p in per_year]); ye = np.array([p["empirical_err_ppm"] for p in per_year])
    w = 1 / ye ** 2
    chi2 = np.sum(w * (yd - np.sum(w * yd) / np.sum(w)) ** 2)
    results["per_year_chi2_const"] = {"chi2": float(chi2), "dof": 3}

    # individual events
    eps = transit_epochs(t2, P2, T02, dur2)
    events = []
    for (n, tm, nin, nexp, nloc) in eps:
        d, e, _ = per_event_depth(t2, f2, fe2, tm, dur2)
        qq = q2[np.argmin(np.abs(t2 - tm))] if len(t2) else -1
        # distance to nearest data gap (>0.5 d) edge
        gaps = [(t2[a], t2[b - 1]) for (a, b) in segments(t2)]
        near_edge = min(min(abs(tm - g0), abs(tm - g1)) for (g0, g1) in gaps)
        events.append({"epoch": n, "t_mid_bkjd": tm, "n_in_box": nin, "n_expected": round(nexp, 1),
                       "coverage": round(nin / nexp, 2), "depth_ppm": float(d * 1e6) if np.isfinite(d) else None,
                       "err_ppm": float(e * 1e6) if np.isfinite(e) else None,
                       "quarter": int(qq), "days_to_nearest_segment_edge": float(near_edge)})
    n_tr = sum(1 for ev in events if ev["coverage"] >= 0.5)
    n_tr_any = sum(1 for ev in events if ev["n_in_box"] > 0)
    results["events"] = events
    results["n_transits_cov50"] = int(n_tr)
    results["n_transits_any_data"] = int(n_tr_any)
    results["n_epochs_in_baseline"] = int(np.floor((t2.max() - T02) / P2) - np.ceil((t2.min() - T02) / P2) + 1)
    print(f"transits with >=50% coverage: {n_tr}; with any data: {n_tr_any}; epochs in baseline: {results['n_epochs_in_baseline']}", flush=True)

    # single-event dominance: SNR without the strongest event
    ev_ok = [ev for ev in events if ev["depth_ppm"] is not None and ev["coverage"] >= 0.5]
    if ev_ok:
        ev_sorted = sorted(ev_ok, key=lambda e: -e["depth_ppm"])
        strongest = ev_sorted[0]
        m_drop = np.abs(t2 - strongest["t_mid_bkjd"]) > 3 * dur2
        d_w, _, _ = box_depth(t2[m_drop], f2[m_drop], fe2[m_drop], P2, T02, dur2)
        s_w, _, _ = empirical_depth_error(t2[m_drop], f2[m_drop], fe2[m_drop], P2, T02, dur2, n=200)
        results["drop_strongest_event"] = {"dropped_epoch": strongest["epoch"], "dropped_depth_ppm": strongest["depth_ppm"],
                                           "depth_without_ppm": float(d_w * 1e6), "err_without_ppm": float(s_w * 1e6),
                                           "snr_without": float(d_w / s_w)}
        print(f"without strongest event (epoch {strongest['epoch']}): {d_w*1e6:.0f}+-{s_w*1e6:.0f} ppm", flush=True)

    # ------------------------------------------------------------------ inverted-flux null (quick, same grid)
    print("inverted-flux null BLS ...", flush=True)
    _, periods_i, res_i, dt_i = run_bls(t1, 2.0 - f1, fe1)
    pw_i = np.asarray(res_i.power)
    sde_i = sde_of(pw_i)
    pk_i = top_peaks(periods_i, pw_i, n=1)[0]
    results["inverted_null"] = {"sde_top": float(sde_i), "period_d": float(res_i.period[pk_i[2]]),
                                "depth_ppm": float(res_i.depth[pk_i[2]] * 1e6)}
    print(f"  inverted top: P={results['inverted_null']['period_d']:.4f} SDE={sde_i:.2f}", flush=True)

    # ------------------------------------------------------------------ SAP flux dilution check
    ts, fs, fes, qs = load("SAP_FLUX")
    phs = phase_fold(ts, P2, T02)
    inmask_s = np.abs(phs * P2) < 0.75 * dur2
    fds, feds, keeps, mads = detrend(ts, fs, fes, window=1.0, mask=inmask_s)
    keeps = keeps | inmask_s
    ts2, fs2, fes2 = ts[keeps], fds[keeps], feds[keeps]
    d_sap, e_sap_formal, n_sap = box_depth(ts2, fs2, fes2, P2, T02, dur2)
    s_sap, _, _ = empirical_depth_error(ts2, fs2, fes2, P2, T02, dur2, n=300)
    ratio_sap_psf = d_sap / dep
    ratio_err = ratio_sap_psf * np.hypot(s_sap / d_sap, emp_sig / dep)
    results["depth_sap"] = {"depth_ppm": float(d_sap * 1e6), "empirical_err_ppm": float(s_sap * 1e6),
                            "formal_err_ppm": float(e_sap_formal * 1e6), "n_in_box": int(n_sap),
                            "snr": float(d_sap / s_sap), "point_to_point_mad_ppm": float(mads * 1e6),
                            "ratio_sap_over_psf": float(ratio_sap_psf), "ratio_err": float(ratio_err)}
    # SAP on the raw per-quarter-normalised flux, same for PSF, for the fold figure
    print(f"SAP depth = {d_sap*1e6:.0f} +- {s_sap*1e6:.0f} ppm ; SAP/PSF = {ratio_sap_psf:.2f} +- {ratio_err:.2f}", flush=True)

    # ------------------------------------------------------------------ noise summary
    cad = np.median(np.diff(t2))
    # 3.4-h binned scatter (approx CDPP)
    nb = int(round(dur2 / cad))
    binned = []
    for (a, b) in segments(t2):
        seg = f2[a:b] - 1.0
        for i in range(0, len(seg) - nb + 1, nb):
            binned.append(np.mean(seg[i:i + nb]))
    results["noise"] = {"ptp_mad_ppm": float(mad2 * 1e6), "binned_dur_scatter_ppm": float(np.std(binned) * 1e6),
                        "n_bins": len(binned), "cadence_d": float(cad)}

    # ------------------------------------------------------------------ comparison with pipeline (numbers quoted by task)
    results["pipeline_quoted"] = {"period_d": PIPELINE_P, "depth_ppm": PIPELINE_DEPTH_PPM, "duration_h": PIPELINE_DUR_H}
    results["comparison"] = {
        "period_frac_diff": float((P2 - PIPELINE_P) / PIPELINE_P),
        "depth_ratio_mine_over_pipeline": float(dep * 1e6 / PIPELINE_DEPTH_PPM),
        "depth_diff_sigma_mine": float((dep * 1e6 - PIPELINE_DEPTH_PPM) / (emp_sig * 1e6)),
        "duration_diff_h": float(dur2 * 24 - PIPELINE_DUR_H),
    }

    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax = axes[0, 0]
    ax.plot(periods, power, lw=0.6, color="0.3")
    for k, pk in enumerate(results["bls"]["top3"]):
        ax.axvline(pk["period_d"], color=["C3", "C1", "C2"][k], ls="--", lw=0.8, alpha=0.8,
                   label=f"peak {k+1}: {pk['period_d']:.4f} d (SDE {pk['sde']:.1f})")
    ax.axvline(PIPELINE_P, color="k", ls=":", lw=0.8, label=f"pipeline P={PIPELINE_P}")
    ax.set_xscale("log"); ax.set_xlabel("period [d]"); ax.set_ylabel("BLS likelihood power")
    ax.set_title(f"{FNAME}: blind BLS 1-50 d (my grid, {len(periods)} periods)")
    ax.legend(fontsize=7)

    ax = axes[0, 1]
    phh = phase_fold(t2, P2, T02) * P2 * 24.0
    s = np.abs(phh) < 24.0
    ax.plot(phh[s], (f2[s] - 1) * 1e6, ".", ms=2, color="0.6", alpha=0.5)
    # binned
    bw = 0.5
    bins = np.arange(-24, 24.01, bw)
    idx = np.digitize(phh[s], bins)
    bx = []; by = []; be = []
    for i in range(1, len(bins)):
        m = idx == i
        if m.sum() > 3:
            bx.append(0.5 * (bins[i - 1] + bins[i])); by.append(np.mean(f2[s][m] - 1) * 1e6)
            be.append(np.std(f2[s][m]) / np.sqrt(m.sum()) * 1e6)
    ax.errorbar(bx, by, be, fmt="o", ms=4, color="C3", label="30-min bins")
    ax.axvspan(-dur2 * 12, dur2 * 12, color="C0", alpha=0.1, label=f"box {dur2*24:.2f} h")
    ax.axhline(-dep * 1e6, color="C0", lw=1, ls="--", label=f"depth {dep*1e6:.0f} ppm")
    ax.set_ylim(-4000, 3000)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("PSF flux - 1 [ppm]")
    ax.set_title(f"P={P2:.5f} d  T0={T02:.4f} BKJD  SNR_emp={dep/emp_sig:.1f}")
    ax.legend(fontsize=7)

    ax = axes[1, 0]
    # odd / even + secondary folds
    for lab, m, c in (("odd", odd, "C1"), ("even", ~odd, "C2")):
        phx = phase_fold(t2[m], P2, T02) * P2 * 24.0
        sx = np.abs(phx) < 24.0
        idx = np.digitize(phx[sx], bins)
        bx = []; by = []
        for i in range(1, len(bins)):
            mm = idx == i
            if mm.sum() > 3:
                bx.append(0.5 * (bins[i - 1] + bins[i])); by.append(np.mean(f2[m][sx][mm] - 1) * 1e6)
        ax.plot(bx, by, "o-", ms=3, lw=0.8, color=c, label=f"{lab}: {([d_odd, d_even][lab=='even'])*1e6:.0f} ppm")
    phx = (phase_fold(t2, P2, T02) - 0.5 + 0.5) % 1.0 - 0.5
    phx = phx * P2 * 24.0
    sx = np.abs(phx) < 24.0
    idx = np.digitize(phx[sx], bins)
    bx = []; by = []
    for i in range(1, len(bins)):
        mm = idx == i
        if mm.sum() > 3:
            bx.append(0.5 * (bins[i - 1] + bins[i])); by.append(np.mean(f2[sx][mm] - 1) * 1e6)
    ax.plot(bx, by, "s-", ms=3, lw=0.8, color="0.4", label=f"phase 0.5: {d_sec*1e6:.0f} ppm ({d_sec/emp_sig:.1f} sig)")
    ax.axvspan(-dur2 * 12, dur2 * 12, color="C0", alpha=0.1)
    ax.set_ylim(-2000, 1500); ax.set_xlabel("hours from centre"); ax.set_ylabel("binned flux - 1 [ppm]")
    ax.set_title(f"odd/even ({oe_sig:.1f} sigma) and secondary")
    ax.legend(fontsize=7)

    ax = axes[1, 1]
    xs = [ev["epoch"] for ev in events if ev["depth_ppm"] is not None]
    ys = [ev["depth_ppm"] for ev in events if ev["depth_ppm"] is not None]
    es = [ev["err_ppm"] for ev in events if ev["depth_ppm"] is not None]
    cov = [ev["coverage"] for ev in events if ev["depth_ppm"] is not None]
    sc = ax.scatter(xs, ys, c=cov, cmap="viridis", vmin=0, vmax=1, zorder=3, s=25)
    ax.errorbar(xs, ys, es, fmt="none", ecolor="0.6", lw=0.7)
    ax.axhline(dep * 1e6, color="C0", ls="--", label=f"mean PSF depth {dep*1e6:.0f} ppm")
    ax.axhline(d_sap * 1e6, color="C4", ls=":", label=f"SAP depth {d_sap*1e6:.0f} ppm")
    ax.axhline(0, color="k", lw=0.5)
    plt.colorbar(sc, ax=ax, label="box coverage")
    ax.set_xlabel("epoch"); ax.set_ylabel("per-event depth [ppm]")
    ax.set_title(f"per-event depths; per-year: " + ", ".join(f"{p['depth_ppm']:.0f}" for p in per_year) + " ppm")
    ax.legend(fontsize=7)
    fig.tight_layout()
    figpath = os.path.join(OUT_DIR, "fold_rederive.png")
    fig.savefig(figpath, dpi=130)
    results["figure"] = figpath

    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        return o

    jpath = os.path.join(OUT_DIR, "rederive_results.json")
    with open(jpath, "w") as fh:
        json.dump(_clean(results), fh, indent=2)
    print("wrote", jpath, "and", figpath)


if __name__ == "__main__":
    main()
