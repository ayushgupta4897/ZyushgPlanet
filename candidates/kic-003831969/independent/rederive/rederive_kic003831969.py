#!/usr/bin/env python3
"""
Independent re-derivation of the transit-like signal in KBONUS-BKG kic-003831969.

Written from scratch with numpy / scipy / astropy only.  Deliberately imports
NOTHING from /home/user/ZyushgPlanet/pipeline/ and uses no wotan / TLS / lightkurve
so that the result is an independent check on the pipeline numbers
(pipeline claim: P = 15.3134 d, depth 997 ppm, duration 3.8 h).

Steps
  1. load LIGHTCURVE_STITCHED, SAP_QUALITY == 0, normalise per quarter,
     detrend with a time-aware 1-day running median (iterative upward clipping)
  2. blind BLS 1-50 d (astropy BoxLeastSquares, own duration-matched frequency grid)
  3. refine the top peak; masked re-detrend; depth with empirical (random-phase-box)
     error, odd/even, secondary at phase 0.5, per-year depths, per-event depths,
     transit count, box duration + trapezoid fit
  4. same depth measurement on SAP_FLUX (dilution check)
  5. inverted-flux null (2 - f) through the identical blind BLS (optional, --null)
Outputs: results.json, events.csv, fold.png in the same directory as this script.
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
        "hlsp_kbonus-bkg_kepler_kepler_kic-003831969_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-003831969/rederive"
BKJD0 = 2454833.0
P_PIPE, DEPTH_PIPE_PPM, DUR_PIPE_H = 15.31340, 997.0, 3.8
RNG = np.random.default_rng(20261017)

# search settings (my own, not the pipeline's)
PMIN, PMAX = 1.0, 50.0
DURS_H = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0])      # blind-search box durations
DUR_MIN_GRID_H = 2.0                                          # sets the frequency step
GRID_OS = float(os.environ.get("GRID_OS", "2.0"))             # grid oversampling
DETREND_WINDOW_D = 1.0
RUN_NULL = "--null" in sys.argv

# figure palette (fixed-order categorical hues, muted text ink)
C_DATA, C_ODD, C_EVEN, C_MODEL, C_INK, C_MUTED = (
    "#9a9a96", "#2a78d6", "#eb6834", "#0b0b0b", "#0b0b0b", "#52514e")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# --------------------------------------------------------------------------- load
def load(col):
    with fits.open(FITS) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d[col], float)
        fe = np.asarray(d[col + "_ERR"], float)
        q = np.asarray(d["QUARTER"], int)
        qual = np.asarray(d["SAP_QUALITY"], int)
    n_all = len(t)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (f > 0)
    t, f, fe, q = t[m], f[m], fe[m], q[m]
    o = np.argsort(t)
    t, f, fe, q = t[o], f[o], fe[o], q[o]
    keep = np.r_[True, np.diff(t) > 0]
    t, f, fe, q = t[keep], f[keep], fe[keep], q[keep]
    for qq in np.unique(q):
        s = q == qq
        med = np.median(f[s])
        f[s] /= med
        fe[s] /= med
    return t, f, fe, q, n_all


# ----------------------------------------------------------------------- detrend
def running_median(t, f, window, use=None, min_pts=8):
    """Time-aware running median of width `window` (days); trend evaluated at every t
    but estimated only from points with use==True."""
    if use is None:
        use = np.ones(len(t), bool)
    tu, fu = t[use], f[use]
    lo = np.searchsorted(tu, t - window / 2.0)
    hi = np.searchsorted(tu, t + window / 2.0)
    tr = np.full(len(t), np.nan)
    for i in range(len(t)):
        if hi[i] - lo[i] >= min_pts:
            tr[i] = np.median(fu[lo[i]:hi[i]])
    return tr


def detrend(t, f, window=DETREND_WINDOW_D, mask=None, nclip=2):
    """Divide by a running median.  mask==True marks points excluded from the trend
    estimate (in-transit).  Upward outliers (>4 sigma, flares/cosmic rays) and
    isolated catastrophic drops (<-10 sigma) are removed iteratively."""
    use = np.ones(len(t), bool) if mask is None else ~mask
    good = np.ones(len(t), bool)
    for it in range(nclip + 1):
        tr = running_median(t, f, window, use=use & good)
        r = f / tr - 1.0
        sig = 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r)))
        bad = ~np.isfinite(tr) | (r > 4 * sig) | (r < -10 * sig)
        good = good & ~bad
    fd = f / tr
    return fd, tr, good, sig


# --------------------------------------------------------------------------- BLS
def period_grid(pmin, pmax, baseline, dur_min, os_=2.0):
    """Frequency grid whose step keeps the phase drift over the baseline below
    dur_min/os_ at every trial period: f_{i+1} = f_i (1 - dur_min/(T*os_))."""
    qq = dur_min / (baseline * os_)
    n = int(np.ceil(np.log(pmax / pmin) / -np.log(1.0 - qq))) + 2
    fr = (1.0 / pmin) * (1.0 - qq) ** np.arange(n)
    fr = fr[fr >= 1.0 / pmax]
    return np.sort(1.0 / fr)


def sde_stats(power):
    p = np.asarray(power, float)
    sde_raw = (p.max() - p.mean()) / p.std()
    # remove a slow baseline (running median over 2% of the grid) before the SDE
    w = max(51, int(0.02 * len(p)) | 1)
    from scipy.ndimage import median_filter
    base = median_filter(p, size=w, mode="nearest")
    pr = p - base
    sde_detr = (pr.max() - pr.mean()) / pr.std()
    return float(sde_raw), float(sde_detr), pr


def top_peaks(periods, power, n=3, sep=0.03):
    """Local maxima, greedily deduplicated so chosen periods differ by > sep."""
    p = np.asarray(power)
    loc = np.r_[False, (p[1:-1] > p[:-2]) & (p[1:-1] >= p[2:]), False]
    idx = np.where(loc)[0]
    idx = idx[np.argsort(p[idx])[::-1]]
    chosen = []
    for i in idx:
        if all(abs(periods[i] / periods[j] - 1) > sep for j in chosen):
            chosen.append(i)
        if len(chosen) == n:
            break
    return chosen


def alias_of(p, pref=P_PIPE, tol=0.005):
    for name, r in [("1", 1), ("2x", 2), ("3x", 3), ("1/2", 0.5), ("1/3", 1 / 3),
                    ("3/2", 1.5), ("2/3", 2 / 3), ("4x", 4), ("1/4", 0.25)]:
        if abs(p / (pref * r) - 1) < tol:
            return name
    return None


def blind_bls(t, f, sig, label):
    base = t.max() - t.min()
    periods = period_grid(PMIN, PMAX, base, DUR_MIN_GRID_H / 24.0, GRID_OS)
    log(f"[{label}] blind BLS: {len(periods)} periods, {len(DURS_H)} durations, "
        f"{len(t)} points")
    bls = BoxLeastSquares(t, f, dy=sig)
    t1 = time.time()
    r = bls.power(periods, DURS_H / 24.0, objective="likelihood", oversample=10)
    log(f"[{label}] done in {time.time()-t1:.1f} s")
    sde_raw, sde_detr, pr = sde_stats(r.power)
    peaks = top_peaks(periods, r.power, 3)
    out = []
    for k, i in enumerate(peaks):
        out.append(dict(rank=k + 1, period_d=float(periods[i]), power=float(r.power[i]),
                        sde_raw=float((r.power[i] - r.power.mean()) / r.power.std()),
                        sde_detr=float((pr[i] - pr.mean()) / pr.std()),
                        depth_ppm=float(r.depth[i] * 1e6),
                        depth_snr_formal=float(r.depth_snr[i]),
                        duration_h=float(r.duration[i] * 24),
                        t0_bkjd=float(r.transit_time[i]),
                        alias_of_pipeline_period=alias_of(periods[i])))
    return dict(n_periods=int(len(periods)), grid_oversample=GRID_OS,
                durations_h=DURS_H.tolist(), sde_raw=sde_raw, sde_detr=sde_detr,
                peaks=out), (periods, r.power, pr)


# ------------------------------------------------------------------ measurements
def phase_dt(t, P, t0):
    return (t - t0 + 0.5 * P) % P - 0.5 * P


def box_depth(t, f, P, t0, dur, ph=0.0):
    inn = np.abs(phase_dt(t, P, t0 + ph * P)) < dur / 2.0
    n = int(inn.sum())
    return (1.0 - f[inn].mean() if n else np.nan), n


def empirical_err(t, f, P, t0, dur, n=2000, exclude=(0.0,), excl_w=None, min_n=10):
    """std of box 'depths' measured at random phases away from the excluded ones."""
    if excl_w is None:
        excl_w = 1.5 * dur / P
    ds = []
    tries = 0
    while len(ds) < n and tries < 20 * n:
        tries += 1
        ph = RNG.uniform(0, 1)
        if any(min(abs(ph - e), 1 - abs(ph - e)) < excl_w for e in exclude):
            continue
        d, nin = box_depth(t, f, P, t0, dur, ph)
        if nin >= min_n:
            ds.append(d)
    ds = np.array(ds)
    return float(np.std(ds)), ds


def depth_block(t, f, P, t0, dur, n_rand=2000, exclude=(0.0,)):
    d, nin = box_depth(t, f, P, t0, dur)
    e_emp, ds = empirical_err(t, f, P, t0, dur, n=n_rand, exclude=exclude)
    return dict(depth_ppm=float(d * 1e6), err_emp_ppm=float(e_emp * 1e6),
                snr_emp=float(d / e_emp) if e_emp > 0 else np.nan, n_in=nin,
                n_random_boxes=int(len(ds)),
                random_box_depth_median_ppm=float(np.median(ds) * 1e6),
                random_box_depth_p99_ppm=float(np.percentile(ds, 99) * 1e6))


def trapezoid(dt, tc, depth, T14, T23):
    """Symmetric trapezoid in flux units (1 = out of transit)."""
    T23 = np.clip(T23, 0, T14)
    x = np.abs(dt - tc)
    m = np.ones_like(dt)
    full = x <= T23 / 2
    ing = (x > T23 / 2) & (x < T14 / 2)
    m[full] = 1 - depth
    m[ing] = 1 - depth * (T14 / 2 - x[ing]) / ((T14 - T23) / 2 + 1e-12)
    return m


def fit_trapezoid(dt, f, sig, dur0, depth0):
    sel = np.abs(dt) < 3 * dur0
    x, y = dt[sel], f[sel]

    def resid(p):
        return (y - trapezoid(x, p[0], p[1], p[2], p[3])) / sig

    p0 = [0.0, depth0, dur0, 0.5 * dur0]
    lb = [-dur0, 0.0, 0.25 * dur0, 0.0]
    ub = [dur0, 0.05, 4 * dur0, 4 * dur0]
    r = least_squares(resid, p0, bounds=(lb, ub))
    chi2_trap = float(np.sum(r.fun ** 2))
    chi2_flat = float(np.sum(((y - 1) / sig) ** 2))
    # box with same duration/depth for comparison
    chi2_box = float(np.sum(((y - trapezoid(x, 0, depth0, dur0, dur0)) / sig) ** 2))
    tc, dep, T14, T23 = r.x
    T23 = min(T23, T14)
    return dict(tc_shift_h=float(tc * 24), depth_ppm=float(dep * 1e6), T14_h=float(T14 * 24),
                T23_h=float(T23 * 24), T23_over_T14=float(T23 / T14),
                chi2_trap=chi2_trap, chi2_box=chi2_box, chi2_flat=chi2_flat, n=int(sel.sum()))


def per_event(t, f, P, t0, dur, sig, local_w=0.5):
    """Per-epoch local box depth: mean(local out-of-transit) - mean(in)."""
    n = np.round((t - t0) / P).astype(int)
    dt = t - (t0 + n * P)
    rows = []
    for k in np.unique(n):
        s = n == k
        inn = s & (np.abs(dt) < dur / 2)
        oot = s & (np.abs(dt) > dur / 2) & (np.abs(dt) < dur / 2 + local_w)
        nin, noot = int(inn.sum()), int(oot.sum())
        if nin == 0:
            continue
        d = (f[oot].mean() if noot >= 5 else 1.0) - f[inn].mean()
        rows.append(dict(epoch=int(k), tmid_bkjd=float(t0 + k * P), n_in=nin, n_oot=noot,
                         depth_ppm=float(d * 1e6), err_ppm=float(sig / np.sqrt(nin) * 1e6),
                         gap_near=bool(noot < 0.6 * (2 * local_w / (29.4 / 1440)))))
    return rows


def year_of(q):
    return np.where(q <= 4, 1, np.where(q <= 8, 2, np.where(q <= 12, 3, 4)))


# --------------------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    res = dict(target="kic-003831969", fits=os.path.basename(FITS), date="2026-10-17",
               method="independent re-derivation; numpy/scipy/astropy only; no pipeline/ code",
               settings=dict(pmin_d=PMIN, pmax_d=PMAX, durations_h=DURS_H.tolist(),
                             detrend="time-aware running median, window 1.0 d, "
                                     "iterative 4-sigma upward clip, -10 sigma drop clip",
                             quality="SAP_QUALITY == 0 only", normalisation="per-quarter median"))

    # 1. load + detrend (PSF flux)
    t, f, fe, q, n_all = load("FLUX")
    fd, tr, good, sig = detrend(t, f)
    log(f"loaded {n_all} cadences, {len(t)} with quality==0, {good.sum()} after clipping; "
        f"robust sigma = {sig*1e6:.0f} ppm")
    res["data"] = dict(n_cadences_file=int(n_all), n_quality0=int(len(t)),
                       n_used=int(good.sum()), baseline_d=float(t.max() - t.min()),
                       t_start_bkjd=float(t.min()), t_end_bkjd=float(t.max()),
                       robust_sigma_ppm_psf=float(sig * 1e6),
                       median_formal_err_ppm_psf=float(np.median(fe) * 1e6))
    tg, fg, qg = t[good], fd[good], q[good]

    # 2. blind BLS
    blind, (periods, power, pr) = blind_bls(tg, fg, sig, "PSF")
    res["blind_bls"] = blind
    for pk in blind["peaks"]:
        log(f"  peak {pk['rank']}: P={pk['period_d']:.5f} d  SDE={pk['sde_detr']:.1f} "
            f"(raw {pk['sde_raw']:.1f})  depth={pk['depth_ppm']:.0f} ppm  dur={pk['duration_h']:.1f} h  "
            f"alias={pk['alias_of_pipeline_period']}")
    best = blind["peaks"][0]
    res["recovered_blind"] = best["alias_of_pipeline_period"] == "1"
    res["pipeline_period_in_top3_or_alias"] = any(
        pk["alias_of_pipeline_period"] is not None for pk in blind["peaks"])

    # 3. refine the top peak (if the pipeline period is an alias, refine the pipeline
    #    period itself as well so the comparison is apples to apples)
    P_ref = best["period_d"]
    if best["alias_of_pipeline_period"] not in (None, "1"):
        log("top peak is an alias of the pipeline period; refining the pipeline period")
        P_ref = P_PIPE
    durs_fine = np.arange(1.0, 8.01, 0.25) / 24.0

    def refine(tt, ff, P0, width=0.003, n=3001):
        bls = BoxLeastSquares(tt, ff, dy=sig)
        per = np.linspace(P0 * (1 - width), P0 * (1 + width), n)
        r = bls.power(per, durs_fine, objective="likelihood", oversample=20)
        i = int(np.argmax(r.power))
        ok = r.power >= r.power[i] - 0.5          # delta loglike = 0.5 -> 1 sigma
        return dict(P=float(per[i]), t0=float(r.transit_time[i]), dur=float(r.duration[i]),
                    depth=float(r.depth[i]), depth_err=float(r.depth_err[i]),
                    P_lo=float(per[ok].min()), P_hi=float(per[ok].max()), power=float(r.power[i]))

    r1 = refine(tg, fg, P_ref)
    # masked re-detrend
    mask = np.abs(phase_dt(t, r1["P"], r1["t0"])) < 1.0 * r1["dur"]
    fd2, tr2, good2, sig2 = detrend(t, f, mask=mask)
    tg, fg, qg = t[good2], fd2[good2], q[good2]
    r2 = refine(tg, fg, r1["P"], width=0.001)
    P, t0, dur = r2["P"], r2["t0"], r2["dur"]
    # move t0 to the transit nearest the data midpoint (minimises P-t0 covariance)
    tmid = 0.5 * (t.min() + t.max())
    t0_mid = t0 + np.round((tmid - t0) / P) * P
    # local t0 scan at fixed P, dur (chi2 of a box)
    shifts = np.linspace(-0.25 * dur, 0.25 * dur, 201)
    chi = np.array([np.sum((fg[np.abs(phase_dt(tg, P, t0 + s)) < dur / 2] - (1 - r2["depth"])) ** 2)
                    for s in shifts])
    log(f"refined: P={P:.6f} d (1sig {r2['P_lo']:.6f}-{r2['P_hi']:.6f}), t0={t0:.5f} BKJD, "
        f"dur={dur*24:.2f} h, box depth={r2['depth']*1e6:.0f} ppm; sigma after masked "
        f"re-detrend={sig2*1e6:.0f} ppm")
    res["ephemeris"] = dict(period_d=P, period_1sig_lo_d=r2["P_lo"], period_1sig_hi_d=r2["P_hi"],
                            period_err_d=float(0.5 * (r2["P_hi"] - r2["P_lo"])),
                            t0_bkjd=float(t0), t0_bjd=float(t0 + BKJD0),
                            t0_mid_bkjd=float(t0_mid), t0_mid_bjd=float(t0_mid + BKJD0),
                            duration_box_h=float(dur * 24), bls_depth_ppm=float(r2["depth"] * 1e6),
                            bls_depth_err_formal_ppm=float(r2["depth_err"] * 1e6),
                            refine_before_mask=r1, robust_sigma_after_masked_detrend_ppm=float(sig2 * 1e6))

    # depth with empirical error, secondary, odd/even, per year
    prim = depth_block(tg, fg, P, t0, dur, n_rand=3000)
    sec = depth_block(tg, fg, P, t0 + 0.5 * P, dur, n_rand=1500, exclude=(0.0, 0.5))
    # significance of the secondary relative to the primary's random-box scatter
    sec["sigma_vs_primary_scatter"] = float(sec["depth_ppm"] / prim["err_emp_ppm"])
    # phase scan: deepest box anywhere except the primary
    phs = np.arange(0.0, 1.0, 0.0025)
    scan = np.array([box_depth(tg, fg, P, t0, dur, ph)[0] for ph in phs])
    away = np.minimum(phs, 1 - phs) > 1.5 * dur / P
    imax = int(np.argmax(np.where(away, scan, -np.inf)))
    nn = np.round((tg - t0) / P).astype(int)
    odd_m, even_m = (nn % 2 == 1), (nn % 2 == 0)
    odd = depth_block(tg[odd_m], fg[odd_m], P, t0, dur, n_rand=1500)
    even = depth_block(tg[even_m], fg[even_m], P, t0, dur, n_rand=1500)
    oe_sigma = abs(odd["depth_ppm"] - even["depth_ppm"]) / np.hypot(odd["err_emp_ppm"], even["err_emp_ppm"])
    years = {}
    yr = year_of(qg)
    for y in (1, 2, 3, 4):
        s = yr == y
        if s.sum() < 100:
            continue
        yb = depth_block(tg[s], fg[s], P, t0, dur, n_rand=1000)
        yb["quarters"] = sorted(int(x) for x in np.unique(qg[s]))
        years[f"year{y}"] = yb
    yd = np.array([years[k]["depth_ppm"] for k in years])
    ye = np.array([years[k]["err_emp_ppm"] for k in years])
    wmean = np.sum(yd / ye ** 2) / np.sum(1 / ye ** 2)
    chi2_years = float(np.sum(((yd - wmean) / ye) ** 2))
    events = per_event(tg, fg, P, t0, dur, sig2)
    exp_in = dur / (29.4 / 1440.0)
    n_tr = sum(1 for e in events if e["n_in"] >= 0.5 * exp_in)
    n_tr_any = len(events)
    ed = np.array([e["depth_ppm"] for e in events if e["n_in"] >= 0.5 * exp_in])
    ee = np.array([e["err_ppm"] for e in events if e["n_in"] >= 0.5 * exp_in])
    ll = (ed / ee) ** 2 * np.sign(ed)                     # per-event chi2 contribution (signed)
    dominance = float(np.max(ll) / np.sum(ll[ll > 0])) if np.any(ll > 0) else np.nan
    chi2_events = float(np.sum(((ed - np.sum(ed / ee**2) / np.sum(1 / ee**2)) / ee) ** 2))
    trap = fit_trapezoid(phase_dt(tg, P, t0), fg, sig2, dur, prim["depth_ppm"] * 1e-6)
    log(f"depth {prim['depth_ppm']:.0f} +/- {prim['err_emp_ppm']:.0f} ppm (SNR_emp {prim['snr_emp']:.1f}); "
        f"secondary {sec['depth_ppm']:.0f} +/- {sec['err_emp_ppm']:.0f}; odd {odd['depth_ppm']:.0f} "
        f"even {even['depth_ppm']:.0f} ({oe_sigma:.1f} sig); transits {n_tr}/{n_tr_any}; "
        f"T14 {trap['T14_h']:.2f} h, T23/T14 {trap['T23_over_T14']:.2f}")
    res["primary"] = prim
    res["secondary_phase05"] = sec
    res["phase_scan"] = dict(deepest_non_primary_phase=float(phs[imax]), deepest_non_primary_ppm=float(scan[imax] * 1e6),
                             primary_over_deepest_other=float(prim["depth_ppm"] / (scan[imax] * 1e6)) if scan[imax] > 0 else None)
    res["odd_even"] = dict(odd=odd, even=even, sigma=float(oe_sigma))
    res["per_year"] = years
    res["per_year_consistency"] = dict(weighted_mean_ppm=float(wmean), chi2=chi2_years, dof=int(len(yd) - 1))
    res["events"] = dict(n_transits_with_ge_half_coverage=int(n_tr), n_epochs_with_any_data=int(n_tr_any),
                         expected_in_transit_points=float(exp_in), n_epochs_in_baseline=int(np.ceil((t.max() - t.min()) / P)),
                         single_event_dominance=dominance, chi2_const_depth=chi2_events, dof=int(len(ed) - 1),
                         max_event_depth_ppm=float(ed.max()), max_event_epoch=int([e for e in events if e["n_in"] >= 0.5 * exp_in][int(np.argmax(ed))]["epoch"]))
    res["trapezoid"] = trap
    res["t0_scan"] = dict(best_shift_h=float(shifts[int(np.argmin(chi))] * 24),
                          delta_chi2_over_scan=float((chi.max() - chi.min()) / sig2 ** 2))
    with open(os.path.join(OUT, "events.csv"), "w") as fh:
        fh.write("epoch,tmid_bkjd,n_in,n_oot,depth_ppm,err_ppm,gap_near\n")
        for e in events:
            fh.write(f"{e['epoch']},{e['tmid_bkjd']:.5f},{e['n_in']},{e['n_oot']},{e['depth_ppm']:.1f},{e['err_ppm']:.1f},{int(e['gap_near'])}\n")

    # 4. SAP flux depth (dilution check)
    ts, fs, fes, qs, _ = load("SAP_FLUX")
    masks = np.abs(phase_dt(ts, P, t0)) < 1.0 * dur
    fsd, _, goods, sigs = detrend(ts, fs, mask=masks)
    sap = depth_block(ts[goods], fsd[goods], P, t0, dur, n_rand=2000)
    sap["robust_sigma_ppm"] = float(sigs * 1e6)
    sap["median_sap_over_psf_flux"] = float(np.median(np.asarray(fits.getdata(FITS, "LIGHTCURVE_STITCHED")["SAP_FLUX"], float)) /
                                            np.median(np.asarray(fits.getdata(FITS, "LIGHTCURVE_STITCHED")["FLUX"], float)))
    ratio = sap["depth_ppm"] / prim["depth_ppm"]
    ratio_err = ratio * np.hypot(sap["err_emp_ppm"] / sap["depth_ppm"], prim["err_emp_ppm"] / prim["depth_ppm"])
    res["sap"] = sap
    res["sap_vs_psf"] = dict(depth_ratio_sap_over_psf=float(ratio), depth_ratio_err=float(ratio_err),
                             crowdsap_catalog=0.762627, flfrcsap_catalog=0.733899, psffrac_catalog=0.754219,
                             expected_ratio_if_on_target_approx="~crowdsap (0.76) if PSF flux is pure target flux; "
                                                                 ">1 would indicate the dip is stronger in the "
                                                                 "contaminated aperture, i.e. likely from a neighbour")
    log(f"SAP depth {sap['depth_ppm']:.0f} +/- {sap['err_emp_ppm']:.0f} ppm; ratio SAP/PSF = {ratio:.2f} +/- {ratio_err:.2f}")

    # 5. inverted-flux null
    if RUN_NULL:
        null, _ = blind_bls(tg, 2.0 - fg, sig2, "NULL 2-f")
        res["inverted_null"] = null
        log(f"null top SDE {null['peaks'][0]['sde_detr']:.1f} at P={null['peaks'][0]['period_d']:.4f} d")

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    a = ax[0, 0]
    a.plot(periods, pr, color=C_DATA, lw=0.6)
    for pk in blind["peaks"]:
        a.axvline(pk["period_d"], color=C_ODD, lw=0.8, alpha=0.6)
        a.text(pk["period_d"], pr.max() * (1.02 - 0.07 * (pk["rank"] - 1)), f"#{pk['rank']} {pk['period_d']:.3f} d",
               color=C_INK, fontsize=8, ha="left")
    a.axvline(P_PIPE, color=C_EVEN, lw=1.0, ls="--", label=f"pipeline P = {P_PIPE} d")
    a.set_xscale("log"); a.set_xlabel("period (d)"); a.set_ylabel("BLS power (baseline-subtracted)")
    a.set_title(f"blind BLS 1-50 d, PSF flux: SDE = {blind['sde_detr']:.1f}", fontsize=10)
    a.legend(fontsize=8, frameon=False)

    a = ax[0, 1]
    ph = ((tg - t0) / P + 0.5) % 1 - 0.5
    nb = 200
    be = np.linspace(-0.5, 0.5, nb + 1)
    bi = np.digitize(ph, be) - 1
    bm = np.array([np.mean(fg[bi == i]) if np.any(bi == i) else np.nan for i in range(nb)])
    a.plot(ph, (fg - 1) * 1e6, ".", color=C_DATA, ms=1, alpha=0.25, rasterized=True)
    a.plot(0.5 * (be[1:] + be[:-1]), (bm - 1) * 1e6, "-", color=C_ODD, lw=1.5, label=f"binned ({nb} bins)")
    a.axvspan(-0.5 * dur / P, 0.5 * dur / P, color=C_EVEN, alpha=0.12, lw=0)
    a.axvline(0.5, color=C_MUTED, lw=0.8, ls=":"); a.axvline(-0.5, color=C_MUTED, lw=0.8, ls=":")
    a.set_ylim(-3000, 3000); a.set_xlabel("phase"); a.set_ylabel("relative flux (ppm)")
    a.set_title(f"full fold P = {P:.5f} d; secondary at 0.5: {sec['depth_ppm']:.0f} +/- {sec['err_emp_ppm']:.0f} ppm", fontsize=10)
    a.legend(fontsize=8, frameon=False, loc="lower left")

    a = ax[1, 0]
    dt_h = phase_dt(tg, P, t0) * 24
    w = 4 * dur * 24
    s = np.abs(dt_h) < w
    a.plot(dt_h[s], (fg[s] - 1) * 1e6, ".", color=C_DATA, ms=2, alpha=0.3, rasterized=True)
    for m_, c_, lab in ((odd_m, C_ODD, f"odd  {odd['depth_ppm']:.0f} ppm"), (even_m, C_EVEN, f"even {even['depth_ppm']:.0f} ppm")):
        be2 = np.arange(-w, w + 1e-9, 0.5)
        ss = s & m_
        bi2 = np.digitize(dt_h[ss], be2) - 1
        bm2 = np.array([np.mean(fg[ss][bi2 == i]) if np.sum(bi2 == i) >= 3 else np.nan for i in range(len(be2) - 1)])
        a.plot(0.5 * (be2[1:] + be2[:-1]), (bm2 - 1) * 1e6, "o-", color=c_, ms=3.5, lw=1.2, label=lab)
    xx = np.linspace(-w, w, 800) / 24
    a.plot(xx * 24, (trapezoid(xx, trap["tc_shift_h"] / 24, trap["depth_ppm"] * 1e-6, trap["T14_h"] / 24, trap["T23_h"] / 24) - 1) * 1e6,
           color=C_MODEL, lw=1.2, label=f"trapezoid T14 {trap['T14_h']:.1f} h")
    a.set_ylim(-2500, 1500); a.set_xlabel("time from mid-transit (h)"); a.set_ylabel("relative flux (ppm)")
    a.set_title(f"depth {prim['depth_ppm']:.0f} +/- {prim['err_emp_ppm']:.0f} ppm (empirical), odd/even {oe_sigma:.1f} sigma", fontsize=10)
    a.legend(fontsize=8, frameon=False, loc="lower left")

    a = ax[1, 1]
    et = np.array([e["tmid_bkjd"] for e in events]); edp = np.array([e["depth_ppm"] for e in events])
    eer = np.array([e["err_ppm"] for e in events]); full = np.array([e["n_in"] >= 0.5 * exp_in for e in events])
    a.errorbar(et[full], edp[full], eer[full], fmt="o", color=C_ODD, ms=4, lw=1, capsize=2, label="per-transit local depth")
    if np.any(~full):
        a.errorbar(et[~full], edp[~full], eer[~full], fmt="o", mfc="none", color=C_MUTED, ms=4, lw=1, capsize=2, label="partial coverage")
    a.axhline(prim["depth_ppm"], color=C_MODEL, lw=1, label="global depth")
    a.axhline(0, color=C_MUTED, lw=0.8, ls=":")
    for k, (yk, yv) in enumerate(years.items()):
        a.text(0.02 + 0.24 * k, 0.95, f"{yk}: {yv['depth_ppm']:.0f}+/-{yv['err_emp_ppm']:.0f}", transform=a.transAxes, fontsize=8, color=C_INK)
    a.set_xlabel("mid-transit time (BKJD)"); a.set_ylabel("depth (ppm)")
    a.set_title(f"{n_tr} transits with >=50% coverage; SAP/PSF depth ratio {ratio:.2f} +/- {ratio_err:.2f}", fontsize=10)
    a.legend(fontsize=8, frameon=False, loc="lower left")
    for aa in ax.ravel():
        aa.tick_params(labelsize=8); aa.grid(alpha=0.15)
        for sp in ("top", "right"):
            aa.spines[sp].set_visible(False)
    fig.suptitle("kic-003831969 independent re-derivation (no pipeline code)", fontsize=11, color=C_INK)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fold.png"), dpi=130)

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        return o
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(clean(res), fh, indent=1)
    log("wrote results.json, events.csv, fold.png")


if __name__ == "__main__":
    main()
