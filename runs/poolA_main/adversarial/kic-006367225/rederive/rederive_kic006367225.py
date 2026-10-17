#!/usr/bin/env python
"""
Independent re-derivation of the transit-like signal in KBONUS-BKG kic-006367225.

Written from scratch with numpy / scipy / astropy only.  Nothing under
ZyushgPlanet/pipeline/ is imported or read.  No wotan, no TLS.

Two passes are run and both are stored:

  pass1_raw     per-quarter median normalisation, 1-d time-based running-median
                detrend (segments split only at data gaps > 0.5 d), one-sided
                4-sigma clip of POSITIVE outliers only.  No negative clipping.
  pass2_robust  as pass 1, plus: detrending segments are also split at step
                discontinuities found by a generic jump detector (median of 6
                cadences before vs after differs by > 10 x point-to-point sigma;
                +-0.05 d guard dropped around each step), and a symmetric 5-sigma
                robust clip of the detrended flux.  Motivated by inspection of
                Q2 (two +6% / +12% one-cadence jumps at BKJD 200.33 and 246.2
                in both PSF and SAP flux); it is period-agnostic.  A 3000-ppm
                transit (~1 sigma per point) cannot trigger either rule.

For each pass:
  blind BLS 1-50 d (astropy BoxLeastSquares, own log-spaced grid) on the flux
  and on the inverted flux 2-f; refinement at the top peak; box depth with a
  random-phase empirical error; odd/even; secondary at phase 0.5; per-year,
  per-season and per-transit depths; trapezoid shape; the same statistics
  forced at the pipeline period 12.9715 d; SAP_FLUX depth at both ephemerides.

Run:  python -I rederive_kic006367225.py            (everything)
      python -I rederive_kic006367225.py --plot-only (figure from cached npz)
"""
import json
import os
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.optimize import least_squares

FNAME = "kic-006367225"
FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-006367225_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006367225/rederive"
SCR = ("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/"
       "scratchpad/rederive-kic-006367225")
os.makedirs(OUT, exist_ok=True)
os.makedirs(SCR, exist_ok=True)

BKJD0 = 2454833.0
P_PIPE = 12.97150          # pipeline value we are asked to test against
CROWDSAP = 0.8019833545    # from targets/poolA_kbonus_bg_dwarfs_ranked.csv
PSFFRAC = 0.5182924833978823
RSTAR = 0.8787632417678832
CADENCE_D = 29.4244 / 1440.0

RNG = np.random.default_rng(20261017)


# ----------------------------------------------------------------------------
# 1. loading and detrending
# ----------------------------------------------------------------------------
def load(col="FLUX"):
    d = fits.open(FITS)["LIGHTCURVE_STITCHED"].data
    t = np.asarray(d["TIME"], float)
    f = np.asarray(d[col], float)
    q = np.asarray(d["QUARTER"], int)
    qual = np.asarray(d["SAP_QUALITY"], int)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(f) & (f > 0)
    o = np.argsort(t[m])
    return t[m][o], f[m][o], q[m][o]


def normalise_quarters(f, q):
    fn = np.empty_like(f)
    for qq in np.unique(q):
        s = q == qq
        fn[s] = f[s] / np.median(f[s])
    return fn


def find_steps(t, fn, q, nside=6, thresh=10.0, gap=0.5):
    """Generic step detector.  At each cadence boundary inside a quarter compare
    the median of `nside` cadences before and after; a step is flagged when the
    jump exceeds `thresh` x sigma_ptp (point-to-point scatter, immune to slow
    excursions) and both sides are internally quiet.  Steps within 0.3 d are
    merged (largest kept).  Returns list of (quarter, t_before, t_after, jump/sig)."""
    steps = []
    for Q in np.unique(q):
        s = np.where(q == Q)[0]
        tq, fq = t[s], fn[s]
        sig = 1.4826 * np.median(np.abs(np.diff(fq))) / np.sqrt(2)
        for k in range(nside, len(s) - nside):
            if tq[k] - tq[k - 1] > gap:
                continue
            a, b = fq[k - nside:k], fq[k:k + nside]
            jump = np.median(b) - np.median(a)
            if abs(jump) > thresh * sig and np.std(a) < 4 * sig and np.std(b) < 4 * sig:
                steps.append((int(Q), float(tq[k - 1]), float(tq[k]), float(jump / sig)))
    merged = []
    for st in sorted(steps, key=lambda x: x[1]):
        if merged and st[1] - merged[-1][2] < 0.3:
            if abs(st[3]) > abs(merged[-1][3]):
                merged[-1] = st
        else:
            merged.append(st)
    return merged


def running_median(t, y, win=1.0, gap=0.5, mask=None, breaks=()):
    """Time-based running median of width `win` days, evaluated at every t.
    Segments are split at gaps > `gap` d and at the times in `breaks`.  If
    `mask` (bool, True = use) is given, only masked-in points enter the median
    (used to exclude transits)."""
    out = np.full_like(y, np.nan)
    cut = np.diff(t) > gap
    for b in breaks:
        k = np.searchsorted(t, b)
        if 0 < k < len(t):
            cut[k - 1] = True
    brk = np.where(cut)[0] + 1
    for seg in np.split(np.arange(len(t)), brk):
        ts, ys = t[seg], y[seg]
        use = np.ones(len(seg), bool) if mask is None else mask[seg]
        tu, yu = ts[use], ys[use]
        if len(tu) < 5:
            continue
        lo = np.searchsorted(tu, ts - win / 2)
        hi = np.searchsorted(tu, ts + win / 2, side="right")
        for k in range(len(seg)):
            if hi[k] - lo[k] >= 3:
                out[seg[k]] = np.median(yu[lo[k]:hi[k]])
    return out


def detrend(t, f, q, win=1.0, transit_mask=None, clip=4.0, robust=False,
            step_thresh=10.0, sym_clip=5.0, guard=0.05):
    """Returns t, y (flat, ~1), q, dy (per-quarter robust scatter), info."""
    fn = normalise_quarters(f, q)
    info = dict(steps=[])
    breaks = ()
    keep0 = np.ones(len(t), bool)
    if robust:
        steps = find_steps(t, fn, q, thresh=step_thresh)
        info["steps"] = steps
        breaks = [0.5 * (s[1] + s[2]) for s in steps]
        for b in breaks:
            keep0 &= np.abs(t - b) > guard
    t, fn, q = t[keep0], fn[keep0], q[keep0]
    if transit_mask is not None:
        transit_mask = transit_mask[keep0]
    tr = running_median(t, fn, win=win, mask=transit_mask, breaks=breaks)
    ok = np.isfinite(tr)
    t, fn, q, tr = t[ok], fn[ok], q[ok], tr[ok]
    y = fn / tr
    mad = 1.4826 * np.median(np.abs(y - np.median(y)))
    keep = y < 1 + clip * mad                     # one-sided: positive outliers only
    if robust and sym_clip:
        keep &= np.abs(y - 1) < sym_clip * mad     # symmetric 5-sigma (|dy| > ~1.3 %)
    info["n_clipped"] = int((~keep).sum())
    info["mad_ppm"] = float(mad * 1e6)
    t, y, q = t[keep], y[keep], q[keep]
    dy = np.empty_like(y)
    for qq in np.unique(q):
        s = q == qq
        dy[s] = 1.4826 * np.median(np.abs(y[s] - np.median(y[s])))
    return t, y, q, dy, info


# ----------------------------------------------------------------------------
# 2. blind BLS
# ----------------------------------------------------------------------------
def period_grid(T, pmin=1.0, pmax=50.0, dur_ref=0.08, oversample=2):
    """Log-spaced grid: successive periods differ by a factor 1 + dur/(T*os),
    so the accumulated phase drift across the baseline is < dur/os."""
    fac = 1.0 + dur_ref / (T * oversample)
    n = int(np.log(pmax / pmin) / np.log(fac)) + 1
    return pmin * fac ** np.arange(n)


DURS = np.array([0.05, 0.08, 0.12, 0.18, 0.25])  # days (1.2 h .. 6 h)


def sde_of(power):
    return float((power.max() - power.mean()) / power.std())


def robust_sde_of(power):
    med = np.median(power)
    return float((power.max() - med) / (1.4826 * np.median(np.abs(power - med))))


def top_peaks(periods, power, n=3, sep=0.02):
    idx = np.argsort(power)[::-1]
    picked = []
    for i in idx:
        if all(abs(periods[i] / periods[j] - 1) > sep for j in picked):
            picked.append(i)
        if len(picked) == n:
            break
    return picked


def alias_of(p, pref=P_PIPE, tol=0.01):
    for num, den in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2),
                     (1, 4), (4, 1), (3, 4), (4, 3)]:
        if abs(p / (pref * num / den) - 1) < tol:
            return f"{num}/{den} x P_pipe" if (num, den) != (1, 1) else "= P_pipe"
    return None


def blind_bls(t, y, dy, label):
    T = t.max() - t.min()
    periods = period_grid(T)
    bls = BoxLeastSquares(t, y, dy)
    t0 = time.time()
    r = bls.power(periods, DURS, objective="likelihood")
    dt = time.time() - t0
    pw = np.asarray(r.power)
    pk = top_peaks(periods, pw)
    peaks = []
    for i in pk:
        peaks.append(dict(period_d=float(periods[i]), power=float(pw[i]),
                          depth_ppm=float(r.depth[i] * 1e6),
                          duration_h=float(r.duration[i] * 24),
                          t0_bkjd=float(r.transit_time[i]),
                          sde_local=float((pw[i] - pw.mean()) / pw.std()),
                          alias_of_pipeline=alias_of(periods[i])))
    # the pipeline period in this periodogram
    w = np.abs(periods / P_PIPE - 1) < 0.003
    j = np.where(w)[0][np.argmax(pw[w])]
    allpk = top_peaks(periods, pw, n=60, sep=0.02)
    rank = next((k + 1 for k, i in enumerate(allpk) if abs(periods[i] / P_PIPE - 1) < 0.01), None)
    at_pipe = dict(period_d=float(periods[j]), power=float(pw[j]),
                   sde_local=float((pw[j] - pw.mean()) / pw.std()),
                   robust_sde_local=float((pw[j] - np.median(pw)) / (1.4826 * np.median(np.abs(pw - np.median(pw))))),
                   depth_ppm=float(r.depth[j] * 1e6), duration_h=float(r.duration[j] * 24),
                   rank_among_distinct_peaks=rank, power_ratio_to_top=float(pw[j] / pw.max()))
    print(f"[{label}] BLS {len(periods)} periods x {len(DURS)} durations in {dt:.1f}s; "
          f"SDE={sde_of(pw):.2f} robustSDE={robust_sde_of(pw):.2f}")
    for p in peaks:
        print("   ", p)
    print("    at P_pipe:", at_pipe)
    np.savez(f"{SCR}/bls_{label}.npz", periods=periods, power=pw,
             depth=r.depth, duration=r.duration, t0=r.transit_time)
    return dict(n_periods=int(len(periods)), durations_d=DURS.tolist(),
                pmin_d=1.0, pmax_d=50.0, sde=sde_of(pw), robust_sde=robust_sde_of(pw),
                top_peaks=peaks, at_pipeline_period=at_pipe, runtime_s=dt)


# ----------------------------------------------------------------------------
# 3. ephemeris refinement and depth statistics
# ----------------------------------------------------------------------------
def refine(t, y, dy, p0, frac=0.003, nper=3001):
    periods = np.linspace(p0 * (1 - frac), p0 * (1 + frac), nper)
    durs = np.arange(0.02, 0.36, 0.01)
    bls = BoxLeastSquares(t, y, dy)
    r = bls.power(periods, durs, objective="likelihood")
    i = int(np.argmax(r.power))
    P, dur, tt = float(periods[i]), float(r.duration[i]), float(r.transit_time[i])
    tmid = 0.5 * (t.min() + t.max())
    tt = float(tt + np.round((tmid - tt) / P) * P)
    st = bls.compute_stats(P, dur, tt)
    half = r.power[i] - 0.5 * (r.power[i] - np.median(r.power))
    above = np.where(r.power >= half)[0]
    j0, j1 = above[above <= i].min(), above[above >= i].max()
    Perr = float(0.5 * (periods[j1] - periods[j0]))
    return dict(P=P, dur=dur, t0=tt, P_err_hwhm=Perr,
                depth_bls=float(r.depth[i]), depth_bls_err=float(r.depth_err[i]),
                stats=st), periods, np.asarray(r.power)


def phase_of(t, P, t0):
    return ((t - t0) / P + 0.5) % 1.0 - 0.5     # in [-0.5, 0.5)


def box_depth(t, y, P, t0, dur, phase_c=0.0, oot_excl=1.5):
    """Depth = mean(out-of-transit) - mean(in-box) for a box of width `dur`
    centred at phase `phase_c`.  Out-of-transit excludes +-oot_excl*dur
    around BOTH the box and the primary (phase 0)."""
    ph = phase_of(t, P, t0)
    dph = ((ph - phase_c + 0.5) % 1.0) - 0.5
    half = 0.5 * dur / P
    inb = np.abs(dph) < half
    oot = (np.abs(dph) > oot_excl * half) & (np.abs(ph) > oot_excl * half)
    if inb.sum() < 3:
        return np.nan, int(inb.sum())
    return float(np.mean(y[oot]) - np.mean(y[inb])), int(inb.sum())


def random_phase_depths(t, y, P, t0, dur, n=3000, avoid=(0.0, 0.5), avoid_w=2.0):
    """Box depth at n random phases (same P, same width) avoiding the primary
    and secondary windows.  Its robust std is the empirical depth error."""
    half = 0.5 * dur / P
    out = []
    while len(out) < n:
        pc = RNG.uniform(-0.5, 0.5)
        if any(abs(((pc - a + 0.5) % 1.0) - 0.5) < avoid_w * half * 2 for a in avoid):
            continue
        d, nin = box_depth(t, y, P, t0, dur, phase_c=pc)
        if np.isfinite(d):
            out.append(d)
    return np.array(out)


def robust_std(a):
    a = a[np.isfinite(a)]
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def per_transit(t, y, q, P, t0, dur, oot_win=1.0, oot_excl=1.5):
    """Per-epoch local depth: local OOT = points within +-oot_win d of the
    transit centre but outside +-oot_excl*dur.  Requires >= 50% of the expected
    in-transit cadences (and >= 2) to be present."""
    epochs = np.arange(np.floor((t.min() - t0) / P), np.ceil((t.max() - t0) / P) + 1)
    exp_n = dur / CADENCE_D
    rows = []
    for n in epochs:
        tc = t0 + n * P
        din = np.abs(t - tc) < dur / 2
        doot = (np.abs(t - tc) < oot_win) & (np.abs(t - tc) > oot_excl * dur / 2)
        nin = int(din.sum())
        if nin < max(2, 0.5 * exp_n) or doot.sum() < 20:
            continue
        d = float(np.mean(y[doot]) - np.mean(y[din]))
        e = float(np.std(y[doot]) / np.sqrt(nin))
        tt = t[(np.abs(t - tc) < 3.0)]
        gaps = np.where(np.diff(tt) > 0.3)[0]
        gap_dist = float(min([min(abs(tt[g] - tc), abs(tt[g + 1] - tc)) for g in gaps] + [3.0]))
        rows.append(dict(epoch=int(n), tc_bkjd=float(tc), tc_bjd=float(tc + BKJD0),
                         quarter=int(np.median(q[din])), n_in=nin, depth_ppm=d * 1e6,
                         err_ppm=e * 1e6, gap_dist_d=gap_dist))
    return rows


def trapezoid(params, x):
    depth, t14, t23 = params
    t14 = abs(t14)
    t23 = min(abs(t23), t14)
    ax = np.abs(x)
    y = np.ones_like(x)
    full = ax <= t23 / 2
    ing = (ax > t23 / 2) & (ax < t14 / 2)
    y[full] = 1 - depth
    y[ing] = 1 - depth * (t14 / 2 - ax[ing]) / max((t14 - t23) / 2, 1e-6)
    return y


def fit_trapezoid(t, y, P, t0, dur):
    ph = phase_of(t, P, t0) * P          # days from mid-transit
    w = np.abs(ph) < 3 * dur
    x, yy = ph[w], y[w]

    def resid(p):
        return trapezoid(p, x) - yy
    best = None
    for t23f in (0.0, 0.5, 0.9):
        r = least_squares(resid, [3e-3, dur, t23f * dur], method="lm")
        if best is None or r.cost < best.cost:
            best = r
    depth, t14, t23 = best.x
    t14 = abs(t14)
    t23 = min(abs(t23), t14)
    return dict(depth_ppm=float(depth * 1e6), T14_h=float(t14 * 24), T23_h=float(t23 * 24),
                shape_T23_over_T14=float(t23 / t14 if t14 > 0 else np.nan),
                n_points=int(w.sum()), rms_resid_ppm=float(np.std(resid(best.x)) * 1e6))


def subset_depth(t, y, P, t0, dur, sel, nrand=1500):
    if sel.sum() < 50:
        return dict(depth_ppm=np.nan, err_ppm=np.nan, n_in=0, n_points=int(sel.sum()))
    d, nin = box_depth(t[sel], y[sel], P, t0, dur)
    if nin < 3:
        return dict(depth_ppm=np.nan, err_ppm=np.nan, n_in=nin, n_points=int(sel.sum()))
    rd = random_phase_depths(t[sel], y[sel], P, t0, dur, n=nrand)
    return dict(depth_ppm=d * 1e6, err_ppm=robust_std(rd) * 1e6, n_in=nin,
                n_points=int(sel.sum()))


def ephemeris_stats(t_raw, f_raw, q_raw, P, t0, dur, robust, label, ts_raw, fs_raw, qs_raw):
    """All depth statistics at a given ephemeris, on a masked re-detrend
    (in-transit points excluded from the running median)."""
    E = {}
    ph = phase_of(t_raw, P, t0)
    tmask = ~(np.abs(ph) * P < 1.5 * dur)
    t2, y2, q2, dy2, info = detrend(t_raw, f_raw, q_raw, win=1.0, transit_mask=tmask, robust=robust)
    # re-refine locally on the masked data
    ref2, rp2, rpow2 = refine(t2, y2, dy2, P, frac=0.001, nper=2001)
    P, dur, t0 = ref2["P"], ref2["dur"], ref2["t0"]
    st = ref2["stats"]
    E["ephemeris"] = dict(period_d=P, period_err_hwhm_d=ref2["P_err_hwhm"], t0_bkjd=t0,
                          t0_bjd=t0 + BKJD0, duration_h=dur * 24,
                          bls_depth_ppm=ref2["depth_bls"] * 1e6,
                          bls_depth_err_ppm=ref2["depth_bls_err"] * 1e6,
                          bls_depth_snr=float(ref2["depth_bls"] / ref2["depth_bls_err"]),
                          astropy_stats=dict(
                              depth_odd_ppm=float(st["depth_odd"][0] * 1e6),
                              depth_even_ppm=float(st["depth_even"][0] * 1e6),
                              depth_half_ppm=float(st["depth_half"][0] * 1e6),
                              depth_phased_ppm=float(st["depth_phased"][0] * 1e6),
                              depth_phased_err_ppm=float(st["depth_phased"][1] * 1e6),
                              harmonic_delta_loglike=float(st["harmonic_delta_log_likelihood"]),
                              n_transit_times=int(len(st["transit_times"]))))
    d_all, nin = box_depth(t2, y2, P, t0, dur)
    rd = random_phase_depths(t2, y2, P, t0, dur, n=4000)
    err_emp = robust_std(rd)
    E["depth"] = dict(box_depth_ppm=d_all * 1e6, empirical_err_ppm=err_emp * 1e6,
                      empirical_snr=float(d_all / err_emp), n_in_transit_points=nin,
                      random_phase_boxes=int(len(rd)),
                      random_phase_depth_mean_ppm=float(rd.mean() * 1e6),
                      random_phase_std_plain_ppm=float(rd.std() * 1e6),
                      max_random_box_depth_ppm=float(rd.max() * 1e6),
                      n_random_boxes_deeper_than_signal=int((rd >= d_all).sum()),
                      formal_white_err_ppm=float(np.std(y2) / np.sqrt(max(nin, 1)) * 1e6))
    epoch = np.round((t2 - t0) / P).astype(int)
    odd = (epoch % 2) != 0
    do = subset_depth(t2, y2, P, t0, dur, odd)
    de = subset_depth(t2, y2, P, t0, dur, ~odd)
    E["odd_even"] = dict(odd=do, even=de,
                         sigma=float(abs(do["depth_ppm"] - de["depth_ppm"]) / np.hypot(do["err_ppm"], de["err_ppm"])))
    d_sec, nsec = box_depth(t2, y2, P, t0, dur, phase_c=0.5)
    E["secondary"] = dict(depth_ppm=d_sec * 1e6, err_ppm=err_emp * 1e6, sigma=float(d_sec / err_emp), n_in=nsec)
    half = 0.5 * dur / P
    scan = []
    for pc in np.linspace(-0.5, 0.5, 400, endpoint=False):
        if abs(pc) < 3 * half:
            continue
        dd, _ = box_depth(t2, y2, P, t0, dur, phase_c=pc)
        scan.append((pc, dd))
    scan = np.array(scan)
    k = np.nanargmax(scan[:, 1])
    E["secondary"].update(deepest_other_phase=float(scan[k, 0]), deepest_other_depth_ppm=float(scan[k, 1] * 1e6),
                          deepest_other_sigma=float(scan[k, 1] / err_emp))
    years = {"Y1_Q1-4": (q2 >= 1) & (q2 <= 4), "Y2_Q5-8": (q2 >= 5) & (q2 <= 8),
             "Y3_Q9-12": (q2 >= 9) & (q2 <= 12), "Y4_Q13-17": (q2 >= 13)}
    E["per_year"] = {k: subset_depth(t2, y2, P, t0, dur, v) for k, v in years.items()}
    E["per_season"] = {f"season_{s}_Q{s}+4n": subset_depth(t2, y2, P, t0, dur, (q2 % 4) == s) for s in range(4)}
    yd = np.array([v["depth_ppm"] for v in E["per_year"].values()])
    ye = np.array([v["err_ppm"] for v in E["per_year"].values()])
    ok = np.isfinite(yd) & np.isfinite(ye)
    wm = np.sum(yd[ok] / ye[ok]**2) / np.sum(1 / ye[ok]**2)
    E["per_year_chi2"] = float(np.sum((yd[ok] - wm)**2 / ye[ok]**2))
    E["per_year_dof"] = int(ok.sum() - 1)
    rows = per_transit(t2, y2, q2, P, t0, dur)
    E["per_transit"] = rows
    dep = np.array([r["depth_ppm"] for r in rows])
    er = np.array([r["err_ppm"] for r in rows])
    E["n_transits"] = int(len(rows))
    E["n_epochs_in_baseline"] = int(np.floor((t2.max() - t2.min()) / P) + 1)
    if len(rows):
        wmean = np.sum(dep / er**2) / np.sum(1 / er**2)
        imax = dep.argmax()
        E["per_transit_summary"] = dict(
            weighted_mean_ppm=float(wmean), median_ppm=float(np.median(dep)), scatter_ppm=float(np.std(dep)),
            chi2_const=float(np.sum((dep - wmean)**2 / er**2)), n_positive=int((dep > 0).sum()),
            max_single_event_ppm=float(dep.max()), max_single_event_tc_bkjd=float(rows[imax]["tc_bkjd"]),
            single_event_dominance=float(dep.max() / dep.sum()) if dep.sum() > 0 else None,
            depth_without_deepest_ppm=float(np.sum(np.delete(dep, imax) / np.delete(er, imax)**2)
                                            / np.sum(1 / np.delete(er, imax)**2)),
            n_near_gap_lt_0p5d=int(sum(r["gap_dist_d"] < 0.5 for r in rows)))
    E["trapezoid"] = fit_trapezoid(t2, y2, P, t0, dur)
    rp_rs = np.sqrt(max(d_all, 0))
    rho = RSTAR / RSTAR**3   # M/Msun ~ R/Rsun for a G/K dwarf
    E["physical"] = dict(rp_over_rs=float(rp_rs), rp_rearth_if_on_target=float(rp_rs * RSTAR * 109.1),
                         expected_T14_h_central_circular=float(13.0 * (P / 365.25)**(1 / 3) * rho**(-1 / 3)),
                         rstar_used=RSTAR, note="PSF flux is a per-source fit; no dilution correction applied")
    # SAP flux at the same ephemeris
    phs = phase_of(ts_raw, P, t0)
    smask = ~(np.abs(phs) * P < 1.5 * dur)
    ts, ys, qs, dys, sinfo = detrend(ts_raw, fs_raw, qs_raw, win=1.0, transit_mask=smask, robust=robust)
    d_sap, nin_s = box_depth(ts, ys, P, t0, dur)
    rds = random_phase_depths(ts, ys, P, t0, dur, n=4000)
    e_sap = robust_std(rds)
    ratio = d_sap / d_all if d_all != 0 else np.nan
    ratio_err = abs(ratio) * np.hypot(e_sap / d_sap if d_sap != 0 else np.inf, err_emp / d_all)
    ds_sec, _ = box_depth(ts, ys, P, t0, dur, phase_c=0.5)
    E["sap"] = dict(n_used=int(len(ts)), point_scatter_ppm=float(robust_std(ys) * 1e6),
                    box_depth_ppm=d_sap * 1e6, empirical_err_ppm=e_sap * 1e6,
                    empirical_snr=float(d_sap / e_sap), n_in=nin_s,
                    ratio_sap_over_psf=float(ratio), ratio_err=float(ratio_err),
                    crowdsap_catalog=CROWDSAP, psffrac_catalog=PSFFRAC,
                    expected_ratio_if_on_target_undiluted_psf=CROWDSAP,
                    secondary_ppm=ds_sec * 1e6, secondary_sigma=float(ds_sec / e_sap),
                    n_steps_split=len(sinfo["steps"]), n_clipped=sinfo["n_clipped"])
    E["detrend_info"] = dict(steps=info["steps"], n_clipped=info["n_clipped"], mad_ppm=info["mad_ppm"],
                             n_used=int(len(t2)))
    np.savez(f"{SCR}/fold_cache_{label}.npz", t=t2, y=y2, q=q2, ts=ts, ys=ys, P=P, t0=t0, dur=dur,
             rp=rp2, rpow=rpow2, rd=rd, scan=scan)
    return E


# ----------------------------------------------------------------------------
# main analysis
# ----------------------------------------------------------------------------
def analyse(robust, label):
    R = dict(label=label, robust_detrend=robust)
    t_raw, f_raw, q_raw = load("FLUX")
    ts_raw, fs_raw, qs_raw = load("SAP_FLUX")
    t, y, q, dy, info = detrend(t_raw, f_raw, q_raw, win=1.0, robust=robust)
    T = t.max() - t.min()
    R["data"] = dict(n_good_quality0=int(len(t_raw)), n_used=int(len(t)), baseline_d=float(T),
                     t_start_bkjd=float(t.min()), t_end_bkjd=float(t.max()),
                     quarters=[int(x) for x in np.unique(q)],
                     point_scatter_ppm=float(robust_std(y) * 1e6),
                     steps_split=info["steps"], n_clipped=info["n_clipped"],
                     detrend=("per-quarter median normalisation; 1.0-d time-based running median "
                              "(segments split at gaps > 0.5 d" + (" and at detected steps; +-0.05 d guard" if robust else "")
                              + "); one-sided 4-sigma clip of positive outliers"
                              + ("; symmetric 5-sigma robust clip" if robust else " only (no negative clipping)")))
    print(f"[{label}] stage 1:", R["data"])
    R["bls_blind"] = blind_bls(t, y, dy, f"{label}_psf")
    R["bls_blind_inverted"] = blind_bls(t, 2.0 - y, dy, f"{label}_psf_inverted")
    best = R["bls_blind"]["top_peaks"][0]
    R["recovered_blind"] = bool(abs(best["period_d"] / P_PIPE - 1) < 0.01)
    R["top_peak_alias_of_pipeline"] = best["alias_of_pipeline"]

    # statistics at the blind top peak
    ref, _, _ = refine(t, y, dy, best["period_d"])
    print(f"[{label}] refined top peak:", ref["P"], ref["dur"] * 24, ref["t0"])
    R["at_top_peak"] = ephemeris_stats(t_raw, f_raw, q_raw, ref["P"], ref["t0"], ref["dur"], robust,
                                       f"{label}_top", ts_raw, fs_raw, qs_raw)
    print(f"[{label}] top-peak depth:", R["at_top_peak"]["depth"])
    # statistics forced at the pipeline period (tight window around 12.9715 d)
    refp, _, _ = refine(t, y, dy, P_PIPE, frac=0.002, nper=2001)
    print(f"[{label}] refined near P_pipe:", refp["P"], refp["dur"] * 24, refp["t0"])
    R["at_pipeline_period"] = ephemeris_stats(t_raw, f_raw, q_raw, refp["P"], refp["t0"], refp["dur"], robust,
                                              f"{label}_pipe", ts_raw, fs_raw, qs_raw)
    print(f"[{label}] P_pipe depth:", R["at_pipeline_period"]["depth"])
    print(f"[{label}] P_pipe per-transit:", R["at_pipeline_period"].get("per_transit_summary"))
    return R


# ----------------------------------------------------------------------------
# 5. figure
# ----------------------------------------------------------------------------
def make_figure(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prim = R["pass2_robust"]
    E = prim["at_top_peak"]
    c = np.load(f"{SCR}/fold_cache_pass2_robust_top.npz")
    b = np.load(f"{SCR}/bls_pass2_robust_psf.npz")
    bi = np.load(f"{SCR}/bls_pass2_robust_psf_inverted.npz")
    b1 = np.load(f"{SCR}/bls_pass1_raw_psf.npz")
    t, y, q, ts, ys = c["t"], c["y"], c["q"], c["ts"], c["ys"]
    P, t0, dur = float(c["P"]), float(c["t0"]), float(c["dur"])

    # reference palette (dataviz skill): neutral ink, categorical slots 1 (blue) and 2 (orange)
    INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#d6d5d0"
    C_PSF, C_SAP, C_ODD, C_EVEN, C_NULL = "#2a78d6", "#eb6834", "#2a78d6", "#eb6834", "#8a8984"
    SURF = "#fcfcfb"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                         "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": SURF, "axes.facecolor": SURF})

    def binned(ph, yy, nb, lo, hi):
        edges = np.linspace(lo, hi, nb + 1)
        idx = np.digitize(ph, edges) - 1
        xb, yb, eb = [], [], []
        for i in range(nb):
            s = idx == i
            if s.sum() >= 3:
                xb.append(0.5 * (edges[i] + edges[i + 1]))
                yb.append(np.mean(yy[s]))
                eb.append(np.std(yy[s]) / np.sqrt(s.sum()))
        return np.array(xb), np.array(yb), np.array(eb)

    fig = plt.figure(figsize=(11, 11))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.1, 1.0, 1.0], hspace=0.55, wspace=0.28)

    # (a) blind periodograms
    ax = fig.add_subplot(gs[0, :])
    ax.plot(b1["periods"], b1["power"] / b1["power"].max() * b["power"].max(), color=C_SAP, lw=0.5, alpha=0.55,
            label="pass 1 (raw detrend; scaled to pass-2 maximum)")
    ax.plot(bi["periods"], bi["power"], color=C_NULL, lw=0.6, label="pass 2, inverted flux (2 - f) null")
    ax.plot(b["periods"], b["power"], color=C_PSF, lw=0.7, label="pass 2, PSF flux")
    ax.set_xscale("log")
    ax.set_xlabel("period (d)")
    ax.set_ylabel("BLS power (delta log L)")
    for k, p in enumerate(prim["bls_blind"]["top_peaks"]):
        ax.annotate(f"#{k+1}  {p['period_d']:.4f} d", (p["period_d"], p["power"]),
                    xytext=(4, 4 - 11 * k), textcoords="offset points", fontsize=8, color=INK2)
    ax.axvline(P_PIPE, color=INK2, lw=0.6, ls=":", label=f"pipeline P = {P_PIPE} d")
    ax.set_title(f"{FNAME}: blind BLS 1-50 d, own grid ({prim['bls_blind']['n_periods']} periods); pass-2 SDE real = "
                 f"{prim['bls_blind']['sde']:.1f}, inverted = {prim['bls_blind_inverted']['sde']:.1f}; "
                 f"pass-1 top peak {R['pass1_raw']['bls_blind']['top_peaks'][0]['period_d']:.3f} d",
                 loc="left", fontsize=9.5, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")

    # (b) fold PSF
    ph = phase_of(t, P, t0) * P * 24
    ax = fig.add_subplot(gs[1, 0])
    w = np.abs(ph) < 12
    ax.plot(ph[w], y[w], ".", ms=1.5, color=GRID, alpha=0.8, rasterized=True)
    xb, yb, eb = binned(ph[w], y[w], 48, -12, 12)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3, color=C_PSF, ecolor=C_PSF, elinewidth=0.8, capsize=0,
                label="PSF flux, 30-min bins")
    tr = E["trapezoid"]
    xx = np.linspace(-12, 12, 600)
    ax.plot(xx, trapezoid([tr["depth_ppm"] * 1e-6, tr["T14_h"] / 24, tr["T23_h"] / 24], xx / 24),
            color=INK, lw=1.2, label=f"trapezoid T14={tr['T14_h']:.2f} h, T23/T14={tr['shape_T23_over_T14']:.2f}")
    ax.axvspan(-dur * 12, dur * 12, color=C_PSF, alpha=0.06, lw=0)
    ax.set_ylim(0.988, 1.008)
    ax.set_xlabel("hours from mid-transit")
    ax.set_ylabel("relative flux")
    d = E["depth"]
    ax.set_title(f"PSF fold at pass-2 top peak: P={P:.5f} d, box depth {d['box_depth_ppm']:.0f} +/- {d['empirical_err_ppm']:.0f} ppm "
                 f"(emp. SNR {d['empirical_snr']:.1f}), {E['n_transits']} transits", loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")

    # (c) SAP fold
    phs = phase_of(ts, P, t0) * P * 24
    ax = fig.add_subplot(gs[1, 1])
    ws = np.abs(phs) < 12
    ax.plot(phs[ws], ys[ws], ".", ms=1.5, color=GRID, alpha=0.8, rasterized=True)
    xb, yb, eb = binned(phs[ws], ys[ws], 48, -12, 12)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3, color=C_SAP, ecolor=C_SAP, elinewidth=0.8, capsize=0,
                label="SAP flux, 30-min bins")
    xb2, yb2, _ = binned(ph[w], y[w], 48, -12, 12)
    ax.plot(xb2, yb2, "-", color=C_PSF, lw=1.0, alpha=0.9, label="PSF binned (for comparison)")
    ax.set_ylim(0.988, 1.008)
    ax.set_xlabel("hours from mid-transit")
    ax.set_ylabel("relative flux")
    s = E["sap"]
    ax.set_title(f"SAP fold: depth {s['box_depth_ppm']:.0f} +/- {s['empirical_err_ppm']:.0f} ppm; "
                 f"SAP/PSF = {s['ratio_sap_over_psf']:.2f} +/- {s['ratio_err']:.2f} (CROWDSAP {CROWDSAP:.2f})",
                 loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")

    # (d) odd / even
    ax = fig.add_subplot(gs[2, 0])
    epoch = np.round((t - t0) / P).astype(int)
    for sel, col, lab, off in (((epoch % 2) != 0, C_ODD, "odd", -0.08), ((epoch % 2) == 0, C_EVEN, "even", 0.08)):
        ww = sel & (np.abs(ph) < 12)
        xb, yb, eb = binned(ph[ww], y[ww], 24, -12, 12)
        ax.errorbar(xb + off, yb, eb, fmt="o", ms=3, color=col, ecolor=col, elinewidth=0.8, capsize=0, label=lab)
    oe = E["odd_even"]
    ax.set_ylim(0.988, 1.008)
    ax.set_xlabel("hours from mid-transit")
    ax.set_ylabel("relative flux")
    ax.set_title(f"odd {oe['odd']['depth_ppm']:.0f}+/-{oe['odd']['err_ppm']:.0f} ppm vs even "
                 f"{oe['even']['depth_ppm']:.0f}+/-{oe['even']['err_ppm']:.0f} ppm  ({oe['sigma']:.1f} sigma)",
                 loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    # (e) secondary
    ax = fig.add_subplot(gs[2, 1])
    ph5 = (((phase_of(t, P, t0) - 0.5 + 0.5) % 1.0) - 0.5) * P * 24
    w5 = np.abs(ph5) < 12
    ax.plot(ph5[w5], y[w5], ".", ms=1.5, color=GRID, alpha=0.8, rasterized=True)
    xb, yb, eb = binned(ph5[w5], y[w5], 48, -12, 12)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3, color=C_PSF, ecolor=C_PSF, elinewidth=0.8, capsize=0)
    ax.axvspan(-dur * 12, dur * 12, color=INK2, alpha=0.08, lw=0)
    ax.set_ylim(0.988, 1.008)
    ax.set_xlabel("hours from phase 0.5")
    ax.set_ylabel("relative flux")
    sec = E["secondary"]
    ax.set_title(f"secondary at phase 0.5: {sec['depth_ppm']:.0f} +/- {sec['err_ppm']:.0f} ppm ({sec['sigma']:.1f} sigma); "
                 f"deepest other box {sec['deepest_other_depth_ppm']:.0f} ppm at phase {sec['deepest_other_phase']:.2f}",
                 loc="left", fontsize=8, color=INK)

    # (f) per-transit depths
    ax = fig.add_subplot(gs[3, 0])
    rows = E["per_transit"]
    tc = np.array([r["tc_bkjd"] for r in rows])
    dp = np.array([r["depth_ppm"] for r in rows])
    ep = np.array([r["err_ppm"] for r in rows])
    ax.errorbar(tc, dp, ep, fmt="o", ms=3, color=C_PSF, ecolor=C_PSF, elinewidth=0.6, capsize=0, label="per-transit")
    ax.axhline(E["depth"]["box_depth_ppm"], color=INK, lw=1.0, label="global box depth")
    ax.axhline(0, color=GRID, lw=0.8)
    for k, (name, v) in enumerate(E["per_year"].items()):
        sel = {"Y1_Q1-4": (q >= 1) & (q <= 4), "Y2_Q5-8": (q >= 5) & (q <= 8),
               "Y3_Q9-12": (q >= 9) & (q <= 12), "Y4_Q13-17": q >= 13}[name]
        if sel.sum() and v["depth_ppm"] is not None and np.isfinite(v["depth_ppm"]):
            x0, x1 = t[sel].min(), t[sel].max()
            ax.fill_between([x0, x1], v["depth_ppm"] - v["err_ppm"], v["depth_ppm"] + v["err_ppm"],
                            color=C_SAP, alpha=0.18, lw=0)
            ax.plot([x0, x1], [v["depth_ppm"]] * 2, color=C_SAP, lw=1.4,
                    label="per-year box depth +/- emp. err" if k == 0 else None)
    ax.set_xlabel("time (BKJD)")
    ax.set_ylabel("depth (ppm)")
    ps = E.get("per_transit_summary", {})
    ax.set_title(f"per-transit local depths: {ps.get('n_positive')}/{E['n_transits']} positive, deepest/sum = "
                 f"{(ps.get('single_event_dominance') or float('nan')):.2f}, chi2(const) = {ps.get('chi2_const', float('nan')):.0f}/{E['n_transits']-1}",
                 loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")

    # (g) random-phase distribution
    ax = fig.add_subplot(gs[3, 1])
    rd = c["rd"] * 1e6
    ax.hist(rd, bins=60, color=GRID, edgecolor=SURF, lw=0.3, label="random-phase boxes (same P, width)")
    ax.axvline(E["depth"]["box_depth_ppm"], color=C_PSF, lw=1.5, label="signal box depth")
    ax.axvline(E["secondary"]["depth_ppm"], color=C_SAP, lw=1.2, ls="--", label="phase-0.5 box")
    ax.set_xlabel("box depth (ppm)")
    ax.set_ylabel("count")
    ax.set_title(f"empirical error = {E['depth']['empirical_err_ppm']:.0f} ppm; "
                 f"{E['depth']['n_random_boxes_deeper_than_signal']}/{E['depth']['random_phase_boxes']} random boxes "
                 f"deeper than signal", loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")

    fig.suptitle(f"Independent re-derivation, {FNAME} (KBONUS-BKG PSF flux, quality==0, pass-2 robust detrend); "
                 f"candidate-language only", x=0.01, ha="left", fontsize=10.5, color=INK, y=0.995)
    out = f"{OUT}/fold_rederive.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print("figure:", out)
    return out


def to_jsonable(o):
    if isinstance(o, dict):
        return {k: to_jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return to_jsonable(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def strip_stats(R):
    for k in ("at_top_peak", "at_pipeline_period"):
        if k in R and "stats" in R[k].get("ephemeris", {}):
            R[k]["ephemeris"].pop("stats")
    return R


if __name__ == "__main__":
    jpath = f"{OUT}/rederive_results.json"
    if "--plot-only" in sys.argv:
        R = json.load(open(jpath))
    else:
        R = dict(fname=FNAME, fits=FITS, date="2026-10-17",
                 method="independent numpy/scipy/astropy re-derivation; no pipeline/ code imported or read",
                 pipeline_reference=dict(period_d=P_PIPE, depth_ppm=3034, duration_h=2.2),
                 primary_pass="pass2_robust")
        R["pass1_raw"] = strip_stats(analyse(robust=False, label="pass1_raw"))
        with open(jpath, "w") as fh:
            json.dump(to_jsonable(R), fh, indent=2)
        R["pass2_robust"] = strip_stats(analyse(robust=True, label="pass2_robust"))
        with open(jpath, "w") as fh:
            json.dump(to_jsonable(R), fh, indent=2)
        print("wrote", jpath)
    if "--no-plot" not in sys.argv:
        R["figure"] = make_figure(R)
        with open(jpath, "w") as fh:
            json.dump(to_jsonable(R), fh, indent=2)
