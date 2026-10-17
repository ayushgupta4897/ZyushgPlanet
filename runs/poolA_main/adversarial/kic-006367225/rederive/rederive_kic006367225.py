#!/usr/bin/env python
"""
Independent re-derivation of the transit-like signal in KBONUS-BKG kic-006367225.

Written from scratch with numpy / scipy / astropy only.  Nothing under
ZyushgPlanet/pipeline/ is imported or read.  No wotan, no TLS.

Stages
  1. load LIGHTCURVE_STITCHED, SAP_QUALITY == 0, per-quarter normalisation,
     1-d time-based running-median detrend (own code), one-sided 4-sigma clip of
     positive outliers only (transits are never clipped).
  2. blind BLS 1-50 d (astropy BoxLeastSquares, own log-spaced grid), plus the
     same search on the inverted flux 2-f as a null.
  3. refine at the best peak; box depth with random-phase empirical error,
     odd/even, secondary at phase 0.5, per-year and per-season depths,
     per-transit depths, transit count, duration, trapezoid shape.
  4. repeat the depth measurement on SAP_FLUX (dilution check).
  5. write results JSON, fold figure.

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


def running_median(t, y, win=1.0, gap=0.5, mask=None):
    """Time-based running median of width `win` days, evaluated at every t.
    Segments are split at gaps > `gap` d.  If `mask` (bool, True = use) is
    given, only masked-in points enter the median (used to exclude transits)."""
    out = np.full_like(y, np.nan)
    brk = np.where(np.diff(t) > gap)[0] + 1
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


def detrend(t, f, q, win=1.0, transit_mask=None, clip=4.0):
    """Returns t, y (flat, ~1), q, dy (per-quarter robust scatter)."""
    fn = normalise_quarters(f, q)
    tr = running_median(t, fn, win=win, mask=transit_mask)
    ok = np.isfinite(tr)
    t, fn, q, tr = t[ok], fn[ok], q[ok], tr[ok]
    y = fn / tr
    # one-sided clip: only remove positive outliers (never clip dips)
    mad = 1.4826 * np.median(np.abs(y - np.median(y)))
    keep = y < 1 + clip * mad
    t, y, q = t[keep], y[keep], q[keep]
    dy = np.empty_like(y)
    for qq in np.unique(q):
        s = q == qq
        dy[s] = 1.4826 * np.median(np.abs(y[s] - np.median(y[s])))
    return t, y, q, dy


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
    print(f"[{label}] BLS {len(periods)} periods x {len(DURS)} durations in {dt:.1f}s; "
          f"SDE={sde_of(pw):.2f}")
    for p in peaks:
        print("   ", p)
    np.savez(f"{SCR}/bls_{label}.npz", periods=periods, power=pw,
             depth=r.depth, duration=r.duration, t0=r.transit_time)
    return dict(n_periods=int(len(periods)), durations_d=DURS.tolist(),
                pmin_d=1.0, pmax_d=50.0, sde=sde_of(pw), top_peaks=peaks,
                runtime_s=dt)


# ----------------------------------------------------------------------------
# 3. ephemeris refinement and depth statistics
# ----------------------------------------------------------------------------
def refine(t, y, dy, p0, frac=0.003, nper=3001):
    periods = np.linspace(p0 * (1 - frac), p0 * (1 + frac), nper)
    durs = np.arange(0.03, 0.31, 0.01)
    bls = BoxLeastSquares(t, y, dy)
    r = bls.power(periods, durs, objective="likelihood")
    i = int(np.argmax(r.power))
    P, dur, tt = float(periods[i]), float(r.duration[i]), float(r.transit_time[i])
    # move reference epoch to the transit nearest the data mid-point
    tmid = 0.5 * (t.min() + t.max())
    tt = tt + np.round((tmid - tt) / P) * P
    st = bls.compute_stats(P, dur, tt)
    # half-width-at-half-depth of the power peak as a period uncertainty proxy
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
    and secondary windows.  Returns the array of depths; its robust std is the
    empirical depth error (includes red noise at this timescale)."""
    half = 0.5 * dur / P
    out = []
    while len(out) < n:
        pc = RNG.uniform(-0.5, 0.5)
        if any(abs(((pc - a + 0.5) % 1.0) - 0.5) < avoid_w * half * 2 for a in avoid):
            continue
        d, nin = box_depth(t, y, P, t0, dur, phase_c=pc)
        if np.isfinite(d):
            out.append(d)
    out = np.array(out)
    return out


def robust_std(a):
    a = a[np.isfinite(a)]
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def per_transit(t, y, q, P, t0, dur, oot_win=1.0, oot_excl=1.5):
    """Per-epoch local depth: local OOT = points within +-oot_win d of the
    transit centre but outside +-oot_excl*dur.  Requires >= 50% of the expected
    in-transit cadences to be present."""
    epochs = np.arange(np.floor((t.min() - t0) / P), np.ceil((t.max() - t0) / P) + 1)
    exp_n = dur / (29.4244 / 1440.0)
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
        # gap proximity: distance from transit centre to the nearest data gap > 0.3 d
        tt = t[(np.abs(t - tc) < 3.0)]
        gaps = np.where(np.diff(tt) > 0.3)[0]
        gap_dist = float(min([min(abs(tt[g] - tc), abs(tt[g + 1] - tc)) for g in gaps] + [3.0]))
        rows.append(dict(epoch=int(n), tc_bkjd=float(tc), tc_bjd=float(tc + BKJD0),
                         quarter=int(np.median(q[din])), n_in=nin, depth_ppm=d * 1e6,
                         err_ppm=e * 1e6, gap_dist_d=gap_dist))
    return rows


def trapezoid(params, x):
    depth, t14, t23 = params
    t14 = abs(t14); t23 = min(abs(t23), t14)
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
    t14 = abs(t14); t23 = min(abs(t23), t14)
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


# ----------------------------------------------------------------------------
# main analysis
# ----------------------------------------------------------------------------
def analyse():
    R = dict(fname=FNAME, fits=FITS, date="2026-10-17",
             method="independent numpy/astropy re-derivation; no pipeline/ code")

    # ---- stage 1: PSF flux
    t_raw, f_raw, q_raw = load("FLUX")
    t, y, q, dy = detrend(t_raw, f_raw, q_raw, win=1.0)
    T = t.max() - t.min()
    R["data"] = dict(n_good=int(len(t_raw)), n_used=int(len(t)), baseline_d=float(T),
                     t_start_bkjd=float(t.min()), t_end_bkjd=float(t.max()),
                     quarters=[int(x) for x in np.unique(q)],
                     point_scatter_ppm=float(robust_std(y) * 1e6),
                     detrend="per-quarter median normalisation; 1.0-d time-based running median "
                             "(segments split at gaps > 0.5 d); one-sided 4-sigma clip of positive "
                             "outliers only")
    print("stage 1:", R["data"])

    # ---- stage 2: blind BLS real + inverted
    R["bls_blind"] = blind_bls(t, y, dy, "psf")
    R["bls_blind_inverted"] = blind_bls(t, 2.0 - y, dy, "psf_inverted")
    best = R["bls_blind"]["top_peaks"][0]
    R["recovered_blind"] = bool(abs(best["period_d"] / P_PIPE - 1) < 0.01)
    R["top_peak_alias_of_pipeline"] = best["alias_of_pipeline"]

    # ---- stage 3: refine at the best peak (use the pipeline-like peak if an alias
    #      of it is the top peak; otherwise refine the true top peak and ALSO P_pipe)
    p_start = best["period_d"]
    ref, rp, rpow = refine(t, y, dy, p_start)
    P, dur, t0 = ref["P"], ref["dur"], ref["t0"]
    print("refined:", P, dur * 24, t0)

    # masked re-detrend so the running median is not biased by the transit
    ph = phase_of(t_raw, P, t0)
    tmask = ~(np.abs(ph) * P < 1.5 * dur)        # True = out of transit
    t2, y2, q2, dy2 = detrend(t_raw, f_raw, q_raw, win=1.0, transit_mask=tmask)
    # re-refine on the masked-detrended data (small shift in dur/depth expected)
    ref2, rp2, rpow2 = refine(t2, y2, dy2, P, frac=0.001, nper=2001)
    P, dur, t0 = ref2["P"], ref2["dur"], ref2["t0"]
    R["ephemeris"] = dict(period_d=P, period_err_hwhm_d=ref2["P_err_hwhm"],
                          t0_bkjd=t0, t0_bjd=t0 + BKJD0, duration_h=dur * 24,
                          bls_depth_ppm=ref2["depth_bls"] * 1e6,
                          bls_depth_err_ppm=ref2["depth_bls_err"] * 1e6,
                          bls_depth_snr=float(ref2["depth_bls"] / ref2["depth_bls_err"]),
                          bls_stats_depth_odd_ppm=float(ref2["stats"]["depth_odd"][0] * 1e6),
                          bls_stats_depth_even_ppm=float(ref2["stats"]["depth_even"][0] * 1e6),
                          bls_stats_depth_half_ppm=float(ref2["stats"]["depth_half"][0] * 1e6),
                          bls_stats_depth_phased_ppm=float(ref2["stats"]["depth_phased"][0] * 1e6),
                          bls_stats_depth_phased_err_ppm=float(ref2["stats"]["depth_phased"][1] * 1e6),
                          bls_stats_harmonic_delta_loglike=float(ref2["stats"]["harmonic_delta_log_likelihood"]),
                          n_transits_bls_stats=int(len(ref2["stats"]["transit_times"])))
    print("ephemeris:", R["ephemeris"])

    # box depth + empirical error (random-phase boxes)
    d_all, nin = box_depth(t2, y2, P, t0, dur)
    rd = random_phase_depths(t2, y2, P, t0, dur, n=4000)
    err_emp = robust_std(rd)
    R["depth"] = dict(box_depth_ppm=d_all * 1e6, empirical_err_ppm=err_emp * 1e6,
                      empirical_snr=float(d_all / err_emp), n_in_transit_points=nin,
                      random_phase_boxes=int(len(rd)),
                      random_phase_depth_mean_ppm=float(rd.mean() * 1e6),
                      random_phase_std_plain_ppm=float(rd.std() * 1e6),
                      max_random_box_depth_ppm=float(rd.max() * 1e6),
                      n_random_boxes_deeper_than_signal=int((rd >= d_all).sum()),
                      formal_white_err_ppm=float(np.std(y2) / np.sqrt(nin) * 1e6))
    print("depth:", R["depth"])

    # odd / even
    ph2 = phase_of(t2, P, t0)
    epoch = np.round((t2 - t0) / P).astype(int)
    odd = (epoch % 2) != 0
    even = ~odd
    do = subset_depth(t2, y2, P, t0, dur, odd)
    de = subset_depth(t2, y2, P, t0, dur, even)
    oe_sig = abs(do["depth_ppm"] - de["depth_ppm"]) / np.hypot(do["err_ppm"], de["err_ppm"])
    R["odd_even"] = dict(odd=do, even=de, sigma=float(oe_sig))
    print("odd/even:", R["odd_even"])

    # secondary at phase 0.5 and deepest box anywhere else
    d_sec, nsec = box_depth(t2, y2, P, t0, dur, phase_c=0.5)
    R["secondary"] = dict(depth_ppm=d_sec * 1e6, err_ppm=err_emp * 1e6,
                          sigma=float(d_sec / err_emp), n_in=nsec)
    # scan phases for the deepest non-primary box
    pcs = np.linspace(-0.5, 0.5, 400, endpoint=False)
    half = 0.5 * dur / P
    scan = []
    for pc in pcs:
        if abs(pc) < 3 * half:
            continue
        dd, _ = box_depth(t2, y2, P, t0, dur, phase_c=pc)
        scan.append((pc, dd))
    scan = np.array(scan)
    k = np.nanargmax(scan[:, 1])
    R["secondary"]["deepest_other_phase"] = float(scan[k, 0])
    R["secondary"]["deepest_other_depth_ppm"] = float(scan[k, 1] * 1e6)
    R["secondary"]["deepest_other_sigma"] = float(scan[k, 1] / err_emp)
    print("secondary:", R["secondary"])

    # per-year (quarter groups) and per-season (Q mod 4 = same CCD channel)
    years = {"Y1_Q1-4": (q2 >= 1) & (q2 <= 4), "Y2_Q5-8": (q2 >= 5) & (q2 <= 8),
             "Y3_Q9-12": (q2 >= 9) & (q2 <= 12), "Y4_Q13-17": (q2 >= 13)}
    R["per_year"] = {k: subset_depth(t2, y2, P, t0, dur, v) for k, v in years.items()}
    R["per_season"] = {f"season_{s}_Q{s}+4n": subset_depth(t2, y2, P, t0, dur, (q2 % 4) == s)
                       for s in range(4)}
    yd = np.array([v["depth_ppm"] for v in R["per_year"].values()])
    ye = np.array([v["err_ppm"] for v in R["per_year"].values()])
    wm = np.nansum(yd / ye**2) / np.nansum(1 / ye**2)
    R["per_year_chi2"] = float(np.nansum((yd - wm)**2 / ye**2))
    R["per_year_dof"] = int(np.isfinite(yd).sum() - 1)
    print("per year:", R["per_year"], R["per_year_chi2"])
    print("per season:", R["per_season"])

    # per-transit depths
    rows = per_transit(t2, y2, q2, P, t0, dur)
    R["per_transit"] = rows
    dep = np.array([r["depth_ppm"] for r in rows]); er = np.array([r["err_ppm"] for r in rows])
    R["n_transits"] = int(len(rows))
    R["n_epochs_in_baseline"] = int(np.floor((t2.max() - t2.min()) / P) + 1)
    R["per_transit_summary"] = dict(
        weighted_mean_ppm=float(np.sum(dep / er**2) / np.sum(1 / er**2)),
        median_ppm=float(np.median(dep)), scatter_ppm=float(np.std(dep)),
        chi2_const=float(np.sum((dep - np.sum(dep / er**2) / np.sum(1 / er**2))**2 / er**2)),
        n_positive=int((dep > 0).sum()),
        max_single_event_ppm=float(dep.max()),
        single_event_dominance=float(dep.max() / dep.sum()) if dep.sum() > 0 else np.nan,
        depth_without_deepest_ppm=float(np.sum(np.delete(dep, dep.argmax()) / np.delete(er, dep.argmax())**2)
                                        / np.sum(1 / np.delete(er, dep.argmax())**2)),
        n_near_gap_lt_0p5d=int(sum(r["gap_dist_d"] < 0.5 for r in rows)))
    print("per transit:", R["per_transit_summary"])

    # trapezoid shape
    R["trapezoid"] = fit_trapezoid(t2, y2, P, t0, dur)
    print("trapezoid:", R["trapezoid"])

    # physical sanity
    rp_rs = np.sqrt(max(d_all, 0))
    R["physical"] = dict(rp_over_rs=float(rp_rs),
                         rp_rearth_if_on_target=float(rp_rs * RSTAR * 109.1),
                         rp_rearth_corrected_for_psffrac_note="PSF flux is a per-source fit; no dilution correction applied",
                         expected_T14_h_central_circular=float(
                             24 * 13.0 / 24 * (P / 365.25)**(1 / 3) * RSTAR / RSTAR**(1 / 3) * 1.0),
                         rstar_used=RSTAR)
    # better estimate of expected duration: T ~ 13 h * (P/yr)^(1/3) * (rho_sun/rho)^(1/3)
    rho = 1.0 / RSTAR**3 * 1.0  # M~R assumed for a G/K dwarf (M/Msun ~ R/Rsun)
    rho = (RSTAR / RSTAR**3)
    R["physical"]["expected_T14_h_central_circular"] = float(13.0 * (P / 365.25)**(1 / 3) * rho**(-1 / 3))

    # ---- stage 4: SAP flux
    ts_raw, fs_raw, qs_raw = load("SAP_FLUX")
    phs = phase_of(ts_raw, P, t0)
    smask = ~(np.abs(phs) * P < 1.5 * dur)
    ts, ys, qs, dys = detrend(ts_raw, fs_raw, qs_raw, win=1.0, transit_mask=smask)
    d_sap, nin_s = box_depth(ts, ys, P, t0, dur)
    rds = random_phase_depths(ts, ys, P, t0, dur, n=4000)
    e_sap = robust_std(rds)
    ratio = d_sap / d_all
    ratio_err = ratio * np.hypot(e_sap / d_sap, err_emp / d_all)
    R["sap"] = dict(n_used=int(len(ts)), point_scatter_ppm=float(robust_std(ys) * 1e6),
                    box_depth_ppm=d_sap * 1e6, empirical_err_ppm=e_sap * 1e6,
                    empirical_snr=float(d_sap / e_sap), n_in=nin_s,
                    ratio_sap_over_psf=float(ratio), ratio_err=float(ratio_err),
                    crowdsap_catalog=CROWDSAP, psffrac_catalog=PSFFRAC,
                    expected_ratio_if_on_target_undiluted_psf=CROWDSAP,
                    note=("If the dip is on the target and the PSF flux is fully de-blended, "
                          "SAP depth should be ~CROWDSAP x PSF depth (SAP aperture contains "
                          "(1-CROWDSAP) foreign flux). A ratio >> 1 would mean the dip is "
                          "stronger in the aperture sum than in the target's PSF fit, i.e. a "
                          "neighbour; a ratio ~ CROWDSAP..1 is consistent with the target."))
    # also SAP odd/even + secondary for completeness
    ds_sec, _ = box_depth(ts, ys, P, t0, dur, phase_c=0.5)
    R["sap"]["secondary_ppm"] = ds_sec * 1e6
    R["sap"]["secondary_sigma"] = float(ds_sec / e_sap)
    print("sap:", R["sap"])

    # ---- cache for plotting
    np.savez(f"{SCR}/fold_cache.npz", t=t2, y=y2, q=q2, ts=ts, ys=ys, P=P, t0=t0, dur=dur,
             rp=rp2, rpow=rpow2, rd=rd, scan=scan)
    return R


# ----------------------------------------------------------------------------
# 5. figure
# ----------------------------------------------------------------------------
def make_figure(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = np.load(f"{SCR}/fold_cache.npz")
    b = np.load(f"{SCR}/bls_psf.npz")
    bi = np.load(f"{SCR}/bls_psf_inverted.npz")
    t, y, q, ts, ys = c["t"], c["y"], c["q"], c["ts"], c["ys"]
    P, t0, dur = float(c["P"]), float(c["t0"]), float(c["dur"])

    # palette: neutral ink + a small fixed categorical set (validated, see notes)
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

    fig = plt.figure(figsize=(11, 10.5))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.1, 1.0, 1.0], hspace=0.5, wspace=0.28)

    # (a) blind periodogram, real vs inverted
    ax = fig.add_subplot(gs[0, :])
    ax.plot(bi["periods"], bi["power"], color=C_NULL, lw=0.6, label="inverted flux (2 - f) null")
    ax.plot(b["periods"], b["power"], color=C_PSF, lw=0.7, label="PSF flux")
    ax.set_xscale("log")
    ax.set_xlabel("period (d)"); ax.set_ylabel("BLS power (delta log L)")
    tp = R["bls_blind"]["top_peaks"]
    for k, p in enumerate(tp):
        ax.annotate(f"#{k+1}  {p['period_d']:.4f} d", (p["period_d"], p["power"]),
                    xytext=(4, 4 - 11 * k), textcoords="offset points", fontsize=8, color=INK2)
    ax.axvline(P_PIPE, color=INK2, lw=0.6, ls=":", label=f"pipeline P = {P_PIPE} d")
    ax.set_title(f"{FNAME}: blind BLS 1-50 d, own grid ({R['bls_blind']['n_periods']} periods); "
                 f"SDE real = {R['bls_blind']['sde']:.1f}, inverted = {R['bls_blind_inverted']['sde']:.1f}",
                 loc="left", fontsize=9.5, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    # (b) full fold PSF with binned points, trapezoid overlay
    ph = phase_of(t, P, t0) * P * 24   # hours
    ax = fig.add_subplot(gs[1, 0])
    w = np.abs(ph) < 12
    ax.plot(ph[w], y[w], ".", ms=1.5, color=GRID, alpha=0.8, rasterized=True)
    xb, yb, eb = binned(ph[w], y[w], 48, -12, 12)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3, color=C_PSF, ecolor=C_PSF, elinewidth=0.8, capsize=0,
                label="PSF flux, 30-min bins")
    tr = R["trapezoid"]
    xx = np.linspace(-12, 12, 600)
    ax.plot(xx, trapezoid([tr["depth_ppm"] * 1e-6, tr["T14_h"] / 24, tr["T23_h"] / 24], xx / 24),
            color=INK, lw=1.2, label=f"trapezoid T14={tr['T14_h']:.2f} h, T23/T14={tr['shape_T23_over_T14']:.2f}")
    ax.axvspan(-dur * 12, dur * 12, color=C_PSF, alpha=0.06, lw=0)
    ax.set_ylim(1 - 6 * R["depth"]["box_depth_ppm"] * 1e-6 * 0 - 0.012, 1.008)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("relative flux")
    d = R["depth"]
    ax.set_title(f"PSF fold: P={P:.5f} d, box depth {d['box_depth_ppm']:.0f} +/- {d['empirical_err_ppm']:.0f} ppm "
                 f"(emp. SNR {d['empirical_snr']:.1f}), {R['n_transits']} transits", loc="left", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")

    # (c) SAP fold on the same axes scale
    phs = phase_of(ts, P, t0) * P * 24
    ax = fig.add_subplot(gs[1, 1])
    ws = np.abs(phs) < 12
    ax.plot(phs[ws], ys[ws], ".", ms=1.5, color=GRID, alpha=0.8, rasterized=True)
    xb, yb, eb = binned(phs[ws], ys[ws], 48, -12, 12)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3, color=C_SAP, ecolor=C_SAP, elinewidth=0.8, capsize=0,
                label="SAP flux, 30-min bins")
    xb2, yb2, _ = binned(ph[w], y[w], 48, -12, 12)
    ax.plot(xb2, yb2, "-", color=C_PSF, lw=1.0, alpha=0.9, label="PSF binned (for comparison)")
    ax.set_ylim(1 - 0.012, 1.008)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("relative flux")
    s = R["sap"]
    ax.set_title(f"SAP fold: depth {s['box_depth_ppm']:.0f} +/- {s['empirical_err_ppm']:.0f} ppm; "
                 f"SAP/PSF = {s['ratio_sap_over_psf']:.2f} +/- {s['ratio_err']:.2f} (CROWDSAP {CROWDSAP:.2f})",
                 loc="left", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")

    # (d) odd vs even binned folds
    ax = fig.add_subplot(gs[2, 0])
    epoch = np.round((t - t0) / P).astype(int)
    for sel, col, lab, off in (((epoch % 2) != 0, C_ODD, "odd", -0.08), ((epoch % 2) == 0, C_EVEN, "even", 0.08)):
        ww = sel & (np.abs(ph) < 12)
        xb, yb, eb = binned(ph[ww], y[ww], 24, -12, 12)
        ax.errorbar(xb + off, yb, eb, fmt="o", ms=3, color=col, ecolor=col, elinewidth=0.8, capsize=0, label=lab)
    oe = R["odd_even"]
    ax.set_ylim(1 - 0.012, 1.008)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("relative flux")
    ax.set_title(f"odd {oe['odd']['depth_ppm']:.0f}+/-{oe['odd']['err_ppm']:.0f} ppm vs even "
                 f"{oe['even']['depth_ppm']:.0f}+/-{oe['even']['err_ppm']:.0f} ppm  ({oe['sigma']:.1f} sigma)",
                 loc="left", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    # (e) secondary: fold around phase 0.5
    ax = fig.add_subplot(gs[2, 1])
    ph5 = (((phase_of(t, P, t0) - 0.5 + 0.5) % 1.0) - 0.5) * P * 24
    w5 = np.abs(ph5) < 12
    ax.plot(ph5[w5], y[w5], ".", ms=1.5, color=GRID, alpha=0.8, rasterized=True)
    xb, yb, eb = binned(ph5[w5], y[w5], 48, -12, 12)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3, color=C_PSF, ecolor=C_PSF, elinewidth=0.8, capsize=0)
    ax.axvspan(-dur * 12, dur * 12, color=INK2, alpha=0.08, lw=0)
    ax.set_ylim(1 - 0.012, 1.008)
    ax.set_xlabel("hours from phase 0.5"); ax.set_ylabel("relative flux")
    sec = R["secondary"]
    ax.set_title(f"secondary at phase 0.5: {sec['depth_ppm']:.0f} +/- {sec['err_ppm']:.0f} ppm ({sec['sigma']:.1f} sigma); "
                 f"deepest other box {sec['deepest_other_depth_ppm']:.0f} ppm at phase {sec['deepest_other_phase']:.2f}",
                 loc="left", fontsize=8.5, color=INK)

    # (f) per-transit depths vs time, per-year means
    ax = fig.add_subplot(gs[3, 0])
    rows = R["per_transit"]
    tc = np.array([r["tc_bkjd"] for r in rows]); dp = np.array([r["depth_ppm"] for r in rows])
    ep = np.array([r["err_ppm"] for r in rows])
    ax.errorbar(tc, dp, ep, fmt="o", ms=3, color=C_PSF, ecolor=C_PSF, elinewidth=0.6, capsize=0, label="per-transit")
    ax.axhline(R["depth"]["box_depth_ppm"], color=INK, lw=1.0, label="global box depth")
    ax.axhline(0, color=GRID, lw=0.8)
    for k, (name, v) in enumerate(R["per_year"].items()):
        sel = {"Y1_Q1-4": (q >= 1) & (q <= 4), "Y2_Q5-8": (q >= 5) & (q <= 8),
               "Y3_Q9-12": (q >= 9) & (q <= 12), "Y4_Q13-17": q >= 13}[name]
        if sel.sum() and np.isfinite(v["depth_ppm"]):
            x0, x1 = t[sel].min(), t[sel].max()
            ax.fill_between([x0, x1], v["depth_ppm"] - v["err_ppm"], v["depth_ppm"] + v["err_ppm"],
                            color=C_SAP, alpha=0.18, lw=0)
            ax.plot([x0, x1], [v["depth_ppm"]] * 2, color=C_SAP, lw=1.4,
                    label="per-year box depth +/- emp. err" if k == 0 else None)
    ax.set_xlabel("time (BKJD)"); ax.set_ylabel("depth (ppm)")
    ps = R["per_transit_summary"]
    ax.set_title(f"per-transit local depths: {ps['n_positive']}/{R['n_transits']} positive, deepest/sum = "
                 f"{ps['single_event_dominance']:.2f}, chi2(const) = {ps['chi2_const']:.0f}/{R['n_transits']-1}",
                 loc="left", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")

    # (g) random-phase depth distribution
    ax = fig.add_subplot(gs[3, 1])
    rd = c["rd"] * 1e6
    ax.hist(rd, bins=60, color=GRID, edgecolor=SURF, lw=0.3, label="random-phase boxes (same P, width)")
    ax.axvline(R["depth"]["box_depth_ppm"], color=C_PSF, lw=1.5, label="signal box depth")
    ax.axvline(R["secondary"]["depth_ppm"], color=C_SAP, lw=1.2, ls="--", label="phase-0.5 box")
    ax.set_xlabel("box depth (ppm)"); ax.set_ylabel("count")
    ax.set_title(f"empirical error = {R['depth']['empirical_err_ppm']:.0f} ppm; "
                 f"{R['depth']['n_random_boxes_deeper_than_signal']}/{R['depth']['random_phase_boxes']} random boxes "
                 f"deeper than signal", loc="left", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")

    fig.suptitle(f"Independent re-derivation, {FNAME} (KBONUS-BKG PSF flux, quality==0); candidate-language only",
                 x=0.01, ha="left", fontsize=10.5, color=INK, y=0.995)
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


if __name__ == "__main__":
    jpath = f"{OUT}/rederive_results.json"
    if "--plot-only" in sys.argv:
        R = json.load(open(jpath))
    else:
        R = analyse()
        R["ephemeris"].pop("stats", None)
        with open(jpath, "w") as fh:
            json.dump(to_jsonable(R), fh, indent=2)
        print("wrote", jpath)
    if "--no-plot" not in sys.argv:
        R["figure"] = make_figure(R)
        with open(jpath, "w") as fh:
            json.dump(to_jsonable(R), fh, indent=2)
