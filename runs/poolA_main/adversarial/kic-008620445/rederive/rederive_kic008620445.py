#!/usr/bin/env python
"""
Independent re-derivation of the kic-008620445 first-screen signal.

Written from scratch for the adversarial / second-derivation step (CLAUDE.md
"Before anything goes public" item 4).  It does NOT import or read anything under
pipeline/.  Only numpy, scipy, astropy (BoxLeastSquares) are used.

Steps
 1. load KBONUS-BKG stitched light curve, SAP_QUALITY==0, finite, flux>0
 2. per-quarter median normalisation; split into gap/quarter segments; detrend each
    segment with a time-aware 1.0 d running median (two passes, 4-MAD clipping);
    drop upper (>4 MAD) outliers only
 3. blind BLS 1-50 d on PSF flux (own geometric period grid), same on inverted flux 2-f
 4. refine the top peak; box depth with random-phase empirical error, odd/even,
    secondary at phase 0.5, per-year depths, number of transits, duration
 5. repeat the depth measurement on SAP_FLUX at the PSF ephemeris (dilution check)
Outputs: results.json, events.csv, bls_cache.npz (scratch), figure made by plot script.
"""
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
import scipy.optimize
import astropy
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares

FNAME = "kic-008620445"
FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-008620445_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-008620445/rederive"
SCRATCH = ("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/"
           "scratchpad/rederive-kic008620445")
os.makedirs(OUT, exist_ok=True)
os.makedirs(SCRATCH, exist_ok=True)

RNG = np.random.default_rng(20261017)
PIPELINE_P = 28.07170          # value quoted in the task text, used only for alias check
WINDOW_D = 1.0                 # running-median window
QMIN = 0.0375                  # 0.9 h, sets the period-grid density
OVERSAMPLE = 2                 # period grid oversampling relative to the shift<QMIN rule
PMIN, PMAX = 1.0, 50.0
DURS = np.array([0.03, 0.045, 0.065, 0.09, 0.13, 0.18, 0.25])
METADATA = dict(crowdsap=0.7049495437, flfrcsap=0.7907097573, psffrac=0.5317729187,
                Rstar_rsun=0.8485, Teff=5211.5, dist_pc=1777.1, G=16.0727,
                tpf_host_kic=8620443)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- loading
def load(col):
    with fits.open(FITS) as f:
        d = f["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float)
        y = np.asarray(d[col], float)
        e = np.asarray(d[col + "_ERR"], float)
        q = np.asarray(d["SAP_QUALITY"])
        Q = np.asarray(d["QUARTER"]).astype(int)
    m = (q == 0) & np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (y > 0)
    o = np.argsort(t[m])
    return t[m][o], y[m][o], e[m][o], Q[m][o]


def running_median(tt, yy, w):
    lo = np.searchsorted(tt, tt - w / 2.0)
    hi = np.searchsorted(tt, tt + w / 2.0, side="right")
    return np.array([np.median(yy[a:b]) for a, b in zip(lo, hi)])


def prepare(col):
    """Return detrended, outlier-cleaned (t, f, fe, Q) for flux column `col`."""
    t, y, e, Q = load(col)
    yn = np.empty_like(y)
    en = np.empty_like(e)
    for qq in np.unique(Q):
        s = Q == qq
        med = np.median(y[s])
        yn[s] = y[s] / med
        en[s] = e[s] / med
    breaks = np.where((np.diff(t) > 0.5) | (np.diff(Q) != 0))[0] + 1
    trend = np.empty_like(yn)
    for s in np.split(np.arange(len(t)), breaks):
        tt, yy = t[s], yn[s]
        tr = running_median(tt, yy, WINDOW_D)
        r = yy - tr
        mad = 1.4826 * np.median(np.abs(r - np.median(r)))
        yy2 = yy.copy()
        bad = np.abs(r) > 4 * mad
        yy2[bad] = tr[bad]
        trend[s] = running_median(tt, yy2, WINDOW_D)
    f = yn / trend
    fe = en / trend
    r = f - 1.0
    mad = 1.4826 * np.median(np.abs(r))
    keep = r < 4 * mad                      # upper outliers only; keep dips
    info = dict(n_quality0=int(len(t)), n_upper_outliers_removed=int((~keep).sum()),
                mad_ppm=float(mad * 1e6), std_ppm=float(np.std(r[keep]) * 1e6),
                n_segments=int(len(breaks) + 1))
    return t[keep], f[keep], fe[keep], Q[keep], info


# --------------------------------------------------------------------------- BLS
def period_grid(T):
    ratio = 1.0 + QMIN / (T * OVERSAMPLE)
    return np.exp(np.arange(np.log(PMIN), np.log(PMAX), np.log(ratio)))


def sde_stats(P, power):
    """Return raw SDE and a periodogram-trend-corrected SDE plus the corrected power."""
    raw_sde = (power.max() - power.mean()) / power.std()
    # trend: running median in index space (grid is geometric so this is ~ log P)
    k = 2001
    pad = np.pad(power, k // 2, mode="edge")
    # cheap running median via strided percentile on chunks
    med = np.empty_like(power)
    step = 200
    for i in range(0, len(power), step):
        j = min(i + step, len(power))
        c = pad[i: j + k - 1]
        med[i:j] = np.median(c[: k])   # block-wise approximation, smooth enough
    dpow = power - med
    sde_dt = (dpow.max() - dpow.mean()) / dpow.std()
    return float(raw_sde), float(sde_dt), dpow


def top_peaks(P, power, n=3, excl=0.02):
    """Greedy top-n peaks with +-excl fractional exclusion around each."""
    pw = power.copy()
    peaks = []
    for _ in range(n):
        i = int(np.argmax(pw))
        peaks.append(i)
        m = np.abs(P / P[i] - 1) < excl
        pw[m] = -np.inf
    return peaks


def alias_of(P, Pref):
    r = P / Pref
    for name, val in [("1:1", 1), ("2:1", 2), ("1:2", 0.5), ("3:1", 3), ("1:3", 1 / 3),
                      ("3:2", 1.5), ("2:3", 2 / 3), ("4:1", 4), ("1:4", 0.25)]:
        if abs(r / val - 1) < 0.01:
            return name
    return None


def run_bls(t, f, fe, label, cache):
    if os.path.exists(cache):
        z = np.load(cache)
        return {k: z[k] for k in z.files}
    T = t.max() - t.min()
    P = period_grid(T)
    bls = BoxLeastSquares(t, f, dy=fe)
    t0 = time.time()
    res = bls.power(P, DURS, oversample=5)
    out = dict(period=res.period, power=res.power, depth=res.depth, depth_err=res.depth_err,
               duration=res.duration, transit_time=res.transit_time,
               depth_snr=res.depth_snr, runtime_s=np.array(time.time() - t0))
    np.savez(cache, **out)
    print(f"[{label}] BLS {len(P)} periods x {len(DURS)} durations in {time.time()-t0:.1f} s",
          flush=True)
    return out


# --------------------------------------------------------------------------- box statistics
def phase_time(t, P, t0):
    return ((t - t0 + 0.5 * P) % P) - 0.5 * P


def epoch_index(t, P, t0):
    return np.round((t - t0) / P).astype(int)


def box_depth(t, f, P, t0, dur, sel=None):
    """Depth = <f_out> - <f_in> (ppm) for a box of width dur centred on t0 mod P.
    sel: optional boolean mask of points allowed (e.g. odd epochs, one year)."""
    if sel is None:
        sel = np.ones(len(t), bool)
    dt = np.abs(phase_time(t, P, t0))
    ins = sel & (dt < dur / 2)
    out = sel & (dt > 1.5 * dur)
    n = int(ins.sum())
    if n == 0:
        return np.nan, 0, np.nan
    d = (f[out].mean() - f[ins].mean()) * 1e6
    formal = f[out].std() / np.sqrt(n) * 1e6
    return float(d), n, float(formal)


def empirical_error(t, f, P, t0, dur, sel=None, nrand=3000, avoid=None):
    """Std of box depths at random phases, avoiding +-2 dur around the given phases."""
    if avoid is None:
        avoid = [0.0, 0.5]
    ds, ns = [], []
    while len(ds) < nrand:
        ph = RNG.uniform(0, 1)
        if any(min(abs(ph - a), 1 - abs(ph - a)) * P < 2 * dur for a in avoid):
            continue
        d, n, _ = box_depth(t, f, P, t0 + ph * P, dur, sel)
        if n > 0:
            ds.append(d)
            ns.append(n)
    ds = np.array(ds)
    return float(ds.std()), float(ds.mean()), float(np.median(ns)), ds


def per_event(t, f, Q, P, t0, dur):
    ep = epoch_index(t, P, t0)
    dt = phase_time(t, P, t0)
    rows = []
    for n in range(ep.min(), ep.max() + 1):
        tc = t0 + n * P
        ins = (ep == n) & (np.abs(dt) < dur / 2)
        loc = (np.abs(t - tc) < 1.0) & (np.abs(dt) > 1.5 * dur)
        nin = int(ins.sum())
        if nin == 0:
            # is there data at all within +-1 d?  (gap diagnostic)
            rows.append(dict(epoch=n, tc_bkjd=tc, n_in=0, depth_ppm=np.nan, err_ppm=np.nan,
                             n_local=int(loc.sum()), quarter=int(Q[loc][0]) if loc.any() else -1,
                             nearest_gap_d=float(np.min(np.abs(t - tc)))))
            continue
        base = f[loc].mean() if loc.sum() > 10 else 1.0
        sig = f[loc].std() if loc.sum() > 10 else np.nan
        rows.append(dict(epoch=n, tc_bkjd=tc, n_in=nin, depth_ppm=(base - f[ins].mean()) * 1e6,
                         err_ppm=sig / np.sqrt(nin) * 1e6, n_local=int(loc.sum()),
                         quarter=int(Q[ins][0]), nearest_gap_d=0.0))
    return rows


def trapezoid_fit(t, f, fe, P, t0, dur, depth_guess):
    """Least-squares trapezoid on the fold within +-0.4 d, 30-min cadence integrated."""
    dt = phase_time(t, P, t0)
    m = np.abs(dt) < 0.4
    x, y, w = dt[m], f[m], 1.0 / fe[m]
    cad = 0.020434
    sub = np.linspace(-0.5, 0.5, 7) * cad

    def model(p):
        dt0, dep, t14, frac = p
        t23 = t14 * np.clip(frac, 0, 1)
        xx = np.abs(x[:, None] - dt0 + sub[None, :])
        mm = np.ones_like(xx)
        full = xx <= t23 / 2
        ing = (xx > t23 / 2) & (xx < t14 / 2)
        mm[full] = 1 - dep
        slope = (xx[ing] - t23 / 2) / max((t14 - t23) / 2, 1e-6)
        mm[ing] = 1 - dep * (1 - slope)
        return mm.mean(axis=1)

    def resid(p):
        return (y - model(p)) * w

    p0 = [0.0, depth_guess, dur, 0.5]
    lo = [-0.1, 0.0, 0.01, 0.0]
    hi = [0.1, 0.05, 0.5, 1.0]
    r = scipy.optimize.least_squares(resid, p0, bounds=(lo, hi))
    chi2 = float(np.sum(r.fun ** 2))
    chi2_flat = float(np.sum(((y - 1) * w) ** 2))
    try:
        J = r.jac
        cov = np.linalg.inv(J.T @ J) * chi2 / max(len(y) - 4, 1)
        perr = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        perr = [np.nan] * 4
    return dict(dt0_d=float(r.x[0]), depth_ppm=float(r.x[1] * 1e6), depth_err_ppm=float(perr[1] * 1e6),
                t14_d=float(r.x[2]), t14_err_d=float(perr[2]), t23_over_t14=float(r.x[3]),
                chi2=chi2, chi2_flat=chi2_flat, n=int(m.sum()))


# --------------------------------------------------------------------------- main
def main():
    t_start = time.time()
    results = dict(fname=FNAME, date="2026-10-17", role="independent re-derivation (no pipeline code)",
                   fits=FITS, fits_sha256=sha256(FITS), script=os.path.abspath(__file__),
                   script_sha256=sha256(os.path.abspath(__file__)),
                   versions=dict(python=platform.python_version(), numpy=np.__version__,
                                 scipy=scipy.__version__, astropy=astropy.__version__),
                   method=dict(quality="SAP_QUALITY==0 & finite & flux>0",
                               normalisation="per-quarter median",
                               detrend=f"time-aware running median, window {WINDOW_D} d, 2 passes, "
                                       "4-MAD clip in trend estimate, per gap(>0.5 d)/quarter segment",
                               outliers="upper >4 MAD removed after flattening; dips kept",
                               bls=f"astropy BoxLeastSquares, geometric period grid {PMIN}-{PMAX} d, "
                                   f"step dlnP = {QMIN}/(T*{OVERSAMPLE}), durations {DURS.tolist()} d, "
                                   "objective=likelihood, oversample=5",
                               sde="(max-mean)/std of BLS power over the grid; also trend-corrected",
                               depth_error="std of box depths at 3000 random phases (>2 dur from "
                                           "phase 0 and 0.5), same epochs/points selection",
                               metadata_used=METADATA))

    # ---------------- PSF flux
    t, f, fe, Q, info = prepare("FLUX")
    results["psf_prep"] = info
    T = float(t.max() - t.min())
    results["baseline_d"] = T
    print("PSF prepared:", info, flush=True)

    bls = run_bls(t, f, fe, "PSF", os.path.join(SCRATCH, "bls_psf.npz"))
    P, pw = bls["period"], bls["power"]
    sde_raw, sde_dt, dpow = sde_stats(P, pw)
    peaks_raw = top_peaks(P, pw)
    peaks_dt = top_peaks(P, dpow)
    i0 = peaks_raw[0]
    results["bls_psf"] = dict(
        n_periods=int(len(P)), runtime_s=float(bls["runtime_s"]),
        sde_raw=sde_raw, sde_trend_corrected=sde_dt,
        top3_raw=[dict(period_d=float(P[i]), power=float(pw[i]), depth_ppm=float(bls["depth"][i] * 1e6),
                       duration_d=float(bls["duration"][i]), t0_bkjd=float(bls["transit_time"][i]),
                       sde_raw=float((pw[i] - pw.mean()) / pw.std()),
                       alias_of_28_0717=alias_of(P[i], PIPELINE_P)) for i in peaks_raw],
        top3_trend_corrected=[dict(period_d=float(P[i]), dpower=float(dpow[i]),
                                   sde_tc=float((dpow[i] - dpow.mean()) / dpow.std()),
                                   alias_of_28_0717=alias_of(P[i], PIPELINE_P)) for i in peaks_dt],
    )
    recovered_blind = alias_of(P[i0], PIPELINE_P) == "1:1"
    results["recovered_blind_top_peak_is_28_07"] = bool(recovered_blind)
    results["recovered_blind_any_alias"] = alias_of(P[i0], PIPELINE_P) is not None
    print("PSF top3:", json.dumps(results["bls_psf"]["top3_raw"], indent=1), flush=True)

    # ---------------- inverted null (2 - f), same grid
    blsi = run_bls(t, 2.0 - f, fe, "INV", os.path.join(SCRATCH, "bls_inv.npz"))
    Pi, pwi = blsi["period"], blsi["power"]
    sde_i_raw, sde_i_dt, dpowi = sde_stats(Pi, pwi)
    ii = int(np.argmax(pwi))
    results["bls_inverted_null"] = dict(sde_raw=sde_i_raw, sde_trend_corrected=sde_i_dt,
                                        top_period_d=float(Pi[ii]), top_depth_ppm=float(blsi["depth"][ii] * 1e6),
                                        top_duration_d=float(blsi["duration"][ii]))
    print("INV top:", results["bls_inverted_null"], flush=True)

    # ---------------- refine around the top peak (or around 28.0717 if top is an alias/other)
    Pc = float(P[i0])
    fine_P = np.arange(Pc - 0.03, Pc + 0.03, 2e-5)
    fine_D = np.arange(0.02, 0.305, 0.005)
    b = BoxLeastSquares(t, f, dy=fe)
    rf = b.power(fine_P, fine_D, oversample=10)
    j = int(np.argmax(rf.power))
    Pb, Db, t0b = float(rf.period[j]), float(rf.duration[j]), float(rf.transit_time[j])
    results["refined"] = dict(period_d=Pb, duration_d=Db, duration_h=Db * 24, t0_bkjd=t0b,
                              t0_bjd=t0b + 2454833.0, bls_depth_ppm=float(rf.depth[j] * 1e6),
                              bls_depth_err_ppm=float(rf.depth_err[j] * 1e6),
                              bls_depth_snr=float(rf.depth_snr[j]),
                              duration_profile=dict(duration_d=fine_D.tolist(), power=None))
    # duration profile at fixed best period (one BLS call per duration)
    prof = []
    for dd in fine_D:
        rr = b.power(np.array([Pb]), np.array([dd]), oversample=10)
        prof.append(float(np.atleast_1d(rr.power)[0]))
    results["refined"]["duration_profile"]["power"] = prof
    # half-max width of duration profile as a crude duration uncertainty
    prof = np.array(prof)
    ok = prof >= prof.max() - 0.5 * (prof.max() - np.median(prof))
    results["refined"]["duration_range_halfmax_d"] = [float(fine_D[ok].min()), float(fine_D[ok].max())]
    print("refined:", Pb, Db, t0b, flush=True)

    # period uncertainty from the fine BLS: range where power > max - 1 (delta loglike = 1 ~ 1 sigma)
    perP = np.unique(rf.period)
    powP = np.array([rf.power[rf.period == p].max() for p in perP])
    okP = powP > powP.max() - 1.0
    results["refined"]["period_range_dloglike1_d"] = [float(perP[okP].min()), float(perP[okP].max())]

    # ---------------- statistics at the best ephemeris
    dep, nin, formal = box_depth(t, f, Pb, t0b, Db)
    err, rmean, nmed, rand_ds = empirical_error(t, f, Pb, t0b, Db)
    results["primary"] = dict(depth_ppm=dep, n_in_transit_points=nin, formal_err_ppm=formal,
                              empirical_err_ppm=err, empirical_snr=dep / err,
                              random_box_mean_ppm=rmean, random_box_median_npts=nmed,
                              frac_random_boxes_deeper=float(np.mean(rand_ds >= dep)))
    print("primary:", results["primary"], flush=True)

    # odd / even
    ep = epoch_index(t, Pb, t0b)
    odd, even = (ep % 2 == 1), (ep % 2 == 0)
    d_o, n_o, f_o = box_depth(t, f, Pb, t0b, Db, odd)
    d_e, n_e, f_e = box_depth(t, f, Pb, t0b, Db, even)
    e_o = empirical_error(t, f, Pb, t0b, Db, odd, nrand=1500)[0]
    e_e = empirical_error(t, f, Pb, t0b, Db, even, nrand=1500)[0]
    results["odd_even"] = dict(odd_depth_ppm=d_o, odd_n=n_o, odd_err_ppm=e_o, odd_formal_err_ppm=f_o,
                               even_depth_ppm=d_e, even_n=n_e, even_err_ppm=e_e, even_formal_err_ppm=f_e,
                               diff_ppm=d_o - d_e, sigma=abs(d_o - d_e) / np.hypot(e_o, e_e))
    print("odd/even:", results["odd_even"], flush=True)

    # secondary at phase 0.5
    d_s, n_s, f_s = box_depth(t, f, Pb, t0b + Pb / 2, Db)
    results["secondary"] = dict(depth_ppm=d_s, n=n_s, formal_err_ppm=f_s, empirical_err_ppm=err,
                                sigma=d_s / err)
    # also scan the deepest secondary anywhere in phase 0.2-0.8 (post-hoc diagnostic)
    phs = np.linspace(0.2, 0.8, 601)
    sec_scan = np.array([box_depth(t, f, Pb, t0b + ph * Pb, Db)[0] for ph in phs])
    k = int(np.nanargmax(sec_scan))
    results["secondary"]["deepest_in_phase_0.2_0.8"] = dict(phase=float(phs[k]), depth_ppm=float(sec_scan[k]),
                                                             sigma=float(sec_scan[k] / err))
    print("secondary:", results["secondary"], flush=True)

    # per-year depths (Kepler years by quarter)
    years = {"Y1_Q1-4": (1, 4), "Y2_Q5-8": (5, 8), "Y3_Q9-12": (9, 12), "Y4_Q13-17": (13, 17)}
    py = {}
    for name, (qa, qb) in years.items():
        sel = (Q >= qa) & (Q <= qb)
        d, n, fo = box_depth(t, f, Pb, t0b, Db, sel)
        e = empirical_error(t, f, Pb, t0b, Db, sel, nrand=1000)[0] if n > 0 else np.nan
        py[name] = dict(depth_ppm=d, n=n, empirical_err_ppm=e, formal_err_ppm=fo,
                        sigma=d / e if n > 0 else np.nan)
    results["per_year"] = py
    yd = np.array([v["depth_ppm"] for v in py.values()])
    ye = np.array([v["empirical_err_ppm"] for v in py.values()])
    wmean = np.nansum(yd / ye ** 2) / np.nansum(1 / ye ** 2)
    results["per_year_chi2_const"] = float(np.nansum(((yd - wmean) / ye) ** 2))
    results["per_year_dof"] = int(np.sum(np.isfinite(yd)) - 1)
    print("per year:", py, flush=True)

    # per-quarter depths (diagnostic)
    pq = {}
    for qq in np.unique(Q):
        sel = Q == qq
        d, n, fo = box_depth(t, f, Pb, t0b, Db, sel)
        pq[int(qq)] = dict(depth_ppm=d, n=n, formal_err_ppm=fo)
    results["per_quarter"] = pq

    # events
    rows = per_event(t, f, Q, Pb, t0b, Db)
    n_expected = len(rows)
    n_obs = sum(1 for r in rows if r["n_in"] >= 1)
    n_obs2 = sum(1 for r in rows if r["n_in"] >= 2)
    results["transits"] = dict(n_epochs_in_baseline=n_expected, n_with_ge1_point=n_obs,
                               n_with_ge2_points=n_obs2,
                               points_per_transit_expected=float(Db / 0.020434))
    ev_d = np.array([r["depth_ppm"] for r in rows if r["n_in"] >= 1])
    ev_e = np.array([r["err_ppm"] for r in rows if r["n_in"] >= 1])
    good = np.isfinite(ev_d) & np.isfinite(ev_e)
    if good.sum() > 2:
        w = 1 / ev_e[good] ** 2
        results["transits"]["event_weighted_mean_depth_ppm"] = float(np.sum(w * ev_d[good]) / np.sum(w))
        results["transits"]["event_chi2_const"] = float(np.sum(w * (ev_d[good] - np.sum(w * ev_d[good]) / np.sum(w)) ** 2))
        results["transits"]["event_dof"] = int(good.sum() - 1)
        # single-event dominance: largest single-event contribution to total (depth/err)^2 sum
        contrib = (ev_d[good] / ev_e[good]) ** 2
        contrib = contrib[ev_d[good] > 0]
        results["transits"]["max_single_event_snr2_fraction"] = float(contrib.max() / contrib.sum()) if contrib.size else np.nan
        # leave-one-out: depth without the deepest event
        kmax = int(np.nanargmax(ev_d))
        ep_list = [r["epoch"] for r in rows if r["n_in"] >= 1]
        sel = ep != ep_list[kmax]
        d_loo, n_loo, _ = box_depth(t, f, Pb, t0b, Db, sel)
        e_loo = empirical_error(t, f, Pb, t0b, Db, sel, nrand=1000)[0]
        results["transits"]["leave_out_deepest_event"] = dict(epoch=ep_list[kmax], depth_ppm=d_loo,
                                                              err_ppm=e_loo, snr=d_loo / e_loo)
    with open(os.path.join(OUT, "events.csv"), "w") as fh:
        fh.write("epoch,tc_bkjd,quarter,n_in,depth_ppm,err_ppm,n_local,nearest_gap_d\n")
        for r in rows:
            fh.write(f"{r['epoch']},{r['tc_bkjd']:.5f},{r['quarter']},{r['n_in']},"
                     f"{r['depth_ppm']:.1f},{r['err_ppm']:.1f},{r['n_local']},{r['nearest_gap_d']:.3f}\n")
    print("transits:", results["transits"], flush=True)

    # trapezoid duration fit
    results["trapezoid_fit"] = trapezoid_fit(t, f, fe, Pb, t0b, Db, dep * 1e-6)
    print("trapezoid:", results["trapezoid_fit"], flush=True)

    # expected duration for a central transit on the catalogued star (sanity, not a cut)
    # T14 ~ (Rstar/a) * P/pi for b=0, circular; a from Kepler's third law with M~R for a dwarf
    R, M = METADATA["Rstar_rsun"], METADATA["Rstar_rsun"]   # crude M ~ R for ~5200 K dwarf
    a_au = (M * (Pb / 365.25) ** 2) ** (1 / 3)
    a_rs = a_au * 215.032 / R
    t14_central_h = Pb / np.pi / a_rs * 24
    results["expected_central_duration_h"] = float(t14_central_h)
    results["duration_ratio_to_central"] = float(Db * 24 / t14_central_h)

    # ---------------- SAP flux at the PSF ephemeris
    ts, fs, fes, Qs, infos = prepare("SAP_FLUX")
    results["sap_prep"] = infos
    d_sap, n_sap, f_sap = box_depth(ts, fs, Pb, t0b, Db)
    e_sap = empirical_error(ts, fs, Pb, t0b, Db)[0]
    # narrow BLS on SAP around the PSF period to see if it is independently detected
    bs = BoxLeastSquares(ts, fs, dy=fes)
    rs = bs.power(fine_P, fine_D, oversample=10)
    js = int(np.argmax(rs.power))
    # local SDE proxy on the narrow window is meaningless; instead compare power at Pb vs the
    # distribution of powers in the same narrow window
    results["sap"] = dict(depth_ppm=d_sap, n=n_sap, formal_err_ppm=f_sap, empirical_err_ppm=e_sap,
                          snr=d_sap / e_sap,
                          narrow_bls_best_period_d=float(rs.period[js]), narrow_bls_best_duration_d=float(rs.duration[js]),
                          narrow_bls_best_depth_ppm=float(rs.depth[js] * 1e6),
                          narrow_bls_peak_within_1e3_d_of_psf_period=bool(abs(rs.period[js] - Pb) < 1e-3))
    ratio = d_sap / dep
    ratio_err = ratio * np.hypot(e_sap / d_sap, err / dep) if d_sap != 0 else np.nan
    # expected ratio if the signal is on the target: target fraction of SAP aperture flux = CROWDSAP
    exp_on_target = METADATA["crowdsap"]
    results["sap_vs_psf"] = dict(ratio=ratio, ratio_err=ratio_err,
                                 expected_ratio_if_on_target_crowdsap=exp_on_target,
                                 psf_median_e_s=None, sap_median_e_s=None)
    with fits.open(FITS) as fh:
        dd = fh["LIGHTCURVE_STITCHED"].data
        results["sap_vs_psf"]["psf_median_e_s"] = float(np.nanmedian(dd["FLUX"]))
        results["sap_vs_psf"]["sap_median_e_s"] = float(np.nanmedian(dd["SAP_FLUX"]))
    print("sap:", results["sap"], results["sap_vs_psf"], flush=True)

    # save fold data for the plot script
    np.savez(os.path.join(SCRATCH, "fold_data.npz"), t=t, f=f, fe=fe, Q=Q, ts=ts, fs=fs,
             P=P, power=pw, dpower=dpow, Pi=Pi, poweri=pwi, Pb=Pb, Db=Db, t0b=t0b,
             rand_ds=rand_ds, fine_D=fine_D, prof=prof)

    results["runtime_s"] = time.time() - t_start
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(results, fh, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
    print("DONE", results["runtime_s"], flush=True)


if __name__ == "__main__":
    main()
