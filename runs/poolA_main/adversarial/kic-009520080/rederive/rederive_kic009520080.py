#!/usr/bin/env python
"""
Independent re-derivation for kic-009520080 (KBONUS-BKG PSF light curve).

Written from scratch with numpy / scipy / astropy only.  Deliberately does NOT
import anything from the repo's pipeline/ directory and does not use wotan or
TLS.  Steps:

  1. load the stitched light curve, keep SAP_QUALITY == 0, normalise each
     quarter by its median, detrend with an iterative time-based running median
     (window 1.0 d, negative excursions protected from the trend estimate)
  2. blind BLS, 1-50 d, own log-uniform period grid, astropy BoxLeastSquares
  3. refine the best peak; measure box depth with an empirical random-phase
     error, odd/even, secondary at phase 0.5, per-year depths, number of
     transits, duration (fine BLS duration scan + trapezoid least squares)
  4. same measurement on SAP_FLUX for a dilution check
  5. dump results JSON + one fold figure.
"""
import json
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
        "hlsp_kbonus-bkg_kepler_kepler_kic-009520080_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-009520080/rederive"
FNAME = "kic-009520080"
PIPE_P = 2.37700          # pipeline value we are testing against (period only)
BKJD0 = 2454833.0
RNG = np.random.default_rng(20261017)

# --------------------------------------------------------------------------
# 1. loading + detrending
# --------------------------------------------------------------------------

def load(col, errcol):
    with fits.open(FITS) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d[col], float)
        e = np.asarray(d[errcol], float)
        q = np.asarray(d["SAP_QUALITY"], int)
        Q = np.asarray(d["QUARTER"], int)
    m = (q == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (f > 0)
    o = np.argsort(t[m])
    return t[m][o], f[m][o], e[m][o], Q[m][o]


def normalise_per_quarter(f, e, Q):
    fn = f.copy()
    en = e.copy()
    for qq in np.unique(Q):
        s = Q == qq
        med = np.nanmedian(f[s])
        fn[s] = f[s] / med
        en[s] = e[s] / med
    return fn, en


def running_median(t, f, good, window=1.0):
    """Time-based running median using only `good` points; evaluated at all t."""
    tg = t[good]
    fg = f[good]
    half = window / 2.0
    lo = np.searchsorted(tg, t - half, side="left")
    hi = np.searchsorted(tg, t + half, side="right")
    out = np.empty_like(f)
    for i in range(len(t)):
        a, b = lo[i], hi[i]
        if b - a < 5:           # pathological edge: widen
            a = max(0, a - 5)
            b = min(len(tg), b + 5)
        out[i] = np.median(fg[a:b])
    return out


def detrend(t, f, window=1.0, niter=3):
    """Iterative running-median detrend.  Points more than 3 sigma below the
    trend are excluded from the trend estimate (protects transits), as are
    points more than 3 sigma above it (flares / cosmic rays)."""
    good = np.ones(len(t), bool)
    trend = running_median(t, f, good, window)
    for _ in range(niter):
        r = f / trend - 1.0
        sig = 1.4826 * np.median(np.abs(r - np.median(r)))
        good = (r > -3 * sig) & (r < 3 * sig)
        trend = running_median(t, f, good, window)
    flat = f / trend
    r = flat - 1.0
    sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    # final clipping of the *search* series: kill positive outliers and gross
    # negative ones (>15 sigma; a 1500 ppm transit at 1600 ppm noise is ~1 sigma)
    keep = (r < 4 * sig) & (r > -15 * sig)
    return flat, keep, sig


# --------------------------------------------------------------------------
# 2. BLS helpers
# --------------------------------------------------------------------------

def period_grid(pmin, pmax, baseline, dmin, oversample=3.0):
    dlnp = dmin / (oversample * baseline)
    n = int(np.ceil(np.log(pmax / pmin) / dlnp))
    return pmin * np.exp(dlnp * np.arange(n + 1))


def top_peaks(periods, power, n=3, exclude_frac=0.02):
    """Return n highest, mutually distinct peaks (distinct = >2 % apart in P)."""
    idx = np.argsort(power)[::-1]
    out = []
    for i in idx:
        p = periods[i]
        if all(abs(p / q - 1) > exclude_frac for q, _ in out):
            out.append((p, power[i]))
        if len(out) >= n:
            break
    return out


def alias_name(p, pref, tol=0.002):
    for k, name in [(1, "same"), (2, "2x"), (0.5, "1/2"), (3, "3x"), (1 / 3, "1/3"),
                    (1.5, "3/2"), (2 / 3, "2/3"), (4, "4x"), (0.25, "1/4")]:
        if abs(p / (pref * k) - 1) < tol:
            return name
    return None


# --------------------------------------------------------------------------
# 3. measurement helpers
# --------------------------------------------------------------------------

def phase_fold(t, P, t0):
    ph = ((t - t0) / P + 0.5) % 1.0 - 0.5
    return ph


def box_mean(ph, f, P, dur, centre=0.0):
    w = dur / P / 2.0
    d = (ph - centre + 0.5) % 1.0 - 0.5
    s = np.abs(d) < w
    if s.sum() == 0:
        return np.nan, 0
    return f[s].mean(), int(s.sum())


def empirical_box_error(ph, f, P, dur, nsamp=3000, avoid=(0.0, 0.5)):
    """std of box means placed at random phases, avoiding +-2 dur around `avoid`."""
    w = dur / P
    vals = []
    tries = 0
    while len(vals) < nsamp and tries < 20 * nsamp:
        tries += 1
        c = RNG.uniform(-0.5, 0.5)
        if any(abs(((c - a + 0.5) % 1.0) - 0.5) < 2 * w for a in avoid):
            continue
        m, n = box_mean(ph, f, P, dur, c)
        if n > 0:
            vals.append(m)
    vals = np.array(vals)
    return np.std(vals), np.mean(vals), len(vals)


def box_depth(ph, f, P, dur, centre=0.0, nsamp=3000):
    """depth (ppm) = out-of-box level - in-box mean; error from random boxes."""
    m_in, n_in = box_mean(ph, f, P, dur, centre)
    d = (ph - centre + 0.5) % 1.0 - 0.5
    out = np.abs(d) > 1.5 * dur / P
    base = np.median(f[out])
    err, _, nb = empirical_box_error(ph, f, P, dur, nsamp)
    depth = (base - m_in) * 1e6
    return dict(depth_ppm=float(depth), err_ppm=float(err * 1e6),
                snr=float(depth / (err * 1e6)) if err > 0 else np.nan,
                n_in=int(n_in), n_boxes=int(nb))


def trapezoid(params, ph, P):
    depth, t14, tf, dt0 = params   # tf = T23/T14 in (0,1)
    t14 = abs(t14)
    tf = np.clip(tf, 0.0, 1.0)
    x = np.abs((ph * P) - dt0)
    t23 = tf * t14
    model = np.ones_like(ph)
    full = x <= t23 / 2
    ing = (x > t23 / 2) & (x < t14 / 2)
    model[full] = 1 - depth
    model[ing] = 1 - depth * (t14 / 2 - x[ing]) / ((t14 - t23) / 2 + 1e-12)
    return model


def fit_trapezoid(ph, f, P, dur0, depth0):
    s = np.abs(ph * P) < 4 * dur0
    x, y = ph[s], f[s]

    def resid(p):
        return y - trapezoid(p, x, P)

    best = None
    for tf in (0.2, 0.5, 0.8):
        p0 = [depth0, dur0, tf, 0.0]
        r = least_squares(resid, p0,
                          bounds=([0, 0.01, 0.0, -dur0], [0.1, 0.5, 1.0, dur0]))
        if best is None or r.cost < best.cost:
            best = r
    depth, t14, tf, dt0 = best.x
    # crude errors from the Jacobian
    try:
        J = best.jac
        cov = np.linalg.inv(J.T @ J) * np.sum(best.fun ** 2) / (len(y) - 4)
        perr = np.sqrt(np.diag(cov))
    except Exception:
        perr = [np.nan] * 4
    return dict(depth_ppm=float(depth * 1e6), depth_err_ppm=float(perr[0] * 1e6),
                t14_h=float(abs(t14) * 24), t14_err_h=float(perr[1] * 24),
                t23_over_t14=float(np.clip(tf, 0, 1)), dt0_d=float(dt0),
                n_fit=int(s.sum()))


def count_transits(t, P, t0, dur):
    n = np.round((t - t0) / P).astype(int)
    intr = np.abs(t - (t0 + n * P)) < dur / 2
    epochs = np.unique(n[intr])
    # expected number of cadences per transit at 29.4 min
    exp_per = dur / (29.4 / 60 / 24)
    cnt = np.array([np.sum(intr & (n == k)) for k in epochs])
    half = epochs[cnt >= 0.5 * exp_per]
    n_total_windows = int(np.floor((t.max() - t0) / P) - np.ceil((t.min() - t0) / P)) + 1
    return dict(n_epochs_with_data=int(len(epochs)),
                n_epochs_ge_half_coverage=int(len(half)),
                n_windows_in_baseline=n_total_windows,
                expected_cadences_per_transit=float(exp_per),
                total_in_transit_cadences=int(intr.sum()))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    res = {"fname": FNAME, "date": "2026-10-17",
           "method": ("own code: quality==0, per-quarter median normalisation, "
                      "iterative 1.0 d running median (3-sigma protected), "
                      "astropy BoxLeastSquares on own log-uniform grid, "
                      "random-phase empirical errors")}
    T0 = time.time()

    # ---------------- PSF flux ----------------
    t, f, e, Q = load("FLUX", "FLUX_ERR")
    res["n_quality0"] = int(len(t))
    fn, en = normalise_per_quarter(f, e, Q)
    flat, keep, sig = detrend(t, fn, window=1.0)
    res["n_after_clip"] = int(keep.sum())
    res["robust_point_scatter_ppm"] = float(sig * 1e6)
    t_s, f_s, e_s, Q_s = t[keep], flat[keep], en[keep], Q[keep]
    baseline = t_s.max() - t_s.min()
    res["baseline_d"] = float(baseline)
    print(f"loaded {len(t)} q0 points, {keep.sum()} after clip, scatter {sig*1e6:.0f} ppm,"
          f" baseline {baseline:.1f} d  [{time.time()-T0:.1f}s]")

    # ---------------- blind BLS 1-50 d ----------------
    durations = np.array([0.03, 0.045, 0.06, 0.08, 0.10, 0.13, 0.17, 0.22])
    periods = period_grid(1.0, 50.0, baseline, durations.min(), oversample=2.5)
    res["bls_grid"] = dict(pmin=1.0, pmax=50.0, n_periods=int(len(periods)),
                           durations_d=durations.tolist(), oversample=2.5,
                           objective="likelihood")
    print(f"BLS grid: {len(periods)} periods")
    bls = BoxLeastSquares(t_s, f_s, dy=e_s)
    pg = bls.power(periods, durations, objective="likelihood")
    power = np.asarray(pg.power)
    print(f"BLS done [{time.time()-T0:.1f}s]")
    sde_all = (power.max() - power.mean()) / power.std()
    # also an SDE with a running-median-flattened periodogram (removes red trend)
    from scipy.ndimage import median_filter
    pflat = power - median_filter(power, size=2001, mode="nearest")
    sde_flat = pflat.max() / pflat.std()
    peaks = top_peaks(periods, power, n=3)
    res["blind_bls"] = dict(
        top3=[dict(period_d=float(p), power=float(pw),
                   alias_of_2p377=alias_name(p, PIPE_P)) for p, pw in peaks],
        sde_raw=float(sde_all), sde_detrended_periodogram=float(sde_flat),
        power_median=float(np.median(power)), power_std=float(power.std()))
    best_P = peaks[0][0]
    res["recovered_blind"] = bool(abs(best_P / PIPE_P - 1) < 0.002)
    print("top3:", res["blind_bls"]["top3"], "SDE raw", sde_all, "flat", sde_flat)

    # ---------------- inverted-flux null on the same grid ----------------
    bls_inv = BoxLeastSquares(t_s, 2.0 - f_s, dy=e_s)
    pinv = np.asarray(bls_inv.power(periods, durations, objective="likelihood").power)
    sde_inv = (pinv.max() - pinv.mean()) / pinv.std()
    pinv_flat = pinv - median_filter(pinv, size=2001, mode="nearest")
    res["inverted_null"] = dict(best_period_d=float(periods[np.argmax(pinv)]),
                                sde_raw=float(sde_inv),
                                sde_detrended_periodogram=float(pinv_flat.max() / pinv_flat.std()))
    print("inverted null:", res["inverted_null"], f"[{time.time()-T0:.1f}s]")

    # ---------------- refine best peak ----------------
    pfine = best_P * np.exp(np.linspace(-3e-4, 3e-4, 1201))
    dfine = np.arange(0.02, 0.2001, 0.0025)
    pgf = bls.power(pfine, dfine, objective="likelihood")
    i = int(np.argmax(pgf.power))
    P = float(pgf.period[i]); dur = float(pgf.duration[i]); t0 = float(pgf.transit_time[i])
    bls_depth = float(pgf.depth[i]); bls_depth_err = float(pgf.depth_err[i])
    # one more very fine pass in period with the best duration
    pfine2 = P * np.exp(np.linspace(-2e-5, 2e-5, 801))
    pgf2 = bls.power(pfine2, dur, objective="likelihood")
    i = int(np.argmax(pgf2.power))
    P = float(pgf2.period[i]); t0 = float(pgf2.transit_time[i])
    # move t0 near the middle of the baseline to decorrelate from P
    nmid = np.round((0.5 * (t_s.min() + t_s.max()) - t0) / P)
    t0_mid = t0 + nmid * P
    stats = bls.compute_stats(P, dur, t0_mid)
    res["ephemeris"] = dict(period_d=P, t0_bkjd=t0_mid, t0_bjd=t0_mid + BKJD0,
                            duration_bls_d=dur, duration_bls_h=dur * 24,
                            bls_depth_ppm=bls_depth * 1e6, bls_depth_err_ppm=bls_depth_err * 1e6,
                            bls_stats_depth_ppm=float(stats["depth"][0] * 1e6),
                            bls_stats_depth_err_ppm=float(stats["depth"][1] * 1e6),
                            bls_stats_depth_odd_ppm=float(stats["depth_odd"][0] * 1e6),
                            bls_stats_depth_even_ppm=float(stats["depth_even"][0] * 1e6),
                            bls_stats_depth_half_ppm=float(stats["depth_half"][0] * 1e6),
                            bls_stats_depth_half_err_ppm=float(stats["depth_half"][1] * 1e6),
                            bls_stats_harmonic_delta_loglike=float(stats["harmonic_delta_log_likelihood"]))
    print(f"refined: P={P:.6f} t0={t0_mid:.5f} dur={dur*24:.2f} h depth={bls_depth*1e6:.0f} ppm"
          f" [{time.time()-T0:.1f}s]")

    # ---------------- box measurements with empirical errors ----------------
    ph = phase_fold(t_s, P, t0_mid)
    prim = box_depth(ph, f_s, P, dur, 0.0)
    sec = box_depth(ph, f_s, P, dur, 0.5)
    n = np.round((t_s - t0_mid) / P).astype(int)
    odd = n % 2 == 1
    d_odd = box_depth(ph[odd], f_s[odd], P, dur, 0.0, nsamp=1500)
    d_even = box_depth(ph[~odd], f_s[~odd], P, dur, 0.0, nsamp=1500)
    oe_sigma = abs(d_odd["depth_ppm"] - d_even["depth_ppm"]) / np.hypot(d_odd["err_ppm"], d_even["err_ppm"])
    res["primary"] = prim
    res["secondary_phase05"] = sec
    res["odd"] = d_odd
    res["even"] = d_even
    res["odd_even_sigma"] = float(oe_sigma)
    print("primary", prim, "\nsecondary", sec, f"\nodd/even {d_odd['depth_ppm']:.0f} / {d_even['depth_ppm']:.0f} -> {oe_sigma:.2f} sigma")

    # per-year (Kepler years by quarter)
    years = {"Y1_Q1-4": (1, 4), "Y2_Q5-8": (5, 8), "Y3_Q9-12": (9, 12), "Y4_Q13-17": (13, 17)}
    py = {}
    for k, (a, b) in years.items():
        s = (Q_s >= a) & (Q_s <= b)
        if s.sum() > 100:
            py[k] = box_depth(ph[s], f_s[s], P, dur, 0.0, nsamp=1500)
            py[k]["n_points"] = int(s.sum())
    res["per_year"] = py
    print("per year:", {k: (round(v["depth_ppm"]), round(v["err_ppm"])) for k, v in py.items()})

    # per-quarter depths (quick look)
    pq = {}
    for qq in np.unique(Q_s):
        s = Q_s == qq
        m_in, n_in = box_mean(ph[s], f_s[s], P, dur, 0.0)
        out = np.abs(ph[s]) > 1.5 * dur / P
        pq[int(qq)] = dict(depth_ppm=float((np.median(f_s[s][out]) - m_in) * 1e6), n_in=int(n_in))
    res["per_quarter_depth_ppm"] = pq

    # transit count / duration
    res["transits"] = count_transits(t_s, P, t0_mid, dur)
    trap = fit_trapezoid(ph, f_s, P, dur, prim["depth_ppm"] * 1e-6)
    res["trapezoid_fit"] = trap
    print("transits:", res["transits"], "\ntrapezoid:", trap)

    # ---------------- SAP flux dilution check ----------------
    ts, fs, es, Qs = load("SAP_FLUX", "SAP_FLUX_ERR")
    fsn, esn = normalise_per_quarter(fs, es, Qs)
    flat_s, keep_s, sig_s = detrend(ts, fsn, window=1.0)
    ph_s = phase_fold(ts[keep_s], P, t0_mid)
    sap_prim = box_depth(ph_s, flat_s[keep_s], P, dur, 0.0)
    sap_sec = box_depth(ph_s, flat_s[keep_s], P, dur, 0.5)
    ratio = sap_prim["depth_ppm"] / prim["depth_ppm"]
    ratio_err = ratio * np.hypot(sap_prim["err_ppm"] / sap_prim["depth_ppm"], prim["err_ppm"] / prim["depth_ppm"])
    res["sap"] = dict(n_points=int(keep_s.sum()), scatter_ppm=float(sig_s * 1e6),
                      primary=sap_prim, secondary=sap_sec,
                      depth_ratio_sap_over_psf=float(ratio), depth_ratio_err=float(ratio_err))
    # also a free BLS on SAP near the period to see if it is independently detected
    bls_sap = BoxLeastSquares(ts[keep_s], flat_s[keep_s], dy=esn[keep_s])
    pgs = bls_sap.power(pfine, dfine, objective="likelihood")
    j = int(np.argmax(pgs.power))
    res["sap"]["bls_near_peak"] = dict(period_d=float(pgs.period[j]), duration_h=float(pgs.duration[j] * 24),
                                       depth_ppm=float(pgs.depth[j] * 1e6), depth_err_ppm=float(pgs.depth_err[j] * 1e6),
                                       depth_snr=float(pgs.depth_snr[j]))
    print("SAP:", res["sap"], f"[{time.time()-T0:.1f}s]")

    # ---------------- summary ----------------
    res["summary"] = dict(
        period_d=P, t0_bjd=t0_mid + BKJD0, duration_h=dur * 24,
        depth_ppm=prim["depth_ppm"], depth_err_ppm=prim["err_ppm"], empirical_snr=prim["snr"],
        sde=float(sde_all), n_transits=res["transits"]["n_epochs_with_data"],
        odd_even_sigma=float(oe_sigma), secondary_ppm=sec["depth_ppm"], secondary_sigma=sec["snr"],
        per_year_depths_ppm=[py[k]["depth_ppm"] for k in py],
        recovered_blind=res["recovered_blind"], sap_over_psf_depth_ratio=float(ratio))
    res["runtime_s"] = float(time.time() - T0)

    with open(f"{OUT}/rederive_results.json", "w") as fh:
        json.dump(res, fh, indent=2, default=float)

    # ---------------- figure ----------------
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    ax = axes[0, 0]
    ax.plot(periods, power, lw=0.4, color="0.3")
    for k, (pp, pw) in enumerate(peaks):
        ax.axvline(pp, color="C3" if k == 0 else "C1", alpha=0.5, lw=1)
    ax.axvline(PIPE_P, color="C2", ls="--", lw=1, label=f"pipeline P={PIPE_P}")
    ax.set_xscale("log"); ax.set_xlabel("period [d]"); ax.set_ylabel("BLS power (dlnL)")
    ax.set_title(f"blind BLS 1-50 d, {len(periods)} periods; SDE={sde_all:.1f} (inv. null {sde_inv:.1f})")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    xh = ph * P * 24
    s = np.abs(xh) < 6
    ax.plot(xh[s], f_s[s], ".", ms=1.5, alpha=0.25, color="0.5")
    nb = 48
    edges = np.linspace(-6, 6, nb + 1)
    idx = np.digitize(xh[s], edges) - 1
    bm = np.array([f_s[s][idx == k].mean() if np.any(idx == k) else np.nan for k in range(nb)])
    be = np.array([f_s[s][idx == k].std() / np.sqrt(max(1, (idx == k).sum())) for k in range(nb)])
    bc = 0.5 * (edges[1:] + edges[:-1])
    ax.errorbar(bc, bm, be, fmt="o", ms=3, color="C0", label="15-min bins")
    xx = np.linspace(-6, 6, 600) / 24 / P
    ax.plot(xx * P * 24, trapezoid([trap["depth_ppm"] * 1e-6, trap["t14_h"] / 24, trap["t23_over_t14"], trap["dt0_d"]], xx, P),
            color="C3", lw=1.5, label=f"trapezoid T14={trap['t14_h']:.2f} h, {trap['depth_ppm']:.0f} ppm")
    ax.set_ylim(1 - 6e-3, 1 + 4e-3)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("PSF flux (detrended)")
    ax.set_title(f"P={P:.5f} d; box depth {prim['depth_ppm']:.0f}±{prim['err_ppm']:.0f} ppm (SNR {prim['snr']:.1f})")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    s = np.abs(((ph - 0.5 + 0.5) % 1.0) - 0.5) * P * 24 < 6
    xh2 = (((ph - 0.5 + 0.5) % 1.0) - 0.5) * P * 24
    idx = np.digitize(xh2[s], edges) - 1
    bm2 = np.array([f_s[s][idx == k].mean() if np.any(idx == k) else np.nan for k in range(nb)])
    be2 = np.array([f_s[s][idx == k].std() / np.sqrt(max(1, (idx == k).sum())) for k in range(nb)])
    ax.plot(xh2[s], f_s[s], ".", ms=1.5, alpha=0.25, color="0.5")
    ax.errorbar(bc, bm2, be2, fmt="o", ms=3, color="C0")
    ax.axhline(1, color="k", lw=0.5)
    ax.set_ylim(1 - 6e-3, 1 + 4e-3)
    ax.set_xlabel("hours from phase 0.5"); ax.set_ylabel("PSF flux (detrended)")
    ax.set_title(f"secondary window: {sec['depth_ppm']:.0f}±{sec['err_ppm']:.0f} ppm ({sec['snr']:.1f} sigma)")

    ax = axes[1, 1]
    # odd / even and SAP folds, binned
    for lab, mask, col in [("odd", odd, "C1"), ("even", ~odd, "C2")]:
        xo = ph[mask] * P * 24
        so = np.abs(xo) < 6
        idx = np.digitize(xo[so], edges) - 1
        bmo = np.array([f_s[mask][so][idx == k].mean() if np.any(idx == k) else np.nan for k in range(nb)])
        ax.plot(bc, bmo, "o-", ms=3, lw=0.8, color=col, label=f"{lab}: {res[lab]['depth_ppm']:.0f}±{res[lab]['err_ppm']:.0f} ppm")
    xs = ph_s * P * 24
    ss = np.abs(xs) < 6
    idx = np.digitize(xs[ss], edges) - 1
    bms = np.array([flat_s[keep_s][ss][idx == k].mean() if np.any(idx == k) else np.nan for k in range(nb)])
    ax.plot(bc, bms - 0.004, "s-", ms=3, lw=0.8, color="C4",
            label=f"SAP (offset -0.004): {sap_prim['depth_ppm']:.0f}±{sap_prim['err_ppm']:.0f} ppm")
    ax.set_ylim(1 - 9e-3, 1 + 3e-3)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("binned flux")
    ax.set_title(f"odd/even {oe_sigma:.1f} sigma; SAP/PSF depth ratio {ratio:.2f}±{ratio_err:.2f}")
    ax.legend(fontsize=8)

    fig.suptitle(f"{FNAME} independent re-derivation (own code, no pipeline/ imports) 2026-10-17", y=0.995)
    fig.tight_layout()
    fig.savefig(f"{OUT}/rederive_fold.png", dpi=130)
    print("wrote", f"{OUT}/rederive_results.json", f"{OUT}/rederive_fold.png", f"[{time.time()-T0:.1f}s]")


if __name__ == "__main__":
    main()
