#!/usr/bin/env python
"""
Independent re-derivation for KBONUS-BKG source kic-006448637 (KIC 6448637).

Written from scratch with numpy / scipy / astropy / matplotlib only.
Nothing under pipeline/ is imported or read; wotan / TLS are not used.

Steps
  1. load LIGHTCURVE_STITCHED, keep SAP_QUALITY == 0, normalise per quarter,
     detrend with a per-quarter 1.0-d running median (iterative, upward-outlier aware)
  2. blind BLS 1-50 d on my own log-uniform period grid (astropy BoxLeastSquares)
  3. at the best ephemeris: box depth with an empirical random-phase-box error,
     odd/even, secondary at phase 0.5, per-year depths, per-event table, duration
     (box SNR scan + trapezoid fit) -- all on a masked re-detrend
  4. the same depth measurement on SAP_FLUX (dilution check)
  5. inverted-flux null (same blind search on -flat)

Run:  . /home/user/ZyushgPlanet/.venv/bin/activate && python -I rederive_kic006448637.py
"""
import hashlib
import json
import os
import platform
import time

import numpy as np
import astropy
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.optimize import least_squares
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FNAME = "kic-006448637"
FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-006448637_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006448637/rederive"
# Numbers the task statement quotes for the pipeline; used ONLY for alias bookkeeping.
REF_P, REF_DEPTH_PPM, REF_DUR_H = 35.5873, 846.0, 4.5
BKJD0 = 2454833.0
CAD = 29.4244 / 1440.0          # long-cadence spacing, d
WINDOW = 1.0                    # running-median window, d
PMIN, PMAX = 1.0, 50.0
DURS = np.array([0.05, 0.07, 0.10, 0.14, 0.19, 0.25, 0.33, 0.42])  # d
OVERSAMPLE = 3
RNG = np.random.default_rng(20261017)
T_START = time.time()


def log(*a):
    print("[%6.1fs]" % (time.time() - T_START), *a, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------- loading
def load(col):
    """Return time (BKJD), per-quarter-normalised flux, its error, quarter -- quality==0 only."""
    d = fits.getdata(FITS, "LIGHTCURVE_STITCHED")
    t = np.asarray(d["TIME"], float)
    f = np.asarray(d[col], float)
    e = np.asarray(d[col + "_ERR"], float)
    q = np.asarray(d["SAP_QUALITY"])
    qu = np.asarray(d["QUARTER"]).astype(int)
    good = (q == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0)
    t, f, e, qu = t[good], f[good], e[good], qu[good]
    o = np.argsort(t)
    t, f, e, qu = t[o], f[o], e[o], qu[o]
    fn = np.empty_like(f)
    en = np.empty_like(e)
    for k in np.unique(qu):
        m = qu == k
        med = np.median(f[m])
        fn[m] = f[m] / med
        en[m] = e[m] / med
    return t, fn, en, qu


# ----------------------------------------------------------------------------- detrending
def running_median(t, f, window, use):
    """Time-windowed running median evaluated at every t, estimated from points with use==True."""
    tu, fu = t[use], f[use]
    lo = np.searchsorted(tu, t - window / 2.0)
    hi = np.searchsorted(tu, t + window / 2.0)
    tr = np.full(len(t), np.nan)
    for i in range(len(t)):
        if hi[i] - lo[i] >= 5:
            tr[i] = np.median(fu[lo[i]:hi[i]])
    bad = ~np.isfinite(tr)
    if bad.any():
        tr[bad] = np.interp(t[bad], t[~bad], tr[~bad])
    return tr


def detrend(t, f, qu, window=WINDOW, mask=None, n_iter=2, clip_hi=4.0, clip_lo=6.0):
    """Per-quarter running median. Points with mask==True (e.g. in-transit) and iteratively
    identified upward outliers are excluded from the trend estimate; dips are never excluded.
    Returns flat = f/trend - 1, trend, keep (outlier mask for the search), per-point robust sigma."""
    trend = np.empty_like(f)
    use0 = np.ones(len(f), bool) if mask is None else ~mask
    sig_q = np.empty_like(f)
    for k in np.unique(qu):
        idx = np.where(qu == k)[0]
        u = use0[idx].copy()
        for _ in range(n_iter):
            tr = running_median(t[idx], f[idx], window, u)
            r = f[idx] / tr - 1.0
            s = 1.4826 * np.median(np.abs(r[u] - np.median(r[u])))
            u = u & (r < clip_hi * s)
        trend[idx] = tr
        r = f[idx] / tr - 1.0
        sig_q[idx] = 1.4826 * np.median(np.abs(r[use0[idx]] - np.median(r[use0[idx]])))
    flat = f / trend - 1.0
    keep = (flat < clip_hi * sig_q) & (flat > -clip_lo * sig_q)
    return flat, trend, keep, sig_q


# ----------------------------------------------------------------------------- BLS
def blind_bls(t, y, dy, pmin=PMIN, pmax=PMAX):
    T = t.max() - t.min()
    dlnp = DURS.min() / (OVERSAMPLE * T)
    periods = np.exp(np.arange(np.log(pmin), np.log(pmax), dlnp))
    bls = BoxLeastSquares(t, y, dy=dy)
    res = bls.power(periods, DURS, objective="likelihood")
    return res, periods


def spectrum_stats(periods, power):
    """Raw Kovacs-style SDE plus an SDE after dividing out a smooth trend in log-period."""
    sde_raw = (power.max() - power.mean()) / power.std()
    lnp = np.log(periods)
    nb = 60
    edges = np.linspace(lnp.min(), lnp.max(), nb + 1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    med = np.array([np.median(power[(lnp >= a) & (lnp < b)]) for a, b in zip(edges[:-1], edges[1:])])
    trend = np.interp(lnp, centers, med)
    pn = power / trend
    sde_norm = (pn.max() - pn.mean()) / pn.std()
    return float(sde_raw), float(sde_norm), pn


def top_peaks(periods, power, n=3, sep=0.02):
    order = np.argsort(power)[::-1]
    peaks = []
    for i in order:
        p = periods[i]
        if all(abs(np.log(p / pp["period_d"])) > sep for pp in peaks):
            peaks.append({"period_d": float(p), "power": float(power[i]), "index": int(i)})
        if len(peaks) >= n:
            break
    return peaks


def alias_label(p, ref=REF_P, tol=5e-3):
    for num, den in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (4, 1)]:
        if abs(p / (ref * num / den) - 1.0) < tol:
            return "%d/%d x %.4f" % (num, den, ref) if (num, den) != (1, 1) else "= %.4f" % ref
    return "none"


def refine(t, y, dy, p0):
    bls = BoxLeastSquares(t, y, dy=dy)
    periods = p0 * np.exp(np.linspace(-0.003, 0.003, 3001))
    durs = np.arange(0.04, 0.50, 0.01)
    res = bls.power(periods, durs, objective="likelihood")
    i = int(np.argmax(res.power))
    return dict(period=float(res.period[i]), t0=float(res.transit_time[i]),
                duration=float(res.duration[i]), depth=float(res.depth[i]),
                depth_err=float(res.depth_err[i]), power=float(res.power[i])), bls


# ----------------------------------------------------------------------------- box statistics
def fold_phase(t, P, t0, center=0.0):
    return ((t - t0) / P - center + 0.5) % 1.0 - 0.5     # cycles, 0 at `center`


def measure(t, y, P, t0, dur, sel=None, center=0.0, ntrial=4000):
    """Box depth (fractional) = out-of-transit mean - in-box mean, with an EMPIRICAL error from
    boxes of the same duration dropped at random phases that avoid phases 0 and 0.5."""
    if sel is None:
        sel = np.ones(len(t), bool)
    tt, yy = t[sel], y[sel]
    ph0 = fold_phase(tt, P, t0, 0.0)
    phc = fold_phase(tt, P, t0, center)
    inbox = np.abs(phc * P) < dur / 2.0
    oot = (np.abs(ph0 * P) > dur) & (np.abs((np.abs(ph0) - 0.5) * P) > dur)
    base = np.mean(yy[oot])
    n_in = int(inbox.sum())
    depth = base - np.mean(yy[inbox]) if n_in else np.nan
    hd = dur / P
    vals, nn = [], []
    for _ in range(ntrial):
        c = RNG.uniform(-0.5, 0.5)
        d0 = abs((c + 0.5) % 1.0 - 0.5)
        d5 = abs(abs(c) - 0.5)
        if d0 < 1.5 * hd or d5 < 1.5 * hd:
            continue
        m = np.abs(((ph0 - c + 0.5) % 1.0 - 0.5) * P) < dur / 2.0
        if m.sum() < 3:
            continue
        vals.append(base - np.mean(yy[m]))
        nn.append(m.sum())
    vals, nn = np.array(vals), np.array(nn)
    err = float(np.std(vals))
    err_mad = float(1.4826 * np.median(np.abs(vals - np.median(vals))))
    sim = np.abs(nn - n_in) <= 0.3 * max(n_in, 1)
    err_similar_n = float(np.std(vals[sim])) if sim.sum() > 30 else np.nan
    white = float(np.std(yy[oot]) / np.sqrt(max(n_in, 1)))
    return dict(depth_ppm=1e6 * float(depth), err_ppm=1e6 * err, err_mad_ppm=1e6 * err_mad,
                err_similar_coverage_ppm=1e6 * err_similar_n, white_err_ppm=1e6 * white,
                snr=float(depth / err) if err > 0 else np.nan, n_in=n_in,
                n_random_boxes=int(len(vals)), null_mean_ppm=1e6 * float(np.mean(vals)))


def per_event_table(t, y, P, t0, dur, qu):
    """Local depth per transit epoch, using a +-1.5 d local baseline, plus gap proximity."""
    epochs = np.round((t - t0) / P).astype(int)
    rows = []
    # data gaps > 0.5 d
    dt = np.diff(t)
    gap_starts = t[:-1][dt > 0.5]
    gap_ends = t[1:][dt > 0.5]
    edges = np.concatenate([[t.min()], gap_starts, gap_ends, [t.max()]])
    n_expect = dur / CAD
    for n in range(epochs.min(), epochs.max() + 1):
        tc = t0 + n * P
        if tc < t.min() - 1 or tc > t.max() + 1:
            continue
        loc = np.abs(t - tc) < 1.5
        inb = np.abs(t - tc) < dur / 2.0
        oo = loc & (np.abs(t - tc) > dur)
        n_in = int(inb.sum())
        if n_in == 0 and oo.sum() == 0:
            rows.append(dict(epoch=n, tc_bkjd=tc, n_in=0, n_oot_local=int(oo.sum()), depth_ppm=np.nan,
                             err_ppm=np.nan, quarter=-1, gap_dist_d=float(np.min(np.abs(edges - tc))),
                             in_data=False))
            continue
        base = np.mean(y[oo]) if oo.sum() > 5 else np.nan
        d = base - np.mean(y[inb]) if n_in else np.nan
        s = np.std(y[oo]) if oo.sum() > 5 else np.nan
        rows.append(dict(epoch=n, tc_bkjd=float(tc), n_in=n_in, n_oot_local=int(oo.sum()),
                         depth_ppm=1e6 * float(d) if np.isfinite(d) else np.nan,
                         err_ppm=1e6 * float(s / np.sqrt(max(n_in, 1))) if np.isfinite(s) else np.nan,
                         quarter=int(np.median(qu[loc])) if loc.any() else -1,
                         gap_dist_d=float(np.min(np.abs(edges - tc))),
                         in_data=bool(n_in >= 0.5 * n_expect)))
    return rows, n_expect


def trapezoid(phd, depth, T14, q, dt):
    T23 = q * T14
    x = np.abs(phd - dt)
    m = np.zeros_like(x)
    full = x <= T23 / 2.0
    ing = (x > T23 / 2.0) & (x < T14 / 2.0)
    m[full] = depth
    if T14 > T23:
        m[ing] = depth * (T14 / 2.0 - x[ing]) / ((T14 - T23) / 2.0)
    return -m


def fit_trapezoid(phd, y, w, depth0, dur0):
    sel = np.abs(phd) < 1.0
    x, yy, ww = phd[sel], y[sel], w[sel]

    def resid(p):
        return (yy - trapezoid(x, *p)) * ww

    p0 = [max(depth0, 1e-5), dur0, 0.6, 0.0]
    lo = [0.0, 0.02, 0.0, -0.15]
    hi = [0.02, 0.8, 1.0, 0.15]
    r = least_squares(resid, p0, bounds=(lo, hi))
    # box (q=1) fit for comparison of chi2
    rb = least_squares(lambda p: (yy - trapezoid(x, p[0], p[1], 1.0, p[2])) * ww,
                       [max(depth0, 1e-5), dur0, 0.0], bounds=([0, 0.02, -0.15], [0.02, 0.8, 0.15]))
    return r, rb, sel.sum()


def year_of_bkjd(t):
    # BJD 2454832.5 = 2009-01-01 00:00 UTC -> BKJD -0.5 ; 2012 is a leap year
    bounds = [(-0.5, 364.5, 2009), (364.5, 729.5, 2010), (729.5, 1094.5, 2011),
              (1094.5, 1460.5, 2012), (1460.5, 1825.5, 2013)]
    yr = np.zeros(len(t), int)
    for a, b, y in bounds:
        yr[(t >= a) & (t < b)] = y
    return yr


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    return o


# ============================================================================= main
def main():
    os.makedirs(OUT, exist_ok=True)
    R = {"fname": FNAME, "date": "2026-10-17", "script": os.path.abspath(__file__),
         "fits": FITS, "fits_sha256": sha256(FITS),
         "env": {"python": platform.python_version(), "numpy": np.__version__, "astropy": astropy.__version__},
         "method": {
             "cadence_selection": "SAP_QUALITY == 0, finite TIME/FLUX/FLUX_ERR, FLUX>0",
             "normalisation": "divide by per-quarter median",
             "detrend": "per-quarter time-windowed running median, window %.2f d, 2 iterations; points >4 sigma "
                        "above the trend excluded from the trend estimate; dips never excluded" % WINDOW,
             "outlier_clip_for_search": "+4 sigma / -6 sigma (per-quarter robust sigma)",
             "bls": "astropy BoxLeastSquares, objective=likelihood, log-uniform periods %.0f-%.0f d, "
                    "dlnP = min(dur)/(%d*baseline), durations %s d, dy = per-quarter robust scatter"
                    % (PMIN, PMAX, OVERSAMPLE, DURS.tolist()),
             "depth": "out-of-transit mean minus in-box mean on a re-detrend with |phase|<dur masked; empirical "
                      "error = std of same-duration boxes at random phases avoiding phases 0 and 0.5",
         }}

    # ---------------------------------------------------------------- 1. load + detrend (PSF)
    t, f, e, qu = load("FLUX")
    log("loaded PSF FLUX: %d good cadences, %.1f-%.1f BKJD, baseline %.1f d, quarters %s"
        % (len(t), t.min(), t.max(), t.max() - t.min(), np.unique(qu).tolist()))
    flat0, trend0, keep0, sig0 = detrend(t, f, qu)
    R["data"] = {"n_quality0": int(len(t)), "t_min_bkjd": float(t.min()), "t_max_bkjd": float(t.max()),
                 "baseline_d": float(t.max() - t.min()), "quarters": np.unique(qu).tolist(),
                 "n_after_outlier_clip": int(keep0.sum()),
                 "per_quarter_robust_sigma_ppm": {int(k): float(1e6 * np.median(sig0[qu == k])) for k in np.unique(qu)},
                 "global_robust_sigma_ppm": float(1e6 * 1.4826 * np.median(np.abs(flat0[keep0])))}
    log("robust sigma per quarter (ppm):", {k: round(v) for k, v in R["data"]["per_quarter_robust_sigma_ppm"].items()})

    # ---------------------------------------------------------------- 2. blind BLS 1-50 d
    ts, ys, dys = t[keep0], flat0[keep0], sig0[keep0]
    res, periods = blind_bls(ts, ys, dys)
    log("blind BLS done: %d trial periods x %d durations" % (len(periods), len(DURS)))
    sde_raw, sde_norm, pnorm = spectrum_stats(periods, res.power)
    peaks = top_peaks(periods, res.power, n=3)
    for pk in peaks:
        i = pk["index"]
        pk.update(depth_ppm=float(1e6 * res.depth[i]), duration_h=float(24 * res.duration[i]),
                  t0_bkjd=float(res.transit_time[i]), alias_of_ref=alias_label(pk["period_d"]),
                  sde_local=float((pk["power"] - res.power.mean()) / res.power.std()))
    peaks_norm = top_peaks(periods, pnorm, n=3)
    for pk in peaks_norm:
        pk["alias_of_ref"] = alias_label(pk["period_d"])
    best = peaks[0]
    R["blind_bls"] = {"n_periods": int(len(periods)), "pmin_d": PMIN, "pmax_d": PMAX,
                      "sde_raw_spectrum": sde_raw, "sde_trend_normalised_spectrum": sde_norm,
                      "top3_peaks": peaks, "top3_peaks_trend_normalised": peaks_norm,
                      "top_peak_is_ref_period": bool(abs(best["period_d"] / REF_P - 1) < 1e-3),
                      "top_peak_alias_of_ref": best["alias_of_ref"]}
    log("top-3 peaks:", [(round(p["period_d"], 5), round(p["power"], 1), p["alias_of_ref"]) for p in peaks])
    log("SDE raw %.2f, trend-normalised %.2f" % (sde_raw, sde_norm))

    # ---------------------------------------------------------------- inverted-flux null (same search)
    res_inv, _ = blind_bls(ts, -ys, dys)
    sde_inv_raw, sde_inv_norm, _ = spectrum_stats(periods, res_inv.power)
    pk_inv = top_peaks(periods, res_inv.power, n=3)
    for pk in pk_inv:
        i = pk["index"]
        pk.update(depth_ppm=float(1e6 * res_inv.depth[i]), duration_h=float(24 * res_inv.duration[i]))
    R["inverted_null"] = {"sde_raw_spectrum": sde_inv_raw, "sde_trend_normalised_spectrum": sde_inv_norm,
                          "top3_peaks": pk_inv}
    log("inverted null: SDE raw %.2f (norm %.2f), top P %.4f" % (sde_inv_raw, sde_inv_norm, pk_inv[0]["period_d"]))

    # ---------------------------------------------------------------- 3. refine at the top peak
    ref1, _ = refine(ts, ys, dys, best["period_d"])
    P, t0, dur = ref1["period"], ref1["t0"], ref1["duration"]
    log("refined (unmasked detrend): P=%.5f t0=%.4f dur=%.2f h depth=%.0f ppm" % (P, t0, 24 * dur, 1e6 * ref1["depth"]))

    # masked re-detrend so the running median does not eat the transit
    mask = np.abs(fold_phase(t, P, t0) * P) < 1.0 * dur
    flat, trend, keep, sig = detrend(t, f, qu, mask=mask)
    # refine once more on the masked re-detrend
    ref2, _ = refine(t[keep], flat[keep], sig[keep], P)
    P, t0, dur = ref2["period"], ref2["t0"], ref2["duration"]
    mask = np.abs(fold_phase(t, P, t0) * P) < 1.0 * dur
    flat, trend, keep, sig = detrend(t, f, qu, mask=mask)
    tk, yk, qk = t[keep], flat[keep], qu[keep]
    log("refined (masked detrend):   P=%.5f t0=%.4f dur=%.2f h depth=%.0f ppm" % (P, t0, 24 * dur, 1e6 * ref2["depth"]))

    # duration: box-SNR scan at fixed P, t0 (my own statistic), then trapezoid fit
    scan = []
    for d_try in np.arange(0.04, 0.50, 0.01):
        m = measure(tk, yk, P, t0, d_try, ntrial=800)
        scan.append((d_try, m["snr"], m["depth_ppm"]))
    scan = np.array(scan)
    dur_snr = float(scan[np.argmax(scan[:, 1]), 0])
    phd = fold_phase(tk, P, t0) * P
    w = 1.0 / sig[keep]
    tr, trb, n_fit = fit_trapezoid(phd, yk, w, ref2["depth"], dur)
    depth_tr, T14, qtr, dt_tr = tr.x
    T23 = qtr * T14
    chi2_tr = float(np.sum(tr.fun ** 2))
    chi2_box = float(np.sum(trb.fun ** 2))
    t0_fit = t0 + dt_tr
    # move the reference epoch to the transit nearest the centre of the data
    n_mid = np.round((0.5 * (t.min() + t.max()) - t0_fit) / P)
    t0_ref = float(t0_fit + n_mid * P)
    R["ephemeris"] = {"period_d": float(P), "t0_bkjd": t0_ref, "t0_bjd": t0_ref + BKJD0,
                      "bls_duration_d": float(dur), "bls_duration_h": float(24 * dur),
                      "box_snr_scan_best_duration_h": 24 * dur_snr,
                      "trapezoid": {"depth_ppm": float(1e6 * depth_tr), "T14_h": float(24 * T14), "T23_h": float(24 * T23),
                                    "T23_over_T14": float(qtr), "dt0_min": float(1440 * dt_tr), "n_points_fit": int(n_fit),
                                    "chi2": chi2_tr, "chi2_box_model": chi2_box},
                      "duration_h_adopted": float(24 * T14),
                      "bls_refine_depth_ppm": float(1e6 * ref2["depth"]),
                      "bls_refine_depth_err_ppm": float(1e6 * ref2["depth_err"]),
                      "bls_refine_depth_ppm_unmasked_detrend": float(1e6 * ref1["depth"])}
    log("duration: BLS %.2f h, box-SNR scan %.2f h, trapezoid T14 %.2f h (T23/T14=%.2f), dt0 %.1f min"
        % (24 * dur, 24 * dur_snr, 24 * T14, qtr, 1440 * dt_tr))

    # ---------------------------------------------------------------- box depth, odd/even, secondary, per-year
    dur_use = float(dur)              # box duration used for all box statistics
    main_m = measure(tk, yk, P, t0_fit, dur_use)
    main_m_unmasked = measure(t[keep0], flat0[keep0], P, t0_fit, dur_use, ntrial=1500)
    epochs = np.round((tk - t0_fit) / P).astype(int)
    odd_m = measure(tk, yk, P, t0_fit, dur_use, sel=(epochs % 2 == 1))
    even_m = measure(tk, yk, P, t0_fit, dur_use, sel=(epochs % 2 == 0))
    oe_sigma = abs(odd_m["depth_ppm"] - even_m["depth_ppm"]) / np.hypot(odd_m["err_ppm"], even_m["err_ppm"])
    sec_m = measure(tk, yk, P, t0_fit, dur_use, center=0.5)
    # secondary scan: deepest box anywhere in phase 0.2-0.8 (max-of-N, for context)
    sec_scan = []
    for c in np.linspace(0.2, 0.8, 121):
        phc = fold_phase(tk, P, t0_fit, c) * P
        m = np.abs(phc) < dur_use / 2
        if m.sum() >= 3:
            sec_scan.append((c, 1e6 * (np.mean(yk[np.abs(fold_phase(tk, P, t0_fit) * P) > dur_use]) - np.mean(yk[m]))))
    sec_scan = np.array(sec_scan)
    yrs = year_of_bkjd(tk)
    per_year = {}
    for y in sorted(set(yrs.tolist())):
        mm = measure(tk, yk, P, t0_fit, dur_use, sel=(yrs == y), ntrial=2500)
        per_year[int(y)] = mm
    rows, n_expect = per_event_table(tk, yk, P, t0_fit, dur_use, qk)
    n_transits = int(sum(r["in_data"] for r in rows))
    n_any = int(sum(r["n_in"] > 0 for r in rows))
    ev = [r for r in rows if r["in_data"]]
    ev_d = np.array([r["depth_ppm"] for r in ev])
    ev_e = np.array([r["err_ppm"] for r in ev])
    wmean = float(np.sum(ev_d / ev_e ** 2) / np.sum(1 / ev_e ** 2))
    chi2_ev = float(np.sum(((ev_d - wmean) / ev_e) ** 2))
    R["box_statistics_psf"] = {
        "box_duration_h_used": 24 * dur_use,
        "depth": main_m, "depth_on_unmasked_detrend": main_m_unmasked,
        "odd": odd_m, "even": even_m, "odd_even_sigma": float(oe_sigma),
        "secondary_phase_0p5": sec_m, "secondary_sigma": float(sec_m["snr"]),
        "deepest_box_phase_0p2_0p8": {"phase": float(sec_scan[np.argmax(sec_scan[:, 1]), 0]),
                                      "depth_ppm": float(sec_scan[:, 1].max())},
        "per_year": per_year,
        "events": {"n_epochs_in_baseline": len(rows), "n_transits_with_>=50pct_coverage": n_transits,
                   "n_epochs_with_any_in_transit_point": n_any, "expected_points_per_transit": float(n_expect),
                   "single_event_weighted_mean_ppm": wmean, "chi2_events_about_mean": chi2_ev, "dof": int(len(ev) - 1),
                   "max_single_event_depth_ppm": float(np.nanmax(ev_d)), "max_single_event_snr": float(np.nanmax(ev_d / ev_e)),
                   "fraction_of_total_depth_snr2_from_strongest_event": float(np.nanmax((ev_d / ev_e) ** 2) / np.nansum((ev_d / ev_e) ** 2))},
    }
    log("depth %.0f +- %.0f ppm (SNR %.1f, n_in %d); unmasked-detrend depth %.0f ppm"
        % (main_m["depth_ppm"], main_m["err_ppm"], main_m["snr"], main_m["n_in"], main_m_unmasked["depth_ppm"]))
    log("odd %.0f+-%.0f, even %.0f+-%.0f -> %.2f sigma; secondary %.0f+-%.0f ppm (%.2f sigma)"
        % (odd_m["depth_ppm"], odd_m["err_ppm"], even_m["depth_ppm"], even_m["err_ppm"], oe_sigma,
           sec_m["depth_ppm"], sec_m["err_ppm"], sec_m["snr"]))
    log("per-year:", {y: "%.0f+-%.0f" % (m["depth_ppm"], m["err_ppm"]) for y, m in per_year.items()})
    log("transits with >=50%% coverage: %d of %d epochs in baseline; chi2 of events about mean %.1f / %d"
        % (n_transits, len(rows), chi2_ev, len(ev) - 1))

    with open(os.path.join(OUT, "events.csv"), "w") as fh:
        fh.write("epoch,tc_bkjd,tc_bjd,quarter,n_in,n_oot_local,depth_ppm,err_ppm,gap_dist_d,counted\n")
        for r in rows:
            fh.write("%d,%.5f,%.5f,%d,%d,%d,%s,%s,%.3f,%d\n" % (
                r["epoch"], r["tc_bkjd"], r["tc_bkjd"] + BKJD0, r["quarter"], r["n_in"], r["n_oot_local"],
                "%.1f" % r["depth_ppm"] if np.isfinite(r["depth_ppm"]) else "",
                "%.1f" % r["err_ppm"] if np.isfinite(r["err_ppm"]) else "", r["gap_dist_d"], int(r["in_data"])))

    # ---------------------------------------------------------------- 4. SAP_FLUX depth (dilution check)
    ts_, fs_, es_, qus_ = load("SAP_FLUX")
    mask_s = np.abs(fold_phase(ts_, P, t0_fit) * P) < 1.0 * dur_use
    flat_s, trend_s, keep_s, sig_s = detrend(ts_, fs_, qus_, mask=mask_s)
    sap_m = measure(ts_[keep_s], flat_s[keep_s], P, t0_fit, dur_use)
    sap_odd = measure(ts_[keep_s], flat_s[keep_s], P, t0_fit, dur_use,
                      sel=(np.round((ts_[keep_s] - t0_fit) / P).astype(int) % 2 == 1), ntrial=1500)
    sap_even = measure(ts_[keep_s], flat_s[keep_s], P, t0_fit, dur_use,
                       sel=(np.round((ts_[keep_s] - t0_fit) / P).astype(int) % 2 == 0), ntrial=1500)
    sap_sec = measure(ts_[keep_s], flat_s[keep_s], P, t0_fit, dur_use, center=0.5, ntrial=1500)
    ratio = sap_m["depth_ppm"] / main_m["depth_ppm"]
    ratio_err = ratio * np.hypot(sap_m["err_ppm"] / sap_m["depth_ppm"], main_m["err_ppm"] / main_m["depth_ppm"])
    R["sap_flux"] = {"n_quality0": int(len(ts_)), "n_after_clip": int(keep_s.sum()),
                     "global_robust_sigma_ppm": float(1e6 * 1.4826 * np.median(np.abs(flat_s[keep_s]))),
                     "depth": sap_m, "odd": sap_odd, "even": sap_even, "secondary_phase_0p5": sap_sec,
                     "sap_over_psf_depth_ratio": float(ratio), "sap_over_psf_depth_ratio_err": float(ratio_err),
                     "target_csv_crowdsap": 0.9389557537, "target_csv_psffrac": 0.5006045397603304,
                     "median_sap_over_psf_flux_level": float(np.median(fits.getdata(FITS, "LIGHTCURVE_STITCHED")["SAP_FLUX"]) /
                                                            np.median(fits.getdata(FITS, "LIGHTCURVE_STITCHED")["FLUX"]))}
    log("SAP depth %.0f +- %.0f ppm (SNR %.1f); SAP/PSF depth ratio %.2f +- %.2f"
        % (sap_m["depth_ppm"], sap_m["err_ppm"], sap_m["snr"], ratio, ratio_err))

    # ---------------------------------------------------------------- summary block
    R["summary"] = {
        "period_d": float(P), "t0_bjd": t0_ref + BKJD0, "duration_h": float(24 * T14),
        "depth_ppm": main_m["depth_ppm"], "depth_err_ppm": main_m["err_ppm"], "empirical_snr": main_m["snr"],
        "sde": sde_raw, "n_transits": n_transits, "odd_even_sigma": float(oe_sigma),
        "secondary_ppm": sec_m["depth_ppm"], "secondary_sigma": float(sec_m["snr"]),
        "per_year_depths_ppm": {str(y): m["depth_ppm"] for y, m in per_year.items()},
        "recovered_blind": R["blind_bls"]["top_peak_is_ref_period"],
        "sap_depth_ppm": sap_m["depth_ppm"], "sap_over_psf_depth_ratio": float(ratio),
    }
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(jsonable(R), fh, indent=2)
    np.savez_compressed(os.path.join(OUT, "bls_periodogram.npz"), period=periods, power=res.power,
                        power_inverted=res_inv.power, depth=res.depth, duration=res.duration)
    log("wrote results.json, events.csv, bls_periodogram.npz")

    # ---------------------------------------------------------------- figure
    make_figure(periods, res.power, res_inv.power, peaks, tk, yk, P, t0_fit, dur_use, tr.x, odd_m, even_m, sec_m,
                per_year, rows, main_m, sde_raw, oe_sigma, ts_[keep_s], flat_s[keep_s], sap_m)
    log("wrote fold.png")


def binned(x, y, width):
    edges = np.arange(x.min(), x.max() + width, width)
    idx = np.digitize(x, edges) - 1
    xb, yb, eb = [], [], []
    for i in range(len(edges) - 1):
        m = idx == i
        if m.sum() >= 3:
            xb.append(np.mean(x[m])); yb.append(np.mean(y[m])); eb.append(np.std(y[m]) / np.sqrt(m.sum()))
    return np.array(xb), np.array(yb), np.array(eb)


def make_figure(periods, power, power_inv, peaks, t, y, P, t0, dur, trp, odd_m, even_m, sec_m, per_year, rows,
                main_m, sde, oe_sigma, t_sap, y_sap, sap_m):
    BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
    INK, MUTED, GRID, SURF = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.facecolor": SURF, "figure.facecolor": SURF, "legend.frameon": False})
    fig, axs = plt.subplots(3, 2, figsize=(12, 12.5))
    fig.suptitle("%s independent re-derivation (own loader / 1-d running median / own BLS grid), 2026-10-17" % FNAME,
                 color=INK, fontsize=11)

    ax = axs[0, 0]
    ax.plot(periods, power_inv, color=ORANGE, lw=0.5, alpha=0.6, label="inverted flux (null)")
    ax.plot(periods, power, color=BLUE, lw=0.6, label="PSF flux")
    for k, pk in enumerate(peaks):
        ax.annotate("%d: %.4f d" % (k + 1, pk["period_d"]), (pk["period_d"], pk["power"]), textcoords="offset points",
                    xytext=(4, 2 + 9 * (2 - k)), color=INK, fontsize=8)
    ax.axvline(REF_P, color=MUTED, ls=":", lw=0.8)
    ax.set_xscale("log"); ax.set_xlabel("period (d)"); ax.set_ylabel("BLS log-likelihood power")
    ax.set_title("blind BLS 1-50 d, raw-spectrum SDE %.1f (dotted: pipeline 35.5873 d)" % sde, fontsize=9, color=INK)
    ax.legend(loc="upper left", fontsize=8)

    ax = axs[0, 1]
    phd = fold_phase(t, P, t0) * P * 24
    m = np.abs(phd) < 24
    ax.plot(phd[m], 1e6 * y[m], ".", ms=2, color=MUTED, alpha=0.35, label="cadences")
    xb, yb, eb = binned(phd[m], 1e6 * y[m], 0.5)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3.5, color=BLUE, lw=1, label="30-min bins")
    xx = np.linspace(-24, 24, 2000)
    ax.plot(xx, 1e6 * trapezoid(xx / 24, trp[0], trp[1], trp[2], trp[3]), color=ORANGE, lw=1.6, label="trapezoid fit")
    ax.set_xlim(-24, 24); ax.set_ylim(-3000, 2000)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("relative flux (ppm)")
    ax.set_title("fold at P=%.5f d: box depth %.0f +- %.0f ppm (empirical SNR %.1f)"
                 % (P, main_m["depth_ppm"], main_m["err_ppm"], main_m["snr"]), fontsize=9, color=INK)
    ax.legend(loc="lower left", fontsize=8)

    ax = axs[1, 0]
    ep = np.round((t - t0) / P).astype(int)
    for sel, col, lab, mm in [(ep % 2 == 1, BLUE, "odd", odd_m), (ep % 2 == 0, ORANGE, "even", even_m)]:
        mm2 = sel & m
        xb, yb, eb = binned(phd[mm2], 1e6 * y[mm2], 1.0)
        ax.errorbar(xb, yb, eb, fmt="o", ms=3.5, color=col, lw=1, label="%s: %.0f +- %.0f ppm" % (lab, mm["depth_ppm"], mm["err_ppm"]))
    ax.set_xlim(-24, 24); ax.set_ylim(-2500, 1500)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("relative flux (ppm), 1-h bins")
    ax.set_title("odd vs even transits: %.2f sigma difference" % oe_sigma, fontsize=9, color=INK)
    ax.legend(loc="lower left", fontsize=8)

    ax = axs[1, 1]
    phs = fold_phase(t, P, t0, 0.5) * P * 24
    ms = np.abs(phs) < 24
    ax.plot(phs[ms], 1e6 * y[ms], ".", ms=2, color=MUTED, alpha=0.35)
    xb, yb, eb = binned(phs[ms], 1e6 * y[ms], 1.0)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3.5, color=BLUE, lw=1, label="1-h bins")
    ax.axvspan(-12 * dur, 12 * dur, color=GRID, alpha=0.6, lw=0, label="box used")
    ax.set_xlim(-24, 24); ax.set_ylim(-2500, 1500)
    ax.set_xlabel("hours from phase 0.5"); ax.set_ylabel("relative flux (ppm)")
    ax.set_title("secondary at phase 0.5: %.0f +- %.0f ppm (%.2f sigma)" % (sec_m["depth_ppm"], sec_m["err_ppm"], sec_m["snr"]),
                 fontsize=9, color=INK)
    ax.legend(loc="lower left", fontsize=8)

    ax = axs[2, 0]
    yrs = sorted(per_year)
    ax.axhline(main_m["depth_ppm"], color=MUTED, ls="--", lw=0.8)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.errorbar(yrs, [per_year[k]["depth_ppm"] for k in yrs], [per_year[k]["err_ppm"] for k in yrs], fmt="o", ms=5,
                color=BLUE, lw=1.2, label="PSF flux, per year")
    ax.errorbar([max(yrs) + 0.6], [sap_m["depth_ppm"]], [sap_m["err_ppm"]], fmt="s", ms=5, color=ORANGE, lw=1.2,
                label="SAP flux, all years")
    ax.set_xticks(yrs + [max(yrs) + 0.6]); ax.set_xticklabels([str(k) for k in yrs] + ["SAP"])
    ax.set_ylabel("box depth (ppm)"); ax.set_title("per-year depths (dashed: global PSF depth)", fontsize=9, color=INK)
    ax.legend(loc="upper left", fontsize=8)

    ax = axs[2, 1]
    ev = [r for r in rows if r["n_in"] > 0 and np.isfinite(r["depth_ppm"])]
    tc = np.array([r["tc_bkjd"] for r in ev]); d = np.array([r["depth_ppm"] for r in ev]); e = np.array([r["err_ppm"] for r in ev])
    full = np.array([r["in_data"] for r in ev])
    ax.axhline(main_m["depth_ppm"], color=MUTED, ls="--", lw=0.8); ax.axhline(0, color=MUTED, lw=0.6)
    ax.errorbar(tc[full], d[full], e[full], fmt="o", ms=4, color=BLUE, lw=1, label=">=50%% coverage (%d)" % full.sum())
    if (~full).any():
        ax.errorbar(tc[~full], d[~full], e[~full], fmt="o", ms=4, mfc="none", color=ORANGE, lw=1, label="partial (%d)" % (~full).sum())
    ax.set_xlabel("mid-transit time (BKJD)"); ax.set_ylabel("local box depth (ppm)")
    ax.set_title("per-event local depths (+-1.5 d baseline)", fontsize=9, color=INK)
    ax.legend(loc="upper left", fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(OUT, "fold.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
