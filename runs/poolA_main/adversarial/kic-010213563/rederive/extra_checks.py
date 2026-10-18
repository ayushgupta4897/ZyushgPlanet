#!/usr/bin/env python
"""Supplementary checks for the kic-010213563 re-derivation (own code, no pipeline imports):
 A. inverted-flux peak at the same period: phase offset and dependence on transit masking
 B. trapezoid fit (T14, T23, depth) and V-shape comparison
 C. per-event depth chi2
Writes extra.json next to results.json."""
import json, os
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.optimize import least_squares

FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-010213563_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-010213563/rederive"
WIN_D, GAP_D, CADENCE = 1.0, 0.5, 0.0204340
R = json.load(open(os.path.join(OUT, "results.json")))
P, T0, DUR = R["ephemeris"]["period_d"], R["ephemeris"]["t0_bkjd"], R["ephemeris"]["duration_h_bls"] / 24.0


def load(column):
    with fits.open(FITS) as f:
        d = f["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float); y = np.asarray(d[column], float)
        e = np.asarray(d[column + "_ERR"], float); q = np.asarray(d["QUARTER"], int)
        qual = np.asarray(d["SAP_QUALITY"], int)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (y > 0)
    o = np.argsort(t[m])
    return t[m][o], y[m][o], e[m][o], q[m][o]


def normalise(t, y, e, q):
    yn, en = np.empty_like(y), np.empty_like(e)
    for qq in np.unique(q):
        s = q == qq; med = np.median(y[s]); yn[s] = y[s] / med; en[s] = e[s] / med
    return yn, en


def running_median(t, y, win, mask, min_pts=8):
    tm, ym = t[mask], y[mask]
    lo = np.searchsorted(tm, t - win / 2); hi = np.searchsorted(tm, t + win / 2)
    tr = np.full_like(y, np.nan)
    for i in range(len(t)):
        if hi[i] - lo[i] >= min_pts:
            tr[i] = np.median(ym[lo[i]:hi[i]])
    return tr


def segments(t, gap):
    return np.split(np.arange(len(t)), np.where(np.diff(t) > gap)[0] + 1)


def detrend(t, y, transit_mask=None, n_iter=3):
    keep = np.ones_like(y, bool)
    use = np.ones_like(y, bool) if transit_mask is None else ~transit_mask
    trend = np.full_like(y, np.nan)
    for seg in segments(t, GAP_D):
        ts, ys = t[seg], y[seg]; good = use[seg].copy()
        for _ in range(n_iter):
            tr = running_median(ts, ys, WIN_D, good); res = ys - tr; ok = np.isfinite(res)
            mad = 1.4826 * np.median(np.abs(res[ok & good] - np.median(res[ok & good])))
            if not np.isfinite(mad) or mad <= 0: break
            ng = use[seg] & ok & (res < 3.5 * mad) & (res > -8 * mad)
            if np.array_equal(ng, good): break
            good = ng
        tr = running_median(ts, ys, WIN_D, good); trend[seg] = tr
        res = ys - tr; ok = np.isfinite(res)
        mad = 1.4826 * np.median(np.abs(res[ok] - np.median(res[ok])))
        keep[seg] = ok & (res < 3.5 * mad) & (res > -8 * mad)
    return y / trend, keep


def phase_days(t, P, t0):
    return ((t - t0 + 0.5 * P) % P) - 0.5 * P


out = {}
t, y, e, q = load("FLUX")
yn, en = normalise(t, y, e, q)
in_tr = np.abs(phase_days(t, P, T0)) < 0.75 * DUR
flat_u, keep_u = detrend(t, yn)                       # unmasked (what the blind search used)
flat_m, keep_m = detrend(t, yn, transit_mask=in_tr)   # transit-masked
mu = keep_u & np.isfinite(flat_u)
mm = (keep_m & np.isfinite(flat_m)) | (in_tr & np.isfinite(flat_m))

# ---------------- A. inverted peak
durs = np.arange(0.04, 0.42, 0.02)
Pgrid = np.linspace(P * 0.997, P * 1.003, 1201)
A = {}
for label, msk, flat in (("unmasked", mu, flat_u), ("masked", mm, flat_m)):
    for sign in ("real", "inverted"):
        f = flat[msk] if sign == "real" else 2.0 - flat[msk]
        r = BoxLeastSquares(t[msk], f, dy=en[msk]).power(Pgrid, durs, objective="likelihood")
        j = int(np.argmax(r.power))
        dphi = phase_days(np.array([float(r.transit_time[j])]), P, T0)[0]
        A[f"{label}_{sign}"] = dict(period_d=float(r.period[j]), power=float(r.power[j]),
                                    depth_ppm=float(r.depth[j] * 1e6), duration_h=float(r.duration[j] * 24),
                                    t0_bkjd=float(r.transit_time[j]), phase_offset_from_transit_h=float(dphi * 24))
        print(f"[A] {label:8s} {sign:8s} P={r.period[j]:.5f} power={r.power[j]:.1f} depth={r.depth[j]*1e6:.0f} ppm "
              f"dur={r.duration[j]*24:.2f} h  offset from transit = {dphi*24:+.2f} h")
# folded profile of the inverted/unmasked flux near transit (shoulders?)
ph = phase_days(t[mu], P, T0) * 24
bins = np.arange(-24, 24.01, 1.0)
idx = np.digitize(ph, bins)
prof = [float(np.mean(flat_u[mu][idx == i + 1]) - 1) * 1e6 for i in range(len(bins) - 1)]
A["unmasked_profile_1h_bins_ppm"] = dict(zip([f"{0.5*(bins[i]+bins[i+1]):+.1f}h" for i in range(len(bins) - 1)], prof))
sh = [(k, v) for k, v in A["unmasked_profile_1h_bins_ppm"].items() if 3 <= abs(float(k[:-1])) <= 12]
print("[A] unmasked folded profile 3-12 h from mid-transit (ppm):", ", ".join(f"{k}:{v:+.0f}" for k, v in sh))
out["inverted_peak"] = A

# ---------------- B. trapezoid fit on masked-detrended flux, |phase| < 1 d
tt, ff, ee = t[mm], flat_m[mm], en[mm]
ph = phase_days(tt, P, T0)
w = np.abs(ph) < 1.0
x, fy, fe = ph[w], ff[w], ee[w]
sig = 1.4826 * np.median(np.abs(fy - np.median(fy)))


def trap(p, x):
    dt0, dep, t14, frac = p
    t23 = t14 * np.clip(frac, 0, 1)
    a = np.abs(x - dt0)
    m = np.ones_like(x)
    half14, half23 = t14 / 2, t23 / 2
    inside = a <= half23
    ramp = (a > half23) & (a < half14)
    m[inside] = 1 - dep
    if half14 > half23:
        m[ramp] = 1 - dep * (half14 - a[ramp]) / (half14 - half23)
    return m


def resid(p, frac_fixed=None):
    pp = list(p) + ([frac_fixed] if frac_fixed is not None else [])
    return (fy - trap(pp, x)) / sig


fit = least_squares(resid, [0.0, 0.0022, 0.25, 0.7], bounds=([-0.1, 0, 0.05, 0], [0.1, 0.02, 0.8, 1]))
fitV = least_squares(lambda p: resid(p, 0.0), [0.0, 0.0022, 0.25], bounds=([-0.1, 0, 0.05], [0.1, 0.02, 0.8]))
fitB = least_squares(lambda p: resid(p, 1.0), [0.0, 0.0022, 0.25], bounds=([-0.1, 0, 0.05], [0.1, 0.02, 0.8]))
chi_t, chi_v, chi_b = [float(np.sum(f.fun ** 2)) for f in (fit, fitV, fitB)]
chi_flat = float(np.sum(((fy - 1) / sig) ** 2))
# crude errors from the Jacobian
try:
    J = fit.jac; cov = np.linalg.inv(J.T @ J) * (chi_t / (len(x) - 4)); perr = np.sqrt(np.diag(cov))
except Exception:
    perr = [np.nan] * 4
dt0, dep, t14, frac = fit.x
B = dict(n_points=int(len(x)), depth_ppm=float(dep * 1e6), depth_err_ppm=float(perr[1] * 1e6),
         T14_h=float(t14 * 24), T14_err_h=float(perr[2] * 24), T23_h=float(t14 * frac * 24), T23_over_T14=float(frac),
         t0_shift_min=float(dt0 * 24 * 60), chi2_trapezoid=chi_t, chi2_pure_V=chi_v, chi2_pure_box=chi_b,
         chi2_flat=chi_flat, delta_chi2_V_minus_trap=chi_v - chi_t, delta_chi2_box_minus_trap=chi_b - chi_t,
         V_depth_ppm=float(fitV.x[1] * 1e6), V_T14_h=float(fitV.x[2] * 24),
         box_depth_ppm=float(fitB.x[1] * 1e6), box_T14_h=float(fitB.x[2] * 24))
print(f"[B] trapezoid: depth={dep*1e6:.0f}+-{perr[1]*1e6:.0f} ppm T14={t14*24:.2f}+-{perr[2]*24:.2f} h "
      f"T23/T14={frac:.2f} (T23={t14*frac*24:.2f} h) t0 shift {dt0*24*60:+.1f} min")
print(f"[B] chi2: flat {chi_flat:.1f}  box {chi_b:.1f}  V {chi_v:.1f}  trapezoid {chi_t:.1f}  (N={len(x)})")
out["trapezoid"] = B

# ---------------- C. per-event chi2
ev = [e for e in R["events"] if e["frac_in"] >= 0.5 and e["depth_ppm"] is not None and e["err_ppm"] is not None]
d = np.array([e["depth_ppm"] for e in ev]); s = np.array([e["err_ppm"] for e in ev])
wm = np.sum(d / s ** 2) / np.sum(1 / s ** 2)
chi2 = float(np.sum(((d - wm) / s) ** 2)); dof = len(d) - 1
z = (d - wm) / s
C = dict(n=len(d), weighted_mean_ppm=float(wm), chi2=chi2, dof=dof, chi2_red=chi2 / dof,
         outliers_gt_2p5sigma=[dict(epoch=e["epoch"], tmid=e["tmid"], depth_ppm=e["depth_ppm"], z=float(zz))
                               for e, zz in zip(ev, z) if abs(zz) > 2.5],
         n_events_consistent_with_zero_depth=int(np.sum(d / s < 2)))
print(f"[C] per-event depths: chi2/dof = {chi2:.1f}/{dof} = {chi2/dof:.2f}; |z|>2.5: "
      + ", ".join(f"epoch {o['epoch']} ({o['depth_ppm']:.0f} ppm, z={o['z']:+.1f})" for o in C["outliers_gt_2p5sigma"]))
out["per_event_chi2"] = C

json.dump(out, open(os.path.join(OUT, "extra.json"), "w"), indent=1)
print("[done] extra.json")
