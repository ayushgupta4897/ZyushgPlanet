#!/usr/bin/env python
"""
Addendum to the independent re-derivation of kic-009520080 (own code, no
pipeline/ imports).  Uses the ephemeris from rederive_results.json and asks
three follow-up questions raised by the comparison with the pipeline:

  A. Is the ~ -190 ppm "secondary" at phase 0.5 part of a smooth orbital
     phase curve (reflection / ellipsoidal / Doppler)?  Fit
     1 + a1 cos(2 pi phi) + b1 sin(2 pi phi) + a2 cos(4 pi phi) + b2 sin(4 pi phi)
     to the out-of-transit flux, and give binned phase-0.5 depth per year.
  B. Per-quarter depths with empirical (random-phase) errors, chi2 against
     the weighted mean, and a Kepler-season (roll-angle) breakdown.
  C. Depth with inverse-variance weights vs plain mean (to explain the
     1490 vs 1568 ppm difference between the two estimators).
"""
import json
import numpy as np
from astropy.io import fits

FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-009520080_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-009520080/rederive"
RNG = np.random.default_rng(20261017)
R = json.load(open(f"{OUT}/rederive_results.json"))
P = R["ephemeris"]["period_d"]; T0 = R["ephemeris"]["t0_bkjd"]; DUR = R["ephemeris"]["duration_bls_d"]


def load(col, errcol):
    with fits.open(FITS) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float); f = np.asarray(d[col], float)
        e = np.asarray(d[errcol], float); q = np.asarray(d["SAP_QUALITY"], int)
        Q = np.asarray(d["QUARTER"], int)
    m = (q == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (f > 0)
    o = np.argsort(t[m])
    return t[m][o], f[m][o], e[m][o], Q[m][o]


def running_median(t, f, good, window=1.0):
    tg, fg = t[good], f[good]
    lo = np.searchsorted(tg, t - window / 2); hi = np.searchsorted(tg, t + window / 2, side="right")
    out = np.empty_like(f)
    for i in range(len(t)):
        a, b = lo[i], hi[i]
        if b - a < 5:
            a = max(0, a - 5); b = min(len(tg), b + 5)
        out[i] = np.median(fg[a:b])
    return out


def detrend(t, f, window=1.0, niter=3):
    good = np.ones(len(t), bool); trend = running_median(t, f, good, window)
    for _ in range(niter):
        r = f / trend - 1; sig = 1.4826 * np.median(np.abs(r - np.median(r)))
        good = (r > -3 * sig) & (r < 3 * sig); trend = running_median(t, f, good, window)
    flat = f / trend; r = flat - 1; sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    return flat, (r < 4 * sig) & (r > -15 * sig), sig


t, f, e, Q = load("FLUX", "FLUX_ERR")
for qq in np.unique(Q):
    s = Q == qq; med = np.nanmedian(f[s]); f[s] /= med; e[s] /= med
flat, keep, sig = detrend(t, f)
t, f, e, Q = t[keep], flat[keep], e[keep], Q[keep]
ph = ((t - T0) / P + 0.5) % 1.0 - 0.5
w = DUR / P / 2
intr = np.abs(ph) < 1.5 * w          # generous transit mask
res = {"fname": "kic-009520080", "ephemeris_used": dict(P=P, t0_bkjd=T0, dur_d=DUR)}

# ---------- A. phase curve ----------
x = ph[~intr]; y = f[~intr] - 1.0
A = np.vstack([np.ones_like(x), np.cos(2 * np.pi * x), np.sin(2 * np.pi * x),
               np.cos(4 * np.pi * x), np.sin(4 * np.pi * x)]).T
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
resid = y - A @ coef
s2 = np.sum(resid ** 2) / (len(y) - 5)
cov = s2 * np.linalg.inv(A.T @ A)
err = np.sqrt(np.diag(cov))
names = ["const", "cos1", "sin1", "cos2", "sin2"]
pc = {n: dict(ppm=float(c * 1e6), err_ppm=float(er * 1e6)) for n, c, er in zip(names, coef, err)}
# model prediction at phase 0.5 relative to phase ~0 baseline: cos1 term flips sign
# note: phase 0 = transit, so a reflection hump would be +cos(2 pi phi) with NEGATIVE
# coefficient (brightest at phi=0.5).  Ellipsoidal = cos(4 pi phi) with negative coefficient.
res["phase_curve_fit"] = pc
# 20-bin phase curve, transit bins excluded in the fit but shown
nb = 20
edges = np.linspace(-0.5, 0.5, nb + 1)
idx = np.digitize(ph, edges) - 1
binned = []
for i in range(nb):
    s = idx == i
    binned.append(dict(phase=float(0.5 * (edges[i] + edges[i + 1])), mean_ppm=float((f[s].mean() - 1) * 1e6),
                       err_ppm=float(f[s].std() / np.sqrt(s.sum()) * 1e6), n=int(s.sum())))
res["phase_curve_20bins"] = binned

# phase-0.5 box depth per year with empirical errors
def box_depth(ph, f, centre, nsamp=1500):
    d = (ph - centre + 0.5) % 1.0 - 0.5
    s = np.abs(d) < w
    out = np.abs(d) > 3 * w
    base = np.median(f[out]); m_in = f[s].mean()
    vals = []
    while len(vals) < nsamp:
        c = RNG.uniform(-0.5, 0.5)
        if min(abs(((c + 0.5) % 1) - 0.5), abs(((c - 0.5 + 0.5) % 1) - 0.5)) < 4 * w:
            continue
        dd = (ph - c + 0.5) % 1.0 - 0.5; ss = np.abs(dd) < w
        if ss.sum() > 0:
            vals.append(f[ss].mean())
    err = np.std(vals)
    return dict(depth_ppm=float((base - m_in) * 1e6), err_ppm=float(err * 1e6),
                sigma=float((base - m_in) / err), n_in=int(s.sum()))

years = {"Y1_Q1-4": (1, 4), "Y2_Q5-8": (5, 8), "Y3_Q9-12": (9, 12), "Y4_Q13-17": (13, 17)}
res["secondary_phase05_per_year"] = {}
for k, (a, b) in years.items():
    s = (Q >= a) & (Q <= b)
    res["secondary_phase05_per_year"][k] = box_depth(ph[s], f[s], 0.5)
# secondary scanned over a +-0.1 phase window: where is the extreme?
scan = []
for c in np.linspace(0.4, 0.6, 41):
    d = (ph - c + 0.5) % 1.0 - 0.5; s = np.abs(d) < w
    scan.append((float(c), float((np.median(f[np.abs(ph) > 3 * w]) - f[s].mean()) * 1e6)))
res["secondary_scan_phase_0p4_0p6"] = scan

# ---------- B. per-quarter depths with empirical errors + seasons ----------
pq = {}
for qq in np.unique(Q):
    s = Q == qq
    r = box_depth(ph[s], f[s], 0.0, nsamp=800)
    r["season"] = int((qq - 1) % 4) if qq >= 1 else -1   # Kepler roll season 0..3 for Q1..Q17 (Q0 excluded)
    pq[int(qq)] = r
d = np.array([v["depth_ppm"] for v in pq.values()]); er = np.array([v["err_ppm"] for v in pq.values()])
wm = np.sum(d / er ** 2) / np.sum(1 / er ** 2)
chi2 = float(np.sum(((d - wm) / er) ** 2))
res["per_quarter"] = dict(depths=pq, weighted_mean_ppm=float(wm), chi2=chi2, dof=int(len(d) - 1),
                          chi2_per_dof=chi2 / (len(d) - 1))
seasons = {}
for sn in range(4):
    qs = [q for q, v in pq.items() if v["season"] == sn]
    s = np.isin(Q, qs)
    seasons[f"season{sn}_Q{qs}"] = box_depth(ph[s], f[s], 0.0, nsamp=800)
ds = np.array([v["depth_ppm"] for v in seasons.values()]); es = np.array([v["err_ppm"] for v in seasons.values()])
wms = np.sum(ds / es ** 2) / np.sum(1 / es ** 2)
res["per_season"] = dict(depths=seasons, weighted_mean_ppm=float(wms),
                         chi2=float(np.sum(((ds - wms) / es) ** 2)), dof=3)

# ---------- C. weighted vs unweighted depth ----------
s_in = np.abs(ph) < w
out = np.abs(ph) > 3 * w
base_med = np.median(f[out]); base_mean = f[out].mean()
wgt = 1 / e ** 2
m_in_w = np.sum(f[s_in] * wgt[s_in]) / np.sum(wgt[s_in])
base_w = np.sum(f[out] * wgt[out]) / np.sum(wgt[out])
res["depth_estimators_ppm"] = dict(
    plain_mean_vs_median_baseline=float((base_med - f[s_in].mean()) * 1e6),
    plain_mean_vs_mean_baseline=float((base_mean - f[s_in].mean()) * 1e6),
    ivar_weighted=float((base_w - m_in_w) * 1e6),
    note="FLUX_ERR-weighted depth is what astropy compute_stats / a weighted BLS return; "
         "the plain mean is what the random-phase box uses")
# per-point error vs scatter: are FLUX_ERR realistic?
res["flux_err_median_ppm"] = float(np.median(e) * 1e6)
res["robust_scatter_ppm"] = float(sig * 1e6)

json.dump(res, open(f"{OUT}/rederive_addendum.json", "w"), indent=1)
print(json.dumps({k: res[k] for k in ["phase_curve_fit", "secondary_phase05_per_year", "per_season",
                                      "depth_estimators_ppm", "flux_err_median_ppm", "robust_scatter_ppm"]}, indent=1))
print("per-quarter chi2 = %.1f / %d dof" % (chi2, len(d) - 1))
for q, v in pq.items():
    print("Q%2d season %d  %5.0f +- %4.0f ppm (n_in=%d)" % (q, v["season"], v["depth_ppm"], v["err_ppm"], v["n_in"]))
print("secondary scan (phase, depth ppm) min/max:", min(scan, key=lambda z: z[1]), max(scan, key=lambda z: z[1]))
print("20-bin phase curve (ppm):", [round(b["mean_ppm"]) for b in binned])
