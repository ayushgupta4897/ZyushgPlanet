#!/usr/bin/env python
"""Supplementary checks for the kic-009655224 re-derivation (same from-scratch code, no pipeline imports):
 - depth robustness vs detrend window (0.5, 1.0, 2.0 d) and vs excluding the first 0.5 d after each gap
 - period uncertainty from the fine-BLS dlogL profile (drop of 0.5)
 - per-year depth chi2 against a constant
 - raw per-quarter SAP_FLUX / FLUX medians from LIGHTCURVE_Qn (is stitched SAP_FLUX rescaled?)
 - per-quarter SAP/PSF depth ratios vs per-quarter raw-flux expectation
"""
import sys, json, os, numpy as np
from astropy.io import fits
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rederive_kic009655224 as R

OUT = R.OUT
res = json.load(open(os.path.join(OUT, "rederive_results.json")))
P = res["ephemeris"]["period_d"]; t0 = res["ephemeris"]["t0_bjd"] - 2454833.0; DUR = res["ephemeris"]["duration_h"] / 24
S = {}

D = R.load(R.FITS)
t_all = D["t"]
# ---- robustness vs detrend window
rob = {}
for win in (0.5, 1.0, 2.0):
    flat, keep, mad = R.detrend(t_all, D["fn"], window=win)
    t, y = t_all[keep], flat[keep]
    d, nin = R.box_depth(t, y, P, t0, DUR)
    rp = R.random_phase_depths(t, y, P, t0, DUR, n=800)
    rob[f"window_{win}d"] = dict(depth_ppm=float(d), err_ppm=float(np.std(rp)), n_in=int(nin), mad_ppm=float(mad * 1e6))
    print(f"window {win} d: depth {d:.0f} +/- {np.std(rp):.0f} ppm")
# exclude 0.5 d after every gap > 0.5 d (thermal settling)
flat, keep, mad = R.detrend(t_all, D["fn"], window=1.0)
t, y = t_all[keep], flat[keep]
gaps = np.where(np.diff(t) > 0.5)[0]
starts = np.concatenate([[t[0]], t[gaps + 1]])
bad = np.zeros_like(t, bool)
for s in starts:
    bad |= (t >= s) & (t < s + 0.5)
d, nin = R.box_depth(t[~bad], y[~bad], P, t0, DUR)
rp = R.random_phase_depths(t[~bad], y[~bad], P, t0, DUR, n=800)
rob["exclude_0p5d_after_gaps"] = dict(depth_ppm=float(d), err_ppm=float(np.std(rp)), n_in=int(nin), n_removed=int(bad.sum()))
print(f"excluding 0.5 d after gaps ({bad.sum()} pts): depth {d:.0f} +/- {np.std(rp):.0f}")
S["robustness"] = rob

# ---- period uncertainty from fine BLS (dlogL drop of 0.5 ~ 1 sigma for a likelihood)
from astropy.timeseries import BoxLeastSquares
dy = np.full_like(y, mad)
pgrid = np.linspace(P * (1 - 0.001), P * (1 + 0.001), 2001)
bls = BoxLeastSquares(t, y, dy)
r = bls.power(pgrid, [DUR], objective="likelihood", oversample=20)
i = int(np.argmax(r.power)); pk = r.power[i]
ok = r.power >= pk - 0.5
P_lo, P_hi = pgrid[ok].min(), pgrid[ok].max()
S["period"] = dict(best_d=float(pgrid[i]), err_lo_d=float(pgrid[i] - P_lo), err_hi_d=float(P_hi - pgrid[i]),
                   pipeline_P_d=R.PIPE_P, pipeline_offset_sigma=float(abs(pgrid[i] - R.PIPE_P) / (0.5 * (P_hi - P_lo))))
print("period", S["period"])
# t0 uncertainty: scan t0 at fixed P
t0grid = t0 + np.linspace(-0.1, 0.1, 801)
ll = []
for tt0 in t0grid:
    ph = R.phase_of(t, P, tt0); inb = np.abs(ph * P) < DUR / 2
    d = np.median(y[~inb]) - np.mean(y[inb]); ll.append(0.5 * d ** 2 * inb.sum() / mad ** 2)
ll = np.array(ll); j = int(np.argmax(ll)); ok = ll >= ll[j] - 0.5
S["t0"] = dict(best_bkjd=float(t0grid[j]), err_d=float(0.5 * (t0grid[ok].max() - t0grid[ok].min())))
print("t0", S["t0"])

# ---- per-year chi2
py = res["per_year"]; dv = np.array([v["depth_ppm"] for v in py.values()]); ev = np.array([v["err_ppm"] for v in py.values()])
w = 1 / ev ** 2; m = np.sum(w * dv) / np.sum(w); chi2 = float(np.sum((dv - m) ** 2 * w))
S["per_year_consistency"] = dict(weighted_mean_ppm=float(m), chi2=chi2, dof=int(len(dv) - 1),
                                 max_pair_sigma=float(max(abs(dv[a] - dv[b]) / np.hypot(ev[a], ev[b]) for a in range(len(dv)) for b in range(a + 1, len(dv)))))
print("per-year", S["per_year_consistency"])
pq = res["per_quarter"]; dv = np.array([v["depth_ppm"] for v in pq.values()]); ev = np.array([v["err_ppm"] for v in pq.values()])
w = 1 / ev ** 2; m = np.sum(w * dv) / np.sum(w)
S["per_quarter_consistency"] = dict(weighted_mean_ppm=float(m), chi2=float(np.sum((dv - m) ** 2 * w)), dof=int(len(dv) - 1))
print("per-quarter", S["per_quarter_consistency"])

# ---- raw per-quarter fluxes from LIGHTCURVE_Qn
raw = {}
with fits.open(R.FITS) as f:
    for i in range(2, len(f), 2):
        h = f[i].header; d = f[i].data; qq = int(h["QUARTER"])
        cols = d.columns.names
        good = (d["SAP_QUALITY"] == 0) if "SAP_QUALITY" in cols else np.ones(len(d), bool)
        ent = dict(columns=cols, n=int(good.sum()))
        for c in ("FLUX", "SAP_FLUX", "PSF_FLAT_FLUX"):
            if c in cols:
                ent["median_" + c] = float(np.nanmedian(d[c][good]))
        for k in ("CROWDSAP", "FLFRCSAP", "PSFFRAC", "NPIXSAP"):
            ent[k] = h.get(k)
        raw[qq] = ent
        print(qq, {k: v for k, v in ent.items() if k != "columns"})
S["raw_per_quarter"] = raw
# stitched medians per quarter (before my normalisation)
st = {}
for qq in np.unique(D["q"]):
    s = D["q"] == qq
    st[int(qq)] = dict(median_FLUX=float(np.median(D["f"][s])), median_SAP=float(np.median(D["sap"][s])))
S["stitched_per_quarter_medians"] = st

# ---- per-quarter SAP/PSF depth ratio vs raw-flux expectation
flat_s, keep_s, mad_s = R.detrend(t_all, D["sapn"], window=1.0)
k = keep & keep_s
tq, yq, ysq, qq_ = t_all[k], R.detrend(t_all, D["fn"], window=1.0)[0][k], flat_s[k], D["q"][k]
pqr = {}
for qq in np.unique(qq_):
    s = qq_ == qq
    dp, n1 = R.box_depth(tq[s], yq[s], P, t0, DUR); ds, n2 = R.box_depth(tq[s], ysq[s], P, t0, DUR)
    rpp = R.random_phase_depths(tq[s], yq[s], P, t0, DUR, n=300); rps = R.random_phase_depths(tq[s], ysq[s], P, t0, DUR, n=300)
    exp = raw[int(qq)].get("median_FLUX", np.nan) / raw[int(qq)].get("median_SAP_FLUX", np.nan) if "median_SAP_FLUX" in raw[int(qq)] else np.nan
    pqr[int(qq)] = dict(psf_depth_ppm=float(dp), psf_err=float(np.std(rpp)), sap_depth_ppm=float(ds), sap_err=float(np.std(rps)),
                        ratio=float(ds / dp) if dp else np.nan, expected_ratio_raw_flux=float(exp))
    print(f"Q{qq}: PSF {dp:.0f}+/-{np.std(rpp):.0f}  SAP {ds:.0f}+/-{np.std(rps):.0f}  ratio {ds/dp:.2f}  expected(raw FLUX/SAP) {exp:.2f}")
S["per_quarter_sap_psf"] = pqr
# season-grouped (same channel) SAP vs PSF
seasons = {}
for qq, ent in raw.items():
    seasons.setdefault(ent.get("NPIXSAP"), [])
json.dump(S, open(os.path.join(OUT, "rederive_supplement.json"), "w"), indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
print("wrote supplement")
