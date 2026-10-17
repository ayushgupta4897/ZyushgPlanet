#!/usr/bin/env python
"""Dilution test by Kepler season (own code). The KBONUS SAP aperture changes with
season (3 px in Q4/8/12/16 with SAP/PSF flux ~1.0; 5 px in Q3/7/11/15 with SAP/PSF
flux ~2.4). If the dip is on the target, the SAP depth must track the target's flux
fraction in the aperture: SAP/PSF depth ratio ~ FLFRCSAP * F_psf / F_sap per season.
Also appends the pipeline comparison block to results.json."""
import json, os
import numpy as np
from astropy.io import fits

OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-010213563/rederive"
FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-010213563_kepler_v1.0_lc.fits")
R = json.load(open(os.path.join(OUT, "results.json")))
X = json.load(open(os.path.join(OUT, "extra.json")))
B = json.load(open(os.path.join(OUT, "aperture_flux_budget.json")))
P, T0, DUR = R["ephemeris"]["period_d"], R["ephemeris"]["t0_bkjd"], R["ephemeris"]["duration_h_bls"] / 24.0
WIN_D, GAP_D = 1.0, 0.5
RNG = np.random.default_rng(7)


def load(column):
    with fits.open(FITS) as f:
        d = f["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float); y = np.asarray(d[column], float)
        e = np.asarray(d[column + "_ERR"], float); q = np.asarray(d["QUARTER"], int)
        qual = np.asarray(d["SAP_QUALITY"], int)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (y > 0)
    o = np.argsort(t[m]); return t[m][o], y[m][o], e[m][o], q[m][o]


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
        if hi[i] - lo[i] >= min_pts: tr[i] = np.median(ym[lo[i]:hi[i]])
    return tr


def segments(t, gap):
    return np.split(np.arange(len(t)), np.where(np.diff(t) > gap)[0] + 1)


def detrend(t, y, transit_mask, n_iter=3):
    keep = np.ones_like(y, bool); use = ~transit_mask; trend = np.full_like(y, np.nan)
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


def box_depth(t, f, P, t0, dur):
    ph = np.abs(phase_days(t, P, t0)); inm = ph < dur / 2; outm = (ph > dur) & (ph < 3 * dur)
    if inm.sum() < 3 or outm.sum() < 10: return np.nan, int(inm.sum())
    return float(np.mean(f[outm]) - np.mean(f[inm])), int(inm.sum())


def rp_err(t, f, P, t0, dur, n=800):
    vals = []
    while len(vals) < n:
        sh = RNG.uniform(0, P); d1 = min(sh, P - sh); d2 = abs(sh - P / 2)
        if d1 < 2 * dur or d2 < 2 * dur: continue
        v, _ = box_depth(t, f, P, t0 + sh, dur)
        if np.isfinite(v): vals.append(v)
    vals = np.array(vals); return float(1.4826 * np.median(np.abs(vals - np.median(vals))))


data = {}
for col in ("FLUX", "SAP_FLUX"):
    t, y, e, q = load(col); yn, en = normalise(t, y, e, q)
    in_tr = np.abs(phase_days(t, P, T0)) < 0.75 * DUR
    flat, keep = detrend(t, yn, in_tr)
    m = (keep & np.isfinite(flat)) | (in_tr & np.isfinite(flat))
    data[col] = (t[m], flat[m], q[m])

groups = {"A_Q4-8-12-16_3px": [4, 8, 12, 16], "B_Q1-5-9-13-17_4px": [1, 5, 9, 13, 17],
          "C_Q2-6-10-14_4px": [2, 6, 10, 14], "D_Q3-7-11-15_5px": [3, 7, 11, 15]}
bq = {b["q"]: b for b in B}
out = {}
print("season          PSF depth        SAP depth        measured SAP/PSF   expected if on target (FLFRCSAP*Fpsf/Fsap)   expected if CROWDSAP keyword were right")
for name, qs in groups.items():
    row = {}
    for col in ("FLUX", "SAP_FLUX"):
        t, f, q = data[col]; s = np.isin(q, qs)
        d, nin = box_depth(t[s], f[s], P, T0, DUR); er = rp_err(t[s], f[s], P, T0, DUR)
        row[col] = dict(depth_ppm=d * 1e6, err_ppm=er * 1e6, n_in=nin)
    ratio = row["SAP_FLUX"]["depth_ppm"] / row["FLUX"]["depth_ppm"]
    rerr = ratio * np.sqrt((row["SAP_FLUX"]["err_ppm"] / row["SAP_FLUX"]["depth_ppm"]) ** 2 + (row["FLUX"]["err_ppm"] / row["FLUX"]["depth_ppm"]) ** 2)
    exp = float(np.mean([bq[qq]["expected_ratio"] for qq in qs]))
    crowd = float(np.mean([bq[qq]["crowdsap"] for qq in qs]))
    fr = float(np.mean([bq[qq]["ext_sap"] / bq[qq]["ext_psf"] for qq in qs]))
    row.update(measured_ratio=ratio, measured_ratio_err=rerr, expected_ratio_on_target=exp, crowdsap_keyword=crowd,
               sap_over_psf_flux=fr, sigma_from_expected=(ratio - exp) / rerr, sigma_from_crowdsap=(ratio - crowd) / rerr)
    out[name] = row
    print(f"{name:18s} {row['FLUX']['depth_ppm']:6.0f}+-{row['FLUX']['err_ppm']:4.0f}   {row['SAP_FLUX']['depth_ppm']:6.0f}+-{row['SAP_FLUX']['err_ppm']:4.0f}   "
          f"{ratio:.2f}+-{rerr:.2f}            {exp:.2f} ({(ratio-exp)/rerr:+.1f} sigma)                  {crowd:.2f} ({(ratio-crowd)/rerr:+.1f} sigma)")
# weighted linear check: measured ratio vs expected ratio across seasons
xs = np.array([v["expected_ratio_on_target"] for v in out.values()]); ys = np.array([v["measured_ratio"] for v in out.values()])
es = np.array([v["measured_ratio_err"] for v in out.values()])
chi2_track = float(np.sum(((ys - xs) / es) ** 2)); chi2_const = float(np.sum(((ys - np.average(ys, weights=1 / es ** 2)) / es) ** 2))
cs = np.array([v["crowdsap_keyword"] for v in out.values()]); chi2_crowd = float(np.sum(((ys - cs) / es) ** 2))
print(f"chi2 (4 seasons): ratio tracks FLFRCSAP*Fpsf/Fsap: {chi2_track:.1f}; constant ratio: {chi2_const:.1f} (3 dof); ratio = CROWDSAP keyword: {chi2_crowd:.1f}")
out["_summary"] = dict(chi2_tracks_aperture_flux_fraction=chi2_track, chi2_constant_ratio=chi2_const, chi2_equals_crowdsap_keyword=chi2_crowd, dof=4)
json.dump(out, open(os.path.join(OUT, "sap_by_season.json"), "w"), indent=1)

# ---------------- comparison block (pipeline numbers copied from vet.json / survivors.csv, read AFTER my numbers)
R["supplementary"] = {"inverted_peak": X["inverted_peak"], "trapezoid": X["trapezoid"], "per_event_chi2": X["per_event_chi2"],
                      "sap_by_season": out, "aperture_flux_budget_file": "aperture_flux_budget.json"}
R["pipeline_comparison"] = {
    "source": "runs/poolA_main/vet/kic-010213563/vet.json and runs/poolA_main/survivors.csv (read after my numbers were final)",
    "period_d": {"mine_blind_top": R["bls_blind"]["best_period_d"], "mine_refined": P, "pipeline_vet": 36.76708318978774, "pipeline_screen": 36.766301710506575},
    "t0_bkjd": {"mine": T0, "pipeline_vet": 164.85416024494916, "diff_min": (T0 - 164.85416024494916) * 1440},
    "depth_ppm": {"mine_box": R["psf_depth"]["depth_ppm"], "mine_box_err": R["psf_depth"]["err_ppm_random_phase_mad"],
                  "mine_per_event_weighted": R["psf_depth"]["per_event_weighted_mean_ppm"], "mine_trapezoid": X["trapezoid"]["depth_ppm"],
                  "pipeline_pass2": 2241.5372292315938, "pipeline_pass2_formal_err": 89.63193502822614, "pipeline_event_weighted": 2204.6336887359084,
                  "pipeline_screen": 1970.2120190003836},
    "duration_h": {"mine_bls_box": R["ephemeris"]["duration_h_bls"], "mine_T14_trapezoid": X["trapezoid"]["T14_h"], "mine_T23_trapezoid": X["trapezoid"]["T23_h"],
                   "pipeline_vet": 5.625, "pipeline_screen": 4.5, "pipeline_t14_max_h_for_rstar": 5.4267340538361735},
    "sde": {"mine_raw": R["bls_blind"]["sde_raw"], "mine_median_detrended": R["bls_blind"]["sde_detrended"], "pipeline_screen": 73.89163585027678},
    "empirical_snr": {"mine": R["psf_depth"]["empirical_snr"], "pipeline_pass2": 20.600198950797104, "pipeline_screen": 22.042642616725022},
    "odd_even_sigma": {"mine_per_event": R["odd_even"]["sigma_per_event"], "mine_box": R["odd_even"]["sigma_box"], "pipeline": 0.7015036801795245},
    "secondary": {"mine_ppm": R["secondary"]["depth_ppm"], "mine_sigma": R["secondary"]["sigma"], "pipeline_ppm": -171.23356267290114, "pipeline_sigma": -1.5736724833808662},
    "n_transits": {"mine_ge50pct_coverage": R["n_transits_observed"], "mine_any_coverage": sum(1 for e in R["events"] if e["n_in"] > 0), "pipeline": 36},
    "per_year_vs_pipeline_quartiles": {"mine_years": [y["depth_ppm"] for y in R["per_year"]], "pipeline_quartiles": [2542.0, 2005.1, 2128.0, 2229.2]},
    "event_chi2": {"mine": f"{X['per_event_chi2']['chi2']:.1f}/{X['per_event_chi2']['dof']}", "pipeline": "70.9/35"},
    "inverted_null": {"mine_top_sde_detrended": R["bls_inverted_null"]["top3_detrended_sde"][0], "pipeline_inv_sde": 16.177595971009552, "pipeline_inv_period": 0.8552175060537904},
}
json.dump(R, open(os.path.join(OUT, "results.json"), "w"), indent=1)
print("[done] sap_by_season.json written, results.json updated with supplementary + pipeline_comparison")
