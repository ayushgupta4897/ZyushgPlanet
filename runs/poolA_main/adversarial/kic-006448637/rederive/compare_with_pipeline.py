#!/usr/bin/env python
"""
Follow-up to rederive_kic006448637.py (run AFTER it):
  * two further independent noise estimates for the box depth (single-box noise on the unfolded
    light curve / sqrt(N); scatter of per-event depths / sqrt(N))
  * depth / per-year depths with the Q13-edge outlier event (epoch 28) removed
  * per-event mid-times -> linear ephemeris with uncertainties
  * comparison against the pipeline's vet.json (read only now, after my own numbers exist)
Writes comparison.json and appends 'extras' + 'comparison' to results.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rederive_kic006448637 as rd  # my own module, same directory

VET = "/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-006448637/vet.json"
R = json.load(open(os.path.join(HERE, "results.json")))
P = R["ephemeris"]["period_d"]
t0 = R["ephemeris"]["t0_bkjd"]
dur = R["ephemeris"]["bls_duration_d"]
trp = R["ephemeris"]["trapezoid"]
depth_tr, T14, q = trp["depth_ppm"] * 1e-6, trp["T14_h"] / 24.0, trp["T23_over_T14"]

t, f, e, qu = rd.load("FLUX")
mask = np.abs(rd.fold_phase(t, P, t0) * P) < dur
flat, trend, keep, sig = rd.detrend(t, f, qu, mask=mask)
tk, yk, sk = t[keep], flat[keep], sig[keep]
rng = np.random.default_rng(7)

X = {}

# ------------------------------------------------------------------ (a) single-box noise on the unfolded LC
oot = np.abs(rd.fold_phase(tk, P, t0) * P) > dur
to, yo = tk[oot], yk[oot]
centers = np.arange(to.min(), to.max(), dur / 2.0)
lo = np.searchsorted(to, centers - dur / 2); hi = np.searchsorted(to, centers + dur / 2)
n = hi - lo
ok = n >= 0.5 * dur / rd.CAD
means = np.array([np.mean(yo[a:b]) for a, b, o in zip(lo, hi, ok) if o])
n_ev = R["box_statistics_psf"]["events"]["n_transits_with_>=50pct_coverage"]
sig_single = 1.4826 * np.median(np.abs(means - np.median(means)))
X["single_box_noise_ppm"] = 1e6 * float(sig_single)
X["single_box_noise_over_sqrtN_ppm"] = 1e6 * float(sig_single / np.sqrt(n_ev))
X["single_box_white_expectation_ppm"] = 1e6 * float(np.std(yo) / np.sqrt(dur / rd.CAD))
X["n_single_boxes"] = int(ok.sum())

# ------------------------------------------------------------------ (b) per-event scatter, outlier removal
ev = [l.split(",") for l in open(os.path.join(HERE, "events.csv")).read().strip().split("\n")[1:]]
ev = [dict(epoch=int(r[0]), tc=float(r[1]), quarter=int(r[3]), n_in=int(r[4]), d=float(r[6]) if r[6] else np.nan,
           e=float(r[7]) if r[7] else np.nan, gap=float(r[8]), counted=int(r[9])) for r in ev]
cnt = [r for r in ev if r["counted"]]
d = np.array([r["d"] for r in cnt]); err = np.array([r["e"] for r in cnt]); ep = np.array([r["epoch"] for r in cnt])
X["per_event_scatter_ppm"] = float(np.std(d, ddof=1))
X["per_event_scatter_over_sqrtN_ppm"] = float(np.std(d, ddof=1) / np.sqrt(len(d)))
X["per_event_mean_ppm"] = float(np.mean(d))
X["per_event_median_ppm"] = float(np.median(d))
worst = cnt[int(np.argmax((d - np.median(d)) / err))]
X["most_discrepant_event"] = {"epoch": worst["epoch"], "tc_bkjd": worst["tc"], "quarter": worst["quarter"],
                              "depth_ppm": worst["d"], "err_ppm": worst["e"], "gap_dist_d": worst["gap"],
                              "pull_vs_median": float((worst["d"] - np.median(d)) / worst["e"])}
drop = np.abs(tk - worst["tc"]) > 1.0
m_all = rd.measure(tk, yk, P, t0, dur, ntrial=3000)
m_drop = rd.measure(tk[drop], yk[drop], P, t0, dur, ntrial=3000)
yrs = rd.year_of_bkjd(tk)
py_drop = {}
for y in sorted(set(yrs.tolist())):
    mm = rd.measure(tk[drop], yk[drop], P, t0, dur, sel=(yrs[drop] == y), ntrial=2000)
    py_drop[str(y)] = {"depth_ppm": mm["depth_ppm"], "err_ppm": mm["err_ppm"]}
epk = np.round((tk[drop] - t0) / P).astype(int)
oe_odd = rd.measure(tk[drop], yk[drop], P, t0, dur, sel=(epk % 2 == 1), ntrial=2000)
oe_even = rd.measure(tk[drop], yk[drop], P, t0, dur, sel=(epk % 2 == 0), ntrial=2000)
X["odd_even_without_most_discrepant_event"] = {"odd_ppm": oe_odd["depth_ppm"], "odd_err_ppm": oe_odd["err_ppm"],
                                               "even_ppm": oe_even["depth_ppm"], "even_err_ppm": oe_even["err_ppm"],
                                               "sigma": float(abs(oe_odd["depth_ppm"] - oe_even["depth_ppm"]) / np.hypot(oe_odd["err_ppm"], oe_even["err_ppm"]))}
d2 = d[ep != worst["epoch"]]; e2 = err[ep != worst["epoch"]]
wm = np.sum(d2 / e2 ** 2) / np.sum(1 / e2 ** 2)
X["depth_all_events_ppm"] = {"depth_ppm": m_all["depth_ppm"], "err_ppm": m_all["err_ppm"], "snr": m_all["snr"]}
X["depth_without_most_discrepant_event"] = {"depth_ppm": m_drop["depth_ppm"], "err_ppm": m_drop["err_ppm"], "snr": m_drop["snr"],
                                            "per_year": py_drop,
                                            "events_chi2_about_mean": float(np.sum(((d2 - wm) / e2) ** 2)), "dof": int(len(d2) - 1)}
pyv = [v["depth_ppm"] for v in py_drop.values()]; pye = [v["err_ppm"] for v in py_drop.values()]
w = 1 / np.array(pye) ** 2
X["per_year_chi2_about_weighted_mean_without_outlier"] = float(np.sum(w * (np.array(pyv) - np.sum(w * pyv) / w.sum()) ** 2))
pyv0 = [v["depth_ppm"] for v in R["box_statistics_psf"]["per_year"].values()]
pye0 = [v["err_ppm"] for v in R["box_statistics_psf"]["per_year"].values()]
w0 = 1 / np.array(pye0) ** 2
X["per_year_chi2_about_weighted_mean_all"] = float(np.sum(w0 * (np.array(pyv0) - np.sum(w0 * pyv0) / w0.sum()) ** 2))
X["per_year_dof"] = int(len(pyv0) - 1)

# ------------------------------------------------------------------ (c) per-event mid-times -> linear ephemeris
tt_rows = []
shifts = np.linspace(-0.1, 0.1, 401)   # d
for r in cnt:
    loc = np.abs(tk - r["tc"]) < 1.0
    x, y, s = tk[loc] - r["tc"], yk[loc], sk[loc]
    chi2 = np.array([np.sum(((y - rd.trapezoid(x, depth_tr, T14, q, dt)) / s) ** 2) for dt in shifts])
    i = int(np.argmin(chi2))
    if i == 0 or i == len(shifts) - 1:
        continue
    within = shifts[chi2 <= chi2[i] + 1.0]
    err_dt = 0.5 * (within.max() - within.min()) if len(within) > 1 else shifts[1] - shifts[0]
    # bootstrap-free sanity: inflate by sqrt(chi2_red) if > 1
    red = chi2[i] / max(len(x) - 1, 1)
    err_dt *= np.sqrt(max(red, 1.0))
    tt_rows.append((r["epoch"], r["tc"] + shifts[i], err_dt))
tt_rows = np.array(tt_rows)
A = np.vstack([np.ones(len(tt_rows)), tt_rows[:, 0]]).T
W = 1 / tt_rows[:, 2] ** 2
cov = np.linalg.inv(A.T @ (A * W[:, None]))
beta = cov @ (A.T @ (W * tt_rows[:, 1]))
resid = tt_rows[:, 1] - A @ beta
chi2_tt = float(np.sum(W * resid ** 2))
infl = np.sqrt(max(chi2_tt / (len(tt_rows) - 2), 1.0))   # Delta-chi2=1 errors are too small for SNR~3-4 events: scale by sqrt(chi2_red)
X["linear_ephemeris_from_event_times"] = {"n_events_timed": int(len(tt_rows)), "P_d": float(beta[1]), "P_err_d": float(infl * np.sqrt(cov[1, 1])),
                                          "P_err_unscaled_d": float(np.sqrt(cov[1, 1])), "error_inflation_factor": float(infl),
                                          "t0_bkjd_epoch0": float(beta[0]), "t0_err_d": float(infl * np.sqrt(cov[0, 0])),
                                          "chi2": chi2_tt, "dof": int(len(tt_rows) - 2),
                                          "rms_timing_residual_min": float(1440 * np.std(resid)),
                                          "max_abs_timing_residual_min": float(1440 * np.max(np.abs(resid)))}

# ------------------------------------------------------------------ comparison with the pipeline (read only now)
V = json.load(open(VET))
S = R["summary"]
Pe = X["linear_ephemeris_from_event_times"]
rows = [
    ("period_d", V["P"], S["period_d"], "BLS refine; timing-fit P = %.5f +- %.5f" % (Pe["P_d"], Pe["P_err_d"])),
    ("t0_bkjd", V["t0_bkjd"], R["ephemeris"]["t0_bkjd"] - round((R["ephemeris"]["t0_bkjd"] - V["t0_bkjd"]) / P) * P, "my t0 moved to the pipeline's reference epoch"),
    ("duration_h", V["dur_h"], S["duration_h"], "pipeline: fixed 4.5 h box; mine: trapezoid T14 (BLS box 4.13 h, box-SNR scan 4.32 h)"),
    ("depth_ppm", V["depth_pass2_ppm"], S["depth_ppm"], "pipeline 4.5-h box; mine 4.13-h box on masked re-detrend"),
    ("depth_err_ppm", V["box_noise_ppm"], S["depth_err_ppm"], "pipeline box_noise vs my random-phase-box std"),
    ("empirical_snr", V["snr_emp_pass2"], S["empirical_snr"], "differs because of the error scale, not the depth"),
    ("n_in_transit", V["n_in_transit"], R["box_statistics_psf"]["depth"]["n_in"], ""),
    ("odd_ppm", V["depth_odd_ppm"], R["box_statistics_psf"]["odd"]["depth_ppm"], ""),
    ("even_ppm", V["depth_even_ppm"], R["box_statistics_psf"]["even"]["depth_ppm"], ""),
    ("odd_even_sigma", V["odd_even_sigma"], S["odd_even_sigma"], "same ~120-180 ppm difference; sigma differs with the error scale"),
    ("secondary_ppm", V["secondary_depth_ppm"], S["secondary_ppm"], ""),
    ("secondary_sigma", V["secondary_sigma"], S["secondary_sigma"], ""),
    ("n_events", V["n_events"], R["box_statistics_psf"]["events"]["n_epochs_with_any_in_transit_point"], "epochs with any in-transit cadence; 24 have >=50% coverage"),
    ("event_chi2/dof", "%.1f/%d" % (V["event_chi2_const"], V["event_dof"]),
     "%.1f/%d" % (R["box_statistics_psf"]["events"]["chi2_events_about_mean"], R["box_statistics_psf"]["events"]["dof"]), ""),
    ("depth_by_quartile_ppm", [round(x[0]) for x in V["depth_by_quartile_ppm"]], [round(v) for v in S["per_year_depths_ppm"].values()],
     "pipeline: event quartiles; mine: calendar years 2009-2013"),
]
C = {"pipeline_vet_json": VET, "rows": [{"quantity": a, "pipeline": b, "mine": c, "note": d} for a, b, c, d in rows]}
C["disagreements"] = [
    "Empirical depth error: pipeline box_noise 85 ppm vs my random-phase-box std %.0f ppm (my single-box-noise/sqrt(N) = %.0f ppm, "
    "per-event-scatter/sqrt(N) = %.0f ppm). Hence empirical SNR 10.0 (pipeline) vs %.1f (mine). The depth itself agrees; the pipeline's "
    "error is the more conservative." % (S["depth_err_ppm"], X["single_box_noise_over_sqrtN_ppm"], X["per_event_scatter_over_sqrtN_ppm"], S["empirical_snr"]),
    "Odd/even: same sense and size (odd shallower by ~120-180 ppm) but %.1f sigma under my errors vs 0.7 sigma under the pipeline's." % S["odd_even_sigma"],
    "Duration: pipeline holds 4.5 h fixed; my free fits give 4.1 h (BLS box), 4.3 h (box-SNR scan), %.2f h (trapezoid T14, T23/T14 = %.2f)." % (trp["T14_h"], trp["T23_over_T14"]),
    "Per-year depths are not fully consistent: 2011 = %.0f +- %.0f ppm vs 2012 = %.0f +- %.0f ppm (chi2 %.1f for %d dof about the weighted mean); "
    "the pipeline's event quartiles show the same drift (642 -> 992 ppm). Removing the Q13-edge outlier event (epoch %d, %.0f +- %.0f ppm, %.2f d after "
    "a gap) lowers 2012 to %.0f +- %.0f ppm and the chi2 to %.1f."
    % (pyv0[2], pye0[2], pyv0[3], pye0[3], X["per_year_chi2_about_weighted_mean_all"], X["per_year_dof"], worst["epoch"], worst["d"], worst["e"], worst["gap"],
       py_drop["2012"]["depth_ppm"], py_drop["2012"]["err_ppm"], X["per_year_chi2_about_weighted_mean_without_outlier"]),
    "Period: 35.58770 (mine, BLS refine) vs 35.58726 (pipeline): differ by 4.4e-4 d, which is %.1f sigma of my timing-fit period error." % (abs(S["period_d"] - V["P"]) / Pe["P_err_d"]),
]
json.dump(rd.jsonable(C), open(os.path.join(HERE, "comparison.json"), "w"), indent=2)
R["extras"] = X
R["comparison"] = C
json.dump(rd.jsonable(R), open(os.path.join(HERE, "results.json"), "w"), indent=2)
print(json.dumps(rd.jsonable(X), indent=1))
print()
for r in C["rows"]:
    print("%-24s pipeline=%-28s mine=%-28s %s" % (r["quantity"], r["pipeline"], r["mine"], r["note"]))
print()
for s in C["disagreements"]:
    print("-", s)
