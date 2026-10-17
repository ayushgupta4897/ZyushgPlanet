#!/usr/bin/env python3
"""Merge the fine period scan and free-ephemeris trapezoid fit into results.json and
write comparison.json against the pipeline's vet.json (read only AFTER the
independent numbers were produced).  Own code; no pipeline/ imports."""
import json, os, sys
import numpy as np, pandas as pd
R = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-003831969/rederive"
S = sys.argv[1]            # scratch dir holding perr.json / trapfit.json
res = json.load(open(f"{R}/results.json"))
perr = json.load(open(f"{S}/perr.json")); trap = json.load(open(f"{S}/trapfit.json"))
res["period_fine_scan_bls"] = dict(**perr, note="BLS delta-loglike=0.5 interval; box edges make the surface spiky, so this interval is optimistic")
res["trapezoid_free_ephemeris"] = dict(**trap, note="least-squares trapezoid with P, t_ref, depth, T14, T23 free; errors from the Jacobian scaled by reduced chi2 (white noise)")
res["adopted"] = dict(period_d=perr["P_best"], period_err_d=trap["P_err_d"],
                      t0_first_bkjd=perr["t0_bkjd"], t0_first_bjd=perr["t0_bkjd"] + 2454833.0,
                      t0_err_min_at_tref=trap["tref_err_min"],
                      depth_box_ppm=res["primary"]["depth_ppm"], depth_box_err_emp_ppm=res["primary"]["err_emp_ppm"],
                      depth_trapezoid_ppm=trap["depth_ppm"], depth_trapezoid_err_ppm=trap["depth_err_ppm"],
                      duration_box_h=res["ephemeris"]["duration_box_h"], T14_trapezoid_h=trap["T14_h"],
                      snr_emp=res["primary"]["snr_emp"], sde_blind=res["blind_bls"]["sde_detr"],
                      n_transits=res["events"]["n_transits_with_ge_half_coverage"])
v = json.load(open("/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-003831969/vet.json"))
a = pd.read_csv(f"{R}/events.csv"); b = pd.read_csv("/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-003831969/events.csv")
idx = np.array([int(np.argmin(np.abs(b.tc_bkjd.values - x))) for x in a.tmid_bkjd])
ok = np.abs(b.tc_bkjd.values[idx] - a.tmid_bkjd.values) < 0.1
dm, dp = a.depth_ppm.values[ok], b.depth.values[idx][ok] * 1e6
ev = dict(n_matched=int(ok.sum()), corr=float(np.corrcoef(dm, dp)[0, 1]), mean_diff_mine_minus_pipeline_ppm=float(np.mean(dm - dp)),
          rms_diff_ppm=float(np.std(dm - dp)), frac_positive_mine=float((a.depth_ppm > 0).mean()), frac_positive_pipeline=float(v["frac_events_positive"]))
def row(name, mine, pipe, unit=""):
    return dict(quantity=name, mine=mine, pipeline=pipe, unit=unit)
yq = [res["per_year"][k]["depth_ppm"] for k in sorted(res["per_year"])]
ye = [res["per_year"][k]["err_emp_ppm"] for k in sorted(res["per_year"])]
cmp_ = dict(
    note="pipeline numbers quoted from runs/poolA_main/vet/kic-003831969/vet.json (run id poolA_main); mine from results.json in this directory",
    rows=[row("period (d)", perr["P_best"], v["P"], "d"),
          row("period err (d) [mine: trapezoid covariance; BLS dLL interval = %.1e]" % (0.5 * (perr["P_1sig_hi"] - perr["P_1sig_lo"])), trap["P_err_d"], None, "d"),
          row("t0 first epoch (BKJD)", perr["t0_bkjd"], v["t0_bkjd"], "BKJD"),
          row("t0 difference (min)", (perr["t0_bkjd"] - v["t0_bkjd"]) * 1440, None, "min"),
          row("box duration (h)", res["ephemeris"]["duration_box_h"], v["dur_h"], "h"),
          row("T14 trapezoid (h) vs pipeline t14_max for Rstar=0.83", trap["T14_h"], v["t14_max_h"], "h"),
          row("box depth (ppm)", res["primary"]["depth_ppm"], v["depth_pass2_ppm"], "ppm"),
          row("box depth error (ppm) [mine empirical random-phase; pipeline formal]", res["primary"]["err_emp_ppm"], v["depth_pass2_formal_err_ppm"], "ppm"),
          row("box noise used for SNR (ppm)", res["primary"]["err_emp_ppm"], v["box_noise_ppm"], "ppm"),
          row("empirical SNR", res["primary"]["snr_emp"], v["snr_emp_pass2"], ""),
          row("trapezoid/TLS full depth (ppm)", trap["depth_ppm"], (1 - v["tls"]["depth"]) * 1e6, "ppm"),
          row("n in-transit points", res["primary"]["n_in"], v["n_in_transit"], ""),
          row("odd depth (ppm)", res["odd_even"]["odd"]["depth_ppm"], v["depth_odd_ppm"], "ppm"),
          row("even depth (ppm)", res["odd_even"]["even"]["depth_ppm"], v["depth_even_ppm"], "ppm"),
          row("odd/even sigma", res["odd_even"]["sigma"], v["odd_even_sigma"], ""),
          row("secondary depth at phase 0.5 (ppm)", res["secondary_phase05"]["depth_ppm"], v["secondary_depth_ppm"], "ppm"),
          row("secondary sigma", res["secondary_phase05"]["snr_emp"], v["secondary_sigma"], ""),
          row("n transits (>=50% coverage)", res["events"]["n_transits_with_ge_half_coverage"], v["n_events"], ""),
          row("per-event chi2 about constant depth / dof", f"{res['events']['chi2_const_depth']:.0f}/{res['events']['dof']}", f"{v['event_chi2_const']:.0f}/{v['event_dof']}", ""),
          row("single-event dominance (max event chi2 share)", res["events"]["single_event_dominance"], None, ""),
          row("per-year depths (ppm) [mine: Kepler years Q1-4,5-8,9-12,13-17]", [round(x) for x in yq], None, "ppm"),
          row("per-year errors (ppm)", [round(x) for x in ye], None, "ppm"),
          row("per-quartile depths (ppm) [pipeline: event quartiles]", None, [round(x[0]) for x in v["depth_by_quartile_ppm"]], "ppm"),
          row("per-year/quartile chi2 about mean / dof", f"{res['per_year_consistency']['chi2']:.1f}/3", None, ""),
          row("SAP_FLUX box depth (ppm)", res["sap"]["depth_ppm"], None, "ppm"),
          row("SAP/PSF depth ratio", res["sap_vs_psf"]["depth_ratio_sap_over_psf"], None, ""),
          row("blind BLS SDE (own definition, 1-50 d)", res["blind_bls"]["sde_detr"], None, ""),
          row("inverted-flux null top SDE (own definition)", res["inverted_null"]["sde_detr"], None, ""),
          row("TLS SDE (pipeline only)", None, v["tls"]["SDE"], ""),
          row("stage2_pass (pipeline)", None, v["stage2_pass"], "")],
    per_event=ev)
res["comparison_with_pipeline"] = cmp_
json.dump(cmp_, open(f"{R}/comparison.json", "w"), indent=1)
json.dump(res, open(f"{R}/results.json", "w"), indent=1)
for r in cmp_["rows"]:
    print(f"{r['quantity']:<75} mine={r['mine']!s:<28} pipe={r['pipeline']!s}")
print("per-event:", ev)
