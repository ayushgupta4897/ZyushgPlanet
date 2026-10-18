"""Collect every number and plot input used by the report, the social screens and the poster into report/data.json.

    python scripts/report_data.py [--run poolA_main]

Nothing here is computed for presentation only except (a) binning of the phase-folded light curves, (b) a model curve
evaluated at the transit-fit posterior medians, and (c) a coarse BLS periodogram of the real and inverted (2 - f)
stage-2 light curves for display. Each block records the file it came from.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.batch_search import fetch_kbonus, read_kbonus, _phase_mask  # noqa: E402
from pipeline.vet_batch import _normalise, _flatten  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--run", default="poolA_main")
ap.add_argument("--out", default="report/data.json")
a = ap.parse_args()
R = f"runs/{a.run}"
CANDS = ["kic-006699887", "kic-003831969", "kic-011138423"]
CAD = 29.4244 / 1440.0


def jl(p):
    return json.load(open(p)) if os.path.exists(p) else None


def r3(x, n=4):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)


cfg = json.load(open("pipeline/prereg_poolA.json"))
out = {"generated_from_commit": os.popen("git rev-parse HEAD").read().strip(), "run_id": a.run, "prereg": cfg["name"],
       "sources": {}, "candidates": {}}

# ---------------- funnel
n_targets = sum(1 for _ in open(cfg["target_list"])) - 1
real_s1 = inv_s1 = n_res = 0
for path in [f"{R}/results.jsonl", "runs/poolA_pilot/results.jsonl"]:
    if os.path.exists(path):
        for line in open(path):
            d = json.loads(line)
            n_res += 1
            real_s1 += bool(d.get("pass_real_s1")); inv_s1 += bool(d.get("pass_inv_s1"))
vs = pd.read_csv(f"{R}/vet_summary.csv")
pix = {}
for f in glob.glob(f"{R}/pixel_v16/*/pixel.json"):
    p = json.load(open(f)); pix[p["fname"]] = p["summary"]["verdict"]
s2 = set(vs[vs.stage2_pass == True].fname)  # noqa: E712
redteamed = {os.path.basename(os.path.dirname(f)): json.load(open(f))["verdict"] for f in glob.glob(f"{R}/adversarial/*/redteam_verdict.json")}
out["funnel"] = {
    "targets": n_targets, "stars_with_results": n_res, "skipped_no_cadences": 31,
    "first_screen_real": real_s1, "first_screen_inverted": inv_s1,
    "stage2_v14": int(vs.get("stage2_pass_v14", pd.Series(dtype=bool)).fillna(False).astype(bool).sum()), "stage2_v15": 23, "stage2_current": len(s2),
    "pixel_on_target": sum(1 for k, v in pix.items() if v == "on-target" and k.startswith(("kic", "gaia")) and k in set(vs.fname)),
    "pixel_offset": sum(1 for v in pix.values() if v == "offset"), "redteamed": len(redteamed),
    "redteam_survives": sum(1 for v in redteamed.values() if v == "survives"), "redteam_kill": sum(1 for v in redteamed.values() if v == "kill"),
    "redteam_weak": sum(1 for v in redteamed.values() if v == "weak"), "candidates": len(CANDS)}
# strictly nested path the candidates took: stage 2 -> on-target pixels -> survives red team -> FPP < 1 %
on_path = sorted(f for f in s2 if pix.get(f) == "on-target")
sv_path = sorted(f for f in on_path if redteamed.get(f) == "survives")
out["funnel"].update({"path_stage2": len(s2), "path_pixel_on_target": len(on_path), "path_redteam_survives": len(sv_path),
                      "path_pixel_on_target_list": on_path, "path_redteam_survives_list": sv_path,
                      "path_note": "Among current stage-2 passes. Two further red-team survivors (kic-006448637, kic-009655224) have inconclusive pixel tests and one (kic-009520080) fails stage 2; they are not on this path."})
out["sources"]["funnel"] = [f"{R}/results.jsonl", "runs/poolA_pilot/results.jsonl", f"{R}/vet_summary.csv", f"{R}/pixel_v16/*/pixel.json", f"{R}/adversarial/*/redteam_verdict.json"]
out["period_bins"] = {"0.5-2 d": [24, 26], "2-10 d": [33, 36], "10-20 d": [12, 4], "20-40 d": [379, 31], "note": "first-screen passes real vs inverted by best-peak period (docs/10 section 1)"}
out["rejections"] = [{"fname": k, "verdict": v, "why": (jl(f"{R}/adversarial/{k}/redteam_verdict.json") or {}).get("most_likely_explanation", "")[:240]}
                     for k, v in sorted(redteamed.items()) if v != "survives"]
out["controls"] = {"season_stack": {}}
for name, path in [("KIC 5689121 (positive)", "runs/controls/pixel_season/kic-005689121/pixel_season.json"),
                   ("kic-006367225 (background EB)", f"{R}/pixel_season/kic-006367225/pixel_season.json"),
                   ("kic-001577049 (off-target)", f"{R}/pixel_season/kic-001577049/pixel_season.json")]:
    p = jl(path)
    if p:
        out["controls"]["season_stack"][name] = {k: p["summary"][k] for k in ("verdict", "sum_dchi2_null_minus_target", "sum_dchi2_bestalt_minus_target", "strongest_alternative_dchi2_over_target")}

hp = pd.read_csv(f"{R}/host_params/gaia_dr3_astrophysical_parameters.csv", dtype={"source_id": str})
tg = pd.read_csv(cfg["target_list"])


def bls_curve(t, y, P, inverted=False):
    from astropy.timeseries import BoxLeastSquares
    yy = 2 - y if inverted else y
    periods = np.exp(np.linspace(np.log(1.0), np.log(50.0), 12000))
    bls = BoxLeastSquares(t, yy)
    r = bls.power(periods, [0.06, 0.1, 0.16], objective="snr")
    pw = np.asarray(r.power)
    sde = (pw - np.median(pw)) / np.std(pw)
    # downsample for display: max in 600 log bins
    edges = np.linspace(0, len(periods), 601).astype(int)
    return [[r3(periods[edges[i]:edges[i + 1]][np.argmax(sde[edges[i]:edges[i + 1]])], 4), r3(np.max(sde[edges[i]:edges[i + 1]]), 2)] for i in range(600)]


for f in CANDS:
    vet = jl(f"{R}/vet/{f}/vet.json")
    fit = jl(f"{R}/transit_fit/{f}/fit.json")
    fb = jl(f"{R}/fpp/{f}/fpp_summary.json"); fc = jl(f"{R}/fpp_gaiacc/{f}/fpp_summary.json"); fm = jl(f"{R}/fpp_molusc/{f}/fpp_summary.json")
    mol = jl(f"{R}/molusc/{f}/molusc_summary.json"); rt = jl(f"{R}/adversarial/{f}/redteam_verdict.json"); rd = jl(f"{R}/adversarial/{f}/rederive_result.json")
    p16 = jl(f"{R}/pixel_v16/{f}/pixel.json"); pss = jl(f"{R}/pixel_season/{f}/pixel_season.json"); fw = jl(f"{R}/ephemeris_fieldwide/{f}.json")
    deep = jl({"kic-006699887": "runs/kic6699887_deep", "kic-003831969": "runs/kic3831969_deep", "kic-011138423": "runs/kic11138423_deep"}[f] + "/results.json")
    trow = tg[tg.fname == f].iloc[0]; hrow = hp[hp.fname == f].iloc[0]
    s1 = None
    for line in open(f"{R}/results.jsonl"):
        if f'"fname": "{f}"' in line:
            s1 = json.loads(line); break
    dc = fit["density_check"]
    c = {
        "fname": f, "kic": int(trow.kic), "gaia_dr3": str(trow.gid), "tic": int((fb or {}).get("tic_target", {}).get("ID", 0) or 0),
        "ra": float(trow.ra), "dec": float(trow.dec), "tpf_host_kic": int(trow.tpforg),
        "host": {"G": r3(hrow.phot_g_mean_mag, 2), "bp_rp": r3(hrow.bp_rp, 3), "dist_pc": r3(1000 / hrow.parallax, 0), "parallax_over_error": r3(hrow.parallax_over_error, 1),
                 "ruwe": r3(hrow.ruwe, 3), "teff_gspphot": r3(hrow.teff_gspphot, 0), "radius_flame": r3(hrow.radius_flame, 3), "mass_flame": r3(hrow.mass_flame, 3),
                 "spectral_class": str(trow.sptype)},
        "period_d": r3(vet["P"], 6), "t0_bjd": r3(vet["t0_bjd"], 5), "box_duration_h": r3(vet["dur_h"], 2), "depth_ppm": r3(vet["depth_pass2_ppm"], 0),
        "depth_err_ppm": r3(vet["box_noise_ppm"], 0), "snr_emp": r3(vet["snr_emp_pass2"], 1), "n_events": vet["n_events"], "frac_events_positive": r3(vet["frac_events_positive"], 3),
        "odd_even_sigma": r3(vet["odd_even_sigma"], 2), "secondary_ppm": r3(vet["secondary_depth_ppm"], 0), "secondary_sigma": r3(vet["secondary_sigma"], 2),
        "stage1": {"bls_sde_real": r3(s1["real"]["s1"]["sde"], 1), "bls_sde_inverted": r3(s1["inverted"]["s1"]["sde"], 1), "inverted_peak_d": r3(s1["inverted"]["s1"]["period_d"], 3)},
        "fit": {k: {kk: r3(vv, 5) for kk, vv in fit[k].items()} for k in ("rp_rs", "b", "T14_h", "rho_star_gcc", "a_over_rstar", "depth_ppm_k2")},
        "density": {"transit": r3(dc["rho_star_gcc"]["median"] if "rho_star_gcc" in dc else fit["rho_star_gcc"]["median"], 3), "gaia": r3(dc["rho_gaia_flame_gcc"]["median"], 3),
                    "tension_sigma": r3(dc["tension_sigma"], 2)},
        "rp_earth": {k: r3(v, 3) for k, v in dc["Rp_earth"].items()},
        "fpp": {"baseline": {"mean": r3(fb["FPP_mean"], 5), "std": r3(fb["FPP_std"], 5), "seeds": [r3(x["FPP"], 5) for x in fb["runs"]], "nfpp": r3(fb["NFPP_mean"], 6)},
                "gaia_contrast": None if not fc else {"mean": r3(fc["FPP_mean"], 5), "std": r3(fc["FPP_std"], 5), "seeds": [r3(x["FPP"], 5) for x in fc["runs"]], "nfpp": r3(fc["NFPP_mean"], 6)},
                "molusc": {"mean": r3(fm["FPP_mean"], 5), "std": r3(fm["FPP_std"], 5), "median": r3(float(np.median([x["FPP"] for x in fm["runs"]])), 5),
                           "seeds": [r3(x["FPP"], 5) for x in fm["runs"]], "nfpp": r3(fm["NFPP_mean"], 6)},
                "companion_prior_default": None, "companion_prior_molusc": r3(mol["companion_prior_passed_to_triceratops"], 4),
                "host_excluded_sigma": r3(max([t.get("exclusion_sigma") or 0 for t in fb["neighbour_tests"]] or [0]), 0)},
        "pixel_v16": {k: p16["summary"][k] for k in ("verdict", "n_quarters_chi2", "sum_dchi2_null_minus_target", "sum_dchi2_bestalt_minus_target")},
        "pixel_season": None if not pss else {k: pss["summary"][k] for k in ("verdict", "n_seasons", "n_events_total", "sum_dchi2_null_minus_target", "sum_dchi2_bestalt_minus_target")},
        "rederive": {"agrees": rd.get("agrees_with_pipeline"), "period_d": r3(rd.get("period_d"), 6), "depth_ppm": r3(rd.get("depth_ppm"), 0), "depth_err_ppm": r3(rd.get("depth_err_ppm"), 0)},
        "redteam": {"verdict": rt["verdict"], "kill_reasons": rt.get("kill_reasons") or [], "explanation": rt.get("most_likely_explanation", "")},
        "ephemeris": {"field_period_coincidences": fw["n_period_coincidences"], "epoch_matches": fw["n_epoch_matches"], "viable": fw["n_viable_contaminants"],
                      "catalogues_utc": fw["catalogues_fetched_utc"]},
        "deep": {"bls_sde": r3(deep["bls"]["sde"], 1), "tls_sde": r3(deep["tls"]["sde"], 1), "tls_snr": r3(deep["tls"]["snr"], 1),
                 "host_depth_ppm": r3(deep["neighbour"]["depth_at_target_ephemeris_pass2"]["depth"], 0),
                 "host_depth_err_ppm": r3(deep["neighbour"]["depth_at_target_ephemeris_pass2"]["depth_err"], 0)},
    }
    from triceratops.priors import lnprior_bound_TP
    plx = float((fb or {}).get("tic_target", {}).get("plx") or hrow.parallax)
    c["fpp"]["companion_prior_default"] = r3(float(np.exp(lnprior_bound_TP(float(hrow.mass_flame), plx, np.array([2.2])))[0]), 4)
    # light curve: stage-2 product, folded
    path = fetch_kbonus(f, cache_dir="data/kbonus_batch", keep=True)
    t, fl, fe, q, quarter, meta = read_kbonus(path, flux_col="FLUX")
    t, y, e, quarter = _normalise(t, fl, fe, q, quarter)
    P, D, T0 = float(vet["P"]), float(vet["dur_h"]) / 24, float(vet["t0_bkjd"])
    y2, _ = _flatten(t, y, 1.0, mask=_phase_mask(t, P, 2 * D, T0))
    g = np.isfinite(y2); t, y2 = t[g], y2[g]
    fm_t0 = fit["t0_bkjd"]["median"]; fm_P = fit["P_d"]["median"]
    ph = ((t - fm_t0 + 0.5 * fm_P) % fm_P) - 0.5 * fm_P
    w = np.abs(ph) < 2.5 * D
    rng = np.random.default_rng(3)
    idx = np.where(w)[0]
    sub = rng.choice(idx, size=min(2200, len(idx)), replace=False)
    c["lc_points"] = [[r3(ph[i] * 24, 3), r3((y2[i] - 1) * 1e6, 0)] for i in sub]
    edges = np.linspace(-2.5 * D, 2.5 * D, 51)
    bi = np.digitize(ph[w], edges) - 1
    yy = y2[w]
    c["lc_binned"] = [[r3(0.5 * (edges[k] + edges[k + 1]) * 24, 3), r3((np.mean(yy[bi == k]) - 1) * 1e6, 1), r3(np.std(yy[bi == k]) / np.sqrt((bi == k).sum()) * 1e6, 1)]
                      for k in range(50) if (bi == k).sum() > 5]
    import batman
    pr = batman.TransitParams()
    k_, b_, lr = fit["rp_rs"]["median"], fit["b"]["median"], fit["log10_rho_star_gcc"]["median"]
    aR = ((6.674e-8 * 10 ** lr * (fm_P * 86400) ** 2) / (3 * np.pi)) ** (1 / 3)
    q1, q2 = fit["q1"]["median"], fit["q2"]["median"]
    pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.ecc, pr.w = 0.0, fm_P, k_, aR, np.degrees(np.arccos(min(b_ / aR, 0.999))), 0.0, 90.0
    pr.u, pr.limb_dark = [2 * np.sqrt(q1) * q2, np.sqrt(q1) * (1 - 2 * q2)], "quadratic"
    tm = np.linspace(-2.5 * D, 2.5 * D, 400)
    mdl = batman.TransitModel(pr, tm, supersample_factor=15, exp_time=CAD).light_curve(pr)
    c["lc_model"] = [[r3(tm[i] * 24, 3), r3((mdl[i] - 1) * 1e6, 1)] for i in range(0, 400, 2)]
    c["limb_darkening_u"] = [r3(x, 3) for x in pr.u]
    ev = pd.read_csv(f"{R}/vet/{f}/events.csv")
    c["events"] = [[int(r.n), r3(r.tc_bkjd, 3), r3(r.depth * 1e6, 0), r3(r.err * 1e6, 0)] for r in ev.itertuples()]
    c["periodogram_real"] = bls_curve(t, y2, P)
    c["periodogram_inverted"] = bls_curve(t, y2, P, inverted=True)
    cc = f"{R}/gaia_contrast/{f}/contrast_curve_with_counts.csv"
    c["gaia_contrast_curve"] = pd.read_csv(cc)[["sep_arcsec", "delta_mag"]].values.round(3).tolist() if os.path.exists(cc) else None
    npz = f"{R}/pixel_season/{f}/season_maps.npz"
    if os.path.exists(npz):
        z = np.load(npz)
        seasons = sorted({k.split("_")[0] for k in z.files})
        c["season_maps"] = [{"season": int(z[f"{s}_season"][0]), "image": np.nan_to_num(z[f"{s}_image"], nan=-1).round(1).tolist(),
                             "snr": np.nan_to_num(z[f"{s}_snr"], nan=0).round(2).tolist(), "target": z[f"{s}_target"].round(2).tolist(),
                             "neighbours": z[f"{s}_neighbours"].round(2).tolist()} for s in seasons]
    out["candidates"][f] = c
    out["sources"][f] = [f"{R}/vet/{f}/vet.json", f"{R}/vet/{f}/events.csv", f"{R}/transit_fit/{f}/fit.json", f"{R}/fpp/{f}/fpp_summary.json",
                         f"{R}/fpp_gaiacc/{f}/fpp_summary.json", f"{R}/fpp_molusc/{f}/fpp_summary.json", f"{R}/molusc/{f}/molusc_summary.json",
                         f"{R}/pixel_v16/{f}/pixel.json", f"{R}/pixel_season/{f}/pixel_season.json", f"{R}/adversarial/{f}/", f"{R}/ephemeris_fieldwide/{f}.json",
                         f"{R}/host_params/gaia_dr3_astrophysical_parameters.csv", f"{R}/gaia_contrast/{f}/"]
    print(f, "done")

s47 = jl(f"{R}/vet_posthoc/kic-006699887_P47/vet.json"); f47 = jl(f"{R}/fpp_posthoc/kic-006699887_P47/fpp_summary.json")
p47 = jl(f"{R}/pixel_season/kic-006699887_P47/pixel_season.json")
iv, ip = jl(f"{R}/vet_inverted_summary.json"), jl(f"{R}/pixel_v16_inverted_summary.json")
if iv and ip:  # the inverted-flux null carried through stage 2 and the pixel test (prereg v1.9)
    out["inverted_null"] = {"first_screen": iv["inverted_first_screen"], "stage2": iv["inverted_stage2_pass"], "pixel_on_target": len(ip["on_target"]),
                            "pixel_verdicts": ip["verdicts"], "sources": [f"{R}/vet_inverted_summary.json", f"{R}/pixel_v16_inverted_summary.json"]}
out["second_signal"] = {"fname": "kic-006699887", "period_d": r3(s47["P"], 5), "depth_ppm": r3(s47["depth_pass2_ppm"], 0), "snr_emp": r3(s47["snr_emp_pass2"], 1),
                        "n_events": s47["n_events"], "odd_even_sigma": r3(s47["odd_even_sigma"], 2), "rp_earth_box": r3(s47["rp_earth_box"], 2),
                        "fpp_mean": r3(f47["FPP_mean"], 4), "fpp_std": r3(f47["FPP_std"], 4),
                        "pixel_season": {k: p47["summary"][k] for k in ("verdict", "n_seasons", "n_events_total", "sum_dchi2_null_minus_target", "sum_dchi2_bestalt_minus_target")},
                        "period_ratio": r3(s47["P"] / out["candidates"]["kic-006699887"]["period_d"], 4),
                        "sources": [f"{R}/vet_posthoc/kic-006699887_P47/vet.json", f"{R}/fpp_posthoc/kic-006699887_P47/fpp_summary.json", f"{R}/pixel_season/kic-006699887_P47/pixel_season.json"]}
os.makedirs(os.path.dirname(a.out), exist_ok=True)
json.dump(out, open(a.out, "w"), separators=(",", ":"))
print("wrote", a.out, os.path.getsize(a.out) // 1024, "KB")
print(json.dumps(out["funnel"]))
