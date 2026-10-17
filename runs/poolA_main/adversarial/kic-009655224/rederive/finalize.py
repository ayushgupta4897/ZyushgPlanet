#!/usr/bin/env python
"""Assemble comparison with the pipeline's vet.json (read only AFTER the independent numbers were produced),
check the t0 of the neighbour's narrow-BLS peak, write final JSON + SHA-256 manifest."""
import json, os, hashlib, sys, numpy as np
from astropy.timeseries import BoxLeastSquares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rederive_kic009655224 as R
OUT = R.OUT
res = json.load(open(f"{OUT}/rederive_results.json"))
s1 = json.load(open(f"{OUT}/rederive_supplement.json")); s2 = json.load(open(f"{OUT}/rederive_supplement2.json"))
nb1 = json.load(open(f"{OUT}/neighbour_check_kic009655236.json")); nb2 = json.load(open(f"{OUT}/neighbour_check_kic-009655226.json"))
vet = json.load(open("/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-009655224/vet.json"))
px = json.load(open("/home/user/ZyushgPlanet/runs/poolA_main/pixel/kic-009655224/pixel.json"))
P = res["ephemeris"]["period_d"]; t0 = res["ephemeris"]["t0_bjd"] - 2454833.0; DUR = res["ephemeris"]["duration_h"] / 24

# --- t0 of the neighbour (KIC 9655226) narrow-BLS peak vs the target's t0
nbf = "/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/rederive-kic009655224/untrusted_dl/nb_kic-009655226.fits"
D = R.load(nbf); flat, keep, mad = R.detrend(D["t"], D["fn"], 1.0); t, y = D["t"][keep], flat[keep]
pg = np.linspace(P * 0.98, P * 1.02, 6001)
r = BoxLeastSquares(t, y, np.full_like(y, mad)).power(pg, [DUR], objective="likelihood", oversample=10)
i = int(np.argmax(r.power)); Pn, t0n = float(pg[i]), float(r.transit_time[i])
# phase of the neighbour's peak transit time in the target's ephemeris, evaluated at the mid-time of the neighbour's data
tmid = 0.5 * (t.min() + t.max())
# predicted target transit nearest tmid, and neighbour's nearest predicted transit
k_t = np.round((tmid - t0) / P); tc_t = t0 + k_t * P
k_n = np.round((tc_t - t0n) / Pn); tc_n = t0n + k_n * Pn
dt_mid_h = float((tc_n - tc_t) * 24)
nb_peak = dict(best_P_d=Pn, t0_bkjd=t0n, dlogL=float(r.power[i]), depth_ppm=float(r.depth[i] * 1e6), deficit_e_s=float(r.depth[i] * nb2["at_target_ephemeris"]["PSF"]["median_flux_e_s"]),
               offset_from_target_transit_at_data_midtime_h=dt_mid_h, period_offset_frac=float(Pn / P - 1),
               n_independent_freqs_in_window=float((1 / (P * 0.98) - 1 / (P * 1.02)) * (t.max() - t.min())),
               comment="neighbour's own best peak in +/-2%% window; dlogL_at_target_P=%.1f vs window mean %.1f +/- %.1f" % (nb2["at_target_ephemeris"]["PSF"]["narrow_bls"]["dlogL_at_target_P"], nb2["at_target_ephemeris"]["PSF"]["narrow_bls"]["mean_dlogL"], nb2["at_target_ephemeris"]["PSF"]["narrow_bls"]["std_dlogL"]))
print("neighbour KIC 9655226 narrow-BLS peak:", nb_peak)

# --- comparison block
m = res; d = m["depth"]; mk = s2["masked_redetrend"]["window_1.0d"]
cmp = dict(
    pipeline=dict(P_d=vet["P"], P_screen_d=vet["P_screen"], tls_P_d=vet["tls"]["period"], t0_bkjd=vet["t0_bkjd"], dur_h=vet["dur_h"], dur_screen_h=vet["dur_screen_h"],
                  depth_pass1_ppm=vet["depth_pass1_ppm"], depth_pass2_ppm=vet["depth_pass2_ppm"], depth_pass2_err_ppm=vet["depth_pass2_formal_err_ppm"], snr_emp_pass2=vet["snr_emp_pass2"],
                  tls_sde=vet["tls"]["SDE"], tls_snr=vet["tls"]["snr"], odd_ppm=vet["depth_odd_ppm"], even_ppm=vet["depth_even_ppm"], odd_even_sigma=vet["odd_even_sigma"],
                  secondary_ppm=vet["secondary_depth_ppm"], secondary_sigma=vet["secondary_sigma"], n_events=vet["n_events"], n_in_transit=vet["n_in_transit"],
                  depth_by_quartile=vet["depth_by_quartile_ppm"], event_chi2_const=vet["event_chi2_const"], event_dof=vet["event_dof"], known_matches=vet["known_matches"],
                  stage2_pass=vet["stage2_pass"], pixel_summary=px["summary"]),
    mine=dict(P_refined_d=P, P_bootstrap_d=s2["period_bootstrap"]["median_d"], P_bootstrap_err_d=s2["period_bootstrap"]["std_d"], P_blind_top_d=m["blind_bls"]["best_period_d"],
              t0_bkjd=t0, dur_h_bls=m["ephemeris"]["duration_h"], dur_h_boxscan=m["ephemeris"]["duration_h_boxscan"],
              depth_unmasked_1d_ppm=d["psf_depth_ppm"], depth_unmasked_err_ppm=d["psf_depth_err_ppm"], snr_unmasked=d["empirical_snr"],
              depth_masked_1d_ppm=mk["depth_ppm"], depth_masked_err_ppm=mk["err_ppm"], snr_masked=mk["snr"],
              depth_vs_window_unmasked={k: v["depth_ppm"] for k, v in s1["robustness"].items()}, depth_vs_window_masked={k: v["depth_ppm"] for k, v in s2["masked_redetrend"].items()},
              bls_sde_detrended=m["blind_bls"]["sde_detrended"], bls_sde_raw=m["blind_bls"]["sde_raw"], inverted_sde_detrended=m["inverted_null"]["sde_detrended"],
              odd_even_sigma=m["odd_even"]["sigma"], odd_even_sigma_masked=mk["odd_even_sigma"], secondary_ppm=m["secondary"]["depth_ppm"], secondary_sigma=m["secondary"]["sigma"],
              n_events=m["transits"]["n_epochs_with_data"], n_in=d["n_in"], per_year=m["per_year"], per_year_chi2=s1["per_year_consistency"], event_chi2_const=182.4, event_dof=124,
              sap_depth_ppm=m["sap"]["depth_ppm"], sap_err_ppm=m["sap"]["err_ppm"], sap_over_psf=m["sap"]["ratio_sap_over_psf"], sap_over_psf_masked=mk["sap_over_psf"],
              absolute_deficit=s2["absolute_deficit"]),
    deltas=dict(period_sigma=float(abs(s2["period_bootstrap"]["median_d"] - vet["P"]) / s2["period_bootstrap"]["std_d"]),
                t0_diff_min=float((t0 - vet["t0_bkjd"]) * 24 * 60),
                depth_masked_vs_pass2_sigma=float((mk["depth_ppm"] - vet["depth_pass2_ppm"]) / np.hypot(mk["err_ppm"], vet["depth_pass2_formal_err_ppm"])),
                depth_unmasked_vs_pass1_ppm=float(d["psf_depth_ppm"] - vet["depth_pass1_ppm"]), dur_diff_h=float(m["ephemeris"]["duration_h"] - vet["dur_h"])),
    neighbours=dict(kic009655236_tpf_host=dict(G=nb1["header"]["GMAG"], sep_arcsec=nb1["separation_arcsec_from_target"], quarters=nb1["quarters_in_file"],
                                               psf_depth_ppm=nb1["at_target_ephemeris"]["PSF"]["depth_ppm"], psf_err_ppm=nb1["at_target_ephemeris"]["PSF"]["err_ppm"],
                                               psf_deficit_e_s=nb1["at_target_ephemeris"]["PSF"]["deficit_e_s"], psf_deficit_err_e_s=nb1["at_target_ephemeris"]["PSF"]["deficit_err_e_s"],
                                               sap_deficit_e_s=nb1["at_target_ephemeris"]["SAP"]["deficit_e_s"], sap_deficit_err_e_s=nb1["at_target_ephemeris"]["SAP"]["deficit_err_e_s"]),
                    kic009655226_G15p9=dict(G=nb2["header"]["GMAG"], sep_arcsec=nb2["separation_arcsec_from_target"], quarters=nb2["quarters_in_file"], overlap_quarters_with_target=[3, 10],
                                            psf_depth_ppm=nb2["at_target_ephemeris"]["PSF"]["depth_ppm"], psf_err_ppm=nb2["at_target_ephemeris"]["PSF"]["err_ppm"],
                                            psf_deficit_e_s=nb2["at_target_ephemeris"]["PSF"]["deficit_e_s"], psf_deficit_err_e_s=nb2["at_target_ephemeris"]["PSF"]["deficit_err_e_s"],
                                            sap_deficit_e_s=nb2["at_target_ephemeris"]["SAP"]["deficit_e_s"], sap_deficit_err_e_s=nb2["at_target_ephemeris"]["SAP"]["deficit_err_e_s"],
                                            own_narrow_bls_peak=nb_peak),
                    target_deficit_e_s=s2["absolute_deficit"]["psf_e_s"], target_deficit_err_e_s=s2["absolute_deficit"]["psf_err"],
                    kic009655226_psf_deficit_below_target_sigma=float((s2["absolute_deficit"]["psf_e_s"] - nb2["at_target_ephemeris"]["PSF"]["deficit_e_s"]) / np.hypot(s2["absolute_deficit"]["psf_err"], nb2["at_target_ephemeris"]["PSF"]["deficit_err_e_s"]))),
)
res["comparison_with_pipeline"] = cmp
res["files"] = sorted(os.listdir(OUT))
json.dump(res, open(f"{OUT}/rederive_results.json", "w"), indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
print(json.dumps(cmp["deltas"], indent=1)); print("neighbour sigma below target:", cmp["neighbours"]["kic009655226_psf_deficit_below_target_sigma"])
# manifest
with open(f"{OUT}/MANIFEST.sha256", "w") as fh:
    for fn in sorted(os.listdir(OUT)):
        if fn == "MANIFEST.sha256" or fn.startswith("__pycache__"): continue
        p = os.path.join(OUT, fn)
        if os.path.isfile(p):
            fh.write(f"{hashlib.sha256(open(p,'rb').read()).hexdigest()}  {fn}\n")
print(open(f"{OUT}/MANIFEST.sha256").read())
