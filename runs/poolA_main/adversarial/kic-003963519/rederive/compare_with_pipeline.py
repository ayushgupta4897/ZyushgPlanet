#!/usr/bin/env python
"""Compare the independent re-derivation (rederive_results.json) with the pipeline's vet.json.
Run AFTER rederive_kic003963519.py; this is the only file here that reads pipeline output."""
import json, numpy as np
D = '/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-003963519/rederive/'
R = json.load(open(D + 'rederive_results.json'))
V = json.load(open('/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-003963519/vet.json'))
X = json.load(open('/home/user/ZyushgPlanet/runs/poolA_main/pixel/kic-003963519/pixel.json'))
P, Pe, t0 = R['ephemeris']['period_d'], R['ephemeris']['period_err_d'], R['ephemeris']['t0_bkjd']
Pv, t0v = V['P'], V['t0_bkjd']
n = np.round((t0 - t0v) / P)
dt0_min = ((t0 - n * P) - t0v) * 24 * 60            # my t0 propagated back to the pipeline epoch
dP_sigma = (P - Pv) / Pe
d, de = R['depth']['box_depth_ppm'], R['depth']['depth_err_empirical_ppm']
dv, dve = V['depth_pass2_ppm'], V['box_noise_ppm']
ev = R['transits']['events']
n_ge1 = sum(1 for e in ev if e['n_in'] >= 1); n_ge3 = sum(1 for e in ev if e['n_in'] >= 3)
rstar = 0.7730915021896362                            # targets CSV Rstar (Rsun), same value vet.json used
rp_re = np.sqrt(d * 1e-6) * rstar * 109.08
sap, sape = R['sap']['box_depth_ppm'], R['sap']['depth_err_empirical_ppm']
ratio = sap / d; ratio_err = ratio * np.sqrt((sape / sap) ** 2 + (de / d) ** 2)
py = R['per_year']; pyd = np.array([p['depth_ppm'] for p in py]); pye = np.array([p['err_ppm'] for p in py])
wm = np.sum(pyd / pye ** 2) / np.sum(1 / pye ** 2); chi2_year = float(np.sum(((pyd - wm) / pye) ** 2))
C = dict(
    period=dict(mine_d=P, mine_err_d=Pe, pipeline_d=Pv, diff_d=P - Pv, diff_sigma=float(dP_sigma)),
    t0=dict(mine_bkjd=t0, pipeline_bkjd=t0v, epochs_between=int(n), diff_at_pipeline_epoch_min=float(dt0_min)),
    depth=dict(mine_box_ppm=d, mine_err_empirical_ppm=de, mine_trapezoid_ppm=R['ephemeris']['trapezoid']['depth'] * 1e6,
               pipeline_pass2_ppm=dv, pipeline_box_noise_ppm=dve, pipeline_pass1_ppm=V['depth_pass1_ppm'],
               diff_ppm=d - dv, diff_sigma=float((d - dv) / np.hypot(de, dve))),
    snr=dict(mine_empirical=R['depth']['empirical_snr'], pipeline_snr_emp_pass2=V['snr_emp_pass2']),
    duration=dict(mine_bls_h=R['ephemeris']['duration_bls_h'], mine_trapezoid_t14_h=R['ephemeris']['duration_trapezoid_t14_h'],
                  pipeline_h=V['dur_h'], pipeline_tls_h=V['tls']['duration'] * 24),
    odd_even=dict(mine_sigma=R['odd_even']['sigma'], pipeline_sigma=V['odd_even_sigma'],
                  mine_odd_even_ppm=[R['odd_even']['odd_ppm'], R['odd_even']['even_ppm']],
                  pipeline_odd_even_ppm=[V['depth_odd_ppm'], V['depth_even_ppm']]),
    secondary=dict(mine_ppm=R['secondary']['depth_ppm'], mine_sigma=R['secondary']['sigma'],
                   pipeline_ppm=V['secondary_depth_ppm'], pipeline_sigma=V['secondary_sigma']),
    n_transits=dict(mine_n_in_ge3=n_ge3, mine_n_in_ge1=n_ge1, mine_epochs_in_baseline=R['transits']['n_epochs_in_baseline'],
                    pipeline_n_events=V['n_events'], pipeline_tls_distinct=V['tls']['distinct_transit_count']),
    event_consistency=dict(mine_chi2_const=R['transits']['chi2_constant_depth'], mine_dof=R['transits']['dof'],
                           pipeline_chi2_const=V['event_chi2_const'], pipeline_dof=V['event_dof']),
    per_year_vs_quartile=dict(mine_per_year_ppm=pyd.tolist(), mine_per_year_err_ppm=pye.tolist(),
                              mine_chi2_const_3dof=chi2_year,
                              pipeline_per_quartile_ppm=[q[0] for q in V['depth_by_quartile_ppm']],
                              pipeline_per_quartile_err_ppm=[q[1] for q in V['depth_by_quartile_ppm']]),
    rp=dict(mine_rp_earth_from_box=float(rp_re), rstar_rsun=rstar, pipeline_rp_earth_box=V['rp_earth_box'], pipeline_tls_rp_earth=V['tls']['rp_earth']),
    sap_vs_psf=dict(ratio=float(ratio), ratio_err=float(ratio_err), catalogue_crowdsap=0.6497, catalogue_psffrac=0.7731, npixsap=2,
                    expected_ratio_if_sap_35pct_contaminated_and_psf_clean=0.65,
                    sigma_from_0p65=float((ratio - 0.65) / ratio_err), sigma_from_1=float((ratio - 1) / ratio_err)),
    sde=dict(mine_bls_1_50d_sde_raw=R['bls']['sde_raw'], mine_bls_1_50d_sde_baseline_subtracted=R['bls']['sde_baseline_subtracted'],
             mine_inverted_null_sde_raw=R['inverted_null']['sde_raw'], mine_inverted_null_sde_baseline_subtracted=R['inverted_null']['sde_detrended'],
             pipeline_tls_narrow_window_sde=V['tls']['SDE'], note='different statistics; not directly comparable'),
    pixel_summary_from_pipeline=X['summary'],
)
json.dump(C, open(D + 'comparison_with_pipeline.json', 'w'), indent=2)
print(json.dumps(C, indent=1))
