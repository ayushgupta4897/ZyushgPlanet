#!/usr/bin/env python
"""Compare the independent re-derivation (results.json) with the pipeline's vet.json / events.csv.
Run AFTER rederive_kic006699887.py.  Reads pipeline outputs only (no pipeline code)."""
import json, numpy as np, pandas as pd, os
D = '/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006699887/rederive'
R = json.load(open(os.path.join(D, 'results.json')))
V = json.load(open('/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-006699887/vet.json'))
E = pd.read_csv('/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-006699887/events.csv')
m = R['pass2']; s = R['summary']
cmp = {}
def row(name, mine, pipe, note=''):
    cmp[name] = {'mine': mine, 'pipeline': pipe, 'diff': (mine - pipe) if isinstance(mine, (int, float)) and isinstance(pipe, (int, float)) else None, 'note': note}
row('period_d', s['period_d'], V['P'], f"my dloglike=1 half-width {R['refined_pass2']['period_err_d_dloglike1']:.2e} d")
row('t0_bkjd', s['t0_bkjd'], V['t0_bkjd'])
# mid-transit time at the middle epoch (epoch 48) under each ephemeris
row('tc_epoch48_bkjd', s['t0_bkjd'] + 48 * s['period_d'], V['t0_bkjd'] + 48 * V['P'], 'drift at mid-baseline, days')
row('duration_h_box', s['duration_h_box'], V['dur_h'], 'my box grid 0.25 h; trapezoid T14 = %.3f h, T23/T14 = %.2f' % (m['trapezoid']['T14_h'], m['trapezoid']['T23_over_T14']))
row('depth_ppm_box', s['depth_ppm'], V['depth_pass2_ppm'], 'mine 3.5 h box on masked re-detrend; pipeline 3.75 h box pass2')
# my depth measured with the pipeline's 3.75 h box width (from duration scan)
ds = np.array(m['duration_scan_h_depth_dll']); d375 = float(ds[np.argmin(np.abs(ds[:, 0] - 3.75)), 1])
row('depth_ppm_at_3.75h_box', d375, V['depth_pass2_ppm'], 'my masked-detrend depth evaluated with the pipeline box width')
row('depth_err_ppm', s['depth_err_ppm'], V['box_noise_ppm'], 'mine = std of 3000 random-phase boxes (same P, dur); pipeline box_noise_ppm. Pipeline formal err %.1f' % V['depth_pass2_formal_err_ppm'])
row('empirical_snr', s['empirical_snr'], V['snr_emp_pass2'])
row('odd_even_sigma', s['odd_even_sigma'], V['odd_even_sigma'])
row('secondary_ppm', s['secondary_ppm'], V['secondary_depth_ppm'])
row('secondary_sigma', s['secondary_sigma'], V['secondary_sigma'])
row('n_transits', s['n_transits'], V['n_events'])
row('event_weighted_depth_ppm', m['events_weighted_mean_depth_ppm'], V['event_weighted_depth_ppm'])
row('event_chi2', m['events_chi2'], V['event_chi2_const'], f"dof mine {m['events_dof']} pipeline {V['event_dof']}")
row('trapezoid_or_tls_depth_ppm', m['trapezoid']['depth_ppm'], (1 - V['tls']['depth']) * 1e6, 'my trapezoid full depth vs pipeline TLS depth')
row('tls_period_d', s['period_d'], V['tls']['period'])
row('sde', R['bls_sde_top'], V['tls']['SDE'], 'NOT the same statistic: mine = SDE of BLS log-likelihood power (skewed); pipeline = TLS SDE. My depth-SNR-spectrum SDE = %.1f' % R.get('bls_sde_depthsnr_spectrum', float('nan')))
cmp['per_year_mine_ppm'] = s['per_year_depths_ppm']
cmp['per_quartile_pipeline_ppm'] = V['depth_by_quartile_ppm'][0::3]
cmp['n_in_transit'] = {'mine': m['box']['n_in_transit_points'], 'pipeline': V['n_in_transit'], 'note': 'different box width (3.5 vs 3.75 h)'}
# per-event comparison: match by epoch number
mine_ev = pd.DataFrame(m['events'])
pe = E.copy(); pe['depth_ppm'] = pe['depth'] * 1e6; pe['err_ppm'] = pe['err'] * 1e6
mine_ev['n'] = np.round((mine_ev['tc_bkjd'] - V['t0_bkjd']) / V['P']).astype(int)
j = mine_ev.merge(pe, on='n', suffixes=('_mine', '_pipe'))
cmp['per_event'] = {'n_matched': int(len(j)), 'pearson_r_depth': float(np.corrcoef(j['depth_ppm_mine'], j['depth_ppm_pipe'])[0, 1]),
                    'median_abs_diff_ppm': float(np.median(np.abs(j['depth_ppm_mine'] - j['depth_ppm_pipe']))),
                    'mean_diff_mine_minus_pipe_ppm': float(np.mean(j['depth_ppm_mine'] - j['depth_ppm_pipe'])),
                    'max_abs_tc_diff_min': float(np.max(np.abs(j['tc_bkjd_mine'] - j['tc_bkjd_pipe'])) * 1440)}
# dilution
cmp['sap_vs_psf'] = {'psf_depth_ppm': s['depth_ppm'], 'sap_depth_ppm': s['sap_depth_ppm'], 'ratio': s['sap_over_psf'],
                     'ratio_err': R['sap']['pass2']['ratio_err'], 'crowdsap_catalogue': 0.8849, 'psffrac_catalogue': 0.7328}
json.dump(cmp, open(os.path.join(D, 'compare.json'), 'w'), indent=2)
for k, v in cmp.items(): print(k, v)
