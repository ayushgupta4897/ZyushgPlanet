#!/usr/bin/env python
"""Per-event transit times (trapezoid template, grid search) -> weighted linear ephemeris.
Gives a period/t0 uncertainty that the box-likelihood (jagged in P) cannot, and a TTV/linear-ephemeris chi2.
Uses the same loader/detrend as rederive_kic006699887.py (imported from that file; no pipeline/ code)."""
import json, os, sys, numpy as np
D = '/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006699887/rederive'
sys.path.insert(0, D)
import importlib.util
spec = importlib.util.spec_from_file_location('rd', os.path.join(D, 'rederive_kic006699887.py'))
# only need the helper functions; avoid executing the full script by reading source up to the 'run' marker
src = open(os.path.join(D, 'rederive_kic006699887.py')).read().split('# run\n')[0]
ns = {}; exec(compile(src, 'rd_helpers', 'exec'), ns)
load, detrend, phase_of, trap = ns['load'], ns['detrend'], ns['phase_of'], ns['trap']
R = json.load(open(os.path.join(D, 'results.json')))
P = R['summary']['period_d']; t0 = R['summary']['t0_bkjd']; dur = R['summary']['duration_h_box'] / 24
tp = R['pass2']['trapezoid']; dep, T14, T23 = tp['depth_ppm'] * 1e-6, tp['T14_h'], tp['T23_h']
t, f, fe, q = load('FLUX')
mask = np.abs(phase_of(t, P, t0)) < 1.0 * dur / P
fl, el, good = detrend(t, f, fe, q, trend_mask=mask)
t, fl, q = t[good], fl[good], q[good]
ep = np.round((t - t0) / P).astype(int)
rows = []
shifts = np.arange(-2.0, 2.0001, 1 / 60)  # hours
for e in np.unique(ep):
    tc0 = t0 + e * P
    w = np.abs(t - tc0) < 0.6  # +-14.4 h window
    if (np.abs(t[w] - tc0) < 0.5 * dur).sum() < 3: continue
    x = (t[w] - tc0) * 24; y = fl[w]
    oot = np.abs(x) > 0.75 * T14
    sig = np.std(y[oot]) if oot.sum() > 10 else R['refined_pass2']['point_scatter_ppm_mad'] * 1e-6
    chi2 = np.array([np.sum((y - trap(x - s, dep, T14, T23))**2) / sig**2 for s in shifts])
    k = int(np.argmin(chi2))
    ok = chi2 <= chi2[k] + 1.0
    # 1-sigma from the delta-chi2=1 interval; require it to be a contiguous interval around k
    lo = shifts[ok].min(); hi = shifts[ok].max()
    err_h = max(0.5 * (hi - lo), 1 / 60)
    chi2_flat = np.sum((y - 1)**2) / sig**2
    rows.append({'epoch': int(e), 'tc_bkjd': float(tc0 + shifts[k] / 24), 'shift_min': float(shifts[k] * 60), 'err_min': float(err_h * 60),
                 'dchi2_vs_flat': float(chi2_flat - chi2[k]), 'n_in_window': int(w.sum())})
E = np.array([r['epoch'] for r in rows]); TC = np.array([r['tc_bkjd'] for r in rows]); S = np.array([r['err_min'] for r in rows]) / 1440
# keep events whose template detection is meaningful (dchi2 > 4) for the fit, report both
keep = np.array([r['dchi2_vs_flat'] > 4 for r in rows])
def wfit(E, TC, S):
    W = 1 / S**2; A = np.vstack([np.ones_like(E), E]).T
    cov = np.linalg.inv(A.T @ (A * W[:, None])); b = cov @ (A.T @ (W * TC))
    resid = TC - A @ b; chi2 = float(np.sum(W * resid**2))
    return b, cov, chi2, resid
b, cov, chi2, resid = wfit(E[keep], TC[keep], S[keep])
scale = np.sqrt(max(chi2 / (keep.sum() - 2), 1.0))
out = {'n_events_timed': int(len(rows)), 'n_used_in_fit(dchi2>4)': int(keep.sum()),
       'P_d': float(b[1]), 'P_err_d': float(np.sqrt(cov[1, 1])), 'P_err_d_scaled_by_sqrt_redchi2': float(np.sqrt(cov[1, 1]) * scale),
       't0_bkjd': float(b[0]), 't0_err_d': float(np.sqrt(cov[0, 0])), 'chi2_linear': chi2, 'dof': int(keep.sum() - 2),
       'rms_resid_min': float(np.std(resid) * 1440), 'median_timing_err_min': float(np.median(S[keep]) * 1440),
       'max_abs_resid_over_err': float(np.max(np.abs(resid) / S[keep])),
       'pipeline_P_offset_sigma': float((15.113611512744262 - b[1]) / (np.sqrt(cov[1, 1]) * scale)),
       'bls_P_offset_sigma': float((P - b[1]) / (np.sqrt(cov[1, 1]) * scale)),
       'events': rows}
json.dump(out, open(os.path.join(D, 'ephemeris.json'), 'w'), indent=2)
for k, v in out.items():
    if k != 'events': print(k, v)
