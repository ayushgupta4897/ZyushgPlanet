#!/usr/bin/env python
"""
Independent re-derivation for kic-006699887 (KBONUS-BKG PSF light curve).
Written from scratch with numpy/scipy/pandas/astropy only -- no code from
/pipeline/ is imported or consulted. No wotan, no TLS.

Steps
 1. load LIGHTCURVE_STITCHED, SAP_QUALITY==0, finite; normalise per quarter;
    detrend per quarter with a 1.0-d centred running median (time-based),
    then clip only POSITIVE outliers (>4 sigma robust) so dips survive.
 2. blind BLS 1-50 d on a log-uniform period grid (own grid), durations 1.5-8 h.
    Top-3 peaks reported after suppressing +-2% neighbours of each found peak.
    Identical run on inverted flux (2-f) as a null.
 3. refine P, t0, duration near the top peak.  Measurements are done twice:
    pass1 on the blind (unmasked) detrend, pass2 after a transit-masked
    re-detrend (in-transit cadences excluded from the running median; this
    removes the known running-median bias that produces positive shoulders
    around a transit and under-estimates depth).  Each pass: box depth with
    empirical error from random-phase boxes; odd/even; secondary at phase 0.5;
    per-year depths; per-event depths; number of transits; trapezoid shape.
 4. repeat depth measurement on SAP_FLUX at the same ephemeris (dilution check).
Outputs: results.json, fold.png, bls_periodogram.npz next to this script.
"""
import json, os, sys, time
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RNG = np.random.default_rng(20261017)
FNAME = 'kic-006699887'
FITS = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-006699887_kepler_v1.0_lc.fits'
OUT = '/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006699887/rederive'
os.makedirs(OUT, exist_ok=True)
PIPE_P = 15.11360   # pipeline value, used only to label aliases

# ----------------------------------------------------------------------------
# 1. load + detrend
# ----------------------------------------------------------------------------
def load(fluxcol):
    d = fits.getdata(FITS, 'LIGHTCURVE_STITCHED')
    t = d['TIME'].astype(float)
    f = d[fluxcol].astype(float)
    fe = d[fluxcol + '_ERR'].astype(float)
    q = d['QUARTER'].astype(int)
    qual = d['SAP_QUALITY'].astype(int)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (f > 0) & (fe > 0)
    o = np.argsort(t[m])
    return t[m][o], f[m][o], fe[m][o], q[m][o]

def detrend(t, f, fe, q, window_d=1.0, clip_sigma=4.0, trend_mask=None):
    """per-quarter: divide by median, divide by a 1-d centred running median (time-based).
    trend_mask: boolean array, True = exclude this cadence from the trend estimate
    (used in pass 2 to mask transits).  Positive outliers are masked from the trend
    and dropped.  Returns normalised flux, err, keep-mask."""
    fl = np.full_like(f, np.nan); el = np.full_like(f, np.nan)
    keep = np.ones(len(t), bool)
    if trend_mask is None: trend_mask = np.zeros(len(t), bool)
    for qq in np.unique(q):
        s = np.where(q == qq)[0]
        med = np.median(f[s][~trend_mask[s]])
        ff = f[s] / med; ee = fe[s] / med
        idx = pd.to_timedelta(t[s], unit='D')
        def rmed(y):
            return pd.Series(y, index=idx).rolling(f'{window_d}D', center=True, min_periods=5).median().values
        y1 = ff.copy(); y1[trend_mask[s]] = np.nan
        tr = rmed(y1)
        tr = np.where(np.isfinite(tr), tr, rmed(ff))
        r = ff / tr
        mad = 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r)))
        up = r > 1 + clip_sigma * mad
        y2 = y1.copy(); y2[up] = np.nan
        tr2 = rmed(y2)
        tr2 = np.where(np.isfinite(tr2), tr2, tr)
        fl[s] = ff / tr2; el[s] = ee / tr2
        keep[s[up]] = False
    good = keep & np.isfinite(fl) & np.isfinite(el)
    return fl, el, good

# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def phase_of(t, P, t0):
    return ((t - t0) / P + 0.5) % 1.0 - 0.5

def box_depth(t, f, P, t0, dur, phase_offset=0.0):
    ph = phase_of(t, P, t0) - phase_offset
    ph = (ph + 0.5) % 1.0 - 0.5
    inb = np.abs(ph) < 0.5 * dur / P
    if inb.sum() < 3: return np.nan, 0
    return 1.0 - np.mean(f[inb]), int(inb.sum())

def empirical_depth(t, f, P, t0, dur, nrand=3000, exclude_cycles=0.1):
    """depth and empirical sigma from random-phase boxes with the same P, dur, data;
    phases within +-exclude_cycles of transit and +-exclude_cycles/2 of phase 0.5 avoided."""
    dep, nin = box_depth(t, f, P, t0, dur)
    deps = []
    while len(deps) < nrand:
        po = RNG.uniform(-0.5, 0.5)
        if abs(po) < exclude_cycles or abs(abs(po) - 0.5) < exclude_cycles / 2: continue
        dd, nn = box_depth(t, f, P, t0, dur, phase_offset=po)
        if np.isfinite(dd) and nn >= 0.5 * nin: deps.append(dd)
    deps = np.array(deps)
    return dep, nin, float(np.std(deps)), float(np.mean(deps)), deps

def bls_top_peaks(P, power, nkeep=3, exclude_frac=0.02):
    pw = power.copy(); peaks = []
    for _ in range(nkeep):
        i = int(np.nanargmax(pw)); peaks.append(i)
        pw[np.abs(P / P[i] - 1) < exclude_frac] = -np.inf
    return peaks

def sde(power):
    return float((np.nanmax(power) - np.nanmean(power)) / np.nanstd(power))

def alias_label(P, Pref):
    for n in range(1, 7):
        for m in range(1, 7):
            if abs(P / Pref - n / m) < 0.003:
                return 'same' if n == m else f'{n}/{m}'
    return 'unrelated'

def trap(x, dep, T14, T23):
    y = np.ones_like(x); a = np.abs(x)
    full = a <= T23 / 2; ing = (a > T23 / 2) & (a < T14 / 2)
    y[full] = 1 - dep
    y[ing] = 1 - dep * (T14 / 2 - a[ing]) / ((T14 - T23) / 2 + 1e-9)
    return y

def measure(t, fl, q, P, t0, dur, sig_pt, label):
    """all ephemeris-level measurements on a detrended light curve."""
    out = {'label': label}
    dep, nin, e_emp, mean_rand, rand_deps = empirical_depth(t, fl, P, t0, dur, nrand=3000)
    ph = phase_of(t, P, t0)
    oot = np.abs(ph) > 1.5 * dur / P
    e_formal = float(np.std(fl[oot]) / np.sqrt(nin))
    out['box'] = {'depth_ppm': dep * 1e6, 'n_in_transit_points': nin, 'depth_err_empirical_ppm': e_emp * 1e6,
                  'depth_err_formal_ppm': e_formal * 1e6, 'random_box_mean_ppm': mean_rand * 1e6,
                  'empirical_snr': dep / e_emp, 'formal_snr': dep / e_formal, 'n_random_boxes': int(len(rand_deps)),
                  'frac_random_boxes_deeper': float(np.mean(rand_deps >= dep)),
                  'oot_level_within_12h_ppm': float((np.mean(fl[(np.abs(ph) * P * 24 < 12) & (np.abs(ph) * P * 24 > 0.75 * dur * 24)]) - 1) * 1e6)}
    # duration scan (box) at fixed P,t0
    dscan = []
    for dd in np.arange(1.0, 10.01, 0.25) / 24:
        inb_ = np.abs(ph) < 0.5 * dd / P
        if inb_.sum() < 3: continue
        dpt = 1 - np.mean(fl[inb_])
        dscan.append((dd * 24, dpt * 1e6, 0.5 * dpt**2 * inb_.sum() / sig_pt**2))
    dscan = np.array(dscan)
    out['duration_scan_h_depth_dll'] = dscan.tolist()
    out['duration_best_from_scan_h'] = float(dscan[np.argmax(dscan[:, 2]), 0])
    # odd/even
    ep = np.round((t - t0) / P).astype(int)
    inb = np.abs(ph) < 0.5 * dur / P
    odd = inb & (ep % 2 == 1); even = inb & (ep % 2 == 0)
    d_odd = 1 - np.mean(fl[odd]); d_even = 1 - np.mean(fl[even])
    sig_eff = e_emp * np.sqrt(nin)
    e_odd = sig_eff / np.sqrt(odd.sum()); e_even = sig_eff / np.sqrt(even.sum())
    out['odd_even'] = {'depth_odd_ppm': d_odd * 1e6, 'err_odd_ppm': e_odd * 1e6, 'n_odd': int(odd.sum()),
                       'depth_even_ppm': d_even * 1e6, 'err_even_ppm': e_even * 1e6, 'n_even': int(even.sum()),
                       'odd_even_sigma': float(abs(d_odd - d_even) / np.hypot(e_odd, e_even))}
    # secondary
    d_sec, n_sec = box_depth(t, fl, P, t0, dur, phase_offset=0.5)
    sec_scan = np.array([(po, box_depth(t, fl, P, t0, dur, phase_offset=po)[0] * 1e6) for po in np.arange(0.15, 0.851, 0.005)])
    kk = np.argmax(sec_scan[:, 1])
    out['secondary'] = {'depth_ppm': d_sec * 1e6, 'n_points': n_sec, 'err_ppm': e_emp * 1e6, 'sigma': d_sec / e_emp,
                        'max_secondary_anywhere_phase0.15-0.85': {'phase': float(sec_scan[kk, 0]), 'depth_ppm': float(sec_scan[kk, 1]), 'sigma': float(sec_scan[kk, 1] * 1e-6 / e_emp)}}
    # per-year
    years = {'Y1_Q1-4': (1, 4), 'Y2_Q5-8': (5, 8), 'Y3_Q9-12': (9, 12), 'Y4_Q13-17': (13, 17)}
    py = {}
    for name, (qa, qb) in years.items():
        s = inb & (q >= qa) & (q <= qb)
        if s.sum() < 3: py[name] = {'depth_ppm': None}; continue
        py[name] = {'depth_ppm': (1 - np.mean(fl[s])) * 1e6, 'err_ppm': sig_eff / np.sqrt(s.sum()) * 1e6, 'n_points': int(s.sum()),
                    'n_transits_with_data': int(len(np.unique(ep[s]))), 't_range_bkjd': [float(t[s].min()), float(t[s].max())]}
    out['per_year'] = py
    vals = np.array([v['depth_ppm'] for v in py.values() if v['depth_ppm'] is not None])
    errs = np.array([v['err_ppm'] for v in py.values() if v['depth_ppm'] is not None])
    out['per_year_chi2'] = float(np.sum(((vals - np.average(vals, weights=1 / errs**2)) / errs)**2)); out['per_year_dof'] = int(len(vals) - 1)
    # per-event
    expected_pts = dur / np.median(np.diff(t))
    events = []
    for e in np.unique(ep[inb]):
        s = inb & (ep == e); tc = t0 + e * P
        loc = (np.abs(t - tc) < 1.0) & (np.abs(t - tc) > 0.75 * dur)
        de = 1 - np.mean(fl[s]); ee = (np.std(fl[loc]) if loc.sum() > 10 else sig_pt) / np.sqrt(s.sum())
        # distance to nearest data gap > 0.5 d
        dt = np.diff(t); gaps = np.where(dt > 0.5)[0]
        gap_edges = np.concatenate([[t.min()], t[gaps], t[gaps + 1], [t.max()]]) if len(gaps) else np.array([t.min(), t.max()])
        events.append({'epoch': int(e), 'tc_bkjd': float(tc), 'n_in': int(s.sum()), 'coverage': float(s.sum() / expected_pts),
                       'depth_ppm': float(de * 1e6), 'err_ppm': float(ee * 1e6), 'n_local_oot': int(loc.sum()),
                       'quarter': int(np.median(q[s])), 'nearest_gap_edge_d': float(np.min(np.abs(gap_edges - tc)))})
    out['events'] = events
    cov = np.array([e['coverage'] for e in events])
    out['n_transits_any_data'] = int(len(events)); out['n_transits_cov_ge_0p5'] = int(np.sum(cov >= 0.5))
    out['n_epochs_in_baseline'] = int(np.floor((t.max() - t0) / P)) + 1
    ed = np.array([e['depth_ppm'] for e in events]); ee_ = np.array([e['err_ppm'] for e in events]); w = 1 / ee_**2
    out['events_weighted_mean_depth_ppm'] = float(np.sum(w * ed) / np.sum(w))
    out['events_chi2'] = float(np.sum(w * (ed - out['events_weighted_mean_depth_ppm'])**2)); out['events_dof'] = int(len(ed) - 1)
    contrib = (np.clip(ed, 0, None) / ee_)**2
    out['single_event_dominance_frac'] = float(contrib.max() / contrib.sum())
    out['deepest_event'] = events[int(np.argmax(ed))]
    out['n_events_depth_gt_2sigma'] = int(np.sum(ed / ee_ > 2))
    s_lo = inb & (ep != out['deepest_event']['epoch'])
    out['depth_without_deepest_event_ppm'] = float((1 - np.mean(fl[s_lo])) * 1e6)
    # trapezoid grid fit within +-12 h
    ph_h = ph * P * 24; w_ = np.abs(ph_h) < 12
    xx, yy = ph_h[w_], fl[w_]
    best = (np.inf, None)
    for T14 in np.arange(1.5, 8.01, 0.125):
        for frac in np.arange(0.0, 1.001, 0.05):
            T23 = T14 * frac
            shape = 1 - trap(xx, 1.0, T14, T23)
            ss = np.sum(shape**2)
            dep_hat = np.sum(shape * (1 - yy)) / ss if ss > 0 else 0
            chi2 = np.sum((yy - (1 - dep_hat * shape))**2)
            if chi2 < best[0]: best = (chi2, (dep_hat, T14, T23))
    chi2_flat = np.sum((yy - 1)**2)
    out['trapezoid'] = {'depth_ppm': best[1][0] * 1e6, 'T14_h': best[1][1], 'T23_h': best[1][2], 'T23_over_T14': best[1][2] / best[1][1],
                        'delta_chi2_vs_flat_sigma_units': float((chi2_flat - best[0]) / sig_pt**2)}
    out['_trap_params'] = best[1]
    return out

# ----------------------------------------------------------------------------
# run
# ----------------------------------------------------------------------------
res = {'fname': FNAME, 'date': '2026-10-17', 'script': os.path.abspath(__file__),
       'notes': ('independent re-derivation; no pipeline/ code used; detrend = per-quarter 1.0-d centred running median '
                 '(pandas time-based), positive-only 4-sigma clipping. pass1 = blind detrend; pass2 = re-detrend with '
                 'in-transit cadences (|phase|<1.0*duration) excluded from the trend. Headline numbers are pass2.')}

t, f, fe, q = load('FLUX')
fl, el, good = detrend(t, f, fe, q)
t, fl, el, q, f, fe = t[good], fl[good], el[good], q[good], f[good], fe[good]
res['n_cadences_used_psf'] = int(len(t)); res['baseline_d'] = float(t.max() - t.min())
res['point_scatter_ppm_std'] = float(np.std(fl) * 1e6)
res['point_scatter_ppm_mad'] = float(1.4826 * np.median(np.abs(fl - np.median(fl))) * 1e6)
sig = res['point_scatter_ppm_mad'] * 1e-6
print('N', len(t), 'baseline', res['baseline_d'], 'scatter std/MAD ppm', res['point_scatter_ppm_std'], res['point_scatter_ppm_mad'])

# ---- 2. blind BLS 1-50 d (own log-uniform grid: dlnP = dmin/(3*T))
Pmin, Pmax = 1.0, 50.0; T = res['baseline_d']; dmin = 2.0 / 24; q_os = 3.0
dlnP = dmin / (q_os * T)
P_grid = np.exp(np.arange(np.log(Pmin), np.log(Pmax), dlnP))
durs = np.array([1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0]) / 24
res['bls_grid'] = {'Pmin': Pmin, 'Pmax': Pmax, 'n_periods': int(len(P_grid)), 'dlnP': dlnP, 'durations_h': (durs * 24).tolist(), 'objective': 'likelihood'}
print('BLS grid', len(P_grid), 'periods')
bls = BoxLeastSquares(t, fl, el)
t1 = time.time(); r = bls.power(P_grid, durs, objective='likelihood'); print('BLS done in', time.time() - t1, 's')
power = np.asarray(r.power); snr_spec = np.asarray(r.depth) / np.asarray(r.depth_err)
peaks = bls_top_peaks(P_grid, power, 3)
res['bls_top3_power'] = [{'period_d': float(P_grid[i]), 'power': float(power[i]), 'sde_local': float((power[i] - np.nanmean(power)) / np.nanstd(power)),
                          'depth_ppm': float(r.depth[i] * 1e6), 'depth_snr': float(snr_spec[i]), 'duration_h': float(r.duration[i] * 24),
                          't0_bkjd': float(r.transit_time[i]), 'relation_to_pipeline_P': alias_label(P_grid[i], PIPE_P)} for i in peaks]
res['bls_sde_top'] = sde(power)
peaks_snr = bls_top_peaks(P_grid, snr_spec, 3)
res['bls_top3_depthsnr'] = [{'period_d': float(P_grid[i]), 'depth_snr': float(snr_spec[i]), 'depth_ppm': float(r.depth[i] * 1e6),
                             'relation_to_pipeline_P': alias_label(P_grid[i], PIPE_P)} for i in peaks_snr]
best_i = peaks[0]; P_best0 = float(P_grid[best_i])
res['recovered_blind'] = bool(abs(P_best0 / PIPE_P - 1) < 0.002)
res['recovered_blind_alias'] = alias_label(P_best0, PIPE_P)
for d_ in res['bls_top3_power']: print('peak', d_)
print('SDE', res['bls_sde_top'])

# ---- inverted-flux null on the blind detrend, same grid
rinv = BoxLeastSquares(t, 2.0 - fl, el).power(P_grid, durs, objective='likelihood')
pinv = np.asarray(rinv.power); j = int(np.nanargmax(pinv))
res['inverted_null_pass1'] = {'sde': sde(pinv), 'best_period_d': float(P_grid[j]), 'depth_ppm': float(rinv.depth[j] * 1e6),
                              'depth_snr': float(rinv.depth[j] / rinv.depth_err[j]), 'max_power': float(pinv[j]), 'real_max_power': float(power[best_i]),
                              'relation_to_pipeline_P': alias_label(P_grid[j], PIPE_P)}
print('inverted null (pass1)', res['inverted_null_pass1'])

# ---- 3. refine around the best peak
P_fine = np.linspace(P_best0 * (1 - 0.003), P_best0 * (1 + 0.003), 4001)
durs_fine = np.arange(1.0, 8.01, 0.25) / 24
rf = bls.power(P_fine, durs_fine, objective='likelihood'); k = int(np.nanargmax(rf.power))
P_best = float(P_fine[k])
P_fine2 = np.linspace(P_best * (1 - 2e-4), P_best * (1 + 2e-4), 2001)
rf2 = bls.power(P_fine2, durs_fine, objective='likelihood'); k2 = int(np.nanargmax(rf2.power))
P_best = float(P_fine2[k2]); dur_best = float(rf2.duration[k2]); t0_best = float(rf2.transit_time[k2])
n0 = np.ceil((t.min() - t0_best) / P_best); t0_first = t0_best + n0 * P_best
pf = np.asarray(rf.power); within = P_fine[pf > pf.max() - 1.0]
st = bls.compute_stats(P_best, dur_best, t0_best)
res['refined_pass1'] = {'period_d': P_best, 'period_err_d_dloglike1': float(0.5 * (within.max() - within.min())),
                        't0_bkjd': t0_first, 't0_bjd': t0_first + 2454833.0, 'duration_h_bls_box': dur_best * 24,
                        'astropy_stats_depth_ppm': float(st['depth'][0] * 1e6), 'astropy_stats_depth_err_ppm': float(st['depth'][1] * 1e6),
                        'astropy_depth_odd_ppm': float(st['depth_odd'][0] * 1e6), 'astropy_depth_even_ppm': float(st['depth_even'][0] * 1e6),
                        'astropy_depth_half_ppm': float(st['depth_half'][0] * 1e6), 'astropy_depth_half_err_ppm': float(st['depth_half'][1] * 1e6),
                        'astropy_harmonic_delta_loglike': float(st['harmonic_delta_log_likelihood'])}
print('refined (pass1)', P_best, t0_first, dur_best * 24, 'h')

m1 = measure(t, fl, q, P_best, t0_first, dur_best, sig, 'pass1_blind_detrend')
print('pass1 box', m1['box']['depth_ppm'], '+-', m1['box']['depth_err_empirical_ppm'], 'oot level +-12h', m1['box']['oot_level_within_12h_ppm'])

# ---- pass 2: transit-masked re-detrend on the raw flux, then re-refine and re-measure
tr_mask = np.abs(phase_of(t, P_best, t0_first)) < 1.0 * dur_best / P_best   # +-1 duration around mid-transit
fl2, el2, good2 = detrend(t, f, fe, q, trend_mask=tr_mask)
t2, fl2, el2, q2 = t[good2], fl2[good2], el2[good2], q[good2]
sig2 = 1.4826 * np.median(np.abs(fl2 - np.median(fl2)))
bls2 = BoxLeastSquares(t2, fl2, el2)
P_fine3 = np.linspace(P_best * (1 - 3e-4), P_best * (1 + 3e-4), 3001)
rf3 = bls2.power(P_fine3, durs_fine, objective='likelihood'); k3 = int(np.nanargmax(rf3.power))
P2 = float(P_fine3[k3]); dur2 = float(rf3.duration[k3]); t02 = float(rf3.transit_time[k3])
n0 = np.ceil((t2.min() - t02) / P2); t02 = t02 + n0 * P2
pf3 = np.asarray(rf3.power); within3 = P_fine3[pf3 > pf3.max() - 1.0]
st2 = bls2.compute_stats(P2, dur2, t02)
res['refined_pass2'] = {'period_d': P2, 'period_err_d_dloglike1': float(0.5 * (within3.max() - within3.min())),
                        't0_bkjd': t02, 't0_bjd': t02 + 2454833.0, 'duration_h_bls_box': dur2 * 24,
                        'astropy_stats_depth_ppm': float(st2['depth'][0] * 1e6), 'astropy_stats_depth_err_ppm': float(st2['depth'][1] * 1e6),
                        'astropy_depth_half_ppm': float(st2['depth_half'][0] * 1e6), 'astropy_depth_half_err_ppm': float(st2['depth_half'][1] * 1e6),
                        'n_cadences': int(len(t2)), 'point_scatter_ppm_mad': float(sig2 * 1e6)}
print('refined (pass2)', P2, t02, dur2 * 24, 'h')
m2 = measure(t2, fl2, q2, P2, t02, dur2, sig2, 'pass2_transit_masked_detrend')
print('pass2 box', m2['box']['depth_ppm'], '+-', m2['box']['depth_err_empirical_ppm'], 'oot level +-12h', m2['box']['oot_level_within_12h_ppm'])
print('pass2 odd/even', m2['odd_even']); print('pass2 secondary', m2['secondary']); print('pass2 per-year', m2['per_year'])
print('pass2 trapezoid', m2['trapezoid'], 'n events', m2['n_transits_cov_ge_0p5'], 'dominance', m2['single_event_dominance_frac'])

# inverted null on the masked detrend (post-hoc diagnostic: shows whether the pass-1 inverted peak was a detrending artefact)
rinv2 = BoxLeastSquares(t2, 2.0 - fl2, el2).power(P_grid, durs, objective='likelihood')
pinv2 = np.asarray(rinv2.power); j2 = int(np.nanargmax(pinv2))
res['inverted_null_pass2_posthoc'] = {'sde': sde(pinv2), 'best_period_d': float(P_grid[j2]), 'depth_ppm': float(rinv2.depth[j2] * 1e6),
                                      'depth_snr': float(rinv2.depth[j2] / rinv2.depth_err[j2]), 'relation_to_pipeline_P': alias_label(P_grid[j2], PIPE_P),
                                      'power_at_transit_period': float(np.max(pinv2[np.abs(P_grid / P2 - 1) < 1e-3]))}
print('inverted null (pass2, post-hoc)', res['inverted_null_pass2_posthoc'])

# ---- 4. SAP flux at the same ephemeris, both passes
ts, fs, fes, qs = load('SAP_FLUX')
fls1, els1, g1 = detrend(ts, fs, fes, qs)
tsm = np.abs(phase_of(ts, P2, t02)) < 1.0 * dur2 / P2
fls2, els2, g2 = detrend(ts, fs, fes, qs, trend_mask=tsm)
sap = {}
for lab, (tt, ff_, gg) in {'pass1': (ts, fls1, g1), 'pass2': (ts, fls2, g2)}.items():
    tt, ff_ = tt[gg], ff_[gg]
    dep_s, nin_s, e_s, mr_s, _ = empirical_depth(tt, ff_, P2, t02, dur2, nrand=3000)
    ref = m1 if lab == 'pass1' else m2
    dp, ep_ = ref['box']['depth_ppm'] * 1e-6, ref['box']['depth_err_empirical_ppm'] * 1e-6
    sap[lab] = {'n_cadences_used': int(len(tt)), 'point_scatter_ppm_mad': float(1.4826 * np.median(np.abs(ff_ - np.median(ff_))) * 1e6),
                'depth_ppm': dep_s * 1e6, 'depth_err_empirical_ppm': e_s * 1e6, 'n_in_transit_points': nin_s, 'empirical_snr': dep_s / e_s,
                'ratio_sap_over_psf_depth': dep_s / dp, 'ratio_err': (dep_s / dp) * float(np.hypot(e_s / dep_s, ep_ / dp))}
rs = BoxLeastSquares(ts[g2], fls2[g2], els2[g2]).power(P_fine3, durs_fine, objective='likelihood'); ks = int(np.nanargmax(rs.power))
sap['narrow_bls_pass2'] = {'best_period_d': float(P_fine3[ks]), 'depth_ppm': float(rs.depth[ks] * 1e6), 'depth_snr': float(rs.depth[ks] / rs.depth_err[ks]), 'duration_h': float(rs.duration[ks] * 24)}
sap['median_sap_flux_e_s'] = float(np.median(fs)); sap['median_psf_flux_e_s'] = float(np.median(f))
sap['psf_over_sap_flux'] = float(np.median(f) / np.median(fs))
res['sap'] = sap
print('SAP', sap)

res['pass1'] = {k: v for k, v in m1.items() if not k.startswith('_')}
res['pass2'] = {k: v for k, v in m2.items() if not k.startswith('_')}
res['summary'] = {'period_d': P2, 'period_err_d': res['refined_pass2']['period_err_d_dloglike1'], 't0_bjd': t02 + 2454833.0, 't0_bkjd': t02,
                  'duration_h_box': dur2 * 24, 'duration_h_trapezoid_T14': m2['trapezoid']['T14_h'],
                  'depth_ppm': m2['box']['depth_ppm'], 'depth_err_ppm': m2['box']['depth_err_empirical_ppm'],
                  'empirical_snr': m2['box']['empirical_snr'], 'sde': res['bls_sde_top'], 'n_transits': m2['n_transits_cov_ge_0p5'],
                  'odd_even_sigma': m2['odd_even']['odd_even_sigma'], 'secondary_ppm': m2['secondary']['depth_ppm'], 'secondary_sigma': m2['secondary']['sigma'],
                  'per_year_depths_ppm': [v['depth_ppm'] for v in m2['per_year'].values()],
                  'sap_depth_ppm': sap['pass2']['depth_ppm'], 'sap_over_psf': sap['pass2']['ratio_sap_over_psf_depth'],
                  'pass1_depth_ppm_biased_low': m1['box']['depth_ppm']}

def _clean(o):
    if isinstance(o, dict): return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)): return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    return o
with open(os.path.join(OUT, 'results.json'), 'w') as fh: json.dump(_clean(res), fh, indent=2)
np.savez_compressed(os.path.join(OUT, 'bls_periodogram.npz'), period=P_grid, power=power, depth_snr=snr_spec, power_inverted_pass1=pinv, power_inverted_pass2=pinv2)

# ----------------------------------------------------------------------------
# figure (pass2 folds; pass1 overlaid where instructive)
# ----------------------------------------------------------------------------
def binned(x, y, nb):
    edges = np.linspace(x.min(), x.max(), nb + 1); idx = np.digitize(x, edges) - 1
    xb = np.array([x[idx == i].mean() if np.any(idx == i) else np.nan for i in range(nb)])
    yb = np.array([y[idx == i].mean() if np.any(idx == i) else np.nan for i in range(nb)])
    eb = np.array([y[idx == i].std() / np.sqrt(max((idx == i).sum(), 1)) if np.any(idx == i) else np.nan for i in range(nb)])
    return xb, yb, eb

fig, axs = plt.subplots(3, 2, figsize=(13, 12))
ax = axs[0, 0]
ax.plot(P_grid, power, lw=0.4, color='0.3', label='PSF flux (blind detrend)')
ax.plot(P_grid, pinv, lw=0.4, color='tab:red', alpha=0.5, label='inverted 2-f (blind detrend)')
ax.axvline(P2, color='tab:blue', ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_xlabel('Period [d]'); ax.set_ylabel('BLS log-likelihood power')
ax.set_title(f'Blind BLS 1-50 d: top {P_best0:.5f} d, SDE {res["bls_sde_top"]:.1f}'); ax.legend(fontsize=8)

ax = axs[0, 1]
ph_h2 = phase_of(t2, P2, t02) * P2 * 24; mm2 = np.abs(ph_h2) < 15
ph_h1 = phase_of(t, P2, t02) * P2 * 24; mm1 = np.abs(ph_h1) < 15
ax.plot(ph_h2[mm2], (fl2[mm2] - 1) * 1e6, '.', ms=2, color='0.8', alpha=0.6)
xb, yb, eb = binned(ph_h1[mm1], (fl[mm1] - 1) * 1e6, 60)
ax.errorbar(xb, yb, eb, fmt='s', ms=3, color='0.45', alpha=0.8, label=f'pass1 blind detrend: {m1["box"]["depth_ppm"]:.0f} ppm (shoulders = running-median bias)')
xb, yb, eb = binned(ph_h2[mm2], (fl2[mm2] - 1) * 1e6, 60)
ax.errorbar(xb, yb, eb, fmt='o', ms=4, color='tab:blue', label=f'pass2 masked detrend: {m2["box"]["depth_ppm"]:.0f} +- {m2["box"]["depth_err_empirical_ppm"]:.0f} ppm')
xx_ = np.linspace(-15, 15, 500); tp = m2['_trap_params']
ax.plot(xx_, (trap(xx_, *tp) - 1) * 1e6, color='tab:orange', lw=1.5, label=f'trapezoid T14={tp[1]:.2f} h, T23/T14={tp[2]/tp[1]:.2f}')
ax.axvspan(-dur2 * 12, dur2 * 12, color='tab:blue', alpha=0.08)
ax.set_ylim(-4000, 3000); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('ppm')
ax.set_title(f'PSF fold at P={P2:.5f} d; empirical SNR {m2["box"]["empirical_snr"]:.1f}'); ax.legend(fontsize=7, loc='lower right')

ax = axs[1, 0]
ep2 = np.round((t2 - t02) / P2).astype(int)
for lab, par, col in [('odd', 1, 'tab:red'), ('even', 0, 'tab:green')]:
    s_all = mm2 & ((ep2 % 2) == par)
    xb, yb, eb = binned(ph_h2[s_all], (fl2[s_all] - 1) * 1e6, 40)
    ax.errorbar(xb, yb, eb, fmt='o', ms=3, color=col, label=f'{lab}: {m2["odd_even"]["depth_"+lab+"_ppm"]:.0f} +- {m2["odd_even"]["err_"+lab+"_ppm"]:.0f} ppm')
ax.set_ylim(-3000, 2000); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('ppm')
ax.set_title(f'odd/even (pass2): {m2["odd_even"]["odd_even_sigma"]:.2f} sigma'); ax.legend(fontsize=8)

ax = axs[1, 1]
ph_sec = ((phase_of(t2, P2, t02) - 0.5) + 0.5) % 1.0 - 0.5; msec = np.abs(ph_sec * P2 * 24) < 15
xb, yb, eb = binned(ph_sec[msec] * P2 * 24, (fl2[msec] - 1) * 1e6, 60)
ax.errorbar(xb, yb, eb, fmt='o', ms=4, color='tab:purple')
ax.axvspan(-dur2 * 12, dur2 * 12, color='tab:purple', alpha=0.08)
ax.set_ylim(-2000, 2000); ax.set_xlabel('hours from phase 0.5'); ax.set_ylabel('ppm')
ax.set_title(f'secondary at phase 0.5 (pass2): {m2["secondary"]["depth_ppm"]:.0f} +- {m2["secondary"]["err_ppm"]:.0f} ppm ({m2["secondary"]["sigma"]:.1f} sigma)')

ax = axs[2, 0]
ev = m2['events']
ax.errorbar([e['tc_bkjd'] for e in ev], [e['depth_ppm'] for e in ev], [e['err_ppm'] for e in ev], fmt='.', ms=3, color='0.6', alpha=0.8, label='per event')
for name, v in m2['per_year'].items():
    if v['depth_ppm'] is None: continue
    ax.plot(v['t_range_bkjd'], [v['depth_ppm']] * 2, color='tab:blue', lw=3)
    ax.fill_between(v['t_range_bkjd'], v['depth_ppm'] - v['err_ppm'], v['depth_ppm'] + v['err_ppm'], color='tab:blue', alpha=0.2)
    ax.text(np.mean(v['t_range_bkjd']), 5200, f'{name}\n{v["depth_ppm"]:.0f}+-{v["err_ppm"]:.0f}', ha='center', fontsize=7)
ax.axhline(m2['box']['depth_ppm'], color='k', ls='--', lw=0.8); ax.axhline(0, color='k', lw=0.5)
ax.set_ylim(-3000, 6000); ax.set_xlabel('BKJD'); ax.set_ylabel('depth [ppm]')
ax.set_title(f'per-event + per-year depths (pass2); per-year chi2={m2["per_year_chi2"]:.1f}/{m2["per_year_dof"]}; events chi2={m2["events_chi2"]:.0f}/{m2["events_dof"]}'); ax.legend(fontsize=8, loc='lower right')

ax = axs[2, 1]
tts, ffs = ts[g2], fls2[g2]
ph_s = phase_of(tts, P2, t02) * P2 * 24; ms_ = np.abs(ph_s) < 15
xb, yb, eb = binned(ph_s[ms_], (ffs[ms_] - 1) * 1e6, 60)
ax.errorbar(xb, yb, eb, fmt='o', ms=4, color='tab:brown', label='SAP_FLUX (masked detrend)')
xb, yb, eb = binned(ph_h2[mm2], (fl2[mm2] - 1) * 1e6, 60)
ax.errorbar(xb, yb, eb, fmt='.', ms=3, color='tab:blue', alpha=0.5, label='PSF FLUX (masked detrend)')
ax.set_ylim(-4000, 3000); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('ppm')
ax.set_title(f'SAP depth {sap["pass2"]["depth_ppm"]:.0f} +- {sap["pass2"]["depth_err_empirical_ppm"]:.0f} ppm; SAP/PSF = {sap["pass2"]["ratio_sap_over_psf_depth"]:.2f} +- {sap["pass2"]["ratio_err"]:.2f}'); ax.legend(fontsize=8)
fig.suptitle(f'{FNAME} independent re-derivation (no pipeline code): P={P2:.5f} d, t0={t02:.4f} BKJD, box dur={dur2*24:.2f} h', fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fold.png'), dpi=130)
print('wrote', OUT)
