#!/usr/bin/env python
"""Red-team summary figure for kic-010213563 (own code)."""
import json, sys, numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

R = '/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-010213563/redteam/'
P, T0, DUR = 36.76708318978774, 164.85416024494916, 5.625 / 24
h = fits.open('/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-010213563_kepler_v1.0_lc.fits')['LIGHTCURVE_STITCHED'].data
t = np.array(h['TIME'], float); f = np.array(h['FLUX'], float); q = np.array(h['QUARTER'], int); qual = np.array(h['SAP_QUALITY'], int)
m = np.isfinite(t) & np.isfinite(f) & (qual == 0) & (f > 0); t, f, q = t[m], f[m], q[m]
for qq in np.unique(q): s = q == qq; f[s] /= np.median(f[s])
ph = ((t - T0) / P + 0.5) % 1 - 0.5; intr = np.abs(ph * P) < 0.75 * DUR
o = np.argsort(t); t, f, ph, intr = t[o], f[o], ph[o], intr[o]
tr = np.full_like(f, np.nan); gaps = np.where(np.diff(t) > 0.5)[0]; ed = np.r_[0, gaps + 1, len(t)]
for a, b in zip(ed[:-1], ed[1:]):
    tt, ff, mm = t[a:b], f[a:b], intr[a:b]
    for i in range(len(tt)):
        s = (~mm) & (np.abs(tt - tt[i]) < 0.5)
        if s.sum() < 5: s = (~mm) & (np.abs(tt - tt[i]) < 1.0)
        tr[a + i] = np.median(ff[s]) if s.sum() >= 3 else np.nan
ff = f / tr; ok = np.isfinite(ff); t, ff, ph = t[ok], ff[ok], ph[ok]
mad = 1.4826 * np.median(np.abs(ff - 1)); keep = ~((ff - 1 > 3.5 * mad)); t, ff, ph = t[keep], ff[keep], ph[keep]

lc = json.load(open(R + 'clean_scan_events.json')); sh = json.load(open(R + 'shape_fit.json')); pf = json.load(open(R + 'pixel_free_fit.json'))
fig, ax = plt.subplots(2, 3, figsize=(17, 9))
# 1 fold + batman best model
hrs = ph * P * 24; w = np.abs(hrs) < 24
ax[0, 0].plot(hrs[w], (ff[w] - 1) * 1e6, '.', ms=2, color='0.6', alpha=.5)
bins = np.arange(-24, 24.1, 0.75); idx = np.digitize(hrs[w], bins)
bm = [np.mean((ff[w] - 1)[idx == i]) * 1e6 for i in range(1, len(bins))]
ax[0, 0].plot(0.5 * (bins[1:] + bins[:-1]), bm, 'o', color='C1', ms=4)
try:
    import batman
    pr = batman.TransitParams(); bf = sh['best_fit']; pr.t0 = 0; pr.per = P; pr.rp = bf['rp']; pr.a = np.exp(bf['aR']); pr.inc = np.degrees(np.arccos(bf['b'] / np.exp(bf['aR']))); pr.ecc = 0; pr.w = 90; pr.limb_dark = 'quadratic'; pr.u = [0.5, 0.2]
    tm = np.linspace(-1, 1, 2000); mod = batman.TransitModel(pr, tm, supersample_factor=15, exp_time=1765.5 / 86400).light_curve(pr)
    ax[0, 0].plot(tm * 24, (mod - 1) * 1e6, 'k-', lw=1, label='batman best (b=%.2f, rho*=%.2f g/cc)' % (bf['b'], sh['rho_star_gcc'][1]))
except Exception as e:
    print('batman plot failed', e)
ax[0, 0].set_ylim(-6000, 4000); ax[0, 0].set_xlabel('hours from mid-transit'); ax[0, 0].set_ylabel('ppm'); ax[0, 0].legend(fontsize=8)
ax[0, 0].set_title('kic-010213563 PSF flux folded at P=36.7671 d (masked 1-d running median)', fontsize=9)
# 2 rho posterior summary vs references
r = sh['rho_star_gcc']
ax[0, 1].errorbar([0], [r[1]], yerr=[[r[1] - r[0]], [r[2] - r[1]]], fmt='o', color='k', capsize=4, label='fit (16/50/84%)')
ax[0, 1].axhline(sh['rho_star_gcc_99pct_upper'], ls=':', color='k', label='99% upper')
for i, (k, v) in enumerate(sh['reference_rho_gcc'].items()):
    ax[0, 1].axhline(v, color='C%d' % (i + 1), lw=1, label='%s %.2f' % (k, v))
ax[0, 1].set_ylabel('stellar density (g/cc), circular orbit'); ax[0, 1].set_xticks([]); ax[0, 1].legend(fontsize=7); ax[0, 1].set_title('density from transit shape (batman, 30-min supersampled)', fontsize=9)
# 3 per-event depths coloured by season (lc_systematics events with flat baseline for event0)
ev = json.load(open(R + 'lc_systematics.json'))['events']
cols = {0: 'C0', 1: 'C1', 2: 'C2', 3: 'C3'}
for e in ev:
    d = e['depth_ppm'] if e['n'] != 0 else lc['psf_events']['event0_depth'][0]
    ax[0, 2].errorbar(e['tc_bkjd'], d, yerr=e['err_ppm'], fmt='o', ms=4, color=cols[e['season']])
for s_, c in cols.items(): ax[0, 2].plot([], [], 'o', color=c, label='season %d' % s_)
ax[0, 2].axhline(lc['psf_events']['wmean'], color='k', lw=1); ax[0, 2].legend(fontsize=8); ax[0, 2].set_xlabel('BKJD'); ax[0, 2].set_ylabel('event depth (ppm)'); ax[0, 2].set_title('per-event PSF depth (flat local baseline), by Kepler season', fontsize=9)
# 4 season depths PSF & SAP
ss = sorted(lc['psf_season']); x = np.arange(4)
ax[1, 0].errorbar(x - 0.1, [lc['psf_season'][s][0] for s in ss], yerr=[lc['psf_season'][s][1] for s in ss], fmt='o', label='PSF depth')
ax[1, 0].errorbar(x + 0.1, [lc['sap_season'][s][0] for s in ss], yerr=[lc['sap_season'][s][1] for s in ss], fmt='s', label='SAP depth')
ax[1, 0].set_xticks(x); ax[1, 0].set_xticklabels(['s0 Q2,6,10,14\nch47', 's1 Q3,7,11,15\nch23', 's2 Q4,8,12,16\nch39', 's3 Q1,5,9,13,17\nch63'], fontsize=8); ax[1, 0].legend(fontsize=8); ax[1, 0].set_ylabel('ppm')
ax[1, 0].set_title('depth by season: PSF chi2=%.1f/3' % lc['psf_season_chi2'][0], fontsize=9)
# 5 column profile of the deficit in the target row per season (from pixel_free_fit stacks are not stored; use wcs-calibrated if present)
try:
    wc = json.load(open(R + 'pixel_wcs_calibrated.json'))
    for i, (k, s) in enumerate(sorted(wc['seasons'].items())):
        ct = s['column_test']; cs = ct['rows_summed_deficit_by_column']; ce = ct['rows_summed_err']
        ax[1, 1].errorbar(np.arange(len(cs)) + 0.08 * i, cs, yerr=ce, fmt='o-', ms=4, label='%s target x=%.2f nb x=%.2f' % (k, s['corrected_px']['target'][0], s['corrected_px']['neighbour'][0]))
    ax[1, 1].set_xlabel('stamp column'); ax[1, 1].set_ylabel('deficit e/s (3 rows around target)'); ax[1, 1].legend(fontsize=7); ax[1, 1].set_title('in-transit missing light by column (WCS-calibrated)', fontsize=9)
except Exception as e:
    ax[1, 1].text(0.1, 0.5, 'wcs-calibrated file missing: %s' % e, fontsize=7)
# 6 text summary of pixel localisation
txt = []
try:
    for k, s in sorted(wc['seasons'].items()):
        f0 = s['fits']['fitted']
        txt.append('%s (sig %.2f px, dx %.2f dy %.2f): chi2 target %.0f | neighbour %.0f | free %.0f at (%.2f,%.2f); 3sig x [%.1f,%.1f]' % (k, s['wcs_calibration']['sigma_px'], s['wcs_calibration']['dx_px'], s['wcs_calibration']['dy_px'], f0['target']['chi2'], f0['neighbour']['chi2'], f0['free_any']['chi2'], f0['free_any']['px'][0], f0['free_any']['px'][1], f0['free_any']['x_3sig'][0], f0['free_any']['x_3sig'][1]))
except Exception:
    for k, s in sorted(pf['summary'].items()):
        txt.append('%s: chi2 target %.0f | neighbour %.0f | free %.0f at (%.1f,%.1f)' % (k, s['chi2_target'], s['chi2_neighbour'], s['chi2_free'], *s['free_px']))
ax[1, 2].axis('off'); ax[1, 2].text(0, 0.95, 'Pixel localisation (single-Gaussian deficit fit)\n' + '\n'.join(txt), fontsize=7, va='top', family='monospace')
plt.tight_layout(); plt.savefig(R + 'redteam_summary.png', dpi=110); print('wrote', R + 'redteam_summary.png')
