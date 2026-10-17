"""Red-team per-pixel check for kic-003963519 on the cached TPFs of KIC 3963530 (own code).
Season-stacked, event-differenced in-transit deficit maps; empirical PRF template taken from the bright
TPF host (KIC 3963530, G=14.86) and placed at the target, at every Gaia source within 45 arcsec (in or out
of the stamp), and on a grid extending 3 px beyond the stamp. Reports chi2 and the implied depth on each source."""
import glob, json, numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import shift as ndshift

P = 6.815624034742515; T0 = 138.28256264993922; DUR = 2.4375 / 24.0
TARGET = 2052383520062380544; HOST = 2052383524365642112
OUTD = '/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic003963519/out/'
gaia = pd.read_csv(OUTD + 'gaia_45arcsec.csv')
files = sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr003963530_lc_Q*/kplr003963530-*_lpd-targ.fits.gz'))
res = {'quarters': [], 'seasons': {}}
G_t = float(gaia.loc[gaia.source_id == TARGET, 'phot_g_mean_mag'].iloc[0]); G_h = float(gaia.loc[gaia.source_id == HOST, 'phot_g_mean_mag'].iloc[0])

stack = {}
for fn in files:
    h = fits.open(fn)
    hd = h[0].header; q = hd['QUARTER']; season = int(hd['SEASON']); ch = hd['CHANNEL']
    d = h[1].data
    t = d['TIME']; fl = d['FLUX']; ql = d['QUALITY']
    ok = np.isfinite(t) & (ql == 0)
    t = t[ok]; fl = fl[ok]
    w = WCS(h[2].header)
    xs, ys = w.all_world2pix(gaia.ra.values, gaia.dec.values, 0)
    pos = {int(s): (float(x), float(y)) for s, x, y in zip(gaia.source_id, xs, ys)}
    ny, nx = fl.shape[1:]
    ph = ((t - T0 + P / 2) % P) - P / 2
    intr = np.abs(ph) < DUR / 2
    oot = (np.abs(ph) > 1.0 * DUR)
    deficit = np.zeros((ny, nx)); var = np.zeros((ny, nx)); nev = 0
    for n in np.unique(np.round((t[intr] - T0) / P).astype(int)):
        tc = T0 + n * P
        a = intr & (np.abs(t - tc) < DUR)
        b = oot & (np.abs(t - tc) < 1.0) & (np.abs(t - tc) > 0.75 * DUR)
        if a.sum() < 2 or b.sum() < 15:
            continue
        img_in = np.nanmean(fl[a], axis=0); img_out = np.nanmedian(fl[b], axis=0)
        sig = np.nanstd(fl[b], axis=0) / np.sqrt(a.sum())
        dmap = img_out - img_in
        good = np.isfinite(dmap) & np.isfinite(sig) & (sig > 0)
        deficit[good] += dmap[good]; var[good] += sig[good]**2; nev += 1
    if nev < 3:
        continue
    img_oot = np.nanmedian(fl[oot], axis=0)
    key = season
    if key not in stack:
        stack[key] = {'def': np.zeros((ny, nx)), 'var': np.zeros((ny, nx)), 'oot': np.zeros((ny, nx)), 'n': 0, 'pos': [], 'shape': (ny, nx), 'quarters': [], 'channel': int(ch)}
    s = stack[key]
    if s['shape'] != (ny, nx):
        res['quarters'].append({'quarter': int(q), 'season': key, 'note': 'shape differs, not stacked', 'shape': [ny, nx]})
        continue
    s['def'] += deficit; s['var'] += var; s['oot'] += img_oot * nev; s['n'] += nev; s['pos'].append(pos); s['quarters'].append(int(q))
    snr = deficit / np.sqrt(var + 1e-30)
    xt, yt = pos[TARGET]
    res['quarters'].append({'quarter': int(q), 'season': key, 'channel': int(ch), 'n_events': nev, 'shape': [ny, nx], 'target_px': [xt, yt], 'host_px': list(pos[HOST]),
                            'max_snr': float(np.nanmax(snr)), 'max_snr_px': [int(np.argmax(np.nanmax(snr, axis=0))), int(np.argmax(np.nanmax(snr, axis=1)))]})

for key, s in stack.items():
    ny, nx = s['shape']; n = s['n']
    dmap = s['def'] / n; sig = np.sqrt(s['var']) / n; oot = s['oot'] / n
    snr = dmap / sig
    pos = {k: (np.mean([p[k][0] for p in s['pos']]), np.mean([p[k][1] for p in s['pos']])) for k in s['pos'][0]}
    xt, yt = pos[TARGET]; xh, yh = pos[HOST]
    tot = float(np.nansum(oot))
    def template(x0, y0):
        return ndshift(np.nan_to_num(oot), (y0 - yh, x0 - xh), order=1, mode='constant', cval=0.0)
    good = np.isfinite(dmap) & np.isfinite(sig) & (sig > 0)
    def fit_at(x0, y0):
        m = template(x0, y0)
        den = np.sum(m[good]**2 / sig[good]**2)
        if den <= 0: return np.inf, 0.0, np.inf, 0.0
        A = np.sum(m[good] * dmap[good] / sig[good]**2) / den
        c2 = np.sum(((dmap[good] - A * m[good]) / sig[good])**2)
        return float(c2), float(A * tot), float(tot / np.sqrt(den)), float(np.sum(m[good]) / tot)
    chi2_null = float(np.sum((dmap[good] / sig[good])**2))
    c2_t, amp_t, amperr_t, frac_t = fit_at(xt, yt)
    # host flux in stamp -> flux of any source of magnitude G (full PRF) = tot_host_full * 10^(-0.4 (G-G_h)); tot is host in-stamp (approx full)
    def full_flux(G): return tot * 10 ** (-0.4 * (G - G_h))
    F_t = full_flux(G_t)
    sources = []
    for _, row in gaia.iterrows():
        sid = int(row.source_id); x0, y0 = pos[sid]
        if x0 < -3 or x0 > nx + 2 or y0 < -3 or y0 > ny + 2:
            continue
        c2, amp, ae, frac = fit_at(x0, y0)
        if not np.isfinite(c2) or frac < 0.02:
            continue
        Fs = full_flux(float(row.phot_g_mean_mag))
        sources.append({'source_id': sid, 'G': float(row.phot_g_mean_mag), 'sep_arcsec': float(row.sep_arcsec), 'px': [x0, y0], 'in_stamp': bool(-0.5 <= x0 <= nx - 0.5 and -0.5 <= y0 <= ny - 0.5),
                        'prf_frac_in_stamp': frac, 'chi2': c2, 'dchi2_vs_null': chi2_null - c2, 'amp_e_s': amp, 'amp_err_e_s': ae,
                        'implied_depth_ppm': amp / Fs * 1e6, 'implied_depth_err_ppm': ae / Fs * 1e6})
    sources.sort(key=lambda r: r['chi2'])
    # grid search including 3 px beyond the stamp
    best = (np.inf, None, None, 0); grid = []
    for gx in np.arange(-3, nx + 3, 0.25):
        for gy in np.arange(-3, ny + 3, 0.25):
            c2, amp, ae, frac = fit_at(gx, gy)
            if frac < 0.02: continue
            grid.append((gx, gy, c2, amp))
            if amp > 0 and c2 < best[0]:
                best = (c2, gx, gy, amp)
    grid = np.array(grid)
    region = grid[(grid[:, 3] > 0) & (grid[:, 2] - best[0] < 2.3)]
    region3 = grid[(grid[:, 3] > 0) & (grid[:, 2] - best[0] < 11.8)]   # 3 sigma for 2 params
    dist_best_target = float(np.hypot(best[1] - xt, best[2] - yt))
    tgt_in_1sig = bool(np.any(np.hypot(region[:, 0] - xt, region[:, 1] - yt) < 0.3)) if len(region) else False
    tgt_in_3sig = bool(np.any(np.hypot(region3[:, 0] - xt, region3[:, 1] - yt) < 0.3)) if len(region3) else False
    # expected on-target deficit (e/s, full-PRF equivalent) for 1060 ppm
    exp_amp = 1.06e-3 * F_t
    res['seasons'][str(key)] = {
        'quarters': s['quarters'], 'channel': s['channel'], 'n_events_stacked': int(n), 'shape': [ny, nx], 'target_px': [xt, yt], 'host_px': [xh, yh],
        'host_in_stamp_e_s': tot, 'target_full_flux_est_e_s': F_t, 'target_prf_frac_in_stamp': frac_t,
        'snr_map': np.round(snr, 1).tolist(), 'deficit_map_e_s': np.round(dmap, 3).tolist(), 'sig_map_e_s': np.round(sig, 3).tolist(),
        'chi2_null': chi2_null, 'chi2_prf_at_target': c2_t, 'amp_at_target_e_s': amp_t, 'amp_err_at_target_e_s': amperr_t,
        'implied_depth_on_target_ppm': amp_t / F_t * 1e6, 'implied_depth_on_target_err_ppm': amperr_t / F_t * 1e6, 'expected_amp_for_1060ppm_e_s': exp_amp,
        'sources': sources,
        'best_fit_px': [best[1], best[2]], 'chi2_best': best[0], 'amp_best_e_s': best[3], 'best_to_target_px': dist_best_target, 'best_to_target_arcsec': dist_best_target * 3.98,
        'target_inside_1sigma_region': tgt_in_1sig, 'target_inside_3sigma_region': tgt_in_3sig,
        'onesigma_region_extent_px': [float(region[:, 0].min()), float(region[:, 0].max()), float(region[:, 1].min()), float(region[:, 1].max())] if len(region) else None,
        'threesigma_region_extent_px': [float(region3[:, 0].min()), float(region3[:, 0].max()), float(region3[:, 1].min()), float(region3[:, 1].max())] if len(region3) else None,
        'n_good_px': int(good.sum()),
    }
    print(f"season {key} ch{s['channel']} Q{s['quarters']} n_ev={n} target=({xt:.2f},{yt:.2f}) prf_frac_in_stamp={frac_t:.2f} host=({xh:.2f},{yh:.2f}) host_in_stamp={tot:.0f} e/s")
    print(f"   chi2 null {chi2_null:.1f} | at target {c2_t:.1f} amp {amp_t:.2f}+-{amperr_t:.2f} e/s -> implied depth on target {amp_t / F_t * 1e6:.0f}+-{amperr_t / F_t * 1e6:.0f} ppm (expected amp {exp_amp:.2f} e/s)")
    print(f"   best ({best[1]:.2f},{best[2]:.2f}) chi2 {best[0]:.1f} amp {best[3]:.2f}; {dist_best_target:.2f} px = {dist_best_target * 3.98:.1f}\" from target; target in 1sig {tgt_in_1sig} 3sig {tgt_in_3sig}")
    for r in sources[:6]:
        print(f"     src {r['source_id']} G={r['G']:.2f} sep={r['sep_arcsec']:.1f}\" px=({r['px'][0]:.1f},{r['px'][1]:.1f}) in={r['in_stamp']} frac={r['prf_frac_in_stamp']:.2f} chi2={r['chi2']:.1f} dchi2={r['dchi2_vs_null']:.1f} amp={r['amp_e_s']:.2f}+-{r['amp_err_e_s']:.2f} depth_on_src={r['implied_depth_ppm']:.0f}+-{r['implied_depth_err_ppm']:.0f} ppm")
    print('   snr map:\n', np.round(snr, 1))

# all-season combined chi2 per source
comb = {}
for key, sr in res['seasons'].items():
    for r in sr['sources']:
        c = comb.setdefault(r['source_id'], {'G': r['G'], 'sep_arcsec': r['sep_arcsec'], 'chi2': 0.0, 'dchi2': 0.0, 'n_seasons': 0, 'amp_w': 0.0, 'w': 0.0})
        c['chi2'] += r['chi2']; c['dchi2'] += r['dchi2_vs_null']; c['n_seasons'] += 1
        if r['implied_depth_err_ppm'] > 0:
            c['amp_w'] += r['implied_depth_ppm'] / r['implied_depth_err_ppm']**2; c['w'] += 1 / r['implied_depth_err_ppm']**2
for sid, c in comb.items():
    c['implied_depth_ppm'] = c['amp_w'] / c['w'] if c['w'] > 0 else None; c['implied_depth_err_ppm'] = 1 / np.sqrt(c['w']) if c['w'] > 0 else None
    c.pop('amp_w'); c.pop('w')
res['combined_by_source'] = sorted(comb.values(), key=lambda c: -c['dchi2'] if c['n_seasons'] == len(res['seasons']) else 1e9)
for c in res['combined_by_source'][:10]:
    print('COMB', c)
json.dump(res, open(OUTD + 'pixel_checks.json', 'w'), indent=1, default=float)
print('written')
