#!/usr/bin/env python
"""Red-team pixel test for kic-010213563 (own code; no pipeline/ imports).

Per quarter and per season: event-differenced per-pixel deficit maps on the TPFs of KIC 10213566,
single-Gaussian deficit fit with the position free (allowed beyond the stamp), chi2 at the Gaia
positions of the target, the 0.97" neighbour and the host, a 'scaled image' systematics model,
and the deficit in e/s compared with depth x KBONUS PSF flux per season.
"""
import json, sys, glob, warnings, numpy as np
warnings.filterwarnings('ignore')
import lightkurve as lk
from astropy.io import fits
from scipy.optimize import minimize

OUT = sys.argv[1]
PNG = sys.argv[2] if len(sys.argv) > 2 else None
P, T0, DUR_H = 36.76708318978774, 164.85416024494916, 5.625
DUR = DUR_H / 24
TPFS = sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr010213566_lc_Q011111111111111111/*targ.fits.gz'))
SRC = {'target': (294.1226612214163, 47.20323165191153, 17.1037),
       'neighbour_0p97': (294.12304409069, 47.20315938049, 18.3494),
       'host_7p8': (294.12313328222, 47.20537958153, 14.2550),
       'nb_5p0_G20p4': (294.121382727, 47.20215180122, 20.36)}
TARGET_LC = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-010213563_kepler_v1.0_lc.fits'


def gauss(shape, x0, y0, sig):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    g = np.exp(-0.5 * ((xx - x0)**2 + (yy - y0)**2) / sig**2)
    return g / (2 * np.pi * sig**2)


def fit_amp(model, data, err):
    w = 1 / err**2
    a = np.sum(w * model * data) / np.sum(w * model**2)
    chi2 = np.sum(w * (data - a * model)**2)
    return a, chi2


def analyse(dmap, derr, oot_img, pos, sig, label):
    """dmap: deficit (e/s per pixel, positive = missing light). pos: dict name->(x,y)."""
    ny, nx = dmap.shape
    ok = np.isfinite(dmap) & np.isfinite(derr) & (derr > 0)
    d, e = dmap[ok], derr[ok]
    out = {'label': label, 'n_pix': int(ok.sum())}
    out['chi2_null'] = float(np.sum((d / e)**2))
    # scaled image (every source dims by the same fraction: systematics / focus)
    a, c = fit_amp(oot_img[ok], d, e)
    out['scaled_image'] = {'chi2': float(c), 'frac_ppm': float(a * 1e6)}
    # gaussian at fixed positions
    for name, (x, y) in pos.items():
        g = gauss(dmap.shape, x, y, sig)
        a, c = fit_amp(g[ok], d, e)
        out['at_' + name] = {'chi2': float(c), 'deficit_e_s': float(a), 'px': [float(x), float(y)]}
    # free position on a grid 3 px beyond the stamp
    xs = np.arange(-3, nx + 3, 0.1); ys = np.arange(-3, ny + 3, 0.1)
    chi = np.full((len(ys), len(xs)), np.nan); amp = np.full_like(chi, np.nan)
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            g = gauss(dmap.shape, x, y, sig)
            if g[ok].sum() < 1e-3:
                continue
            a, c = fit_amp(g[ok], d, e)
            chi[j, i] = c; amp[j, i] = a
    jb, ib = np.unravel_index(np.nanargmin(chi), chi.shape)
    cmin = chi[jb, ib]
    reg = chi < cmin + 11.83   # 3 sigma for 2 parameters
    yy, xx = np.meshgrid(ys, xs, indexing='ij')
    out['free'] = {'px': [float(xs[ib]), float(ys[jb])], 'chi2': float(cmin), 'deficit_e_s': float(amp[jb, ib]),
                   'x_3sig': [float(np.nanmin(xx[reg])), float(np.nanmax(xx[reg]))], 'y_3sig': [float(np.nanmin(yy[reg])), float(np.nanmax(yy[reg]))]}
    tx, ty = pos['target']
    out['free']['offset_from_target_px'] = float(np.hypot(xs[ib] - tx, ys[jb] - ty))
    out['free']['offset_from_target_arcsec'] = float(np.hypot(xs[ib] - tx, ys[jb] - ty) * 3.98)
    # which named sources lie inside the 3-sigma region
    inside = []
    for name, (x, y) in pos.items():
        i = np.argmin(np.abs(xs - x)); j = np.argmin(np.abs(ys - y))
        if 0 <= i < len(xs) and 0 <= j < len(ys) and np.isfinite(chi[j, i]):
            inside.append((name, bool(chi[j, i] < cmin + 11.83), float(chi[j, i] - cmin)))
    out['sources_dchi2_vs_free'] = inside
    out['psf_sigma_px'] = float(sig)
    return out, (xs, ys, chi)


results = {'ephemeris': {'P': P, 't0_bkjd': T0, 'dur_h': DUR_H}, 'quarters': [], 'seasons': {}}
season_stack = {}
figs = []
for path in TPFS:
    tpf = lk.KeplerTargetPixelFile(path, quality_bitmask='default')
    q = int(tpf.quarter); season = int(tpf.hdu[0].header['SEASON']); chan = int(tpf.hdu[0].header['CHANNEL'])
    t = tpf.time.value; cube = np.array(tpf.flux.value, float); cerr = np.array(tpf.flux_err.value, float)
    good = np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    t, cube, cerr = t[good], cube[good], cerr[good]
    w = tpf.wcs
    pos = {}
    for name, (ra, dec, G) in SRC.items():
        x, y = w.all_world2pix(ra, dec, 0)
        pos[name] = (float(x), float(y))
    ny, nx = cube.shape[1:]
    oot_img = np.nanmedian(cube, axis=0)
    # PSF width from the out-of-transit image: fit Gaussian at the host position + flat background
    hx, hy = pos['host_7p8']
    def chi_sig(p):
        sg, bg = p
        g = gauss((ny, nx), hx, hy, abs(sg))
        a, c = fit_amp(g, oot_img - bg, np.sqrt(np.abs(oot_img)) + 1)
        return c
    best = min((minimize(chi_sig, [s0, np.nanmin(oot_img)], method='Nelder-Mead') for s0 in (0.5, 0.8, 1.2)), key=lambda r: r.fun)
    sig = float(np.clip(abs(best.x[0]), 0.4, 1.3))
    # events in this quarter
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    epochs = np.round((t - T0) / P).astype(int)
    dmaps, derrs, nev = [], [], 0
    for n_ in np.unique(epochs):
        s = epochs == n_
        dt = (t - (T0 + n_ * P))
        i = s & (np.abs(dt) < 0.4 * DUR)             # core
        o = s & (np.abs(dt) > 0.75 * DUR) & (np.abs(dt) < 1.5)
        if i.sum() < 4 or o.sum() < 20:
            continue
        # per-pixel linear baseline in time from the oot cadences
        A = np.vstack([np.ones(o.sum()), dt[o]]).T
        coef, *_ = np.linalg.lstsq(A, cube[o].reshape(o.sum(), -1), rcond=None)
        base_in = (np.vstack([np.ones(i.sum()), dt[i]]).T @ coef)
        resid = cube[o].reshape(o.sum(), -1) - A @ coef
        pix_err = np.std(resid, axis=0) / np.sqrt(i.sum())
        dmap = (base_in - cube[i].reshape(i.sum(), -1)).mean(axis=0)
        dmaps.append(dmap.reshape(ny, nx)); derrs.append(pix_err.reshape(ny, nx)); nev += 1
    if nev == 0:
        results['quarters'].append({'quarter': q, 'season': season, 'channel': chan, 'error': 'no usable events'}); continue
    W = np.array([1 / e**2 for e in derrs]); D = np.array(dmaps)
    dmap = np.sum(W * D, axis=0) / np.sum(W, axis=0); derr = 1 / np.sqrt(np.sum(W, axis=0))
    res, _ = analyse(dmap, derr, oot_img, pos, sig, 'Q%d' % q)
    res.update({'quarter': q, 'season': season, 'channel': chan, 'n_events': nev, 'stamp': [int(ny), int(nx)], 'positions_px': pos})
    results['quarters'].append(res)
    key = 'season%d' % season
    season_stack.setdefault(key, []).append((dmap, derr, oot_img, pos, sig, q, chan))

# KBONUS PSF flux of the target per quarter (median e/s) for the deficit budget
h = fits.open(TARGET_LC)['LIGHTCURVE_STITCHED'].data
qq = np.array(h['QUARTER'], int); fl = np.array(h['FLUX'], float); ok = np.isfinite(fl) & (np.array(h['SAP_QUALITY'], int) == 0)
psf_flux_q = {int(k): float(np.nanmedian(fl[ok & (qq == k)])) for k in np.unique(qq)}
results['target_psf_flux_e_s_by_quarter'] = psf_flux_q

for key, items in sorted(season_stack.items()):
    # stack quarters of the same season (identical geometry up to sub-pixel shifts; use the first quarter's positions)
    shapes = set(x[0].shape for x in items)
    if len(shapes) != 1:
        results['seasons'][key] = {'error': 'shape mismatch %s' % shapes}; continue
    W = np.array([1 / x[1]**2 for x in items]); D = np.array([x[0] for x in items])
    dmap = np.sum(W * D, axis=0) / np.sum(W, axis=0); derr = 1 / np.sqrt(np.sum(W, axis=0))
    oot = np.mean([x[2] for x in items], axis=0)
    pos = {k: (float(np.mean([x[3][k][0] for x in items])), float(np.mean([x[3][k][1] for x in items]))) for k in items[0][3]}
    sig = float(np.mean([x[4] for x in items]))
    res, (xs, ys, chi) = analyse(dmap, derr, oot, pos, sig, key)
    res['quarters'] = [x[5] for x in items]; res['channel'] = items[0][6]; res['positions_px'] = pos
    res['psf_flux_e_s_mean'] = float(np.mean([psf_flux_q.get(x[5], np.nan) for x in items]))
    results['seasons'][key] = res
    figs.append((key, dmap, derr, oot, pos, xs, ys, chi, res))

# summary across seasons
summ = {}
for key, r in results['seasons'].items():
    if 'free' in r:
        summ[key] = {'free_px': r['free']['px'], 'offset_arcsec': r['free']['offset_from_target_arcsec'], 'x_3sig': r['free']['x_3sig'], 'y_3sig': r['free']['y_3sig'],
                     'deficit_e_s_at_target': r['at_target']['deficit_e_s'], 'deficit_e_s_free': r['free']['deficit_e_s'],
                     'chi2_null': r['chi2_null'], 'chi2_scaled_image': r['scaled_image']['chi2'], 'chi2_target': r['at_target']['chi2'],
                     'chi2_neighbour': r['at_neighbour_0p97']['chi2'], 'chi2_host': r['at_host_7p8']['chi2'], 'chi2_free': r['free']['chi2'],
                     'deficit_at_host_e_s': r['at_host_7p8']['deficit_e_s'], 'psf_flux_e_s': r['psf_flux_e_s_mean']}
results['summary'] = summ
json.dump(results, open(OUT, 'w'), indent=1, default=float)
print(json.dumps(summ, indent=1))

if PNG and figs:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    n = len(figs)
    fig, axes = plt.subplots(3, n, figsize=(4 * n, 11))
    for k, (key, dmap, derr, oot, pos, xs, ys, chi, res) in enumerate(figs):
        ax = axes[0, k]; im = ax.imshow(oot, origin='lower', cmap='gray_r'); ax.set_title('%s Q%s ch%d: out-of-transit' % (key, res['quarters'], res['channel']), fontsize=9)
        for name, (x, y) in pos.items():
            ax.plot(x, y, 'o', mfc='none', ms=8, label=name);
        ax.legend(fontsize=6, loc='upper left')
        ax = axes[1, k]; im = ax.imshow(dmap / derr, origin='lower', cmap='RdBu_r', vmin=-6, vmax=6); ax.set_title('deficit SNR map (red = missing light)', fontsize=9)
        for name, (x, y) in pos.items():
            ax.plot(x, y, 'o', mfc='none', ms=8)
        ax.plot(*res['free']['px'], 'kx', ms=10)
        plt.colorbar(im, ax=ax, fraction=0.04)
        ax = axes[2, k]; cm = np.nanmin(chi)
        im = ax.imshow(np.clip(chi - cm, 0, 50), origin='lower', extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap='viridis')
        ax.contour(xs, ys, chi - cm, levels=[2.3, 6.2, 11.8], colors='w', linewidths=0.8)
        for name, (x, y) in pos.items():
            ax.plot(x, y, 'o', mfc='none', ms=8, mec='r')
        ax.add_patch(plt.Rectangle((-0.5, -0.5), dmap.shape[1], dmap.shape[0], fill=False, ec='k', ls='--'))
        ax.set_title('dchi2 of free Gaussian position (1,2,3 sigma)', fontsize=9)
        plt.colorbar(im, ax=ax, fraction=0.04)
    plt.tight_layout(); plt.savefig(PNG, dpi=110)
    print('wrote', PNG)
