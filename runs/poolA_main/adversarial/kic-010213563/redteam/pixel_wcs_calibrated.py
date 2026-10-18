#!/usr/bin/env python
"""WCS- and PSF-calibrated deficit localisation for kic-010213563 (own code; no pipeline/ imports).

For each season stack: fit the out-of-transit image with a multi-Gaussian source model (host, target,
0.97" neighbour, G=20.4 and G=21.0 sources at their Gaia pixel positions) with a common (dx, dy) shift
and a common sigma -> calibrates the TPF WCS and the PSF width in that season. Then refit the deficit
map (same event-differencing as pixel_free_fit.py) with corrected positions, for several sigmas, and
report chi2 for: target, neighbour, free-inside-stamp, free-anywhere, plus the column ratio test.
"""
import json, sys, glob, warnings, numpy as np
warnings.filterwarnings('ignore')
import lightkurve as lk
from scipy.optimize import minimize

OUT = sys.argv[1]
P, T0, DUR_H = 36.76708318978774, 164.85416024494916, 5.625
DUR = DUR_H / 24
TPFS = sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr010213566_lc_Q011111111111111111/*targ.fits.gz'))
SRC = {'host': (294.12313328222, 47.20537958153, 14.2550), 'target': (294.1226612214163, 47.20323165191153, 17.1037),
       'neighbour': (294.12304409069, 47.20315938049, 18.3494), 'nb20': (294.121382727, 47.20215180122, 20.36),
       'nb21': (294.12778174586, 47.2030559186, 21.02), 'nb18_15as': (294.11725696692, 47.20526939268, 18.48), 'nb17_19as': (294.12638542652, 47.20793099831, 17.79)}


def gauss(shape, x0, y0, sig):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return np.exp(-0.5 * ((xx - x0)**2 + (yy - y0)**2) / sig**2) / (2 * np.pi * sig**2)


def lin_fit(models, data, err):
    """non-negative-ish linear least squares for amplitudes; returns amps, chi2"""
    A = np.array([m.ravel() for m in models]).T / err.ravel()[:, None]
    b = data.ravel() / err.ravel()
    # ridge-regularised normal equations (models of off-stamp sources are ~0 and make lstsq singular)
    ATA = A.T @ A; lam = 1e-9 * np.trace(ATA) / max(len(models), 1)
    coef = np.linalg.solve(ATA + lam * np.eye(len(models)), A.T @ b)
    chi2 = np.sum((b - A @ coef)**2)
    return coef, chi2


def deficit_map(tpf):
    t = tpf.time.value; cube = np.array(tpf.flux.value, float)
    ny, nx = cube.shape[1:]
    valid = np.isfinite(np.nanmedian(cube, axis=0))          # pixels that exist in this stamp (NaN corners are zero-weight)
    cube[:, ~valid] = 0.0
    good = np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1); t, cube = t[good], cube[good]
    epochs = np.round((t - T0) / P).astype(int)
    dmaps, derrs = [], []
    for n_ in np.unique(epochs):
        s = epochs == n_; dt = t - (T0 + n_ * P)
        i = s & (np.abs(dt) < 0.4 * DUR); o = s & (np.abs(dt) > 0.75 * DUR) & (np.abs(dt) < 1.5)
        if i.sum() < 4 or o.sum() < 20: continue
        A = np.vstack([np.ones(o.sum()), dt[o]]).T
        coef, *_ = np.linalg.lstsq(A, cube[o].reshape(o.sum(), -1), rcond=None)
        base_in = np.vstack([np.ones(i.sum()), dt[i]]).T @ coef
        resid = cube[o].reshape(o.sum(), -1) - A @ coef
        dm = (base_in - cube[i].reshape(i.sum(), -1)).mean(axis=0).reshape(ny, nx)
        de = (np.std(resid, axis=0) / np.sqrt(i.sum())).reshape(ny, nx)
        dm[~valid] = 0.0; de[~valid] = np.inf
        dmaps.append(dm); derrs.append(de)
    oot = np.nanmedian(cube, axis=0); oot[~valid] = 0.0
    return dmaps, derrs, oot


results = {'seasons': {}}
stacks = {}
for path in TPFS:
    tpf = lk.KeplerTargetPixelFile(path, quality_bitmask='default')
    season = int(tpf.hdu[0].header['SEASON']); q = int(tpf.quarter)
    dm, de, oot = deficit_map(tpf)
    pos = {k: tuple(map(float, tpf.wcs.all_world2pix(v[0], v[1], 0))) for k, v in SRC.items()}
    stacks.setdefault(season, []).append((q, dm, de, oot, pos, int(tpf.hdu[0].header['CHANNEL'])))

for season, items in sorted(stacks.items()):
    shapes = set(x[3].shape for x in items)
    if len(shapes) != 1: continue
    ny, nx = items[0][3].shape
    oot = np.mean([x[3] for x in items], axis=0)
    pos = {k: (float(np.mean([x[4][k][0] for x in items])), float(np.mean([x[4][k][1] for x in items]))) for k in SRC}
    ooterr = np.sqrt(np.abs(oot) / 30.) + 2.0; ooterr[oot == 0] = np.inf   # per-pixel scatter of a ~1-month median image; rough, only sets relative weights
    # ---- WCS / PSF calibration on the oot image: common shift + sigma; amplitudes + flat background linear
    names = list(SRC)
    def chi_oot(p):
        dx, dy, sg = p
        if not (0.3 < sg < 1.5 and abs(dx) < 1.5 and abs(dy) < 1.5): return 1e30
        models = [np.ones((ny, nx))] + [gauss((ny, nx), pos[k][0] + dx, pos[k][1] + dy, sg) for k in names]
        coef, c2 = lin_fit(models, oot, ooterr)
        return c2
    best = min((minimize(chi_oot, [x0, y0, s0], method='Nelder-Mead', options={'xatol': 1e-4, 'fatol': 1e-3, 'maxiter': 3000})
                for x0 in (-0.3, 0, 0.3) for y0 in (-0.3, 0, 0.3) for s0 in (0.5, 0.8)), key=lambda r: r.fun)
    dx, dy, sg = best.x
    models = [np.ones((ny, nx))] + [gauss((ny, nx), pos[k][0] + dx, pos[k][1] + dy, sg) for k in names]
    coef, c2 = lin_fit(models, oot, ooterr)
    amps = dict(zip(names, map(float, coef[1:])))
    # chi2 profile for dx to get its uncertainty (1-sigma where dchi2=1, after rescaling chi2/dof to 1)
    dof = oot.size - len(models) - 3
    scale = c2 / dof
    prof = []
    for ddx in np.linspace(-0.6, 0.6, 61):
        r = minimize(lambda p: chi_oot([dx + ddx, p[0], p[1]]), [dy, sg], method='Nelder-Mead', options={'xatol': 1e-4, 'fatol': 1e-3})
        prof.append((ddx, (r.fun - c2) / scale))
    prof = np.array(prof); within = prof[prof[:, 1] < 1, 0]
    dx_err = float(max(within.max(), -within.min())) if len(within) else None
    cal = {'dx_px': float(dx), 'dy_px': float(dy), 'sigma_px': float(sg), 'chi2_red_oot': float(scale), 'dx_1sig_px': dx_err,
           'fitted_amps_e_s': amps, 'flux_ratio_target_plus_nb_over_host': float((amps['target'] + amps['neighbour']) / amps['host']),
           'gaia_ratio_target_plus_nb_over_host': float(10**(-0.4 * (17.1037 - 14.255)) + 10**(-0.4 * (18.3494 - 14.255)))}
    # ---- deficit stack
    W = np.array([1 / e**2 for x in items for e in x[2]]); D = np.array([d for x in items for d in x[1]])
    dmap = np.sum(W * D, axis=0) / np.sum(W, axis=0); derr = 1 / np.sqrt(np.sum(W, axis=0))
    nev = sum(len(x[1]) for x in items)
    if nev == 0:
        results['seasons']['season%d' % season] = {'error': 'no usable events', 'quarters': [x[0] for x in items]}; continue
    derr[~np.isfinite(derr)] = np.inf; dmap[~np.isfinite(dmap)] = 0.0
    cpos = {k: (pos[k][0] + dx, pos[k][1] + dy) for k in SRC}
    out = {'quarters': [x[0] for x in items], 'channel': items[0][5], 'n_events': nev, 'stamp': [ny, nx], 'wcs_calibration': cal,
           'gaia_px': pos, 'corrected_px': cpos, 'chi2_null': float(np.sum((dmap / derr)**2)), 'fits': {}}
    for sig_lab, sig in [('fitted', sg), ('0.45', 0.45), ('0.60', 0.60), ('0.75', 0.75), ('0.90', 0.90), ('1.10', 1.10)]:
        f = {}
        for k in ['target', 'neighbour', 'host', 'nb20', 'nb21']:
            g = gauss((ny, nx), cpos[k][0], cpos[k][1], sig)
            a, c = lin_fit([g], dmap, derr)
            f[k] = {'chi2': float(c), 'deficit_e_s': float(a[0])}
        # free position: inside stamp and anywhere (up to 2.5 px beyond)
        xs = np.arange(-2.5, nx + 2.5, 0.05); ys = np.arange(-2.5, ny + 2.5, 0.05)
        chi = np.full((len(ys), len(xs)), np.inf); amp = np.zeros_like(chi)
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                g = gauss((ny, nx), x, y, sig)
                if g.sum() < 0.02: continue
                a, c = lin_fit([g], dmap, derr); chi[j, i] = c; amp[j, i] = a[0]
        jb, ib = np.unravel_index(np.argmin(chi), chi.shape)
        inside = (xs[None, :] > -0.5) & (xs[None, :] < nx - 0.5) & (ys[:, None] > -0.5) & (ys[:, None] < ny - 0.5)
        chi_in = np.where(inside, chi, np.inf); ji, ii = np.unravel_index(np.argmin(chi_in), chi.shape)
        reg = chi < chi[jb, ib] + 11.83
        yy, xx = np.meshgrid(ys, xs, indexing='ij')
        f['free_any'] = {'px': [float(xs[ib]), float(ys[jb])], 'chi2': float(chi[jb, ib]), 'deficit_e_s': float(amp[jb, ib]),
                         'x_3sig': [float(xx[reg].min()), float(xx[reg].max())], 'y_3sig': [float(yy[reg].min()), float(yy[reg].max())]}
        f['free_inside'] = {'px': [float(xs[ii]), float(ys[ji])], 'chi2': float(chi_in[ji, ii]), 'deficit_e_s': float(amp[ji, ii])}
        # 1-sigma x range along the best row (dchi2 < 1 for one parameter, at fixed y=best)
        row = chi[jb]; ok = row < chi[jb, ib] + 1
        f['free_any']['x_1sig_at_best_y'] = [float(xs[ok].min()), float(xs[ok].max())]
        f['dchi2_target_minus_free'] = f['target']['chi2'] - f['free_any']['chi2']
        f['dchi2_neighbour_minus_free'] = f['neighbour']['chi2'] - f['free_any']['chi2']
        f['dchi2_target_minus_neighbour'] = f['target']['chi2'] - f['neighbour']['chi2']
        # offset of the target from the free position in arcsec, and the minimum offset consistent at 3 sigma
        tx, ty = cpos['target']
        dist = np.hypot(xx - tx, yy - ty)
        f['min_offset_target_to_3sig_region_arcsec'] = float(dist[reg].min() * 3.98)
        f['offset_target_to_best_arcsec'] = float(np.hypot(xs[ib] - tx, ys[jb] - ty) * 3.98)
        out['fits'][sig_lab] = f
    # ---- column ratio test in the target's row: observed deficit in last two columns vs expected for target/neighbour
    r_t = int(round(cpos['target'][1]))
    rows = [r for r in (r_t - 1, r_t, r_t + 1) if 0 <= r < ny]
    obs = {int(c): [float(dmap[r_t, c]), float(derr[r_t, c])] for c in range(nx)}
    ratio_obs = float(dmap[r_t, nx - 2] / dmap[r_t, nx - 1]) if dmap[r_t, nx - 1] != 0 else None
    exp = {}
    for k in ['target', 'neighbour']:
        for sig in (0.45, 0.6, 0.75, 0.9):
            g = gauss((ny, nx), cpos[k][0], cpos[k][1], sig)
            exp['%s_sig%.2f' % (k, sig)] = float(g[r_t, nx - 2] / g[r_t, nx - 1])
    # sum over rows r_t-1..r_t+1 for robustness
    colsum = dmap[rows].sum(axis=0); colerr = np.sqrt((derr[rows]**2).sum(axis=0))
    out['column_test'] = {'target_row': r_t, 'row_deficit_by_column_e_s': obs, 'ratio_col_nx-2_over_nx-1_observed': ratio_obs,
                          'ratio_err': float(np.hypot(derr[r_t, nx - 2] / dmap[r_t, nx - 1], dmap[r_t, nx - 2] * derr[r_t, nx - 1] / dmap[r_t, nx - 1]**2)) if dmap[r_t, nx - 1] != 0 else None,
                          'expected_ratio': exp, 'rows_summed_deficit_by_column': [float(v) for v in colsum], 'rows_summed_err': [float(v) for v in colerr],
                          'ratio_rows_summed': float(colsum[nx - 2] / colsum[nx - 1]), 'ratio_rows_summed_err': float(np.hypot(colerr[nx - 2] / colsum[nx - 1], colsum[nx - 2] * colerr[nx - 1] / colsum[nx - 1]**2))}
    results['seasons']['season%d' % season] = out
    print('season', season, 'Q', out['quarters'], 'cal', json.dumps(cal), '\n  fitted-sigma fits:', json.dumps(out['fits']['fitted']), '\n  column test:', json.dumps(out['column_test']))
json.dump(results, open(OUT, 'w'), indent=1, default=float)
