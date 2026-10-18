#!/usr/bin/env python
"""Red-team transit-shape fit for kic-010213563 (own code; batman + emcee, 30-min supersampling).

Fits rp/rs, b, a/rs (-> stellar density), t0, P on the transit-masked flattened KBONUS PSF flux within
+-1 d of each transit. Quadratic limb darkening fixed for a ~5000 K dwarf. Also a fit with the depth
undiluted by the unmodelled 0.97" neighbour. Reports the posterior on rho_star and compares with a
0.79 Rsun K dwarf (2.3 g/cc), a 0.9 Rsun G dwarf (1.6 g/cc), and what the 'dip on neighbour' case needs.
"""
import json, sys, numpy as np
from astropy.io import fits
import batman, emcee

OUT = sys.argv[1]
TARGET = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-010213563_kepler_v1.0_lc.fits'
P0, T00, DUR_H = 36.76708318978774, 164.85416024494916, 5.625
DUR = DUR_H / 24
DROP_EVENT0 = '--drop0' in sys.argv

h = fits.open(TARGET)['LIGHTCURVE_STITCHED'].data
t = np.array(h['TIME'], float); f = np.array(h['FLUX'], float); q = np.array(h['QUARTER'], int); qual = np.array(h['SAP_QUALITY'], int)
m = np.isfinite(t) & np.isfinite(f) & (qual == 0) & (f > 0)
t, f, q = t[m], f[m], q[m]
for qq in np.unique(q):
    s = q == qq; f[s] /= np.nanmedian(f[s])
ph = ((t - T00) / P0 + 0.5) % 1 - 0.5
intr = np.abs(ph * P0) < 0.75 * DUR
# running median 1 d with transits masked, segments split at gaps
order = np.argsort(t); t, f, intr = t[order], f[order], intr[order]
trend = np.full_like(f, np.nan)
gaps = np.where(np.diff(t) > 0.5)[0]; edges = np.r_[0, gaps + 1, len(t)]
for a, b in zip(edges[:-1], edges[1:]):
    tt, ff, mm = t[a:b], f[a:b], intr[a:b]
    for i in range(len(tt)):
        s = (~mm) & (np.abs(tt - tt[i]) < 0.5)
        if s.sum() < 5: s = (~mm) & (np.abs(tt - tt[i]) < 1.0)
        trend[a + i] = np.median(ff[s]) if s.sum() >= 3 else np.nan
ff = f / trend
ok = np.isfinite(ff)
t, ff, ph = t[ok], ff[ok], ph[ok]
# clip upward outliers out of transit (cosmic rays etc.)
oot = ~(np.abs(ph * P0) < 0.75 * DUR)
mad = 1.4826 * np.median(np.abs(ff[oot] - 1))
keep = ~(oot & (np.abs(ff - 1) > 5 * mad))
t, ff, ph = t[keep], ff[keep], ph[keep]
win = np.abs(ph * P0) < 1.0
if DROP_EVENT0:
    win &= np.round((t - T00) / P0) != 0
tw, fw = t[win], ff[win]
err = mad
print('points in fit window', win.sum(), 'per-point scatter ppm', mad * 1e6)

params = batman.TransitParams()
params.t0 = T00; params.per = P0; params.rp = 0.047; params.a = 40.; params.inc = 89.; params.ecc = 0.; params.w = 90.
params.limb_dark = 'quadratic'; params.u = [0.50, 0.20]   # ~5000 K dwarf, Kepler band
model = batman.TransitModel(params, tw, supersample_factor=15, exp_time=1765.5 / 86400)

G = 6.674e-11; day = 86400.
def rho_from_aR(aR, per):
    return 3 * np.pi / (G * (per * day)**2) * aR**3 / 1000   # g/cc

def lnprob(p):
    t0, per, rp, b, lnaR, lnf = p
    aR = np.exp(lnaR)
    if not (0 < rp < 0.3 and 0 <= b < 1 + rp and 2 < aR < 400 and -2 < lnf < 2): return -np.inf
    inc = np.degrees(np.arccos(b / aR))
    params.t0 = t0; params.per = per; params.rp = rp; params.a = aR; params.inc = inc
    mod = model.light_curve(params)
    s2 = (err * np.exp(lnf))**2
    return -0.5 * np.sum((fw - mod)**2 / s2 + np.log(2 * np.pi * s2))

ndim, nw = 6, 32
p0 = np.array([T00, P0, 0.047, 0.3, np.log(40.), 0.0])
pos = p0 + np.array([1e-3, 1e-5, 2e-3, 0.1, 0.1, 0.02]) * np.random.default_rng(1).standard_normal((nw, ndim))
pos[:, 3] = np.abs(pos[:, 3])
sampler = emcee.EnsembleSampler(nw, ndim, lnprob)
sampler.run_mcmc(pos, 3000, progress=False)
chain = sampler.get_chain(discard=1000, flat=True)
lp = sampler.get_log_prob(discard=1000, flat=True)
t0s, pers, rps, bs, lnaRs, lnfs = chain.T
aRs = np.exp(lnaRs)
rhos = rho_from_aR(aRs, pers)
# durations
inc = np.arccos(bs / aRs)
T14 = pers / np.pi * np.arcsin(np.sqrt(np.clip((1 + rps)**2 - bs**2, 0, None)) / (aRs * np.sin(inc))) * 24
T23 = pers / np.pi * np.arcsin(np.sqrt(np.clip((1 - rps)**2 - bs**2, 0, None)) / (aRs * np.sin(inc))) * 24
def pct(x): return [float(np.percentile(x, 16)), float(np.percentile(x, 50)), float(np.percentile(x, 84))]
best = chain[np.argmax(lp)]
res = {'drop_event0': DROP_EVENT0, 'n_points': int(win.sum()), 'scatter_ppm': float(mad * 1e6), 'limb_darkening_u': params.u,
       't0_bkjd': pct(t0s), 'period_d': pct(pers), 'rp_rs': pct(rps), 'depth_rp2_ppm': pct(rps**2 * 1e6), 'b': pct(bs), 'a_rs': pct(aRs),
       'rho_star_gcc': pct(rhos), 'rho_star_gcc_95pct_upper': float(np.percentile(rhos, 95)), 'rho_star_gcc_99pct_upper': float(np.percentile(rhos, 99)),
       'T14_h': pct(T14), 'T23_h': pct(T23), 'T23_over_T14': pct(T23 / T14), 'noise_scale': pct(np.exp(lnfs)),
       'best_fit': dict(zip(['t0', 'per', 'rp', 'b', 'aR', 'lnf'], map(float, best))), 'max_lnprob': float(lp.max()),
       'frac_posterior_rho_gt_2p0': float(np.mean(rhos > 2.0)), 'frac_posterior_rho_gt_1p5': float(np.mean(rhos > 1.5)), 'frac_posterior_b_gt_0p7': float(np.mean(bs > 0.7)),
       'reference_rho_gcc': {'K_dwarf_0p79Rsun_0p80Msun': 2.29, 'G_dwarf_0p90Rsun_0p90Msun': 1.74, 'G_dwarf_0p95Rsun_0p95Msun': 1.56, 'sun': 1.41}}
# fixed-b comparisons: chi2 at b=0 vs best
def chi2_at(t0, per, rp, b, aR):
    params.t0 = t0; params.per = per; params.rp = rp; params.a = aR; params.inc = np.degrees(np.arccos(b / aR))
    return float(np.sum((fw - model.light_curve(params))**2 / err**2))
res['chi2_best'] = chi2_at(*best[:5].tolist()[:4], np.exp(best[4]))
# best model constrained to K-dwarf density (aR from rho=2.29) : profile over rp,b,t0 quickly via emcee-lite grid
aR_K = (2.29 * 1000 * G * (P0 * day)**2 / (3 * np.pi))**(1 / 3)
aR_G = (1.74 * 1000 * G * (P0 * day)**2 / (3 * np.pi))**(1 / 3)
res['aR_for_K_dwarf'] = float(aR_K); res['aR_for_G_dwarf_0p9'] = float(aR_G)
from scipy.optimize import minimize
for lab, aRfix in [('K_dwarf_2p29', aR_K), ('G_dwarf_1p74', aR_G)]:
    def c2(p):
        t0, rp, b = p
        if not (0 < rp < 0.3 and 0 <= b < 1): return 1e30
        return chi2_at(t0, P0, rp, b, aRfix)
    r = min((minimize(c2, [T00, 0.05, b0], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-3, 'maxiter': 4000}) for b0 in (0.0, 0.3, 0.6, 0.85)), key=lambda r: r.fun)
    res['fixed_density_' + lab] = {'chi2': float(r.fun), 'dchi2_vs_free': float(r.fun - res['chi2_best']), 't0': float(r.x[0]), 'rp_rs': float(r.x[1]), 'b': float(r.x[2])}
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
