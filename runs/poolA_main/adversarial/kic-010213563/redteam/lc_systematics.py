#!/usr/bin/env python
"""Red-team light-curve checks for kic-010213563 (own code, no pipeline/ imports).

Target KBONUS LC (PSF + SAP), the TPF host's KBONUS LC (kic-010213566), per-event depths vs gaps,
per-quarter/season depths, secondary scan at all phases, Lomb-Scargle, inverted flux at P with a
transit-masked detrend, odd/even, and the duration/density arithmetic.
"""
import json, sys, numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares, LombScargle

OUT = sys.argv[1]
TARGET = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-010213563_kepler_v1.0_lc.fits'
HOST = '/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic010213563/dl/hlsp_kbonus-bkg_kepler_kepler_kic-010213566_kepler_v1.0_lc.fits'
P, T0, DUR_H = 36.76708318978774, 164.85416024494916, 5.625   # pipeline ephemeris (vet.json)
DUR = DUR_H / 24.0
res = {'ephemeris': {'P': P, 't0_bkjd': T0, 'dur_h': DUR_H, 'source': 'vet.json poolA_main'}}


def load(path, col='FLUX'):
    h = fits.open(path)
    d = h['LIGHTCURVE_STITCHED'].data
    t = np.array(d['TIME'], float); f = np.array(d[col], float); q = np.array(d['QUARTER'], int)
    qual = np.array(d['SAP_QUALITY'], int)
    m = np.isfinite(t) & np.isfinite(f) & (qual == 0) & (f > 0)
    t, f, q = t[m], f[m], q[m]
    # per-quarter median normalise
    fn = f.copy()
    for qq in np.unique(q):
        s = q == qq
        fn[s] = f[s] / np.nanmedian(f[s])
    return t, fn, q, f


def running_median(t, f, mask, win=1.0):
    """time-windowed running median, computed from unmasked points only, segments split at gaps>0.5 d."""
    trend = np.full_like(f, np.nan)
    order = np.argsort(t)
    t, f, mask = t[order], f[order], mask[order]
    gaps = np.where(np.diff(t) > 0.5)[0]
    edges = np.r_[0, gaps + 1, len(t)]
    for a, b in zip(edges[:-1], edges[1:]):
        tt, ff, mm = t[a:b], f[a:b], mask[a:b]
        good = ~mm
        for i in range(len(tt)):
            s = good & (np.abs(tt - tt[i]) < win / 2)
            if s.sum() >= 5:
                trend[a + i] = np.median(ff[s])
            else:
                s = good & (np.abs(tt - tt[i]) < win)
                trend[a + i] = np.median(ff[s]) if s.sum() >= 3 else np.nan
    out = np.full_like(trend, np.nan); out[order] = trend
    return out


def phase(t, P=P, t0=T0):
    return ((t - t0) / P + 0.5) % 1 - 0.5


def flatten(t, f, P=P, t0=T0, dur=DUR):
    ph = phase(t, P, t0)
    mask = np.abs(ph * P) < 0.75 * dur
    tr = running_median(t, f, mask)
    ff = f / tr
    ok = np.isfinite(ff)
    return t[ok], ff[ok], mask[ok]


def box_depth(t, ff, P=P, t0=T0, dur=DUR, nrand=2000, seed=1):
    ph = phase(t, P, t0) * P
    i = np.abs(ph) < dur / 2
    o = (np.abs(ph) > dur) & (np.abs(ph) < 3 * dur)
    dep = np.mean(ff[o]) - np.mean(ff[i])
    rng = np.random.default_rng(seed)
    fake = []
    for _ in range(nrand):
        dt = rng.uniform(2 * dur, P - 2 * dur)
        ph2 = ((t - t0 - dt) / P + 0.5) % 1 - 0.5
        ph2 *= P
        if abs(dt - P / 2) < 1.5 * dur:
            continue
        i2 = np.abs(ph2) < dur / 2; o2 = (np.abs(ph2) > dur) & (np.abs(ph2) < 3 * dur)
        if i2.sum() > 20 and o2.sum() > 40:
            fake.append(np.mean(ff[o2]) - np.mean(ff[i2]))
    fake = np.array(fake)
    err = 1.4826 * np.median(np.abs(fake - np.median(fake)))
    return dep, err, i.sum(), float(np.std(fake))


# ---------------------------------------------------------------- target PSF and SAP
t, fpsf, q, fraw_psf = load(TARGET, 'FLUX')
t2, fsap, q2, fraw_sap = load(TARGET, 'SAP_FLUX')
res['target'] = {'n_cadences': int(len(t)), 'median_psf_flux_e_s': float(np.median(fraw_psf)), 'median_sap_flux_e_s': float(np.median(fraw_sap))}
tt, ff, msk = flatten(t, fpsf)
d, e, n, sd = box_depth(tt, ff)
res['target']['psf_box_depth_ppm'] = [d * 1e6, e * 1e6, int(n)]
res['target']['psf_scatter_ppm'] = float(1.4826 * np.median(np.abs(ff[~msk] - 1)) * 1e6)
ts, fs, msks = flatten(t2, fsap)
d2, e2, n2, _ = box_depth(ts, fs)
res['target']['sap_box_depth_ppm'] = [d2 * 1e6, e2 * 1e6, int(n2)]
res['target']['sap_over_psf'] = d2 / d

# per-event depths (PSF), with local linear baseline, +-3 dur window
ph = phase(tt) * P
epochs = np.round((tt - T0) / P).astype(int)
ev = []
for n_ in np.unique(epochs):
    s = epochs == n_
    w = s & (np.abs(ph) < 3 * DUR)
    i = s & (np.abs(ph) < DUR / 2)
    o = s & (np.abs(ph) > DUR) & (np.abs(ph) < 3 * DUR)
    if i.sum() < 4 or o.sum() < 12:
        continue
    # local linear baseline from out-of-transit points
    c = np.polyfit(tt[o] - (T0 + n_ * P), ff[o], 1)
    base = np.polyval(c, tt[i] - (T0 + n_ * P))
    dep = np.mean(base - ff[i])
    err = np.std(ff[o] - np.polyval(c, tt[o] - (T0 + n_ * P))) / np.sqrt(i.sum())
    tc = T0 + n_ * P
    # distance to nearest gap >0.5 d in the full time series
    allt = np.sort(t)
    gi = np.where(np.diff(allt) > 0.5)[0]
    gap_edges = np.r_[allt[0], allt[gi], allt[gi + 1], allt[-1]]
    gd = float(np.min(np.abs(gap_edges - tc)))
    qq = int(np.median(q[np.abs(t - tc) < 1]))
    # position inside quarter (fraction)
    tq = t[q == qq]
    frac_in_q = float((tc - tq.min()) / (tq.max() - tq.min()))
    ev.append({'n': int(n_), 'tc_bkjd': float(tc), 'depth_ppm': float(dep * 1e6), 'err_ppm': float(err * 1e6), 'n_in': int(i.sum()),
               'n_out': int(o.sum()), 'gap_dist_d': gd, 'quarter': qq, 'season': (qq - 2) % 4 if qq > 0 else 3, 'frac_in_quarter': frac_in_q})
res['events'] = ev
dep = np.array([x['depth_ppm'] for x in ev]); err = np.array([x['err_ppm'] for x in ev])
w = 1 / err**2
wm = np.sum(w * dep) / np.sum(w)
res['events_summary'] = {'n': len(ev), 'weighted_mean_ppm': float(wm), 'weighted_err_ppm': float(1 / np.sqrt(np.sum(w))),
                         'chi2_const': float(np.sum((dep - wm)**2 / err**2)), 'dof': len(ev) - 1,
                         'n_positive': int((dep > 0).sum()), 'n_gt2sig': int((dep / err > 2).sum()),
                         'max_single_event_share': float(np.max(w * dep) / np.sum(w * dep)),
                         'leave_one_out_range_ppm': [float(min((np.sum(w * dep) - w[i] * dep[i]) / (np.sum(w) - w[i]) for i in range(len(ev)))),
                                                     float(max((np.sum(w * dep) - w[i] * dep[i]) / (np.sum(w) - w[i]) for i in range(len(ev))))]}
gd = np.array([x['gap_dist_d'] for x in ev])
for lab, s in [('within_1d_of_gap', gd < 1), ('beyond_1d_of_gap', gd >= 1), ('within_2d_of_gap', gd < 2)]:
    if s.sum():
        res['events_summary'][lab] = {'n': int(s.sum()), 'wmean_ppm': float(np.sum(w[s] * dep[s]) / np.sum(w[s])), 'err_ppm': float(1 / np.sqrt(np.sum(w[s])))}
# shallow events
res['events_summary']['shallow_events_lt_half_mean'] = [x for x in ev if x['depth_ppm'] < 0.5 * wm]
# per quarter and per season
pq = {}
for qq in sorted(set(x['quarter'] for x in ev)):
    s = np.array([x['quarter'] == qq for x in ev])
    pq[int(qq)] = [float(np.sum(w[s] * dep[s]) / np.sum(w[s])), float(1 / np.sqrt(np.sum(w[s]))), int(s.sum())]
res['depth_by_quarter_ppm'] = pq
pss = {}
for ss in range(4):
    s = np.array([x['season'] == ss for x in ev])
    if s.sum():
        pss[ss] = [float(np.sum(w[s] * dep[s]) / np.sum(w[s])), float(1 / np.sqrt(np.sum(w[s]))), int(s.sum())]
res['depth_by_season_ppm'] = pss
vals = np.array([v[0] for v in pss.values()]); errs = np.array([v[1] for v in pss.values()])
m = np.sum(vals / errs**2) / np.sum(1 / errs**2)
res['season_chi2'] = [float(np.sum((vals - m)**2 / errs**2)), len(vals) - 1]
# depth vs position in quarter (thermal ramps happen early in each quarter / after downlinks)
fq = np.array([x['frac_in_quarter'] for x in ev])
for lab, s in [('first_10pct_of_quarter', fq < 0.1), ('rest_of_quarter', fq >= 0.1)]:
    if s.sum():
        res['events_summary'][lab] = {'n': int(s.sum()), 'wmean_ppm': float(np.sum(w[s] * dep[s]) / np.sum(w[s])), 'err_ppm': float(1 / np.sqrt(np.sum(w[s])))}

# odd / even
odd = np.array([x['n'] % 2 == 1 for x in ev])
do = np.sum(w[odd] * dep[odd]) / np.sum(w[odd]); de = np.sum(w[~odd] * dep[~odd]) / np.sum(w[~odd])
eo = 1 / np.sqrt(np.sum(w[odd])); ee = 1 / np.sqrt(np.sum(w[~odd]))
res['odd_even'] = {'odd_ppm': float(do), 'odd_err': float(eo), 'even_ppm': float(de), 'even_err': float(ee), 'sigma': float((do - de) / np.hypot(eo, ee))}

# secondary scan at all phases (box of transit duration), both polarities
phs = np.linspace(-0.5, 0.5, 400)
scan = []
oot = np.abs(ph) > 1.5 * DUR
ref = np.mean(ff[oot])
for p_ in phs:
    dph = ((ph / P - p_) + 0.5) % 1 - 0.5
    i = np.abs(dph * P) < DUR / 2
    if i.sum() > 30:
        scan.append((p_, (ref - np.mean(ff[i])) * 1e6, i.sum()))
scan = np.array(scan)
far = np.abs(scan[:, 0]) > 0.02
sig = 1.4826 * np.median(np.abs(scan[far, 1] - np.median(scan[far, 1])))
i5 = np.argmin(np.abs(scan[:, 0] - 0.5)) if np.any(np.abs(scan[:, 0]) > 0.49) else np.argmin(np.abs(np.abs(scan[:, 0]) - 0.5))
res['secondary_scan'] = {'box_noise_ppm_from_scan': float(sig), 'depth_at_phase_0p5_ppm': float(scan[i5, 1]),
                         'max_dip_outside_transit': [float(scan[far][np.argmax(scan[far, 1]), 0]), float(np.max(scan[far, 1]))],
                         'max_brightening_outside_transit': [float(scan[far][np.argmin(scan[far, 1]), 0]), float(np.min(scan[far, 1]))],
                         'z_std': float(np.std(scan[far, 1] / sig)), 'z_max': float(np.max(scan[far, 1] / sig)), 'z_min': float(np.min(scan[far, 1] / sig))}
# wider secondary box (12 h) at phase 0.5
dph = ((ph / P - 0.5) + 0.5) % 1 - 0.5
i = np.abs(dph * P) < 0.25
res['secondary_scan']['depth_12h_box_phase0p5_ppm'] = float((ref - np.mean(ff[i])) * 1e6)
res['secondary_scan']['n_12h'] = int(i.sum())

# Lomb-Scargle on out-of-transit quarter-normalised (not flattened) flux
oo = np.abs(phase(t) * P) > DUR
freq = np.linspace(1 / 60., 1 / 0.1, 200000)
ls = LombScargle(t[oo], fpsf[oo])
pw = ls.power(freq)
top = np.argsort(pw)[::-1]
peaks = []
for i in top:
    per = 1 / freq[i]
    if all(abs(per / p_ - 1) > 0.05 for p_, _, _ in peaks):
        amp = np.sqrt(pw[i]) * np.std(fpsf[oo]) * np.sqrt(2) * 1e6  # rough half-amplitude
        peaks.append((float(per), float(pw[i]), float(amp)))
    if len(peaks) >= 5:
        break
res['lomb_scargle'] = {'top_peaks_period_power_halfamp_ppm': peaks,
                       'power_at_P': float(ls.power(np.array([1 / P]))[0]), 'power_at_2P': float(ls.power(np.array([1 / (2 * P)]))[0]),
                       'power_at_P_over_2': float(ls.power(np.array([2 / P]))[0])}

# inverted flux at P with masked detrend: BLS 35-39 d
inv = 2 - ff
bls = BoxLeastSquares(tt, inv)
periods = np.exp(np.linspace(np.log(35), np.log(39), 4000))
r = bls.power(periods, [0.15, 0.2, 0.25, 0.3], objective='snr')
res['inverted_at_P'] = {'best_period': float(periods[np.argmax(r.power)]), 'best_depth_ppm': float(r.depth[np.argmax(r.power)] * 1e6),
                        'best_snr_stat': float(np.max(r.power)), 'depth_at_P_ppm': float(bls.compute_stats(P, DUR, T0)['depth'][0] * 1e6)}
# real flux same window for comparison
r2 = bls.power(periods, [0.15, 0.2, 0.25, 0.3], objective='snr')
bls2 = BoxLeastSquares(tt, ff)
r2 = bls2.power(periods, [0.15, 0.2, 0.25, 0.3], objective='snr')
res['real_at_P_same_grid'] = {'best_period': float(periods[np.argmax(r2.power)]), 'best_snr_stat': float(np.max(r2.power))}

# ---------------------------------------------------------------- host star KIC 10213566
th, fh, qh, frawh = load(HOST, 'FLUX')
th2, fh2, qh2, frawh2 = load(HOST, 'SAP_FLUX')
res['host'] = {'kic': 10213566, 'median_psf_flux_e_s': float(np.median(frawh)), 'median_sap_flux_e_s': float(np.median(frawh2)), 'n_cadences': int(len(th))}
tth, ffh, mh = flatten(th, fh)
dh, eh, nh, _ = box_depth(tth, ffh)
res['host']['psf_box_depth_ppm'] = [dh * 1e6, eh * 1e6, int(nh)]
res['host']['psf_scatter_ppm'] = float(1.4826 * np.median(np.abs(ffh[~mh] - 1)) * 1e6)
tsh, fsh, _ = flatten(th2, fh2)
dh2, eh2, nh2, _ = box_depth(tsh, fsh)
res['host']['sap_box_depth_ppm'] = [dh2 * 1e6, eh2 * 1e6, int(nh2)]
# if the dip were on the host and leaking: target PSF flux 2195 e/s, host 26513 e/s. The host itself would show
# at least depth_target * F_target_psf / F_host_psf  (if ALL target flux were host light)
res['host']['min_expected_if_on_host_ppm'] = float(d * 1e6 * np.median(fraw_psf) / np.median(frawh))
# host's own BLS 1-50 d quick look (is the host a variable/EB?)
bh = BoxLeastSquares(tth, ffh)
per = np.exp(np.linspace(np.log(1), np.log(50), 60000))
rh = bh.power(per, [0.1, 0.2, 0.3], objective='snr')
res['host']['bls_top'] = [float(per[np.argmax(rh.power)]), float(np.max(rh.power)), float(rh.depth[np.argmax(rh.power)] * 1e6)]
lsh = LombScargle(th, fh); pwh = lsh.power(freq)
res['host']['ls_top_period_power'] = [float(1 / freq[np.argmax(pwh)]), float(np.max(pwh))]

# ---------------------------------------------------------------- duration / density arithmetic
G = 6.674e-11; Msun = 1.989e30; Rsun = 6.957e8; day = 86400.
Pd = P * day
def rho_from_T14_T23(T14_h, T23_h, k):
    T14 = T14_h * 3600; T23 = T23_h * 3600
    s14 = np.sin(np.pi * T14 / Pd)**2; s23 = np.sin(np.pi * T23 / Pd)**2
    b2 = ((1 - k)**2 - (s23 / s14) * (1 + k)**2) / (1 - s23 / s14)
    aR2 = ((1 + k)**2 - b2 * (1 - s14)) / s14
    rho = 3 * np.pi / (G * Pd**2) * aR2**1.5
    return float(rho / 1000), float(np.sqrt(max(b2, 0))), float(np.sqrt(aR2))   # g/cc, b, a/R
def rho_b0(T14_h, k):
    T14 = T14_h * 3600
    aR = (1 + k) / np.sin(np.pi * T14 / Pd)
    return float(3 * np.pi / (G * Pd**2) * aR**3 / 1000)
# dilution: neighbour G=18.35 vs target G=17.10 -> flux ratio
fr = 10**(-0.4 * (18.349396 - 17.103693))
res['dilution'] = {'neighbour_over_target_flux': float(fr), 'depth_undiluted_if_on_target_ppm': float(d * 1e6 * (1 + fr)),
                   'depth_on_neighbour_if_on_neighbour_ppm': float(d * 1e6 * (1 + 1 / fr))}
dens = {}
for lab, dep_ppm in [('as_measured', d * 1e6), ('undiluted_on_target', d * 1e6 * (1 + fr)), ('on_neighbour', d * 1e6 * (1 + 1 / fr))]:
    k = np.sqrt(dep_ppm * 1e-6)
    dens[lab] = {'k': float(k), 'rho_b0_T14_5p85h_gcc': rho_b0(5.85, k), 'rho_b0_T14_6p35h_gcc': rho_b0(6.35, k),
                 'rho_from_trapezoid_T14_5p85_T23_4p43': rho_from_T14_T23(5.85, 4.43, k),
                 'rho_from_trapezoid_T14_6p35_T23_4p43': rho_from_T14_T23(6.35, 4.43, k)}
res['density'] = dens
res['density']['reference'] = {'K_dwarf_0p79Rsun_0p80Msun_gcc': float(0.80 * Msun / (4 / 3 * np.pi * (0.79 * Rsun)**3) / 1000),
                               'G_dwarf_0p95Rsun_0p95Msun_gcc': float(0.95 * Msun / (4 / 3 * np.pi * (0.95 * Rsun)**3) / 1000),
                               'M_dwarf_0p5Rsun_0p5Msun_gcc': float(0.5 * Msun / (4 / 3 * np.pi * (0.5 * Rsun)**3) / 1000),
                               'sun_gcc': float(Msun / (4 / 3 * np.pi * Rsun**3) / 1000)}
json.dump(res, open(OUT, 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in res.items() if k != 'events'}, indent=1, default=float))
