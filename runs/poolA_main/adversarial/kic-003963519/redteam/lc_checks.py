"""Red-team light-curve checks for kic-003963519 (own code, no pipeline imports)."""
import json, numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.optimize import least_squares

FITS = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-003963519_kepler_v1.0_lc.fits'
OUT = '/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic003963519/out/lc_checks.json'
P = 6.815624034742515; T0 = 138.28256264993922; DUR = 2.4375 / 24.0
rng = np.random.default_rng(42)
res = {'P': P, 'T0_bkjd': T0, 'dur_d': DUR}

f = fits.open(FITS)
d = f['LIGHTCURVE_STITCHED'].data
t_all = d['TIME']; q_all = d['QUARTER']; qual = d['SAP_QUALITY']
qext = {f[i].header['QUARTER']: i for i in range(1, len(f)) if f[i].header.get('EXTNAME', '').startswith('LIGHTCURVE_Q')}

def detrend(t, y, mask_transit=True, win=1.0, gap=0.5):
    ph = ((t - T0 + P / 2) % P) - P / 2
    intr = np.abs(ph) < DUR * 0.75 if mask_transit else np.zeros_like(t, bool)
    order = np.argsort(t); t = t[order]; y = y[order]; intr = intr[order]
    trend = np.full_like(y, np.nan)
    br = np.where(np.diff(t) > gap)[0] + 1
    for s in np.split(np.arange(len(t)), br):
        ts, ys, ms = t[s], y[s], intr[s]
        tt, yy = ts[~ms], ys[~ms]
        for j, i in enumerate(s):
            m = np.abs(tt - ts[j]) < win / 2
            if m.sum() > 8:
                trend[i] = np.median(yy[m])
    return t, y / trend - 1

def box_depth(t, r, P=P, t0=T0, dur=DUR, nnull=500):
    ph = ((t - t0 + P / 2) % P) - P / 2
    a = np.abs(ph) < dur / 2
    b = (np.abs(ph) > dur) & (np.abs(ph) < 3 * dur)
    if a.sum() < 5 or b.sum() < 5:
        return np.nan, np.nan, int(a.sum())
    dep = -(np.mean(r[a]) - np.mean(r[b]))
    nulls = []
    for k in range(nnull):
        sh = rng.uniform(0.08 * P, 0.92 * P)
        ph2 = ((t - t0 - sh + P / 2) % P) - P / 2
        aa = np.abs(ph2) < dur / 2; bb = (np.abs(ph2) > dur) & (np.abs(ph2) < 3 * dur)
        if aa.sum() > 5 and bb.sum() > 5:
            nulls.append(-(np.mean(r[aa]) - np.mean(r[bb])))
    return dep * 1e6, (np.std(nulls) * 1e6 if nnull > 1 else np.nan), int(a.sum())

# ---------- 0. header info per quarter ----------
season = {0: 3, 1: 3, 2: 0, 3: 1, 4: 2, 5: 3, 6: 0, 7: 1, 8: 2, 9: 3, 10: 0, 11: 1, 12: 2, 13: 3, 14: 0, 15: 1, 16: 2, 17: 3}
hdr = {}
for q, i in qext.items():
    h = f[i].header
    hdr[int(q)] = dict(PSFFRAC=h['PSFFRAC'], CROWDSAP=h['CROWDSAP'], FLFRCSAP=h['FLFRCSAP'], NPIXSAP=h['NPIXSAP'], CHANNEL=h['CHANNEL'], season=season[q])
res['header_by_quarter'] = hdr

# ---------- 1. per-quarter / per-season depth in each flux column ----------
good = (qual == 0)
cols = {}
for col in ['FLUX', 'SAP_FLUX', 'PSF_FLAT_FLUX']:
    y = d[col].astype(float)
    ok = good & np.isfinite(y) & np.isfinite(t_all) & (y > 0)
    tt = []; rr = []; qq = []
    for q in np.unique(q_all[ok]):
        m = ok & (q_all == q)
        yq = y[m] / np.nanmedian(y[m])
        tq, rq = detrend(t_all[m], yq, mask_transit=(col != 'PSF_FLAT_FLUX'))
        tt.append(tq); rr.append(rq); qq.append(np.full(len(tq), q))
    tt = np.concatenate(tt); rr = np.concatenate(rr); qq = np.concatenate(qq)
    okr = np.isfinite(rr)
    mad = 1.4826 * np.median(np.abs(rr[okr] - np.median(rr[okr])))
    okr &= (rr < 4 * mad) & (rr > -10 * mad)
    tt, rr, qq = tt[okr], rr[okr], qq[okr]
    cols[col] = (tt, rr, qq)
    dep, err, n = box_depth(tt, rr)
    per_q = {}
    for q in np.unique(qq):
        m = qq == q
        dq, eq, nq = box_depth(tt[m], rr[m], nnull=300)
        per_q[int(q)] = [round(dq, 1), round(eq, 1), nq]
    per_s = {}
    for s in [0, 1, 2, 3]:
        m = np.isin(qq, [q for q in season if season[q] == s])
        if m.sum() > 100:
            ds, es, ns = box_depth(tt[m], rr[m], nnull=400); per_s[s] = [round(ds, 1), round(es, 1), ns]
    res[col] = {'depth_all_ppm': round(dep, 1), 'err_ppm': round(err, 1), 'n_in': n, 'rms_ppm': float(np.std(rr) * 1e6),
                'per_quarter': per_q, 'per_season': per_s, 'median_flux_e_s': float(np.nanmedian(y[ok]))}
    print(col, 'all', round(dep), '+-', round(err), 'rms', round(np.std(rr) * 1e6), 'per season', per_s)

# PSF_FLUX_NOVA (per-quarter extensions)
tt = []; rr = []; qq = []
for q, i in qext.items():
    dd = f[i].data
    y = dd['PSF_FLUX_NOVA'].astype(float); t = dd['TIME']; ql = dd['SAP_QUALITY']
    ok = (ql == 0) & np.isfinite(y) & np.isfinite(t) & (y > 0)
    if ok.sum() < 200:
        continue
    tq, rq = detrend(t[ok], y[ok] / np.nanmedian(y[ok]))
    tt.append(tq); rr.append(rq); qq.append(np.full(len(tq), q))
tt = np.concatenate(tt); rr = np.concatenate(rr); qq = np.concatenate(qq)
okr = np.isfinite(rr); mad = 1.4826 * np.median(np.abs(rr[okr] - np.median(rr[okr]))); okr &= (rr < 4 * mad) & (rr > -10 * mad)
tt, rr, qq = tt[okr], rr[okr], qq[okr]
dep, err, n = box_depth(tt, rr)
per_s = {}
for s in [0, 1, 2, 3]:
    m = np.isin(qq, [q for q in season if season[q] == s])
    if m.sum() > 100:
        ds, es, ns = box_depth(tt[m], rr[m], nnull=300); per_s[s] = [round(ds, 1), round(es, 1), ns]
res['PSF_FLUX_NOVA'] = {'depth_all_ppm': round(dep, 1), 'err_ppm': round(err, 1), 'n_in': n, 'per_season': per_s, 'n_points': int(len(tt))}
print('PSF_FLUX_NOVA all', round(dep), '+-', round(err), per_s)

# depth vs PSFFRAC / CROWDSAP / NPIXSAP (PSF FLUX per quarter)
pq = res['FLUX']['per_quarter']
res['depth_vs_header'] = {}
for key in ['PSFFRAC', 'CROWDSAP', 'NPIXSAP', 'FLFRCSAP']:
    x = np.array([hdr[q][key] for q in pq], float); yv = np.array([pq[q][0] for q in pq]); ev = np.array([pq[q][1] for q in pq])
    w = 1 / ev**2
    A = np.vstack([x, np.ones_like(x)]).T
    coef = np.linalg.lstsq(A * np.sqrt(w)[:, None], yv * np.sqrt(w), rcond=None)[0]
    cov = np.linalg.inv((A * w[:, None]).T @ A)
    res['depth_vs_header'][key] = {'slope_ppm_per_unit': float(coef[0]), 'slope_err': float(np.sqrt(cov[0, 0])), 'intercept': float(coef[1])}
    print('depth vs', key, 'slope', round(coef[0]), '+-', round(np.sqrt(cov[0, 0])))
yv = np.array([pq[q][0] for q in pq]); ev = np.array([pq[q][1] for q in pq]); w = 1 / ev**2
res['depth_vs_header']['chi2_const_quarters'] = float(np.sum(w * (yv - np.sum(w * yv) / np.sum(w))**2)); res['depth_vs_header']['dof'] = int(len(yv) - 1)
print('per-quarter chi2 const', res['depth_vs_header']['chi2_const_quarters'], 'dof', len(yv) - 1)

# ---------- 2. per-event depths, gap proximity, random-ephemeris calibration ----------
tt, rr, qq = cols['FLUX']
def per_event(t, r, t0, P=P, dur=DUR):
    n0 = int(np.floor((t.min() - t0) / P)); n1 = int(np.ceil((t.max() - t0) / P))
    out = []
    for n in range(n0, n1 + 1):
        tc = t0 + n * P
        a = np.abs(t - tc) < dur / 2
        b = (np.abs(t - tc) > dur) & (np.abs(t - tc) < 4 * dur)
        if a.sum() >= 2 and b.sum() >= 10:
            dep = -(np.mean(r[a]) - np.mean(r[b]))
            err = np.std(r[b]) * np.sqrt(1 / a.sum() + 1 / b.sum())
            out.append((tc, dep, err, a.sum(), b.sum()))
    return np.array(out)
ev = per_event(tt, rr, T0)
dep_w = np.sum(ev[:, 1] / ev[:, 2]**2) / np.sum(1 / ev[:, 2]**2)
chi2 = np.sum(((ev[:, 1] - dep_w) / ev[:, 2])**2)
res['per_event'] = {'n_events': int(len(ev)), 'weighted_depth_ppm': float(dep_w * 1e6), 'chi2_const': float(chi2), 'dof': int(len(ev) - 1),
                    'frac_positive': float(np.mean(ev[:, 1] > 0)), 'max_single_snr2_frac': float(np.max((ev[:, 1] / ev[:, 2])**2) / np.sum((ev[:, 1] / ev[:, 2])**2))}
cal = []
for k in range(150):
    sh = rng.uniform(0.1 * P, 0.9 * P)
    e2 = per_event(tt, rr, T0 + sh)
    dw = np.sum(e2[:, 1] / e2[:, 2]**2) / np.sum(1 / e2[:, 2]**2)
    cal.append(np.sum(((e2[:, 1] - dw) / e2[:, 2])**2) / (len(e2) - 1))
cal = np.array(cal)
res['per_event']['chi2_dof_random_ephemerides'] = {'median': float(np.median(cal)), 'p16': float(np.percentile(cal, 16)), 'p84': float(np.percentile(cal, 84)), 'max': float(cal.max()),
                                                   'real_chi2_dof': float(chi2 / (len(ev) - 1)), 'frac_random_ge_real': float(np.mean(cal >= chi2 / (len(ev) - 1)))}
print('per-event chi2/dof real', chi2 / (len(ev) - 1), 'random median', np.median(cal), 'p84', np.percentile(cal, 84), 'frac>=real', np.mean(cal >= chi2 / (len(ev) - 1)))
tsorted = np.sort(tt); gaps = np.where(np.diff(tsorted) > 0.5)[0]
edges = np.concatenate([[tsorted[0]], tsorted[gaps], tsorted[gaps + 1], [tsorted[-1]]])
gd = np.array([np.min(np.abs(edges - tc)) for tc in ev[:, 0]])
def wmean(x, e): return np.sum(x / e**2) / np.sum(1 / e**2), np.sqrt(1 / np.sum(1 / e**2))
for lab, near in [('near_gap_1d', gd < 1.0), ('near_gap_2d', gd < 2.0)]:
    dn, en = wmean(ev[near, 1], ev[near, 2]); df_, ef = wmean(ev[~near, 1], ev[~near, 2])
    res['per_event'][lab] = {'n': int(near.sum()), 'depth_ppm': float(dn * 1e6), 'err': float(en * 1e6)}
    res['per_event'][lab.replace('near', 'far')] = {'n': int((~near).sum()), 'depth_ppm': float(df_ * 1e6), 'err': float(ef * 1e6)}
    print(lab, near.sum(), round(dn * 1e6), '+-', round(en * 1e6), ' far', (~near).sum(), round(df_ * 1e6), '+-', round(ef * 1e6))
# top-5 strongest events: are they special?
top = np.argsort(ev[:, 1] / ev[:, 2])[::-1][:6]
res['per_event']['strongest_events'] = [[float(ev[i, 0]), float(ev[i, 1] * 1e6), float(ev[i, 2] * 1e6), float(gd[i])] for i in top]
# drop the 10 strongest events: does the depth survive?
keep = np.ones(len(ev), bool); keep[np.argsort(ev[:, 1] / ev[:, 2])[::-1][:10]] = False
dk, ek = wmean(ev[keep, 1], ev[keep, 2])
res['per_event']['depth_without_10_strongest'] = [float(dk * 1e6), float(ek * 1e6)]
print('depth without 10 strongest events', round(dk * 1e6), '+-', round(ek * 1e6))
# per-event depth vs position within quarter (thermal ramps): first third vs last third of each quarter
qs = np.array([qq[np.argmin(np.abs(tt - tc))] for tc in ev[:, 0]])
frac_in_q = np.zeros(len(ev))
for q in np.unique(qs):
    m = qs == q; tq = tt[qq == q]
    frac_in_q[m] = (ev[m, 0] - tq.min()) / (tq.max() - tq.min() + 1e-9)
for lab, m in [('first_third_of_quarter', frac_in_q < 1 / 3), ('middle_third', (frac_in_q >= 1 / 3) & (frac_in_q < 2 / 3)), ('last_third', frac_in_q >= 2 / 3)]:
    dd_, ee_ = wmean(ev[m, 1], ev[m, 2]); res['per_event'][lab] = [int(m.sum()), float(dd_ * 1e6), float(ee_ * 1e6)]
    print(lab, m.sum(), round(dd_ * 1e6), '+-', round(ee_ * 1e6))
# per-year t0 shift
t0s = []
for (a, b) in [(100, 470), (470, 840), (840, 1210), (1210, 1600)]:
    m = (tt > a) & (tt < b)
    shifts = np.linspace(-0.08, 0.08, 161)
    deps = [box_depth(tt[m], rr[m], t0=T0 + s, nnull=1)[0] for s in shifts]
    best = shifts[int(np.nanargmax(deps))]
    t0s.append([a, b, float(best * 24 * 60), float(np.nanmax(deps))])
res['per_year_best_t0_shift_min'] = t0s
print('per-year t0 shift (min) and depth', t0s)

# ---------- 3. quality flags, background, PSF chi2 and centroids in transit ----------
ph = ((t_all - T0 + P / 2) % P) - P / 2
intr_all = np.abs(ph) < DUR / 2
okt = np.isfinite(t_all)
res['quality'] = {'frac_flagged_in_transit': float(np.mean(qual[okt & intr_all] != 0)), 'frac_flagged_overall': float(np.mean(qual[okt] != 0)),
                  'bits_present': sorted(set(int(b) for b in np.unique(qual[okt]) if b))}
bk_in = []; rc_in = []; cen = {}
for q, i in qext.items():
    dd = f[i].data; t = dd['TIME']; ql = dd['SAP_QUALITY']; ok = (ql == 0) & np.isfinite(t)
    ph = ((t - T0 + P / 2) % P) - P / 2
    a = ok & (np.abs(ph) < DUR / 2); b = ok & (np.abs(ph) > DUR) & (np.abs(ph) < 3 * DUR)
    bkg = dd['SAP_BKG'].astype(float); rc = dd['RED_CHI2'].astype(float)
    if np.isfinite(bkg[b]).sum() > 50:
        bk_in.append(np.nanmedian(bkg[a]) / np.nanmedian(bkg[b]) - 1)
    if np.isfinite(rc[b]).sum() > 50:
        rc_in.append(np.nanmedian(rc[a]) / np.nanmedian(rc[b]) - 1)
    # centroid shift in transit (local: each event vs +-1 d neighbours), in pixels
    cc = dd['CENTROID_COLUMN'].astype(float); cr = dd['CENTROID_ROW'].astype(float)
    dcol = []; drow = []; wts = []
    for n in np.unique(np.round((t[a] - T0) / P).astype(int)):
        tc = T0 + n * P
        aa = ok & (np.abs(t - tc) < DUR / 2); bb = ok & (np.abs(t - tc) > DUR) & (np.abs(t - tc) < 1.0)
        if aa.sum() >= 2 and bb.sum() >= 15 and np.isfinite(cc[bb]).sum() > 10:
            dcol.append(np.nanmean(cc[aa]) - np.nanmedian(cc[bb])); drow.append(np.nanmean(cr[aa]) - np.nanmedian(cr[bb]))
            wts.append(aa.sum() / (np.nanstd(cc[bb])**2 + np.nanstd(cr[bb])**2 + 1e-12))
    if dcol:
        dcol = np.array(dcol); drow = np.array(drow)
        cen[int(q)] = {'n': int(len(dcol)), 'dcol_px': float(np.mean(dcol)), 'dcol_err': float(np.std(dcol) / np.sqrt(len(dcol))), 'drow_px': float(np.mean(drow)), 'drow_err': float(np.std(drow) / np.sqrt(len(drow))),
                       'col_rms_oot_px': float(np.nanstd(cc[b])), 'row_rms_oot_px': float(np.nanstd(cr[b]))}
res['sap_bkg_in_over_out_minus1_per_quarter'] = [float(x) for x in bk_in]
res['red_chi2_in_over_out_minus1_per_quarter'] = [float(x) for x in rc_in]
res['centroid_shift_in_transit'] = cen
print('bkg in/out-1 median', np.median(bk_in), ' redchi2 in/out-1 median', np.median(rc_in))
allc = np.array([[v['dcol_px'], v['dcol_err'], v['drow_px'], v['drow_err']] for v in cen.values()])
if len(allc):
    mc, ec = wmean(allc[:, 0], allc[:, 1]); mr, er = wmean(allc[:, 2], allc[:, 3])
    res['centroid_shift_in_transit_all'] = {'dcol_px': float(mc), 'dcol_err': float(ec), 'drow_px': float(mr), 'drow_err': float(er)}
    print('centroid shift in transit (all quarters, px): col', mc, '+-', ec, ' row', mr, '+-', er)

# ---------- 4. secondary scan over all phases, odd/even, trapezoid shape ----------
tt, rr, qq = cols['FLUX']
phases = np.linspace(0.05, 0.95, 181)
scan = np.array([box_depth(tt, rr, t0=T0 + p_ * P, nnull=1)[0] for p_ in phases])
_, err_box, _ = box_depth(tt, rr, nnull=600)
res['phase_scan'] = {'box_err_ppm': float(err_box), 'max_dip_ppm': float(np.nanmax(scan)), 'max_dip_phase': float(phases[int(np.nanargmax(scan))]),
                     'max_bump_ppm': float(np.nanmin(scan)), 'max_bump_phase': float(phases[int(np.nanargmin(scan))]),
                     'depth_at_phase_0.5_ppm': float(scan[np.argmin(np.abs(phases - 0.5))]),
                     'depth_at_phase_0.5_3xdur': float(box_depth(tt, rr, t0=T0 + 0.5 * P, dur=3 * DUR, nnull=1)[0])}
print('phase scan: max dip', np.nanmax(scan), 'at', phases[np.nanargmax(scan)], 'err', err_box, ' sec@0.5', res['phase_scan']['depth_at_phase_0.5_ppm'])
n_ep = np.round((tt - T0) / P).astype(int)
m_odd = n_ep % 2 == 1
do, eo, _ = box_depth(tt[m_odd], rr[m_odd], nnull=400); de, ee, _ = box_depth(tt[~m_odd], rr[~m_odd], nnull=400)
res['odd_even'] = {'odd': [do, eo], 'even': [de, ee], 'sigma': float(abs(do - de) / np.hypot(eo, ee))}
print('odd/even', do, eo, de, ee)
# trapezoid fit (T14, T23, depth) on binned phase curve, with 29.4-min integration via supersampling
ph = ((tt - T0 + P / 2) % P) - P / 2
m = np.abs(ph) < 0.4
def trap(params, x):
    dep, t14, t23, dt = params
    out = np.zeros_like(x)
    ss = np.linspace(-0.5, 0.5, 7) * (29.4 / 60 / 24)
    for s in ss:
        xx = np.abs(x + s - dt)
        y = np.zeros_like(x)
        full = xx <= t23 / 2; y[full] = dep
        ing = (xx > t23 / 2) & (xx < t14 / 2); y[ing] = dep * (t14 / 2 - xx[ing]) / ((t14 - t23) / 2 + 1e-9)
        out += y
    return -out / len(ss)
def resid(p):
    if p[2] > p[1] or p[2] < 0 or p[1] <= 0: return 1e3 * np.ones(m.sum())
    return (rr[m] - trap(p, ph[m])) * 1e3
best = None
for t23_0 in [0.02, 0.06, 0.09]:
    r0 = least_squares(resid, [1.1e-3, 0.105, t23_0, 0.0], method='lm' if False else 'trf')
    if best is None or r0.cost < best.cost: best = r0
dep, t14, t23, dt = best.x
chi2_trap = np.sum(best.fun**2)
# box of same t14
r_box = least_squares(lambda p: (rr[m] - trap([p[0], p[1], p[1] * 0.999, p[2]], ph[m])) * 1e3, [1.1e-3, 0.1, 0.0])
chi2_box = np.sum(r_box.fun**2)
# V shape (t23 = 0)
r_v = least_squares(lambda p: (rr[m] - trap([p[0], p[1], 0.0, p[2]], ph[m])) * 1e3, [1.5e-3, 0.12, 0.0])
chi2_v = np.sum(r_v.fun**2)
res['trapezoid'] = {'depth_ppm': float(dep * 1e6), 'T14_h': float(t14 * 24), 'T23_h': float(t23 * 24), 'T23_over_T14': float(t23 / t14), 'dt_min': float(dt * 24 * 60),
                    'chi2_trap': float(chi2_trap), 'chi2_box_T14': float(chi2_box), 'box_T14_h': float(r_box.x[1] * 24), 'chi2_V': float(chi2_v), 'V_T14_h': float(r_v.x[1] * 24), 'V_depth_ppm': float(r_v.x[0] * 1e6), 'n': int(m.sum())}
print('trapezoid', res['trapezoid'])

# ---------- 5. Lomb-Scargle on raw normalised flux, and depth after notching the 0.62-d sinusoid ----------
y = d['FLUX'].astype(float); ok = good & np.isfinite(y) & (y > 0)
yn = np.concatenate([y[ok & (q_all == q)] / np.nanmedian(y[ok & (q_all == q)]) for q in np.unique(q_all[ok])])
tn = np.concatenate([t_all[ok & (q_all == q)] for q in np.unique(q_all[ok])])
freq = np.linspace(1 / 60.0, 1 / 0.25, 300000)
pw = LombScargle(tn, yn).power(freq)
top = np.argsort(pw)[::-1][:8]
res['lomb_scargle_raw'] = {'top_periods_d': [float(1 / freq[i]) for i in top], 'top_powers': [float(pw[i]) for i in top],
                           'power_at_P': float(pw[np.argmin(np.abs(1 / freq - P))]), 'power_at_P/2': float(pw[np.argmin(np.abs(1 / freq - P / 2))])}
print('LS top periods', res['lomb_scargle_raw']['top_periods_d'][:5], res['lomb_scargle_raw']['top_powers'][:3], 'power at P', res['lomb_scargle_raw']['power_at_P'])
# LS on the detrended flux (1-d window kills > few d), restricted 0.3-3 d
tt, rr, qq = cols['FLUX']
freq2 = np.linspace(1 / 3.0, 1 / 0.3, 200000)
pw2 = LombScargle(tt, rr).power(freq2)
i2 = int(np.argmax(pw2)); p_ls = 1 / freq2[i2]
ls = LombScargle(tt, rr); model = ls.model(tt, freq2[i2])
amp = (model.max() - model.min()) / 2
res['ls_detrended'] = {'best_period_d': float(p_ls), 'power': float(pw2[i2]), 'half_amp_ppm': float(amp * 1e6), 'P_over_pls': float(P / p_ls),
                       'fap': float(ls.false_alarm_probability(pw2[i2]))}
print('LS detrended best', p_ls, 'power', pw2[i2], 'half amp ppm', amp * 1e6, 'P/pls', P / p_ls)
rr_notch = rr - (model - np.mean(model))
dep_n, err_n, _ = box_depth(tt, rr_notch, nnull=400)
res['ls_detrended']['depth_after_notch_ppm'] = [float(dep_n), float(err_n)]
# is the 0.62-d amplitude different in vs out of transit? (if the sinusoid is on the same star it should be, trivially tiny; just record)
print('depth after notching 0.62-d sinusoid', dep_n, '+-', err_n)
# inverted-flux check with my own box: best box in 2-f around P (should be ~0)
dep_i, err_i, _ = box_depth(tt, -rr, nnull=200)
res['inverted_box_at_P_ppm'] = [float(dep_i), float(err_i)]

json.dump(res, open(OUT, 'w'), indent=1, default=float)
print('written', OUT)
