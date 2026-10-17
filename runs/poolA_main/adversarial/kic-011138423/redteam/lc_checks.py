"""Red-team light-curve checks for kic-011138423 (own code, no pipeline imports)."""
import json, sys, numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle

FITS = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-011138423_kepler_v1.0_lc.fits'
OUT = '/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic011138423/out/lc_checks.json'
P = 4.894451164394656; T0 = 125.33662554571406; DUR = 1.95 / 24.0
rng = np.random.default_rng(42)
res = {}

f = fits.open(FITS)
d = f['LIGHTCURVE_STITCHED'].data
t_all = d['TIME']; q_all = d['QUARTER']; qual = d['SAP_QUALITY']

def detrend(t, y, mask_transit=True, win=1.0, gap=0.5):
    """per-segment running median (transit-masked) -> relative flux."""
    ph = ((t - T0 + P / 2) % P) - P / 2
    intr = np.abs(ph) < DUR * 0.75 if mask_transit else np.zeros_like(t, bool)
    order = np.argsort(t); t = t[order]; y = y[order]; intr = intr[order]
    trend = np.full_like(y, np.nan)
    # segments
    br = np.where(np.diff(t) > gap)[0] + 1
    segs = np.split(np.arange(len(t)), br)
    for s in segs:
        ts, ys, ms = t[s], y[s], intr[s]
        tt, yy = ts[~ms], ys[~ms]
        for j, i in enumerate(s):
            m = np.abs(tt - ts[j]) < win / 2
            if m.sum() > 8:
                trend[i] = np.median(yy[m])
    r = y / trend - 1
    return t, r

def box_depth(t, r, P=P, t0=T0, dur=DUR, nnull=500):
    ph = ((t - t0 + P / 2) % P) - P / 2
    a = np.abs(ph) < dur / 2
    b = (np.abs(ph) > dur) & (np.abs(ph) < 3 * dur)
    dep = -(np.mean(r[a]) - np.mean(r[b]))
    nulls = []
    for k in range(nnull):
        sh = rng.uniform(0.08 * P, 0.92 * P)
        ph2 = ((t - t0 - sh + P / 2) % P) - P / 2
        aa = np.abs(ph2) < dur / 2; bb = (np.abs(ph2) > dur) & (np.abs(ph2) < 3 * dur)
        if aa.sum() > 5 and bb.sum() > 5:
            nulls.append(-(np.mean(r[aa]) - np.mean(r[bb])))
    return dep * 1e6, np.std(nulls) * 1e6, int(a.sum())

# ---------- 1. per-quarter / per-season depth in each flux column ----------
season = {0: 3, 1: 3, 2: 0, 3: 1, 4: 2, 5: 3, 6: 0, 8: 2, 9: 3, 10: 0, 11: 1, 12: 2, 13: 3, 14: 0, 15: 1, 16: 2, 17: 3}
psffrac = {}
for i in range(2, 36, 2):
    h = f[i].header; psffrac[h['QUARTER']] = (h['PSFFRAC'], h['CROWDSAP'], h['FLFRCSAP'], h['NPIXSAP'], h['CHANNEL'])
res['psffrac_by_quarter'] = {int(k): v for k, v in psffrac.items()}

good = (qual == 0)
cols = {}
for col in ['FLUX', 'SAP_FLUX', 'PSF_FLAT_FLUX']:
    y = d[col].astype(float)
    ok = good & np.isfinite(y) & np.isfinite(t_all)
    # per-quarter normalise then detrend
    tt = []; rr = []; qq = []
    for q in np.unique(q_all[ok]):
        m = ok & (q_all == q)
        yq = y[m] / np.nanmedian(y[m])
        tq, rq = detrend(t_all[m], yq, mask_transit=(col != 'PSF_FLAT_FLUX'))
        tt.append(tq); rr.append(rq); qq.append(np.full(len(tq), q))
    tt = np.concatenate(tt); rr = np.concatenate(rr); qq = np.concatenate(qq)
    okr = np.isfinite(rr)
    # clip
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
        ds, es, ns = box_depth(tt[m], rr[m], nnull=400)
        per_s[s] = [round(ds, 1), round(es, 1), ns]
    res[col] = {'depth_all_ppm': round(dep, 1), 'err_ppm': round(err, 1), 'n_in': n,
                'per_quarter': per_q, 'per_season': per_s, 'median_flux_e_s': float(np.nanmedian(y[ok]))}
    print(col, 'all', round(dep), '+-', round(err), 'per season', per_s)

# per-quarter PSF_FLUX_NOVA (per-quarter extensions)
nova = {}
tt = []; rr = []; qq = []
for i in range(2, 36, 2):
    h = f[i].header; q = h['QUARTER']; dd = f[i].data
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

# depth vs PSFFRAC correlation (PSF FLUX, per quarter)
pq = res['FLUX']['per_quarter']
x = np.array([psffrac[q][0] for q in pq]); yv = np.array([pq[q][0] for q in pq]); ev = np.array([pq[q][1] for q in pq])
w = 1 / ev**2
A = np.vstack([x, np.ones_like(x)]).T
coef = np.linalg.lstsq(A * np.sqrt(w)[:, None], yv * np.sqrt(w), rcond=None)[0]
cov = np.linalg.inv((A * w[:, None]).T @ A)
res['depth_vs_psffrac'] = {'slope_ppm_per_unit': float(coef[0]), 'slope_err': float(np.sqrt(cov[0, 0])),
                           'intercept': float(coef[1]), 'chi2_const': float(np.sum(w * (yv - np.sum(w * yv) / np.sum(w))**2)), 'dof_const': int(len(yv) - 1)}
print('depth vs PSFFRAC slope', coef[0], '+-', np.sqrt(cov[0, 0]), 'chi2_const', res['depth_vs_psffrac']['chi2_const'], 'dof', len(yv) - 1)

# ---------- 2. per-event depths with local-noise errors and random-ephemeris calibration ----------
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
            # gap distance: nearest time with a data gap > 0.5 d
            out.append((tc, dep, err, a.sum(), b.sum()))
    return np.array(out)
ev = per_event(tt, rr, T0)
dep_w = np.sum(ev[:, 1] / ev[:, 2]**2) / np.sum(1 / ev[:, 2]**2)
chi2 = np.sum(((ev[:, 1] - dep_w) / ev[:, 2])**2)
res['per_event'] = {'n_events': int(len(ev)), 'weighted_depth_ppm': float(dep_w * 1e6), 'chi2_const': float(chi2), 'dof': int(len(ev) - 1),
                    'frac_positive': float(np.mean(ev[:, 1] > 0)), 'max_single_snr2_frac': float(np.max((ev[:, 1] / ev[:, 2])**2) / np.sum((ev[:, 1] / ev[:, 2])**2))}
# calibration: random t0 (same P) -> chi2/dof distribution of per-event "depths" about their weighted mean
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
# gap proximity: events within 1 d of a gap edge vs others
tsorted = np.sort(tt); gaps = np.where(np.diff(tsorted) > 0.5)[0]
edges = np.concatenate([[tsorted[0]], tsorted[gaps], tsorted[gaps + 1], [tsorted[-1]]])
gd = np.array([np.min(np.abs(edges - tc)) for tc in ev[:, 0]])
near = gd < 1.0
def wmean(x, e): return np.sum(x / e**2) / np.sum(1 / e**2), np.sqrt(1 / np.sum(1 / e**2))
dn, en = wmean(ev[near, 1], ev[near, 2]); df_, ef = wmean(ev[~near, 1], ev[~near, 2])
res['per_event']['near_gap_1d'] = {'n': int(near.sum()), 'depth_ppm': float(dn * 1e6), 'err': float(en * 1e6)}
res['per_event']['far_from_gap'] = {'n': int((~near).sum()), 'depth_ppm': float(df_ * 1e6), 'err': float(ef * 1e6)}
print('near-gap', near.sum(), dn * 1e6, en * 1e6, ' far', (~near).sum(), df_ * 1e6, ef * 1e6)
# depth excluding all events within 1 d of a gap
# first-half vs second-half of each quarter (thermal/pointing ramps)
half = []
for q in np.unique(qq):
    m = qq == q; tm = np.median(tt[m]); half.append((tt[m] < tm))
# per-year linear ephemeris (TTV / period drift) via per-year box t0 scan
tt_years = [(100, 470), (470, 840), (840, 1210), (1210, 1600)]
t0s = []
for (a, b) in tt_years:
    m = (tt > a) & (tt < b)
    shifts = np.linspace(-0.08, 0.08, 161)
    deps = [box_depth(tt[m], rr[m], t0=T0 + s, nnull=1)[0] for s in shifts]
    best = shifts[int(np.argmax(deps))]
    t0s.append([a, b, float(best * 24 * 60), float(max(deps))])
res['per_year_best_t0_shift_min'] = t0s
print('per-year t0 shift (min) and depth', t0s)

# ---------- 3. quality flags and background in transit ----------
ph = ((t_all - T0 + P / 2) % P) - P / 2
intr_all = np.abs(ph) < DUR / 2
okt = np.isfinite(t_all)
res['quality'] = {'frac_flagged_in_transit': float(np.mean(qual[okt & intr_all] != 0)), 'frac_flagged_overall': float(np.mean(qual[okt] != 0)),
                  'bits_present': sorted(set(int(b) for b in np.unique(qual[okt]) if b))}
# background and PSF red chi2 in vs out (per-quarter extensions), as relative difference
bk_in = []; bk_out = []; rc_in = []; rc_out = []
for i in range(2, 36, 2):
    dd = f[i].data; t = dd['TIME']; ql = dd['SAP_QUALITY']; ok = (ql == 0) & np.isfinite(t)
    ph = ((t - T0 + P / 2) % P) - P / 2
    a = ok & (np.abs(ph) < DUR / 2); b = ok & (np.abs(ph) > DUR) & (np.abs(ph) < 3 * DUR)
    bkg = dd['SAP_BKG'].astype(float); rc = dd['RED_CHI2'].astype(float)
    if np.isfinite(bkg[b]).sum() > 50:
        bk_in.append(np.nanmedian(bkg[a]) / np.nanmedian(bkg[b]) - 1)
    if np.isfinite(rc[b]).sum() > 50:
        rc_in.append(np.nanmedian(rc[a]) / np.nanmedian(rc[b]) - 1)
res['sap_bkg_in_over_out_minus1_per_quarter'] = [float(x) for x in bk_in]
res['red_chi2_in_over_out_minus1_per_quarter'] = [float(x) for x in rc_in]
print('bkg in/out-1 median', np.median(bk_in), ' redchi2 in/out-1 median', np.median(rc_in))

# ---------- 4. secondary scan over all phases, odd/even, half-period ----------
tt, rr, qq = cols['FLUX']
phases = np.linspace(0.05, 0.95, 181)
scan = []
for p_ in phases:
    dp, _, _ = box_depth(tt, rr, t0=T0 + p_ * P, nnull=1)
    scan.append(dp)
scan = np.array(scan)
_, err_box, _ = box_depth(tt, rr, nnull=600)
res['phase_scan'] = {'box_err_ppm': float(err_box), 'max_dip_ppm': float(scan.max()), 'max_dip_phase': float(phases[int(np.argmax(scan))]),
                     'max_bump_ppm': float(scan.min()), 'max_bump_phase': float(phases[int(np.argmin(scan))]),
                     'depth_at_phase_0.5_ppm': float(scan[np.argmin(np.abs(phases - 0.5))]),
                     'depth_at_phase_0.5_3xdur': float(box_depth(tt, rr, t0=T0 + 0.5 * P, dur=3 * DUR, nnull=1)[0])}
print('phase scan: max dip', scan.max(), 'at', phases[np.argmax(scan)], 'err', err_box, ' sec@0.5', res['phase_scan']['depth_at_phase_0.5_ppm'])
# odd/even
n_ep = np.round((tt - T0) / P).astype(int)
m_odd = n_ep % 2 == 1
do, eo, _ = box_depth(tt[m_odd], rr[m_odd], nnull=400); de, ee, _ = box_depth(tt[~m_odd], rr[~m_odd], nnull=400)
res['odd_even'] = {'odd': [do, eo], 'even': [de, ee], 'sigma': float(abs(do - de) / np.hypot(eo, ee))}
print('odd/even', do, eo, de, ee)

# ---------- 5. Lomb-Scargle on raw normalised flux (both polarities identical for LS) ----------
y = d['FLUX'].astype(float); ok = good & np.isfinite(y)
yn = np.concatenate([y[ok & (q_all == q)] / np.nanmedian(y[ok & (q_all == q)]) for q in np.unique(q_all[ok])])
tn = np.concatenate([t_all[ok & (q_all == q)] for q in np.unique(q_all[ok])])
freq = np.linspace(1 / 60.0, 1 / 0.25, 200000)
pw = LombScargle(tn, yn).power(freq)
top = np.argsort(pw)[::-1][:8]
res['lomb_scargle_raw'] = {'top_periods_d': [float(1 / freq[i]) for i in top], 'top_powers': [float(pw[i]) for i in top],
                           'power_at_P': float(pw[np.argmin(np.abs(1 / freq - P))]), 'power_at_P/2': float(pw[np.argmin(np.abs(1 / freq - P / 2))]),
                           'amp_at_P_ppm': float(np.sqrt(pw[np.argmin(np.abs(1 / freq - P))]) * np.std(yn) * 1e6 * np.sqrt(2))}
print('LS top periods', res['lomb_scargle_raw']['top_periods_d'][:5], 'power at P', res['lomb_scargle_raw']['power_at_P'])

json.dump(res, open(OUT, 'w'), indent=1, default=float)
print('written', OUT)
