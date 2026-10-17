"""Red-team light-curve tests for kic-005950408 (own code; no pipeline/ imports).
Outputs rt_lc_tests.json and rt_lc_tests.png in the same directory as the script's --out.
"""
import json, sys, os
import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle

FITS = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005950408_kepler_v1.0_lc.fits'
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
P = 38.60671217842687; T0 = 151.4674802687019; DUR_H = 0.75   # pipeline vet.json ephemeris
res = {'ephemeris': {'P_d': P, 't0_bkjd': T0, 'dur_h_pipeline': DUR_H}}

f = fits.open(FITS)
st = f['LIGHTCURVE_STITCHED'].data
t = st['TIME'].astype(float); q = st['QUARTER'].astype(int); cad = st['CADENCENO'].astype(int)
psf = st['FLUX'].astype(float); psfe = st['FLUX_ERR'].astype(float)
sap = st['SAP_FLUX'].astype(float); qual = st['SAP_QUALITY'].astype(int)
fmask = st['FLATTEN_MASK'].astype(int)
good = (qual == 0) & np.isfinite(t) & np.isfinite(psf) & np.isfinite(sap)
# per-quarter extras keyed by cadence number
extra = {}
chan = {}
for qq in range(18):
    try:
        h = f[f'LIGHTCURVE_Q{qq}']
    except KeyError:
        continue
    d = h.data; chan[qq] = dict(season=int(h.header['SEASON']), module=int(h.header['MODULE']), output=int(h.header['OUTPUT']),
                                channel=int(h.header['CHANNEL']), row=float(h.header['ROW']), col=float(h.header['COLUMN']),
                                flfrcsap=float(h.header['FLFRCSAP']), psffrac=float(h.header['PSFFRAC']), npixsap=int(h.header['NPIXSAP']))
    for c, rc, cc, rr, bk, nova in zip(d['CADENCENO'], d['RED_CHI2'], d['CENTROID_COLUMN'], d['CENTROID_ROW'], d['SAP_BKG'], d['PSF_FLUX_NOVA']):
        extra[int(c)] = (float(rc), float(cc), float(rr), float(bk), float(nova))
res['channels'] = chan

t = t[good]; q = q[good]; cad = cad[good]; psf = psf[good]; psfe = psfe[good]; sap = sap[good]; fmask = fmask[good]
o = np.argsort(t); t, q, cad, psf, psfe, sap, fmask = [a[o] for a in (t, q, cad, psf, psfe, sap, fmask)]
rc2 = np.array([extra.get(c, (np.nan,)*5)[0] for c in cad]); ccol = np.array([extra.get(c, (np.nan,)*5)[1] for c in cad])
crow = np.array([extra.get(c, (np.nan,)*5)[2] for c in cad]); bkg = np.array([extra.get(c, (np.nan,)*5)[3] for c in cad])
nova = np.array([extra.get(c, (np.nan,)*5)[4] for c in cad])

# per-quarter flux levels (evidence of contamination by season)
lv = {}
for qq in np.unique(q):
    m = q == qq
    lv[int(qq)] = dict(n=int(m.sum()), med_psf=float(np.median(psf[m])), med_sap=float(np.median(sap[m])), med_nova=float(np.nanmedian(nova[m])),
                       med_bkg=float(np.nanmedian(bkg[m])), med_redchi2=float(np.nanmedian(rc2[m])),
                       raw_mad_psf_ppm=float(1.4826*np.median(np.abs(psf[m]/np.median(psf[m])-1))*1e6),
                       raw_mad_sap_ppm=float(1.4826*np.median(np.abs(sap[m]/np.median(sap[m])-1))*1e6),
                       cen_col_std=float(np.nanstd(ccol[m])), cen_row_std=float(np.nanstd(crow[m])), **chan[int(qq)])
res['quarter_levels'] = lv

def local_rel(y, win=1.0, excl=0.12):
    """relative flux vs a local running median (window +/-win d, excluding +/-excl d around the point)"""
    out = np.full_like(y, np.nan)
    for i in range(len(y)):
        m = (np.abs(t - t[i]) < win) & (np.abs(t - t[i]) > excl) & (q == q[i])
        if m.sum() >= 10:
            out[i] = y[i]/np.median(y[m]) - 1
    return out

# normalise per quarter then local detrend
def qnorm(y):
    z = y.copy()
    for qq in np.unique(q):
        m = q == qq; z[m] = y[m]/np.median(y[m])
    return z
psfn = qnorm(psf); sapn = qnorm(sap); novan = qnorm(np.where(np.isfinite(nova), nova, np.nan))
psf_rel = local_rel(psfn); sap_rel = local_rel(sapn)

# gaps and quarter edges
dt = np.diff(t); gap_idx = np.where(dt > 0.3)[0]
gap_edges = np.sort(np.concatenate([[t[0]], t[gap_idx], t[gap_idx+1], [t[-1]]]))
def gapdist(x): return float(np.min(np.abs(gap_edges - x)))
qedge = {}
for qq in np.unique(q):
    m = q == qq; qedge[int(qq)] = (float(t[m].min()), float(t[m].max()))
def qedgedist(x, qq): a, b = qedge[int(qq)]; return float(min(x-a, b-x))

ph = ((t - T0)/P + 0.5) % 1 - 0.5
epoch = np.round((t - T0)/P).astype(int)
rows = []
for dur_h, tag in ((0.75, 'box075'), (0.9, 'box090')):
    inbox = np.abs(ph*P) < dur_h/48.0
    tab = []
    for i in np.where(inbox)[0]:
        tab.append(dict(epoch=int(epoch[i]), t=float(t[i]), quarter=int(q[i]), channel=chan[int(q[i])]['channel'], cadenceno=int(cad[i]),
                        psf_rel_ppm=float(psf_rel[i]*1e6), sap_rel_ppm=float(sap_rel[i]*1e6), psf_err_ppm=float(psfe[i]/psf[i]*1e6),
                        red_chi2=float(rc2[i]), red_chi2_over_qmed=float(rc2[i]/lv[int(q[i])]['med_redchi2']),
                        dcol_px=float(ccol[i]-np.nanmedian(ccol[q == q[i]])), drow_px=float(crow[i]-np.nanmedian(crow[q == q[i]])),
                        flatten_mask=int(fmask[i]), gap_dist_d=gapdist(t[i]), qedge_dist_d=qedgedist(t[i], q[i])))
    res[tag] = dict(n_in=len(tab), table=tab)
    v = np.array([r['psf_rel_ppm'] for r in tab]); s = np.array([r['sap_rel_ppm'] for r in tab])
    out = ~inbox & np.isfinite(psf_rel)
    res[tag]['depth_mean_psf_ppm'] = float(-np.mean(v)); res[tag]['depth_median_psf_ppm'] = float(-np.median(v))
    res[tag]['depth_mean_sap_ppm'] = float(-np.nanmean(s)); res[tag]['depth_median_sap_ppm'] = float(-np.nanmedian(s))
    keep = v > -10000
    res[tag]['depth_mean_psf_without_outliers_below_-1pct_ppm'] = float(-np.mean(v[keep])); res[tag]['n_outliers_below_-1pct'] = int((~keep).sum())
    res[tag]['in_box_std_ppm'] = float(np.std(v)); res[tag]['formal_err_from_inbox_std_ppm'] = float(np.std(v)/np.sqrt(len(v)))
    res[tag]['out_of_box_mad_ppm'] = float(1.4826*np.median(np.abs(psf_rel[out] - np.median(psf_rel[out])))*1e6)
    # per-epoch
    ep = {}
    for r in tab: ep.setdefault(r['epoch'], []).append(r['psf_rel_ppm'])
    res[tag]['per_epoch_mean_ppm'] = {int(k): [float(np.mean(x)), len(x)] for k, x in ep.items()}
    # per-season
    seas = {}
    for r in tab: seas.setdefault(chan[r['quarter']]['season'], []).append(r['psf_rel_ppm'])
    res[tag]['per_season_depth_ppm'] = {int(k): [float(-np.mean(x)), len(x)] for k, x in seas.items()}

# neighbourhood of the outlier cadences (single-cadence dropout test)
outl = []
for r in res['box075']['table'] + res['box090']['table']:
    if r['psf_rel_ppm'] < -10000 and r['cadenceno'] not in [x['cadenceno'] for x in outl]:
        i = int(np.where(cad == r['cadenceno'])[0][0])
        nb = []
        for j in range(max(0, i-4), min(len(t), i+5)):
            nb.append(dict(dt_min=float((t[j]-t[i])*1440), psf_rel_ppm=float(psf_rel[j]*1e6), sap_rel_ppm=float(sap_rel[j]*1e6),
                           nova_rel_ppm=float((novan[j]/np.nanmedian(novan[np.abs(t-t[j])<1])-1)*1e6) if np.isfinite(novan[j]) else None,
                           red_chi2=float(rc2[j]), bkg=float(bkg[j]), dcol_px=float(ccol[j]-np.nanmedian(ccol[q == q[j]])), drow_px=float(crow[j]-np.nanmedian(crow[q == q[j]]))))
        outl.append(dict(cadenceno=r['cadenceno'], t=r['t'], epoch=r['epoch'], quarter=r['quarter'], neighbours=nb))
res['outlier_neighbourhoods'] = outl

# census of deep single-cadence dropouts per channel
cen = {}
for qq in np.unique(q):
    m = (q == qq) & np.isfinite(psf_rel)
    cen[int(qq)] = {'channel': chan[int(qq)]['channel'], 'n': int(m.sum()), 'n_psf_below_-1pct': int((psf_rel[m] < -0.01).sum()), 'n_psf_below_-3pct': int((psf_rel[m] < -0.03).sum()),
                    'n_sap_below_-1pct': int((sap_rel[m] < -0.01).sum()), 'n_psf_above_+1pct': int((psf_rel[m] > 0.01).sum())}
res['dropout_census'] = cen
bych = {}
for qq, c in cen.items():
    b = bych.setdefault(c['channel'], dict(n=0, psf_below_1pct=0, psf_below_3pct=0, sap_below_1pct=0))
    b['n'] += c['n']; b['psf_below_1pct'] += c['n_psf_below_-1pct']; b['psf_below_3pct'] += c['n_psf_below_-3pct']; b['sap_below_1pct'] += c['n_sap_below_-1pct']
res['dropout_census_by_channel'] = bych

# random-ephemeris null for the box-depth statistic (mean in-box relative flux), same P range as prereg (0.5-40 d), 0.75 h box
rng = np.random.default_rng(5950408)
fin = np.isfinite(psf_rel); tt = t[fin]; yy = psf_rel[fin]
N = 30000; depths = np.full(N, np.nan); nins = np.zeros(N, int)
Ps = np.exp(rng.uniform(np.log(20.0), np.log(40.0), N)); t0s = rng.uniform(0, 1, N)*Ps + tt[0]
for k in range(N):
    phk = ((tt - t0s[k])/Ps[k] + 0.5) % 1 - 0.5
    m = np.abs(phk*Ps[k]) < DUR_H/48.0
    nins[k] = m.sum()
    if m.sum() >= 20: depths[k] = -np.mean(yy[m])*1e6
ok = np.isfinite(depths)
res['random_ephemeris_null'] = dict(N=N, n_valid=int(ok.sum()), P_range_d=[20, 40], dur_h=DUR_H,
                                    frac_depth_ge_4104ppm=float(np.mean(depths[ok] >= 4104)), n_depth_ge_4104ppm=int((depths[ok] >= 4104).sum()),
                                    frac_depth_ge_3000ppm=float(np.mean(depths[ok] >= 3000)), depth_p50=float(np.median(depths[ok])), depth_p99=float(np.percentile(depths[ok], 99)),
                                    depth_p999=float(np.percentile(depths[ok], 99.9)), depth_max=float(np.nanmax(depths)), depth_robust_std_ppm=float(1.4826*np.median(np.abs(depths[ok]-np.median(depths[ok])))),
                                    depth_std_ppm=float(np.std(depths[ok])), median_n_in=float(np.median(nins[ok])))
# same null on dataset without channel 25
fin2 = fin & np.array([chan[int(x)]['channel'] != 25 for x in q]); tt2 = t[fin2]; yy2 = psf_rel[fin2]
d2 = np.full(N, np.nan)
for k in range(N):
    phk = ((tt2 - t0s[k])/Ps[k] + 0.5) % 1 - 0.5
    m = np.abs(phk*Ps[k]) < DUR_H/48.0
    if m.sum() >= 20: d2[k] = -np.mean(yy2[m])*1e6
ok2 = np.isfinite(d2)
res['random_ephemeris_null_no_ch25'] = dict(n_valid=int(ok2.sum()), frac_depth_ge_4104ppm=float(np.mean(d2[ok2] >= 4104)), depth_p99=float(np.percentile(d2[ok2], 99)),
                                            depth_p999=float(np.percentile(d2[ok2], 99.9)), depth_max=float(np.nanmax(d2)), depth_robust_std_ppm=float(1.4826*np.median(np.abs(d2[ok2]-np.median(d2[ok2])))))

# Lomb-Scargle, both polarities, on per-quarter-normalised PSF flux (not locally detrended) and on local-detrended
for nm, y in (('qnorm', psfn - 1), ('local', psf_rel)):
    m = np.isfinite(y) & (np.abs(y) < 0.05)
    freq = np.linspace(1/100., 1/0.1, 200000)
    ls = LombScargle(t[m], y[m]); pw = ls.power(freq)
    i = np.argmax(pw); per = 1/freq[i]
    amp = ls.model_parameters(freq[i])
    res[f'ls_{nm}'] = dict(best_period_d=float(per), power=float(pw[i]), half_amp_ppm=float(np.hypot(amp[1], amp[2])*1e6), fap=float(ls.false_alarm_probability(pw[i])),
                           ratio_P_over_Pls=float(P/per))
    # top 5 peaks
    idx = np.argsort(pw)[::-1][:2000]; peaks = []
    for j in idx:
        if all(abs(1/freq[j] - pp) > 0.02*pp for pp, _ in peaks): peaks.append((float(1/freq[j]), float(pw[j])))
        if len(peaks) >= 5: break
    res[f'ls_{nm}']['top5'] = peaks

# duration / density
Rs = 0.9126939129829408; Ms = 1.0; G = 6.674e-11; Rsun = 6.957e8; Msun = 1.989e30
a = (G*Ms*Msun*(P*86400)**2/(4*np.pi**2))**(1/3)
k = np.sqrt(4104e-6)
T14_b0 = P*24/np.pi*np.arcsin(Rs*Rsun*(1+k)/a)
for dur in (0.75, 0.9):
    b_req = np.sqrt((1+k)**2 - (a/(Rs*Rsun)*np.sin(dur/24*np.pi/P))**2)
    res[f'duration_geometry_{dur}h'] = dict(a_over_Rs=float(a/(Rs*Rsun)), T14_b0_h=float(T14_b0), b_required=float(b_req), rp_earth_box=float(k*Rs*109.1))
res['duration_geometry_note'] = 'Ms=1.0 Msun assumed; Rs from target table; b_required is the impact parameter needed to make T14 equal the quoted duration for a (1+k) chord'

json.dump(res, open(os.path.join(OUT, 'rt_lc_tests.json'), 'w'), indent=1)

# figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(4, 1, figsize=(14, 14))
for qq in np.unique(q):
    m = q == qq; c = 'C3' if chan[int(qq)]['channel'] == 25 else 'C0'
    ax[0].plot(t[m], psfn[m], '.', ms=1, color=c)
ax[0].set_ylim(0.5, 1.3); ax[0].set_ylabel('PSF flux / quarter median'); ax[0].set_title('kic-005950408 PSF flux; red = channel 25 (season 3); dashed = pipeline ephemeris epochs')
for n in range(-1, 40):
    tc = T0 + n*P
    if t[0] < tc < t[-1]: ax[0].axvline(tc, color='k', lw=0.5, ls='--')
ax[1].plot(t, psf_rel*1e6, '.', ms=1, color='0.5'); ax[1].set_ylim(-60000, 20000); ax[1].set_ylabel('local-detrended PSF [ppm]')
ib = np.abs(ph*P) < 0.75/48
ax[1].plot(t[ib], psf_rel[ib]*1e6, 'o', ms=4, color='C1', label='in 0.75-h box at pipeline ephemeris'); ax[1].legend()
for n in range(-1, 40):
    tc = T0 + n*P
    if t[0] < tc < t[-1]: ax[1].axvline(tc, color='k', lw=0.5, ls='--')
# zoom on the three carrying epochs
for k_, ep_ in enumerate((0, 10, 29)):
    tc = T0 + ep_*P; m = np.abs(t - tc) < 1.5
    a2 = ax[2] if k_ == 0 else ax[2].inset_axes([0.35 + 0.33*(k_-1), 0.1, 0.3, 0.8])
    if k_ == 0: a2 = ax[2]
    a2.plot((t[m]-tc)*24, psf_rel[m]*1e6, 'o-', ms=3, color='C0', label='PSF'); a2.plot((t[m]-tc)*24, sap_rel[m]*1e6, 's-', ms=3, color='C2', label='SAP')
    a2.axvspan(-0.375, 0.375, color='C1', alpha=0.3); a2.set_title(f'epoch {ep_}  t={tc:.2f} Q{int(q[m][0]) if m.sum() else -1}')
    a2.set_xlim(-30, 30)
ax[2].set_xlabel('hours from pipeline mid-transit'); ax[2].set_ylabel('ppm'); ax[2].legend(loc='lower left')
ax[3].hist(np.clip(depths[ok], -20000, 20000), bins=200, histtype='step', color='C0', label='all data (box075, P 20-40 d, random phase)')
ax[3].hist(np.clip(d2[ok2], -20000, 20000), bins=200, histtype='step', color='C2', label='without channel-25 quarters')
ax[3].axvline(4104, color='C3', ls='--', label='pipeline depth 4104 ppm'); ax[3].set_yscale('log'); ax[3].set_xlabel('box depth of a random ephemeris [ppm]'); ax[3].legend()
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'rt_lc_tests.png'), dpi=110)
print('done')
