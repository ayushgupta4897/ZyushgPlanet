# Synthesis-stage spot check (2026-10-17), own code, no pipeline imports.
# Question: which cadences carry the 28.07-d box depth of kic-008620445, and are they segment-edge artefacts?
import numpy as np, pandas as pd
from astropy.io import fits
f='data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-008620445_kepler_v1.0_lc.fits'
d=fits.open(f)['LIGHTCURVE_STITCHED'].data
t=d['TIME']; fl=d['FLUX']; sap=d['SAP_FLUX']; q=d['SAP_QUALITY']; qq=d['QUARTER']
ok=np.isfinite(t)&np.isfinite(fl)&(q==0)
fn=fl.copy(); sn=sap.copy()
for Q in np.unique(qq[ok]):
    m=ok&(qq==Q); fn[m]=fl[m]/np.nanmedian(fl[m]); sn[m]=sap[m]/np.nanmedian(sap[m])
P=28.07173237974542; t0=144.14049296937884; dur=0.75/24   # vet.json ephemeris
for lo,hi in [(200.1,200.6),(1126.4,1126.9)]:
    m=ok&(t>lo)&(t<hi)
    print(f'\n--- BKJD {lo}-{hi} (quarter {np.unique(qq[m])}) idx, BKJD, epoch at P_pipe, PSF ppm, SAP ppm, dt from previous good cadence (d)')
    for i in np.where(m)[0]:
        prev=np.where(ok[:i])[0]; dtp=t[i]-t[prev[-1]] if len(prev) else np.nan
        print(f'{i:6d} {t[i]:10.4f} ep={(t[i]-t0)/P:8.4f} psf={1e6*(fn[i]-1):9.0f} sap={1e6*(sn[i]-1):9.0f} dt_prev={dtp:6.3f}')
ph=((t-t0)/P+0.5)%1-0.5; inbox=ok&(np.abs(ph*P)<dur/2); ep=np.round((t-t0)/P).astype(int)
print('\nin-box cadences:',inbox.sum())
rows=[]
for e in np.unique(ep[inbox]):
    mi=inbox&(ep==e); tc=t0+e*P; mo=ok&(np.abs(t-tc)<1.0)&(~inbox)
    if mo.sum()<10: continue
    rows.append((e,tc,mi.sum(),1e6*(np.median(fn[mo])-np.mean(fn[mi])),1e6*(np.median(sn[mo])-np.mean(sn[mi]))))
rows=np.array(rows); w=rows[:,2]
print('epoch  tc  n_in  psf_depth_ppm  sap_depth_ppm   (local depth = median of +-1 d out-of-box minus in-box mean)')
for r in rows: print(f'{int(r[0]):4d} {r[1]:9.3f} {int(r[2]):3d} {r[3]:9.0f} {r[4]:9.0f}')
print('\ncadence-weighted mean PSF depth all epochs: %.0f ppm'%(np.sum(rows[:,3]*w)/w.sum()))
keep=~np.isin(rows[:,0],[2,35])
print('without epochs 2 and 35: PSF %.0f ppm, SAP %.0f ppm'%(np.sum(rows[keep,3]*w[keep])/w[keep].sum(), np.sum(rows[keep,4]*w[keep])/w[keep].sum()))
print('epochs 2+35 share of total depth*n: %.2f'%(np.sum(rows[~keep,3]*w[~keep])/np.sum(rows[:,3]*w)))
print('alias comb: (1126.651-200.284) / 28.07173 = %.4f ; / 40.27715 = %.4f'%((1126.651-200.284)/28.07173,(1126.651-200.284)/40.27715))
e=pd.read_csv('runs/poolA_main/vet/kic-008620445/events.csv'); wt=1/e.err**2; mu=np.sum(wt*e.depth)/wt.sum()
print('pipeline events.csv (17 events with n_in>=2, formal errors): weighted mean %.0f +- %.0f ppm, chi2 %.1f / %d'%(1e6*mu,1e6/np.sqrt(wt.sum()),np.sum(wt*(e.depth-mu)**2),len(e)-1))
