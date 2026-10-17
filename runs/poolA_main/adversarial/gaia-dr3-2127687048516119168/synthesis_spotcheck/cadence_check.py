# Synthesis-stage spot check (2026-10-17), own code, no pipeline imports.
# Question: which cadences carry the 32.87-d box depth of gaia-dr3-2127687048516119168; Q13 share of pixel_v16.
import numpy as np, json
from astropy.io import fits
f='data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2127687048516119168_kepler_v1.0_lc.fits'
d=fits.open(f)['LIGHTCURVE_STITCHED'].data
t=d['TIME']; fl=d['FLUX']; sap=d['SAP_FLUX']; q=d['SAP_QUALITY']; qq=d['QUARTER']
ok=np.isfinite(t)&np.isfinite(fl)&(q==0)
fn=fl.copy(); sn=sap.copy()
for Q in np.unique(qq[ok]):
    m=ok&(qq==Q); fn[m]=fl[m]/np.nanmedian(fl[m]); sn[m]=sap[m]/np.nanmedian(sap[m])
P=32.868695202714555; t0=163.9810557924211; dur=1.5/24   # vet.json ephemeris
for lo,hi in [(1182.6,1183.3),(1215.6,1216.1)]:
    m=ok&(t>lo)&(t<hi)
    print(f'\n--- BKJD {lo}-{hi} (quarter {np.unique(qq[m])}) idx, BKJD, epoch, PSF ppm, SAP ppm, dt_prev')
    for i in np.where(m)[0]:
        prev=np.where(ok[:i])[0]; dtp=t[i]-t[prev[-1]] if len(prev) else np.nan
        print(f'{i:6d} {t[i]:10.4f} ep={(t[i]-t0)/P:8.4f} psf={1e6*(fn[i]-1):8.0f} sap={1e6*(sn[i]-1):8.0f} dt_prev={dtp:6.3f}')
gaps=np.where(np.diff(t[ok])>0.5)[0]; tt=t[ok]
print('\nsegment starts after gaps >0.5 d between BKJD 1170 and 1280:',[round(float(tt[g+1]),3) for g in gaps if 1170<tt[g+1]<1280])
ph=((t-t0)/P+0.5)%1-0.5; inbox=ok&(np.abs(ph*P)<dur/2); ep=np.round((t-t0)/P).astype(int)
rows=[]
for e in np.unique(ep[inbox]):
    mi=inbox&(ep==e); tc=t0+e*P; mo=ok&(np.abs(t-tc)<1.0)&(~inbox)
    if mo.sum()<10: continue
    rows.append((e,tc,mi.sum(),1e6*(np.median(fn[mo])-np.mean(fn[mi])),1e6*(np.median(sn[mo])-np.mean(sn[mi]))))
rows=np.array(rows); w=rows[:,2]
print('\nepoch tc n_in psf_depth sap_depth (|depth|>400 ppm or epochs 31,32)')
for r in rows:
    if abs(r[3])>400 or int(r[0]) in (31,32): print(f'{int(r[0]):4d} {r[1]:9.3f} {int(r[2]):3d} {r[3]:8.0f} {r[4]:8.0f}')
print('n epochs with data:',len(rows),' in-box cadences:',int(w.sum()))
print('cadence-weighted mean PSF depth: %.0f ppm; SAP %.0f'%(np.sum(rows[:,3]*w)/w.sum(),np.sum(rows[:,4]*w)/w.sum()))
keep=~np.isin(rows[:,0],[31,32]); n=keep.sum()
print('without epochs 31,32: PSF %.0f +- %.0f ppm (scatter/sqrt n), SAP %.0f +- %.0f'%(np.sum(rows[keep,3]*w[keep])/w[keep].sum(), np.std(rows[keep,3])/np.sqrt(n), np.sum(rows[keep,4]*w[keep])/w[keep].sum(), np.std(rows[keep,4])/np.sqrt(n)))
print('epochs 31+32 share of depth*n: %.2f'%(np.sum(rows[~keep,3]*w[~keep])/np.sum(rows[:,3]*w)))
print('1215.779 - 1182.911 = %.3f d (pipeline P = %.4f d)'%(1215.779-1182.911,P))
p=json.load(open('runs/poolA_main/pixel_v16/gaia-dr3-2127687048516119168/pixel.json'))
print('\npixel_v16 per quarter:')
for qd in p['quarters']:
    if 'error' in qd: print(' Q%d: %s'%(qd['quarter'],qd['error'])); continue
    keys=[k for k in qd if 'chi2' in k or 'events' in k or 'snr' in k.lower() or 'conclusive' in k or 'depth' in k.lower()]
    print(' Q%d:'%qd['quarter'], {k:(round(qd[k],2) if isinstance(qd[k],float) else qd[k]) for k in keys})
