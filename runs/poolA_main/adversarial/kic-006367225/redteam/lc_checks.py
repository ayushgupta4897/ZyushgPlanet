import numpy as np, json, sys
from astropy.io import fits
from astropy.timeseries import LombScargle
fn='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-006367225_kepler_v1.0_lc.fits'
P=12.971527521795288; T0=142.75300012519583; DUR=2.25/24
f=fits.open(fn)
def detrend(t,y,win=1.0,mask=None,gap=0.5):
    y=y.copy(); out=np.full_like(y,np.nan)
    br=np.where(np.diff(t)>gap)[0]+1; segs=np.split(np.arange(len(t)),br)
    for s in segs:
        ts,ys=t[s],y[s]; ms=np.ones(len(s),bool) if mask is None else ~mask[s]
        for i in range(len(s)):
            w=(np.abs(ts-ts[i])<win/2)&ms
            out[s[i]]=np.nanmedian(ys[w]) if w.sum()>=5 else np.nan
    return y/out
def phase(t): return ((t-T0)/P+0.5)%1-0.5
res={'quarters':{}}
allq=[]
for h in f:
    if not h.name.startswith('LIGHTCURVE_Q'): continue
    q=int(h.name.split('Q')[1]); hd=h.header; d=h.data
    t=d['TIME']; good=(d['SAP_QUALITY']==0)&np.isfinite(t)&np.isfinite(d['FLUX'])&np.isfinite(d['SAP_FLUX'])
    t=t[good]; rec={'module':hd.get('MODULE'),'output':hd.get('OUTPUT'),'channel':hd.get('CHANNEL'),'season':hd.get('SEASON'),'n':int(good.sum())}
    intr=np.abs(phase(t))*P<0.75*DUR; mask=np.abs(phase(t))*P<1.5*DUR
    for col in ('FLUX','SAP_FLUX'):
        y=d[col][good].astype(float); raw=np.nanmedian(y)
        yn=detrend(t,y/raw,mask=mask); ok=np.isfinite(yn)
        oot=(np.abs(phase(t))*P>1.5*DUR)&(np.abs(phase(t))*P<6/24)&ok
        inn=intr&ok
        if inn.sum()<3: rec[col]=None; continue
        sig=1.4826*np.nanmedian(np.abs(yn[ok]-np.nanmedian(yn[ok])))
        dep=np.nanmedian(yn[oot])-np.nanmean(yn[inn]); err=sig*np.sqrt(1/inn.sum()+1/max(oot.sum(),1))
        rec[col]={'raw_median_e_s':float(raw),'depth_ppm':float(dep*1e6),'err_ppm':float(err*1e6),'n_in':int(inn.sum()),'scatter_ppm':float(sig*1e6),'abs_deficit_e_s':float(dep*raw),'abs_err_e_s':float(err*raw)}
    res['quarters'][q]=rec
    allq.append((q,t,d['FLUX'][good].astype(float),d['SAP_FLUX'][good].astype(float),rec['module']))
# group by module
print('Q  mod ch  | PSF depth ppm (err)  raw e/s  abs e/s | SAP depth ppm (err)  raw e/s  abs e/s')
for q,r in sorted(res['quarters'].items()):
    a,b=r['FLUX'],r['SAP_FLUX']
    if a and b: print(f"{q:2d} {r['module']:3} {str(r['channel']):>3} | {a['depth_ppm']:7.0f} ({a['err_ppm']:4.0f}) {a['raw_median_e_s']:7.0f} {a['abs_deficit_e_s']:5.2f}+-{a['abs_err_e_s']:4.2f} | {b['depth_ppm']:7.0f} ({b['err_ppm']:4.0f}) {b['raw_median_e_s']:7.0f} {b['abs_deficit_e_s']:5.2f}+-{b['abs_err_e_s']:4.2f}")
mods={}
for q,r in res['quarters'].items():
    if r['FLUX'] and r['SAP_FLUX']: mods.setdefault(r['module'],[]).append(r)
print('\nper-module inverse-variance means:')
summ={}
for m,L in sorted(mods.items()):
    o={}
    for col in ('FLUX','SAP_FLUX'):
        for key in ('depth_ppm','abs_deficit_e_s'):
            ek='err_ppm' if key=='depth_ppm' else 'abs_err_e_s'
            v=np.array([r[col][key] for r in L]); e=np.array([r[col][ek] for r in L]); w=1/e**2
            o[f'{col}_{key}']=(float((v*w).sum()/w.sum()),float(w.sum()**-0.5))
    summ[m]=o
    print(m,[q for q,r in sorted(res['quarters'].items()) if r['module']==m], {k:(round(v[0],2),round(v[1],2)) for k,v in o.items()})
# chi2 of constancy for abs deficit
for col in ('FLUX','SAP_FLUX'):
    for key in ('depth_ppm','abs_deficit_e_s'):
        v=np.array([summ[m][f'{col}_{key}'][0] for m in summ]); e=np.array([summ[m][f'{col}_{key}'][1] for m in summ]); w=1/e**2; mu=(v*w).sum()/w.sum()
        print(f'chi2 const {col} {key}: {((v-mu)**2*w).sum():.1f} / {len(v)-1} dof ; mean {mu:.3g}')
res['per_module']=summ
# stitched: fold at KOI-1798 ephemeris, LS, both
ds=f['LIGHTCURVE_STITCHED'].data; g=(ds['SAP_QUALITY']==0)&np.isfinite(ds['TIME'])&np.isfinite(ds['FLUX'])
t=ds['TIME'][g]; y=ds['FLUX'][g].astype(float); 
# per-quarter normalise
qq=ds['QUARTER'][g] if 'QUARTER' in ds.names else np.zeros(len(t))
for u in np.unique(qq): y[qq==u]/=np.nanmedian(y[qq==u])
mask=np.abs(phase(t))*P<1.5*DUR
yn=detrend(t,y,mask=mask); ok=np.isfinite(yn)
def boxdepth(Pk,t0k,durk):
    ph=((t-t0k)/Pk+0.5)%1-0.5; inn=(np.abs(ph)*Pk<0.5*durk)&ok; oot=(np.abs(ph)*Pk>durk)&(np.abs(ph)*Pk<4*durk)&ok
    sig=1.4826*np.nanmedian(np.abs(yn[ok]-1)); dep=np.nanmedian(yn[oot])-np.nanmean(yn[inn]); return dep*1e6, sig*np.sqrt(1/inn.sum()+1/oot.sum())*1e6, int(inn.sum())
print('\nown box depth at pipeline ephemeris (stitched PSF):',boxdepth(P,T0,DUR))
print('box depth at KOI-1798.01 ephemeris (P=12.964718,t0=136.143342,dur 1.83h):',boxdepth(12.964718,136.143342,1.8303/24))
print('box depth at KOI-1798.01 period, best phase scan:')
best=None
for ph0 in np.linspace(0,12.964718,400,endpoint=False):
    r=boxdepth(12.964718,136.143342+ph0,1.8303/24)
    if best is None or r[0]>best[0][0]: best=(r,ph0)
print('  deepest:',best)
res['koi1798_fold']={'at_koi_ephemeris':boxdepth(12.964718,136.143342,1.8303/24),'best_phase':best}
# LS
oo=ok&~mask
for lab,yy in (('PSF',yn),):
    freq=np.linspace(1/50,1/0.1,200000); pw=LombScargle(t[oo],yy[oo]).power(freq)
    idx=np.argsort(pw)[::-1][:8]; print('LS top peaks (period d, power, amp ppm):')
    peaks=[]
    for i in idx:
        A=LombScargle(t[oo],yy[oo]).model_parameters(freq[i]); amp=np.hypot(A[1],A[2])*1e6 if len(A)>2 else np.nan
        peaks.append((float(1/freq[i]),float(pw[i]),float(amp))); print(f'  {1/freq[i]:.4f} {pw[i]:.5f} {amp:.0f}')
    for Pt in (P,P/2,2*P,P/3):
        j=np.argmin(np.abs(1/freq-Pt)); print(f'  LS power at {Pt:.3f} d: {pw[j]:.5f} (max {pw.max():.5f})')
    res['LS_peaks']=peaks
# events vs gaps
import csv
ev=list(csv.DictReader(open('/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-006367225/events.csv')))
dep=np.array([float(e['depth']) for e in ev])*1e6; err=np.array([float(e['err']) for e in ev])*1e6; gd=np.array([float(e['gap_dist_d']) for e in ev]); tc=np.array([float(e['tc_bkjd']) for e in ev])
# assign module per event
qedges={q:(tt.min(),tt.max()) for q,tt,*_ in allq}; modof=np.array([next((r['module'] for qq2,(a,b) in qedges.items() for r in [res['quarters'][qq2]] if a-0.5<=x<=b+0.5),-1) for x in tc])
print('\nevents: n',len(ev),' near gap (<1d):',(gd<1).sum(),' mean depth near-gap vs far:',dep[gd<1].mean() if (gd<1).any() else None, dep[gd>=1].mean())
for m in sorted(set(modof)):
    s=modof==m; w=1/err[s]**2; print(f'  module {m}: n={s.sum()} weighted mean event depth {(dep[s]*w).sum()/w.sum():.0f} ppm, n>5000ppm: {(dep[s]>5000).sum()}')
print('  corr(depth, gap_dist):',np.corrcoef(dep,gd)[0,1])
res['events_by_module']={int(m):{'n':int((modof==m).sum()),'wmean_depth_ppm':float((dep[modof==m]/err[modof==m]**2).sum()/(1/err[modof==m]**2).sum())} for m in set(modof)}
json.dump(res,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic006367225/lc_checks.json','w'),indent=1,default=float)
