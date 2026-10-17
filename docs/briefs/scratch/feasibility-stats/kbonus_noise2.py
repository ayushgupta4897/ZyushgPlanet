import numpy as np, lightkurve as lk, warnings, glob, pandas as pd
from wotan import flatten
warnings.filterwarnings('ignore')
fn=glob.glob('lkcache/**/*.fits', recursive=True); print(fn)
lc=lk.read(fn[0])
print("N", len(lc), "time range", lc.time.value.min(), lc.time.value.max(), "quarters", np.unique(lc['quarter'].value))
tt=np.asarray(lc.time.value,float)
def cdpp_like(tt, fn_, hours, window=1.0):
    ok=np.isfinite(fn_)&np.isfinite(tt); tt=tt[ok]; fn_=fn_[ok]
    flat=flatten(tt, fn_, method='biweight', window_length=window)
    g=np.isfinite(flat); flat=flat[g]
    n=max(1,int(round(hours/ (29.4/60))))
    cs=np.convolve(flat-1, np.ones(n)/n, mode='valid')
    return 1.4826*np.nanmedian(np.abs(cs-np.nanmedian(cs)))*1e6, np.nanstd(cs)*1e6, flat.size
rows=[]
for col in ['flux','sap_flux','psf_flat_flux']:
    fl=np.asarray(lc[col].value,float); 
    ok=np.isfinite(fl)&(fl>0)
    # normalise per quarter
    q=np.asarray(lc['quarter'].value)
    fn_=np.full(fl.size,np.nan)
    for qq in np.unique(q[ok]):
        m=ok&(q==qq); fn_[m]=fl[m]/np.nanmedian(fl[m])
    err=np.asarray(lc[col+'_err'].value,float)/np.nanmedian(fl[ok])
    p2p=1.4826*np.nanmedian(np.abs(np.diff(fn_[ok])-np.nanmedian(np.diff(fn_[ok]))))/np.sqrt(2)*1e6
    r=dict(col=col, median_flux_err_ppm=np.nanmedian(err[ok])*1e6, p2p_ppm=p2p)
    for h in [2.3,3,6,12]:
        mad,std,n=cdpp_like(tt[ok],fn_[ok],h)
        r[f'cdpp{h}h_mad_ppm']=mad; r[f'cdpp{h}h_std_ppm']=std
    rows.append(r)
df=pd.DataFrame(rows); pd.set_option('display.width',250); print(df.round(1).to_string(index=False))
df.to_csv('kbonus_target_noise.csv',index=False)
# per-quarter 2.3-h and 6-h CDPP for PSF flux
fl=np.asarray(lc['flux'].value,float); q=np.asarray(lc['quarter'].value); ok=np.isfinite(fl)&(fl>0)
print("quarter  npts  cdpp2.3h  cdpp6h  (PSF 'flux')")
for qq in np.unique(q[ok]):
    m=ok&(q==qq)
    if m.sum()<500: continue
    f_=fl[m]/np.nanmedian(fl[m])
    a,_,_=cdpp_like(tt[m],f_,2.3); b,_,_=cdpp_like(tt[m],f_,6)
    print(qq, m.sum(), round(a), round(b))
# expected SNR of the candidate: depth / cdpp(2.3h) * sqrt(Ntransits) ; simple fold check of depth
P=13.769448; T0=2455630.245431-2454833.0  # BKJD
f_=np.full(fl.size,np.nan)
for qq in np.unique(q[ok]):
    m=ok&(q==qq); f_[m]=fl[m]/np.nanmedian(fl[m])
flat=flatten(tt[ok],f_[ok],method='biweight',window_length=1.0)
ph=((tt[ok]-T0+0.5*P)%P)-0.5*P
intr=np.abs(ph)<(2.3/24/2); out=(np.abs(ph)>0.2)&(np.abs(ph)<1.0)
d=(np.nanmedian(flat[out])-np.nanmedian(flat[intr]))*1e6
ntr=len(np.unique(np.round((tt[ok][intr]-T0)/P)))
sig_in=1.4826*np.nanmedian(np.abs(flat[out]-np.nanmedian(flat[out])))/np.sqrt(intr.sum())*1e6
print(f"fold check (PSF flux, biweight 1d): in-transit pts {intr.sum()}, transits with data {ntr}, box depth {d:.0f} ppm, +/- {sig_in:.0f} ppm (white-noise estimate) -> SNR {d/sig_in:.1f}")
cd=df.loc[df.col=='flux','cdpp2.3h_mad_ppm'].item()
print(f"expected SNR from CDPP: depth 2165 ppm / {cd:.0f} ppm * sqrt({ntr}) = {2165/cd*np.sqrt(ntr):.1f}")
