import numpy as np, pandas as pd, glob, warnings, lightkurve as lk
from wotan import flatten
warnings.filterwarnings('ignore')
sel=pd.read_csv('kbonus_noise_sample.csv'); 
elig=pd.read_csv('kbonus_eligible_background_mdwarfs.csv')
def cdpp_like(tt, fn_, hours, window=1.0):
    ok=np.isfinite(fn_)&np.isfinite(tt); tt=tt[ok]; fn_=fn_[ok]
    flat=flatten(tt, fn_, method='biweight', window_length=window); flat=flat[np.isfinite(flat)]
    n=max(1,int(round(hours/(29.4/60)))); cs=np.convolve(flat-1, np.ones(n)/n, mode='valid')
    return 1.4826*np.nanmedian(np.abs(cs-np.nanmedian(cs)))*1e6
rows=[]
for fn in sorted(glob.glob('lkcache/mastDownload/HLSP/hlsp_kbonus-bkg_*/*_lc.fits')):
    lc=lk.read(fn); kic=int(lc.meta.get('KEPLERID') or 0); G=lc.meta.get('GMAG')
    row=elig[elig.kic==kic] if kic else elig[elig.gaia_designation==lc.meta.get('GAIAID','')]
    tt=np.asarray(lc.time.value,float); q=np.asarray(lc['quarter'].value)
    out=dict(kic=kic, G=round(float(G),2), kepmag=round(float(lc.meta.get('KEPMAG') or np.nan),2), psffrac=round(float(row.psffrac.iloc[0]),2) if len(row) else np.nan, crowdsap=round(float(row.crowdsap.iloc[0]),2) if len(row) else np.nan, nq=len(np.unique(q)), npts=len(lc))
    for col in ['flux','sap_flux']:
        fl=np.asarray(lc[col].value,float); ok=np.isfinite(fl)&(fl>0)
        fn_=np.full(fl.size,np.nan)
        for qq in np.unique(q[ok]):
            m=ok&(q==qq); fn_[m]=fl[m]/np.nanmedian(fl[m])
        out[f'{col}_p2p']=round(1.4826*np.nanmedian(np.abs(np.diff(fn_[ok])-np.nanmedian(np.diff(fn_[ok]))))/np.sqrt(2)*1e6)
        out[f'{col}_cdpp2.3h']=round(cdpp_like(tt[ok],fn_[ok],2.3)); out[f'{col}_cdpp6h']=round(cdpp_like(tt[ok],fn_[ok],6))
    rows.append(out); print(out, flush=True)
df=pd.DataFrame(rows).sort_values('G'); df.to_csv('kbonus_noise_vs_G.csv',index=False)
pd.set_option('display.width',250); print(df.to_string(index=False))
# compare with DR25 Kepler cool-dwarf median CDPP at the same magnitude
c=pd.read_csv('dl/dr25_cool_dwarfs.csv'); c=c[c.rrmscdpp06p0>0]
for _,r in df.iterrows():
    m=c[(c.kepmag>r.G-0.25)&(c.kepmag<r.G+0.25)]
    print(f"G {r.G}: KBONUS PSF 6h {r['flux_cdpp6h']} ppm vs Kepler DR25 cool-dwarf median 6h CDPP at Kp~G: {m.rrmscdpp06p0.median():.0f} ppm (n={len(m)}) -> ratio {r['flux_cdpp6h']/m.rrmscdpp06p0.median():.2f}" if len(m) else f"G {r.G}: no Kepler comparison stars")
