import numpy as np, pandas as pd, lightkurve as lk, astropy.units as u, warnings, time
from astropy.coordinates import SkyCoord
from wotan import flatten
warnings.filterwarnings('ignore')
sel=pd.read_csv('kbonus_noise_sample.csv')
def cdpp_like(tt, fn_, hours, window=1.0):
    ok=np.isfinite(fn_)&np.isfinite(tt); tt=tt[ok]; fn_=fn_[ok]
    flat=flatten(tt, fn_, method='biweight', window_length=window); flat=flat[np.isfinite(flat)]
    n=max(1,int(round(hours/(29.4/60)))); cs=np.convolve(flat-1, np.ones(n)/n, mode='valid')
    return 1.4826*np.nanmedian(np.abs(cs-np.nanmedian(cs)))*1e6
rows=[]
for _,r in sel.iterrows():
    t0=time.time()
    try:
        sr=lk.search_lightcurve(SkyCoord(r.ra,r.dec,unit='deg'), author='KBONUS-BKG', radius=3*u.arcsec)
        sr=sr[[r.gaia_designation in str(n) for n in sr.table['target_name']]]
        if len(sr)==0: print("no product", r.gaia_designation); continue
        lc=sr[0].download(download_dir='lkcache')
        tt=np.asarray(lc.time.value,float); q=np.asarray(lc['quarter'].value)
        out=dict(rank=int(r['rank']), gaia=r.gaia_designation, kic=int(r.kic), G=round(r.phot_g_mean_mag,2), psffrac=round(r.psffrac,2), crowdsap=round(r.crowdsap,2), nq=int(r.n_psf_q), npts=len(lc))
        for col in ['flux','sap_flux']:
            fl=np.asarray(lc[col].value,float); ok=np.isfinite(fl)&(fl>0)
            fn_=np.full(fl.size,np.nan)
            for qq in np.unique(q[ok]):
                m=ok&(q==qq); fn_[m]=fl[m]/np.nanmedian(fl[m])
            out[f'{col}_p2p']=round(1.4826*np.nanmedian(np.abs(np.diff(fn_[ok])-np.nanmedian(np.diff(fn_[ok]))))/np.sqrt(2)*1e6)
            out[f'{col}_cdpp2.3h']=round(cdpp_like(tt[ok],fn_[ok],2.3)); out[f'{col}_cdpp6h']=round(cdpp_like(tt[ok],fn_[ok],6))
        out['sec']=round(time.time()-t0,1)
        rows.append(out); print(out, flush=True)
    except Exception as e:
        print("err", r.gaia_designation, repr(e)[:200], flush=True)
df=pd.DataFrame(rows); df.to_csv('kbonus_noise_vs_G.csv',index=False)
pd.set_option('display.width',250); print(df.to_string(index=False))
