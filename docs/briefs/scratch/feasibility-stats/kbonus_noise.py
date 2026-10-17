import numpy as np, lightkurve as lk, astropy.units as u, time, json, sys, warnings
from astropy.coordinates import SkyCoord
warnings.filterwarnings('ignore')
t0=time.time()
sr = lk.search_lightcurve(SkyCoord(285.907827,40.941988,unit='deg'), author='KBONUS-BKG', radius=20*u.arcsec)
print(sr)
# Keep only products for our Gaia source
sr = sr[[ 'Gaia DR3 2103792461826460288' in str(n) for n in sr.table['target_name']]]
print("N products for target:", len(sr))
lcc = sr.download_all(download_dir='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/feasibility-stats/lkcache')
print("download time %.1f s"%(time.time()-t0))
lc0 = lcc[0]
print("columns:", lc0.columns)
print("meta keys sample:", [k for k in list(lc0.meta.keys())[:60]])
for k in ['FLUX_COLUMN','PSFFRAC','CROWDSAP','KEPMAG','QUARTER','OBJECT','GMAG','RPMAG','BPMAG','KICID','TPFKIC','ORIGIN','NPIX']:
    if k in lc0.meta: print(k, lc0.meta[k])
import pickle
rows=[]
for lc in lcc:
    try:
        q = lc.meta.get('QUARTER')
        for col in ['flux','psf_flux','sap_flux','flux_corr','psf_flux_nova']:
            if col in lc.columns:
                fl = np.asarray(lc[col].value, dtype=float)
                tt = np.asarray(lc.time.value, dtype=float)
                ok = np.isfinite(fl)&np.isfinite(tt)&(fl>0)
                if ok.sum()<200: continue
                fl=fl[ok]; tt=tt[ok]
                med=np.nanmedian(fl); fn=fl/med
                # point-to-point scatter (MAD-based, like lightkurve estimate_cdpp's sigma via std of diffs/sqrt2)
                p2p = 1.4826*np.nanmedian(np.abs(np.diff(fn)-np.nanmedian(np.diff(fn))))/np.sqrt(2)
                # 6-h CDPP-like: biweight detrend (wotan) then std of 13-cadence (6.5h) boxcar means
                from wotan import flatten
                flat, trend = flatten(tt, fn, method='biweight', window_length=1.0, return_trend=True)
                good=np.isfinite(flat)
                flat=flat[good]
                n=13
                if flat.size>n*5:
                    cs=np.convolve(flat-1, np.ones(n)/n, mode='valid')
                    cdpp6 = 1.4826*np.nanmedian(np.abs(cs-np.nanmedian(cs)))*1e6
                else: cdpp6=np.nan
                rows.append(dict(quarter=q,col=col,npts=int(fl.size),p2p_ppm=float(p2p*1e6),cdpp6h_ppm=float(cdpp6)))
    except Exception as e:
        print("err", e)
import pandas as pd
df=pd.DataFrame(rows)
pd.set_option('display.width',200)
print(df.to_string())
print(df.groupby('col')[['p2p_ppm','cdpp6h_ppm']].median())
df.to_csv('kbonus_noise.csv',index=False)
print("total time %.1f s"%(time.time()-t0))
