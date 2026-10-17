import numpy as np, time
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-006448637_kepler_v1.0_lc.fits'
d=fits.getdata(f,'LIGHTCURVE_STITCHED')
t=d['TIME']; fl=d['FLUX']; fe=d['FLUX_ERR']; sap=d['SAP_FLUX']; q=d['SAP_QUALITY']; qu=d['QUARTER']; fm=d['FLATTEN_MASK']
print('N',len(t),'finite time',np.isfinite(t).sum(),'finite flux',np.isfinite(fl).sum(),'finite sap',np.isfinite(sap).sum())
print('quality==0:',(q==0).sum(),' unique quality bits sample:',np.unique(q)[:20])
print('quarters:',np.unique(qu), [ (int(k),(qu==k).sum()) for k in np.unique(qu)])
print('FLATTEN_MASK unique',np.unique(fm), (fm!=0).sum())
g=(q==0)&np.isfinite(t)&np.isfinite(fl)&np.isfinite(fe)
print('good',g.sum(),'tmin',t[g].min(),'tmax',t[g].max(),'baseline',t[g].max()-t[g].min())
dt=np.diff(t[g]); print('median cadence (d)',np.median(dt), 'in min',np.median(dt)*1440)
for k in np.unique(qu):
    m=g&(qu==k)
    print(int(k), m.sum(), 'tspan %.1f-%.1f'%(t[m].min(),t[m].max()), 'median flux %.1f'%np.median(fl[m]), 'sap %.1f'%np.median(sap[m]), 'p2p-std ppm %.0f'%(1e6*np.std(np.diff(fl[m]))/np.sqrt(2)/np.median(fl[m])), 'median err ppm %.0f'%(1e6*np.median(fe[m])/np.median(fl[m])))
# time a small BLS
tt=t[g]; ff=fl[g]/np.median(fl[g])
bls=BoxLeastSquares(tt,ff)
per=np.exp(np.linspace(np.log(1),np.log(50),5000))
t0=time.time(); r=bls.power(per,[0.08,0.12,0.17,0.25,0.35],objective='snr'); print('5000 periods x5 dur: %.1fs'%(time.time()-t0))
