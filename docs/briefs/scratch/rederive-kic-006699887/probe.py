import numpy as np, time, pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
fn='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-006699887_kepler_v1.0_lc.fits'
d=fits.getdata(fn,'LIGHTCURVE_STITCHED')
t=d['TIME'].astype(float); f=d['FLUX'].astype(float); fe=d['FLUX_ERR'].astype(float); q=d['QUARTER'].astype(int); qual=d['SAP_QUALITY'].astype(int)
m=(qual==0)&np.isfinite(t)&np.isfinite(f)&np.isfinite(fe)&(f>0)
print('N total',len(t),'good',m.sum()); print('quarters',np.unique(q[m]))
print('time span',t[m].min(),t[m].max())
print('median cadence (d)',np.median(np.diff(t[m])))
print('FLATTEN_MASK unique',np.unique(d['FLATTEN_MASK']))
print('SAP_QUALITY unique counts',np.unique(qual,return_counts=True))
# per-quarter stats
for qq in np.unique(q[m]):
    s=m&(q==qq); print(qq, s.sum(), np.nanmedian(f[s]), np.nanstd(f[s])/np.nanmedian(f[s])*1e6, 'ppm; SAP med', np.nanmedian(d['SAP_FLUX'][s]))
