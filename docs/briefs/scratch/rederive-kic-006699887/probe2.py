import numpy as np, time, pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
fn='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-006699887_kepler_v1.0_lc.fits'
d=fits.getdata(fn,'LIGHTCURVE_STITCHED')
t=d['TIME'].astype(float); f=d['FLUX'].astype(float); fe=d['FLUX_ERR'].astype(float); q=d['QUARTER'].astype(int); qual=d['SAP_QUALITY'].astype(int)
m=(qual==0)&np.isfinite(t)&np.isfinite(f)&np.isfinite(fe)&(f>0)
t,f,fe,q=t[m],f[m],fe[m],q[m]
# quick normalise+detrend per quarter via time-based rolling median
fl=np.empty_like(f)
for qq in np.unique(q):
    s=q==qq
    med=np.median(f[s]); ff=f[s]/med
    ser=pd.Series(ff,index=pd.to_timedelta(t[s],unit='D'))
    tr=ser.rolling('1D',center=True,min_periods=5).median().values
    fl[s]=ff/tr
fe_n=fe/np.median(f)
print('detrended std ppm',np.std(fl)*1e6, 'MAD ppm',1.4826*np.median(np.abs(fl-np.median(fl)))*1e6)
bls=BoxLeastSquares(t,fl,fe_n)
durs=np.array([1.5,2,3,4,5,6,8])/24
P=np.exp(np.linspace(np.log(1),np.log(50),10000))
t0=time.time(); r=bls.power(P,durs,objective='likelihood'); print('10k periods, 7 durations:',time.time()-t0,'s')
i=np.argmax(r.power); print('best',P[i],r.depth[i]*1e6,r.duration[i]*24)
