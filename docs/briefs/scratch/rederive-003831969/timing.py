import time, numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-003831969_kepler_v1.0_lc.fits'
d=fits.getdata(f,'LIGHTCURVE_STITCHED')
m=d['SAP_QUALITY']==0
t=np.asarray(d['TIME'][m],float); fl=np.asarray(d['FLUX'][m],float)
for q in np.unique(d['QUARTER'][m]):
    s=d['QUARTER'][m]==q; fl[s]/=np.median(fl[s])
T=t.max()-t.min()
q=(2/24)/(T*2.0); n=int(np.ceil(np.log(50)/-np.log(1-q)))
print('baseline',T,'grid size at OS=2, dur_min=2h:',n)
bls=BoxLeastSquares(t,fl,dy=0.001)
per=np.linspace(10,10.5,3000)
t1=time.time(); r=bls.power(per,np.array([2,3,4,5,6,8,10])/24,objective='likelihood',oversample=10); dt=time.time()-t1
print('3000 periods x 7 durations:',dt,'s ->',dt/3000*n/60,'min for full grid')
per=np.linspace(45,46,500)
t1=time.time(); r=bls.power(per,np.array([2,3,4,5,6,8,10])/24,objective='likelihood',oversample=10); dt=time.time()-t1
print('500 long periods:',dt,'s per 500')
