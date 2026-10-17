import numpy as np, time
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-006367225_kepler_v1.0_lc.fits'
d=fits.open(f)['LIGHTCURVE_STITCHED'].data
t=np.array(d['TIME'],float); fl=np.array(d['FLUX'],float); q=np.array(d['QUARTER']); qual=np.array(d['SAP_QUALITY'])
m=(qual==0)&np.isfinite(t)&np.isfinite(fl)&(fl>0)
t,fl,q=t[m],fl[m],q[m]
print('good cadences',len(t))
# per-quarter normalise
fn=np.empty_like(fl)
for qq in np.unique(q):
    s=q==qq; fn[s]=fl[s]/np.nanmedian(fl[s])
# running median 1 d, time-based, per segment
def running_median(t,y,win=1.0,gap=0.5):
    out=np.empty_like(y)
    brk=np.where(np.diff(t)>gap)[0]+1
    segs=np.split(np.arange(len(t)),brk)
    for seg in segs:
        ts=t[seg]; ys=y[seg]
        lo=np.searchsorted(ts,ts-win/2); hi=np.searchsorted(ts,ts+win/2,side='right')
        for k,(a,b) in enumerate(zip(lo,hi)):
            out[seg[k]]=np.median(ys[a:b])
    return out
t0=time.time(); tr=running_median(t,fn); print('runmed s',time.time()-t0)
y=fn/tr
mad=np.median(np.abs(y-np.median(y)))*1.4826
print('point-to-point MAD-std ppm',mad*1e6, 'ptp std/ sqrt2', np.std(np.diff(y))/np.sqrt(2)*1e6)
keep=y<1+4*mad
y=y[keep]; t2=t[keep]
dy=np.full_like(y,mad)
T=t2.max()-t2.min()
# grid: log-spaced with factor 1+dur_ref/(T*os)
dur_ref=0.08; os_=2
fac=1+dur_ref/(T*os_)
n=int(np.log(50/1.0)/np.log(fac))+1
print('n periods full grid',n)
periods=1.0*fac**np.arange(n)
durs=[0.05,0.08,0.12,0.18,0.25]
bls=BoxLeastSquares(t2,y,dy)
sub=periods[::20]
t0=time.time(); r=bls.power(sub,durs,objective='likelihood'); dt=time.time()-t0
print('subgrid',len(sub),'took',dt,'-> full est',dt*20/60,'min')
i=np.argmax(r.power); print('coarse best',r.period[i],r.depth[i],r.duration[i])
