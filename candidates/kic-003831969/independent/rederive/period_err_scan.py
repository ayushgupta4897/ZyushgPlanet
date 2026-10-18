# finer period / t0 uncertainty from the BLS delta-loglike surface (own code)
import numpy as np, json
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-003831969_kepler_v1.0_lc.fits'
d=fits.getdata(f,'LIGHTCURVE_STITCHED')
m=d['SAP_QUALITY']==0
t=np.asarray(d['TIME'][m],float); fl=np.asarray(d['FLUX'][m],float); q=d['QUARTER'][m]
for qq in np.unique(q):
    s=q==qq; fl[s]/=np.median(fl[s])
# simple running median (1 d), transit-masked, same approach as the main script
P0,t00,dur=15.313409,144.28070,3.75/24
def rmed(t,f,use,w=1.0):
    tu,fu=t[use],f[use]; lo=np.searchsorted(tu,t-w/2); hi=np.searchsorted(tu,t+w/2)
    tr=np.full(len(t),np.nan)
    for i in range(len(t)):
        if hi[i]-lo[i]>=8: tr[i]=np.median(fu[lo[i]:hi[i]])
    return tr
dt=(t-t00+0.5*P0)%P0-0.5*P0
use=np.abs(dt)>dur
tr=rmed(t,fl,use); fd=fl/tr; r=fd-1
sig=1.4826*np.nanmedian(np.abs(r-np.nanmedian(r)))
g=np.isfinite(fd)&(r<4*sig)&(r>-10*sig)
t,fd=t[g],fd[g]
bls=BoxLeastSquares(t,fd,dy=sig)
per=np.linspace(P0-0.0005,P0+0.0005,4001)   # step 2.5e-7 d
res=bls.power(per,np.array([3.5,3.75,4.0])/24,objective='likelihood',oversample=40)
i=int(np.argmax(res.power)); ok=res.power>=res.power[i]-0.5
print('P_best',per[i],'1sig range',per[ok].min(),per[ok].max(),'-> sigma_P ~',0.5*(per[ok].max()-per[ok].min()))
# t0 scan at best P: chi2 of a box vs t0 shift
P=per[i]; dep=res.depth[i]; D=res.duration[i]; t0=res.transit_time[i]
sh=np.linspace(-0.15,0.15,601)*D
ll=[]
for s in sh:
    inn=np.abs((t-(t0+s)+0.5*P)%P-0.5*P)<D/2
    ll.append(-0.5*np.sum((fd[inn]-(1-dep))**2)/sig**2 + 0.5*np.sum((fd[inn]-1)**2)/sig**2)
ll=np.array(ll); j=int(np.argmax(ll)); ok=ll>=ll[j]-0.5
print('t0 best',t0+sh[j],'1sig half-width (min)',0.5*(sh[ok].max()-sh[ok].min())*1440)
json.dump(dict(P_best=float(per[i]),P_1sig_lo=float(per[ok].min()) if False else float(per[res.power>=res.power[i]-0.5].min()),
               P_1sig_hi=float(per[res.power>=res.power[i]-0.5].max()),
               t0_bkjd=float(t0+sh[j]),t0_1sig_min=float(0.5*(sh[ok].max()-sh[ok].min())*1440),
               duration_h=float(D*24),depth_ppm=float(dep*1e6)),
          open('perr.json','w'),indent=1)
