# trapezoid fit with P and t0 free -> covariance-based period/t0 errors (own code)
import numpy as np, json
from astropy.io import fits
from scipy.optimize import least_squares
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-003831969_kepler_v1.0_lc.fits'
d=fits.getdata(f,'LIGHTCURVE_STITCHED'); m=d['SAP_QUALITY']==0
t=np.asarray(d['TIME'][m],float); fl=np.asarray(d['FLUX'][m],float); q=d['QUARTER'][m]
for qq in np.unique(q):
    s=q==qq; fl[s]/=np.median(fl[s])
P0,t00,dur=15.31336,144.2824,3.75/24
def rmed(t,f,use,w=1.0):
    tu,fu=t[use],f[use]; lo=np.searchsorted(tu,t-w/2); hi=np.searchsorted(tu,t+w/2)
    tr=np.full(len(t),np.nan)
    for i in range(len(t)):
        if hi[i]-lo[i]>=8: tr[i]=np.median(fu[lo[i]:hi[i]])
    return tr
dt=(t-t00+0.5*P0)%P0-0.5*P0
tr=rmed(t,fl,np.abs(dt)>dur); fd=fl/tr; r=fd-1
sig=1.4826*np.nanmedian(np.abs(r-np.nanmedian(r)))
g=np.isfinite(fd)&(r<4*sig)&(r>-10*sig); t,fd=t[g],fd[g]
tmid=0.5*(t.min()+t.max()); ncyc=np.round((tmid-t00)/P0); tref0=t00+ncyc*P0   # reference epoch near the middle
def model(p):
    P,tref,dep,T14,T23=p
    T23=min(max(T23,0),T14)
    x=np.abs((t-tref+0.5*P)%P-0.5*P)
    mm=np.ones_like(t)
    full=x<=T23/2; ing=(x>T23/2)&(x<T14/2)
    mm[full]=1-dep; mm[ing]=1-dep*(T14/2-x[ing])/((T14-T23)/2+1e-12)
    return mm
res=least_squares(lambda p:(fd-model(p))/sig,[P0,tref0,1.0e-3,4.3/24,2.7/24],
                  x_scale=[1e-4,1e-3,1e-4,1e-2,1e-2],diff_step=[2e-6,5e-5,1e-3,1e-3,1e-3])
J=res.jac; cov=np.linalg.inv(J.T@J); chi2=np.sum(res.fun**2); dof=len(fd)-5
cov*=chi2/dof; err=np.sqrt(np.diag(cov))
P,tref,dep,T14,T23=res.x
out=dict(P_d=float(P),P_err_d=float(err[0]),tref_bkjd=float(tref),tref_err_min=float(err[1]*1440),
         tref_bjd=float(tref+2454833),t0_first_bkjd=float(tref-ncyc*P),
         depth_ppm=float(dep*1e6),depth_err_ppm=float(err[2]*1e6),T14_h=float(T14*24),T14_err_h=float(err[3]*24),
         T23_h=float(T23*24),T23_err_h=float(err[4]*24),T23_over_T14=float(T23/T14),chi2_red=float(chi2/dof),
         rp_rs_from_depth=float(np.sqrt(dep)))
print(json.dumps(out,indent=1)); json.dump(out,open('trapfit.json','w'),indent=1)
