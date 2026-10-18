import numpy as np, json
from astropy.io import fits
from scipy.optimize import least_squares
P=15.31335; T0=144.28486; DUR=3.75/24
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-003831969_kepler_v1.0_lc.fits'
d=fits.open(f)['LIGHTCURVE_STITCHED'].data
t=d['TIME']; fl=d['FLUX']; q=d['QUARTER']; ok=(d['SAP_QUALITY']==0)&np.isfinite(t)&np.isfinite(fl); t,fl,q=t[ok],fl[ok],q[ok]
fn=np.empty_like(fl)
for qq in np.unique(q): m=q==qq; fn[m]=fl[m]/np.nanmedian(fl[m])
ph=((t-T0)/P+0.5)%1-0.5; mask=np.abs(ph)<(DUR/P)*0.8
src=fn.copy(); src[mask]=np.nan; r=np.empty_like(fn)
for i in range(len(t)):
    m=np.abs(t-t[i])<0.5; v=src[m]; v=v[np.isfinite(v)]; r[i]=fn[i]/np.nanmedian(v)-1 if len(v)>5 else np.nan
ok=np.isfinite(r); r=r[ok]; ph=ph[ok]
s=1.4826*np.nanmedian(np.abs(r)); keep=r<4*s; r=r[keep]; ph=ph[keep]
x=ph*P*24  # hours from mid-transit
w=np.abs(x)<12; x=x[w]; r=r[w]
def trap(p,x,ss=15):
    dep,T14,fin=p  # fin = ingress fraction of T14 (0..0.5)
    Ting=max(fin*T14,1e-3)
    # supersample over 29.4-min cadence
    off=np.linspace(-0.245,0.245,ss)
    xx=x[:,None]+off[None,:]
    a=np.abs(xx); m=np.zeros_like(xx)
    m[a<T14/2-Ting]=1
    sl=(a>=T14/2-Ting)&(a<T14/2); m[sl]=(T14/2-a[sl])/Ting
    return -dep*m.mean(1)
def res(p): return (trap(p,x)-r)/s
out={}
best=least_squares(res,[1e-3,4.2,0.2],bounds=([0,2,0.005],[5e-3,8,0.5])); c0=2*best.cost
print('free fit: depth %.0f ppm T14 %.2f h ingress frac %.3f  -> T23/T14=%.2f chi2=%.1f N=%d'%(best.x[0]*1e6,best.x[1],best.x[2],1-2*best.x[2],c0,len(x)))
for fin in [0.005,0.05,0.1,0.15,0.2,0.3,0.4,0.5]:
    def rr(p): return res([p[0],p[1],fin])
    b=least_squares(rr,[1e-3,4.2],bounds=([0,2],[5e-3,8])); print('  fixed ingress frac %.3f (T23/T14=%.2f): depth %.0f T14 %.2f dchi2=%.2f'%(fin,1-2*fin,b.x[0]*1e6,b.x[1],2*b.cost-c0))
# b-k implied: for k=0.033, box-shaped (ingress frac = k/(1+k)~0.032 for b=0)
