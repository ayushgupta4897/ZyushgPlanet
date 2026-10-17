import numpy as np, json
from astropy.io import fits
P=2.377014357204058; t0=132.15453895516694; dur=1.125/24
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-009520080_kepler_v1.0_lc.fits'
d=fits.open(f)['LIGHTCURVE_STITCHED'].data
m=np.isfinite(d['TIME'])&np.isfinite(d['FLUX'])&(d['SAP_QUALITY']==0)&(d['FLUX']>0)&np.isfinite(d['SAP_FLUX'])
t=d['TIME'][m]; q=d['QUARTER'][m]
ph=((t-t0)/P+0.5)%1-0.5; intr=np.abs(ph*P)<0.75*dur
def detrend(flux):
    y=flux.copy()
    for qq in np.unique(q):
        s=q==qq; y[s]=y[s]/np.nanmedian(y[s])
    tr=np.full_like(y,np.nan)
    for i in range(len(t)):
        w=(np.abs(t-t[i])<0.5)&(~intr); tr[i]=np.median(y[w]) if w.sum()>5 else np.nan
    return y/tr-1
R={}
for nm,col in [('PSF','FLUX'),('SAP','SAP_FLUX')]:
    r=detrend(d[col][m]); sig=1.4826*np.nanmedian(np.abs(r-np.nanmedian(r)))
    r[(r>5*sig)|(r<-15*sig)]=np.nan; R[nm]=(r,sig)
r,sig=R['PSF']; ok=np.isfinite(r)
sin=np.abs(ph*P)<dur/2
out={}
for nm,(rr,ss) in R.items():
    o=np.isfinite(rr); print(nm,'depth %.0f ppm (n %d), sigma %.0f'%(-np.mean(rr[o&sin])*1e6,(o&sin).sum(),ss*1e6)); out[nm+'_depth_ppm']=float(-np.mean(rr[o&sin])*1e6)
# phase scan with mean / median / clipped mean; box = transit duration; step 0.005
def box(c,rr):
    s=np.isfinite(rr)&(np.abs(((ph-c+0.5)%1)-0.5)*P<dur/2); x=rr[s]*1e6
    xc=x[np.abs(x-np.median(x))<3*sig*1e6]; return np.array([x.mean(),np.median(x),xc.mean(),len(x)])
centers=np.arange(0.06,0.941,0.005); scan=np.array([box(c,r) for c in centers])
k=np.argmin(np.abs(centers-0.5)); far=np.abs(centers-0.5)>0.03
for j,lab in enumerate(['mean','median','clipped_mean']):
    z=(scan[k,j]-scan[far,j].mean())/scan[far,j].std(); rank=(scan[:,j]>=scan[k,j]).sum()
    print(f'phase-0.5 {lab}: {scan[k,j]:+.0f} ppm; scan std {scan[far,j].std():.0f}; z={z:.2f}; rank (1=most positive/brightest) {rank}/{len(centers)}')
    out['phase05_'+lab]=dict(value_ppm=float(scan[k,j]),scan_std_ppm=float(scan[far,j].std()),z=float(z),rank=int(rank),nboxes=len(centers))
# most extreme boxes anywhere (positive and negative) for the mean
o=np.argsort(scan[:,0]); print('5 most positive boxes (phase, mean):',[(round(centers[i],3),round(scan[i,0])) for i in o[-5:][::-1]]); print('5 most negative:',[(round(centers[i],3),round(scan[i,0])) for i in o[:5]])
# contribution of positive outliers to the phase-0.5 mean
s05=np.isfinite(r)&(np.abs(((ph-0.5+0.5)%1)-0.5)*P<dur/2); x=r[s05]*1e6
print('phase-0.5 box: n %d, mean %.0f, median %.0f, mean of points < +3sigma %.0f, n points > +3 sigma: %d (expected for Gaussian: %.1f)'%(len(x),x.mean(),np.median(x),x[x<3*sig*1e6].mean(),(x>3*sig*1e6).sum(),len(x)*0.00135))
# control: same stats at the in-transit box with inverted sign for reference and typical oot box
sx=np.isfinite(r)&(np.abs(((ph-0.25+0.5)%1)-0.5)*P<dur/2); y=r[sx]*1e6
print('phase-0.25 control box: n %d, mean %.0f, median %.0f, n>3sig %d'%(len(y),y.mean(),np.median(y),(y>3*sig*1e6).sum()))
# fraction of ALL oot points above 3 sigma
oot=np.isfinite(r)&(np.abs(ph)>0.1); print('all-oot fraction >3sigma: %.5f (n %d)'%((r[oot]>3*sig).mean(),oot.sum()))
# odd/even timing (argmin) and shape by width
ev=np.round((t-t0)/P).astype(int)
def bd(sel,off,w=dur):
    s=sel&np.isfinite(r)&(np.abs(ph*P-off)<w/2); return r[s].mean()*1e6
for lab,sel in [('odd',ev%2==1),('even',ev%2==0),('all',ev>-9)]:
    offs=np.arange(-0.05,0.0501,0.002); dd=np.array([bd(sel,o) for o in offs]); i=np.argmin(dd)
    # parabolic refinement
    if 0<i<len(offs)-1:
        a,b,c=dd[i-1],dd[i],dd[i+1]; x0=offs[i]+0.002*0.5*(a-c)/(a-2*b+c)
    else: x0=offs[i]
    print(f'{lab}: best centre offset {x0*1440:+.1f} min, depth {dd[i]:.0f} ppm; depth vs box width (h): '+', '.join(f'{w}:{bd(sel,0,w/24):.0f}' for w in (0.5,0.75,1.0,1.125,1.5,2.0,2.5)))
    out[lab+'_offset_min']=float(x0*1440)
json.dump(out,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic009520080/out/lc_tests2.json','w'),indent=1)
