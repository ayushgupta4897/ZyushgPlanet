import numpy as np, glob
from astropy.io import fits
from astropy.wcs import WCS
P=15.31335; T0=144.28486; DUR=3.75/24
tgt=(285.0768009047113,38.93627480915959)
np.set_printoptions(precision=2, suppress=True, linewidth=150)
for qsel in [[4,8,12,16],[2,6,10,14]]:
    Ds=[];Ws=[];meds=[]
    for f in sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr003831977_lc_Q011111111111111111/*.fits.gz')):
        h=fits.open(f); q=h[0].header['QUARTER']
        if q not in qsel: continue
        d=h[1].data; t=d['TIME']; fl=d['FLUX']; qual=d['QUALITY']
        ok=(qual==0)&np.isfinite(t)&(np.nansum(fl,axis=(1,2))>0); t,fl=t[ok],fl[ok]
        w=DUR/P; ev=np.round((t-T0)/P).astype(int); diffs=[];wts=[]
        for e in np.unique(ev):
            m=ev==e; tt=t[m]; ff=fl[m]; pp=(tt-(T0+e*P))/P
            i=np.abs(pp)<0.4*w; o=(np.abs(pp)>1.2*w)&(np.abs(pp)<4*w)
            if i.sum()<3 or o.sum()<8: continue
            x=tt-(T0+e*P); A=np.vstack([np.ones(o.sum()),x[o]]).T
            coef=np.linalg.lstsq(A,ff[o].reshape(o.sum(),-1),rcond=None)[0]
            base_in=(np.vstack([np.ones(i.sum()),x[i]]).T@coef).reshape(i.sum(),*ff.shape[1:])
            dimg=(base_in-ff[i]).mean(0); resid=ff[o]-(A@coef).reshape(o.sum(),*ff.shape[1:])
            var=resid.var(0)*(1/i.sum()+1/o.sum()); diffs.append(dimg); wts.append(1/np.maximum(var,1e-6))
        diffs=np.array(diffs); wts=np.array(wts); Ds.append((diffs*wts).sum(0)/wts.sum(0)); Ws.append(wts.sum(0)); meds.append(np.nanmedian(fl,0))
        wcs=WCS(h[2].header); tx,ty=wcs.all_world2pix(tgt[0],tgt[1],0)
    W=np.array(Ws); D=(np.array(Ds)*W).sum(0)/W.sum(0); E=1/np.sqrt(W.sum(0))
    print('quarters',qsel,'target px (x,y)=',round(float(tx),2),round(float(ty),2))
    print('median image (e/s):'); print(np.nanmean(meds,0))
    print('deficit (e/s):'); print(D)
    print('deficit SNR:'); print(D/E)
    print('per-pixel chi2 contribution (null):'); print((D/E)**2)
    # scaled-image model: deficit proportional to image (all stars dim together / host dims)
    img=np.nanmean(meds,0); ok=np.isfinite(D)
    a=np.nansum(D[ok]*img[ok]/E[ok]**2)/np.nansum(img[ok]**2/E[ok]**2); print('scaled-image amp',a*1e6,'ppm, chi2',np.nansum(((D-a*img)/E)[ok]**2))
    # expected deficit from target gaussian at (tx,ty), sigma 0.6, total 3.2 e/s*depth
    yy,xx=np.indices(D.shape); g=np.exp(-((xx-tx)**2+(yy-ty)**2)/(2*0.6**2)); g/=g.sum()*1.0
    A=np.nansum(D[ok]*g[ok]/E[ok]**2)/np.nansum(g[ok]**2/E[ok]**2); print('target-gaussian total deficit in stamp',A,'e/s, chi2',np.nansum(((D-A*g)/E)[ok]**2), 'residual SNR map:'); print(((D-A*g)/E))
