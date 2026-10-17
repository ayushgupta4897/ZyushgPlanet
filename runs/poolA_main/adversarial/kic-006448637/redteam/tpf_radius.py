import glob, json, numpy as np
from astropy.io import fits
from astropy.wcs import WCS
P=35.58726; T0=186.68617; DUR=4.5/24
ra,dec=293.6827470811608,41.85317701406968; bra,bdec=293.662990,41.874123; nra,ndec=293.6866113558775,41.85313790492249
def phase(t): return ((t-T0)/P+0.5)%1-0.5
rows=[]; cen=[]; 
for fpath in sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr006448574_lc_Q*/*.fits.gz')):
    q=fits.getheader(fpath,0)['QUARTER']
    if q in (0,1,3,6,7,9,17): continue
    d=fits.getdata(fpath,1); h2=fits.getheader(fpath,2); w=WCS(h2)
    tx,ty=[float(v) for v in w.all_world2pix(ra,dec,0)]; bx,by=[float(v) for v in w.all_world2pix(bra,bdec,0)]; nx,ny=[float(v) for v in w.all_world2pix(nra,ndec,0)]
    t=d['TIME']; fl=d['FLUX']; good=np.isfinite(t)&(d['QUALITY']==0); t=t[good]; fl=fl[good]
    ny_,nx_=fl.shape[1:]; yy,xx=np.mgrid[0:ny_,0:nx_]; valid=np.isfinite(fl).mean(axis=0)>0.9; bleed=np.abs(xx-round(bx))<=0.5
    # per-event local differenced image
    intr=np.abs(phase(t))*P<DUR/2
    dm=np.zeros((ny_,nx_)); vm=np.zeros((ny_,nx_))
    for e in np.unique(np.round((t[intr]-T0)/P)):
        tc=T0+e*P; mi=np.abs(t-tc)<DUR/2; mo=(np.abs(t-tc)<2.0)&(np.abs(t-tc)>0.75*DUR)
        if mi.sum()<3 or mo.sum()<20: continue
        diff=np.nanmean(fl[mo],axis=0)-np.nanmean(fl[mi],axis=0); var=np.nanvar(fl[mo],axis=0)*(1/mi.sum()+1/mo.sum()); var[~np.isfinite(var)|(var<=0)]=np.inf
        dm+=diff/var; vm+=1/var
    dmap=dm/vm; dmap[~np.isfinite(dmap)]=0
    line=f"Q{q:2d}"
    for r in (1.0,1.6,2.5,3.5):
        m=(np.hypot(xx-tx,yy-ty)<=r)&valid&~bleed
        s=dmap[m].sum(); e=np.sqrt((1/vm[m]).sum()); line+=f"  r<={r}: {s:5.1f}+-{e:4.1f} e/s ({m.sum():2d}px)"
    # also annulus 2.5-4.5 px (should be ~0 if on-target)
    m=(np.hypot(xx-tx,yy-ty)>2.5)&(np.hypot(xx-tx,yy-ty)<=4.5)&valid&~bleed; s=dmap[m].sum(); e=np.sqrt((1/vm[m]).sum()); line+=f"  annulus2.5-4.5: {s:5.1f}+-{e:4.1f} ({m.sum()}px)"
    # deficit centroid within r<=2.5 of target (positive deficits only weighted)
    m=(np.hypot(xx-tx,yy-ty)<=2.5)&valid&~bleed&(dmap>0)
    cx=(dmap[m]*xx[m]).sum()/dmap[m].sum(); cy=(dmap[m]*yy[m]).sum()/dmap[m].sum()
    line+=f"  centroid offset from target: dx={cx-tx:+.2f} dy={cy-ty:+.2f} px (toward nbr: dx={nx-tx:+.1f},dy={ny-ty:+.1f}; toward host dy={by-ty:+.0f})"
    print(line); cen.append((cx-tx,cy-ty))
cen=np.array(cen); print('mean centroid offset over quarters: dx=%.2f+-%.2f dy=%.2f+-%.2f px (1 px = 3.98 arcsec)'%(cen[:,0].mean(),cen[:,0].std()/np.sqrt(len(cen)),cen[:,1].mean(),cen[:,1].std()/np.sqrt(len(cen))))
