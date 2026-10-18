# Own pixel test: per-event local differencing on KIC 9655236 TPFs. Unclipped total stamp deficit (e/s),
# per-pixel deficit at the target / G=15.9 neighbour / host pixels, and difference-image centroid with
# bootstrap errors. Quarter frames combined in e/s.
import numpy as np, glob, json
from astropy.io import fits
from astropy.wcs import WCS
P=6.703464464595031; t0=265.183161835447; dur=2.4375/24
srcs={'target':(294.28267299472,46.34150351561),'G15.9_835648':(294.28296961419,46.34006010397),'host_G16.4_837184':(294.2852862259,46.33963206258),'G18.4_646464':(294.28723639435,46.33632135739)}
files=sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr009655236_lc_Q000111011101111000/*.fits.gz'))
res={}; tot=[]; pixsum={k:[] for k in srcs}
for fn in files:
    h=fits.open(fn); Q=h[0].header['QUARTER']; d=h[1].data; ap=h[2].data
    t=d['TIME']; fl=d['FLUX']; qu=d['QUALITY']
    goodpx=np.isfinite(fl).mean(axis=0)>0.9
    ok=np.isfinite(t)&(qu==0)&np.isfinite(fl[:,goodpx]).all(axis=1)
    t,fl=t[ok],fl[ok]; fl=np.where(goodpx[None],fl,np.nan)
    w=WCS(h[2].header)  # pixel wcs of stamp
    pos={}
    for k,(ra,dec) in srcs.items():
        x,y=w.all_world2pix(ra,dec,0); pos[k]=(float(x),float(y))
    ny,nx=fl.shape[1:]
    # events
    n0=int(np.floor((t.min()-t0)/P)); n1=int(np.ceil((t.max()-t0)/P))
    diffs=[]; 
    for n in range(n0,n1+1):
        tc=t0+n*P; dt=t-tc
        ib=np.abs(dt)<dur/2*0.9; ob=(np.abs(dt)>dur/2*1.5)&(np.abs(dt)<dur/2*1.5+0.25)
        if ib.sum()<3 or ob.sum()<8: continue
        # local linear trend per pixel from out-of-transit, evaluated at in-transit times
        A=np.vstack([np.ones(ob.sum()),dt[ob]]).T
        coef=np.linalg.lstsq(A,fl[ob].reshape(ob.sum(),-1),rcond=None)[0]
        model=(np.vstack([np.ones(ib.sum()),dt[ib]]).T@coef).reshape(ib.sum(),ny,nx)
        diffs.append((model-fl[ib]).mean(axis=0))  # positive = light lost
    diffs=np.array(diffs); nev=len(diffs)
    if nev<3: continue
    D=diffs.mean(axis=0); E=diffs.std(axis=0,ddof=1)/np.sqrt(nev)
    valid=np.isfinite(D)&(ap>0) if np.any(ap>0) else np.isfinite(D)
    valid=np.isfinite(D)&goodpx
    total=np.nansum(D[valid]); toterr=np.sqrt(np.nansum(E[valid]**2))
    # per-source pixel (nearest pixel) deficit, and 3x3 sum around it
    out={'Q':int(Q),'n_events':nev,'stamp':[int(ny),int(nx)],'total_deficit_e_s':[float(total),float(toterr)],'pos':pos,'px':{}}
    for k,(x,y) in pos.items():
        xi,yi=int(round(x)),int(round(y))
        if 0<=xi<nx and 0<=yi<ny and np.isfinite(D[yi,xi]):
            out['px'][k]=[float(D[yi,xi]),float(E[yi,xi])]; pixsum[k].append((D[yi,xi],E[yi,xi]))
        else: out['px'][k]=None
    # difference-image centroid (all pixels, unclipped, weights = D) only if total significant
    yy,xx=np.mgrid[0:ny,0:nx]
    if total>0:
        cx=np.nansum((D*xx)[valid])/total; cy=np.nansum((D*yy)[valid])/total
        # bootstrap over events
        bs=[]
        for b in range(300):
            idx=np.random.RandomState(Q).randint(0,nev,nev); Db=diffs[idx].mean(axis=0); tb=np.nansum(Db[valid])
            if tb>0: bs.append((np.nansum((Db*xx)[valid])/tb,np.nansum((Db*yy)[valid])/tb))
        bs=np.array(bs)
        out['diff_centroid_px']=[float(cx),float(cy)]; out['diff_centroid_boot_std']=[float(bs[:,0].std()),float(bs[:,1].std())] if len(bs) else None
        out['dist_to_src_px']={k:float(np.hypot(cx-x,cy-y)) for k,(x,y) in pos.items()}
    # flux-weighted out-of-transit centroid for reference
    img=np.nanmedian(fl,axis=0); out['flux_centroid_px']=[float(np.nansum(img*xx)/np.nansum(img)),float(np.nansum(img*yy)/np.nansum(img))]
    out['target_in_stamp_margin_px']=float(min(pos['target'][0],nx-1-pos['target'][0],pos['target'][1],ny-1-pos['target'][1]))
    tot.append((total,toterr)); res[int(Q)]=out
    print('Q%2d nev=%2d stamp %dx%d total deficit %.2f +/- %.2f e/s | target px %.2f,%.2f (margin %.2f) | per-px: %s | diffcent dist: %s'%(Q,nev,ny,nx,total,toterr,pos['target'][0],pos['target'][1],out['target_in_stamp_margin_px'],{k:(None if v is None else '%.2f+/-%.2f'%tuple(v)) for k,v in out['px'].items()},{k:'%.1f'%v for k,v in out.get('dist_to_src_px',{}).items()}))
tot=np.array(tot); w=1/tot[:,1]**2; print('\nWEIGHTED MEAN total stamp deficit: %.2f +/- %.2f e/s (expected ~1.07 e/s if the eclipse is on the target and fully inside the stamp)'%(np.sum(w*tot[:,0])/np.sum(w),1/np.sqrt(np.sum(w))))
for k,v in pixsum.items():
    v=np.array(v); 
    if len(v): w=1/v[:,1]**2; print('weighted mean deficit in pixel nearest %s: %.3f +/- %.3f e/s (%d quarters)'%(k,np.sum(w*v[:,0])/np.sum(w),1/np.sqrt(np.sum(w)),len(v)))
json.dump(res,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic009655224/work/pixel_check.json','w'),indent=1,default=float)
