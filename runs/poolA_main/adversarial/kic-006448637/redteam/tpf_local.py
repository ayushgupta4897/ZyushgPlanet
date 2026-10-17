# Red-team pixel test: box photometry at the target corner of the saturated-star TPF (KIC 6448574), independent of pipeline/
import glob, json, numpy as np
from astropy.io import fits
from astropy.wcs import WCS
P=35.58726; T0=186.68617; DUR=4.5/24
ra,dec=293.6827470811608,41.85317701406968           # target KIC 6448637
nra,ndec=293.6866113558775,41.85313790492249          # KIC 6448650 (G=16.03) 10.4" E
bra,bdec=293.662990,41.874123                         # saturated host KIC 6448574 Kp 8.75
def runmed(t,f,w=1.0,mask=None):
    out=np.empty_like(f)
    for i in range(len(t)):
        m=(np.abs(t-t[i])<w/2)
        if mask is not None: m&=~mask
        out[i]=np.nanmedian(f[m]) if m.sum()>5 else np.nan
    return out
def phase(t): return ((t-T0)/P+0.5)%1-0.5
def box_depth(t,f,intr):
    # f normalised, detrended; depth = mean(out)-mean(in); err from random-phase boxes
    d=np.nanmean(f[~intr])-np.nanmean(f[intr])
    rng=np.random.default_rng(1); vals=[]
    ph=phase(t)
    for k in range(600):
        c=rng.uniform(0.08,0.92)-0.5
        m=np.abs(ph-c)*P<DUR/2
        if m.sum()>=0.5*intr.sum(): vals.append(np.nanmean(f[~m])-np.nanmean(f[m]))
    return d, np.nanstd(vals), intr.sum()
res={}
maps=[]
for fpath in sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr006448574_lc_Q*/*.fits.gz')):
    h0=fits.getheader(fpath,0); q=h0['QUARTER']
    d=fits.getdata(fpath,1); h2=fits.getheader(fpath,2); ap=(fits.getdata(fpath,2)&2)>0
    w=WCS(h2); tx,ty=[float(v) for v in w.all_world2pix(ra,dec,0)]; nx,ny=[float(v) for v in w.all_world2pix(nra,ndec,0)]; bx,by=[float(v) for v in w.all_world2pix(bra,bdec,0)]
    t=d['TIME']; fl=d['FLUX']; qual=d['QUALITY']
    good=np.isfinite(t)&(qual==0)
    t=t[good]; fl=fl[good]
    ny_,nx_=fl.shape[1:]
    if not (0<=ty<ny_-0.5): res[q]={'note':'target outside stamp','ty':ty}; print(f'Q{q}: target outside stamp (y={ty:.2f})'); continue
    # box pixels around target: radius 1.5 px, exclude pixels within 1 px of the saturated column (bx) when |y-by| large? bleed runs along columns at x~bx
    yy,xx=np.mgrid[0:ny_,0:nx_]
    valid=np.isfinite(fl).mean(axis=0)>0.9
    tb=(np.hypot(xx-tx,yy-ty)<=1.6)&valid
    nb=(np.hypot(xx-nx,yy-ny)<=1.6)&valid
    bleed=(np.abs(xx-round(bx))<=0.5)
    ph=phase(t); intr=np.abs(ph)*P<DUR/2
    nev=len(np.unique(np.round((t[intr]-T0)/P)))
    out={}
    for name,m in [('target_box',tb&~bleed),('nbr_box',nb&~bleed),('host_aperture',ap),('target_box_inclbleed',tb)]:
        if m.sum()==0: continue
        f=np.nansum(fl[:,m],axis=1)
        tr=runmed(t,f,1.0,mask=intr); fn=f/tr
        ok=np.isfinite(fn)
        dep,err,nin=box_depth(t[ok],fn[ok],intr[ok])
        out[name]={'npix':int(m.sum()),'mean_e_s':float(np.nanmedian(f)),'depth_ppm':float(dep*1e6),'err_ppm':float(err*1e6),'depth_e_s':float(dep*np.nanmedian(f)),'n_in':int(nin)}
    # local deficit map: per event, in-transit mean minus local out-of-transit mean (within +-2 d, excluding +-1.5 dur), weighted sum over events
    dm=np.zeros((ny_,nx_)); vm=np.zeros((ny_,nx_))
    for e in np.unique(np.round((t[intr]-T0)/P)):
        tc=T0+e*P
        mi=np.abs(t-tc)<DUR/2; mo=(np.abs(t-tc)<2.0)&(np.abs(t-tc)>0.75*DUR)
        if mi.sum()<3 or mo.sum()<20: continue
        diff=np.nanmean(fl[mo],axis=0)-np.nanmean(fl[mi],axis=0)
        var=np.nanvar(fl[mo],axis=0)*(1/mi.sum()+1/mo.sum()); var[~np.isfinite(var)|(var<=0)]=np.inf
        dm+=diff/var; vm+=1/var
    dmap=dm/vm; smap=dmap*np.sqrt(vm)   # deficit (e/s) and SNR
    # 7x7 cutouts around target
    iy,ix=int(round(ty)),int(round(tx))
    y0,y1=max(0,iy-3),min(ny_,iy+4); x0,x1=max(0,ix-3),min(nx_,ix+4)
    cut_s=smap[y0:y1,x0:x1]; cut_d=dmap[y0:y1,x0:x1]
    # significance of deficit in target box vs neighbour box (sum of deficits / sqrt(sum var))
    def boxsig(m):
        m=m&np.isfinite(dmap)
        return float(dmap[m].sum()), float(np.sqrt((1/vm[m]).sum()))
    out['deficit_target_box_e_s']=boxsig(tb&~bleed); out['deficit_nbr_box_e_s']=boxsig(nb&~bleed)
    out['n_events']=int(nev); out['target_px']=[tx,ty]; out['nbr_px']=[nx,ny]; out['host_px']=[bx,by]; out['stamp']=[ny_,nx_]
    out['snr_cutout_rows_y%d-%d_cols_x%d-%d'%(y0,y1-1,x0,x1-1)]=np.round(cut_s,1).tolist()
    res[q]=out
    print(f"Q{q:2d} nev={nev} target=({tx:.1f},{ty:.1f}) nbr=({nx:.1f},{ny:.1f}) host=({bx:.1f},{by:.1f}) stamp={ny_}x{nx_}")
    for k in ('target_box','target_box_inclbleed','nbr_box','host_aperture'):
        if k in out: o=out[k]; print(f"   {k:22s} npix={o['npix']:3d} flux={o['mean_e_s']:9.0f} e/s depth={o['depth_ppm']:7.0f}+-{o['err_ppm']:5.0f} ppm = {o['depth_e_s']:6.1f} e/s n_in={o['n_in']}")
    print(f"   deficit sum: target box {out['deficit_target_box_e_s'][0]:.2f}+-{out['deficit_target_box_e_s'][1]:.2f} e/s ; nbr box {out['deficit_nbr_box_e_s'][0]:.2f}+-{out['deficit_nbr_box_e_s'][1]:.2f} e/s")
    print('   SNR cutout (rows y%d..%d, cols x%d..%d; target at row %d col %d):'%(y0,y1-1,x0,x1-1,iy,ix)); print(np.array2string(np.round(cut_s,1),max_line_width=200))
json.dump(res,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic006448637/tpf_local.json','w'),indent=1,default=float)
