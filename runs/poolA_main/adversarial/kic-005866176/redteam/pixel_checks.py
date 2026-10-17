"""Red-team pixel-level check for kic-005866176 on the Kepler TPFs of the host KIC 5866192 (own code).
Event-differenced deficit images stacked per season; deficit centroid vs Gaia positions of target and host;
template (out-of-transit image, shifted) fit at target vs host. Also run with post-gap events excluded."""
import json, glob, numpy as np, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import shift as ndshift
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005866176/redteam/'
TPFS=sorted(glob.glob('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam_kic005866176/tpf/mastDownload/Kepler/*/kplr005866192-*_lpd-targ.fits.gz'))
P=29.433016127730333; T0=172.57296158611774; DUR=1.95/24.
RA_T,DEC_T=288.4755283869849,41.10860915391378   # target KIC 5866176 (Gaia DR3 2102115328638366592)
RA_H,DEC_H=288.47974,41.1105                     # host KIC 5866192
POSTGAP_EVENTS=[231.476,937.845]  # the two post-gap deep events (BKJD)
R={'tpfs':[], 'seasons':{}}
stacks={}
for path in TPFS:
    f=fits.open(path); d=f[1].data; h1=f[1].header; h0=f[0].header
    q=h0['QUARTER']; ch=h0['CHANNEL']; season=h0['SEASON']
    t=d['TIME']; flux=d['FLUX']; ferr=d['FLUX_ERR']; qual=d['QUALITY']
    collected=np.isfinite(flux).mean(axis=0)>0.5
    ok=np.isfinite(t)&(qual==0)&np.isfinite(flux[:,collected]).all(axis=1)
    t=t[ok]; flux=np.nan_to_num(flux[ok]); ferr=ferr[ok]
    # WCS for the aperture extension (pixel coords)
    w=WCS(f[2].header)
    xt,yt=[float(v) for v in w.all_world2pix(RA_T,DEC_T,0)]; xh,yh=[float(v) for v in w.all_world2pix(RA_H,DEC_H,0)]
    ny,nx=flux.shape[1:]
    info={'quarter':int(q),'channel':int(ch),'season':int(season),'shape':[int(ny),int(nx)],'target_xy':[float(xt),float(yt)],'host_xy':[float(xh),float(yh)],'sep_px':float(np.hypot(xt-xh,yt-yh))}
    # per-event difference images
    n0=int(np.floor((t.min()-T0)/P)); n1=int(np.ceil((t.max()-T0)/P))
    diffs=[]; wts=[]; outs=[]; tcs=[]
    for n in range(n0,n1+1):
        tc=T0+n*P
        inn=np.abs(t-tc)<DUR/2; out=(np.abs(t-tc)<1.0)&(np.abs(t-tc)>1.5*DUR)
        if inn.sum()<2 or out.sum()<10: continue
        # local linear trend per pixel from out-of-transit, evaluated at in-transit times
        X=np.vstack([np.ones(out.sum()),t[out]-tc]).T
        beta=np.linalg.lstsq(X,flux[out].reshape(out.sum(),-1),rcond=None)[0]
        model_in=(np.vstack([np.ones(inn.sum()),t[inn]-tc]).T@beta).reshape(inn.sum(),ny,nx)
        resid=flux[out].reshape(out.sum(),-1)-(X@beta); sd=resid.std(axis=0).reshape(ny,nx)
        diff=(model_in-flux[inn]).mean(axis=0)   # positive = deficit in transit
        var=sd**2/inn.sum()+sd**2/out.sum()
        wt=np.where(collected&(var>0),1/np.where(var>0,var,1),0.0)
        diffs.append(diff); wts.append(wt); outs.append(np.median(flux[out],axis=0)); tcs.append(tc)
    info['n_events']=len(diffs); info['n_collected_pixels']=int(collected.sum())
    R['tpfs'].append(info)
    key=(season,ch)
    stacks.setdefault(key,{'diffs':[],'wts':[],'outs':[],'tcs':[],'quarters':[],'target_xy':(xt,yt),'host_xy':(xh,yh),'shape':(ny,nx)})
    stacks[key]['diffs']+=diffs; stacks[key]['wts']+=wts; stacks[key]['outs']+=outs; stacks[key]['tcs']+=tcs; stacks[key]['quarters'].append(int(q))
    print('Q',q,'ch',ch,'season',season,'shape',(ny,nx),'target',(round(float(xt),2),round(float(yt),2)),'host',(round(float(xh),2),round(float(yh),2)),'events',len(diffs),flush=True)

def analyse(key,excl_postgap):
    S=stacks[key]; ny,nx=S['shape']; xt,yt=S['target_xy']; xh,yh=S['host_xy']
    sel=[i for i,tc in enumerate(S['tcs']) if not (excl_postgap and any(abs(tc-pg)<0.5 for pg in POSTGAP_EVENTS))]
    D=np.array([S['diffs'][i] for i in sel]); W=np.array([S['wts'][i] for i in sel])
    Wsum=W.sum(0); good=Wsum>0
    stack=np.where(good,(D*W).sum(0)/np.where(good,Wsum,1),0.0); err=np.where(good,1/np.sqrt(np.where(good,Wsum,1)),np.inf); snr=np.where(good,stack/err,0.0)
    outimg=np.median(np.array([S['outs'][i] for i in sel]),axis=0)
    # significance-weighted deficit centroid (positive snr only)
    yy,xx=np.mgrid[0:ny,0:nx]; wpos=np.clip(snr,0,None)**2
    cx=(wpos*xx).sum()/wpos.sum(); cy=(wpos*yy).sum()/wpos.sum()
    # template fit: out-of-transit image (host-dominated PRF) shifted to trial position, amplitude by weighted LSQ
    def fit_at(x,y):
        tmpl=ndshift(outimg,(y-yh,x-xh),order=1,mode='constant',cval=0.0); tmpl=np.clip(tmpl,0,None)
        if tmpl.sum()<=0: return np.nan,np.nan,np.nan
        tmpl=tmpl/tmpl.sum()
        wv=W.sum(0); A=(wv*tmpl*stack).sum()/(wv*tmpl**2).sum(); Ae=1/np.sqrt((wv*tmpl**2).sum())
        chi2=(wv*(stack-A*tmpl)**2).sum(); return A,Ae,chi2
    chi2_null=(W.sum(0)*stack**2).sum()
    At,Aet,c2t=fit_at(xt,yt); Ah,Aeh,c2h=fit_at(xh,yh)
    # grid search
    best=(np.inf,None)
    for gx in np.arange(-1.0,nx+0.01,0.25):
        for gy in np.arange(-1.0,ny+0.01,0.25):
            A,Ae,c2=fit_at(gx,gy)
            if np.isfinite(c2) and A>0 and c2<best[0]: best=(c2,(float(gx),float(gy),float(A),float(Ae)))
    # pixel nearest target and host
    it=(int(round(yt)),int(round(xt))); ih=(int(round(yh)),int(round(xh)))
    inside=lambda p: 0<=p[0]<ny and 0<=p[1]<nx
    res={'quarters':S['quarters'],'n_events':len(sel),'target_xy':[float(xt),float(yt)],'host_xy':[float(xh),float(yh)],'target_inside_stamp':bool(inside(it)),
         'deficit_centroid_xy':[float(cx),float(cy)],'centroid_dist_to_target_px':float(np.hypot(cx-xt,cy-yt)),'centroid_dist_to_host_px':float(np.hypot(cx-xh,cy-yh)),
         'snr_at_target_pixel':float(snr[it]) if inside(it) else None,'snr_at_host_pixel':float(snr[ih]) if inside(ih) else None,
         'deficit_at_target_pixel_e_s':float(stack[it]) if inside(it) else None,'deficit_err_target_pixel_e_s':float(err[it]) if inside(it) else None,
         'max_snr_pixel_yx':[int(np.unravel_index(np.argmax(snr),snr.shape)[0]),int(np.unravel_index(np.argmax(snr),snr.shape)[1])],'max_snr':float(snr.max()),'min_snr':float(snr.min()),
         'chi2_null':float(chi2_null),'fit_at_target':{'A_e_s':float(At),'A_err':float(Aet),'chi2':float(c2t)},'fit_at_host':{'A_e_s':float(Ah),'A_err':float(Aeh),'chi2':float(c2h)},
         'best_grid':{'chi2':float(best[0]),'xy_A_Aerr':best[1]},'ndata':int(ny*nx),
         'out_image_total_e_s':float(outimg.sum()),'snr_map':snr.round(2).tolist(),'stack_e_s':stack.round(3).tolist()}
    return res,stack,snr,outimg
fig,ax=plt.subplots(2,len(stacks),figsize=(4*len(stacks),8))
for j,key in enumerate(sorted(stacks)):
    for i,ex in enumerate([False,True]):
        res,stack,snr,outimg=analyse(key,ex)
        R['seasons'][f'season{key[0]}_ch{key[1]}'+('_noPostGap' if ex else '_all')]=res
        a=ax[i,j]; im=a.imshow(snr,origin='lower',cmap='RdBu_r',vmin=-max(4,abs(snr).max()),vmax=max(4,abs(snr).max()))
        a.plot(res['target_xy'][0],res['target_xy'][1],'g*',ms=14,label='target'); a.plot(res['host_xy'][0],res['host_xy'][1],'k+',ms=14,mew=2,label='host')
        a.plot(res['deficit_centroid_xy'][0],res['deficit_centroid_xy'][1],'mo',ms=8,label='deficit centroid')
        a.set_title(f"S{key[0]} ch{key[1]} {'no post-gap' if ex else 'all'} n={res['n_events']}\nSNR@target {res['snr_at_target_pixel']} max {res['max_snr']:.1f}",fontsize=9)
        plt.colorbar(im,ax=a,fraction=0.04)
        if i==0 and j==0: a.legend(fontsize=7)
plt.tight_layout(); plt.savefig(OUT+'pixel_checks.png',dpi=110)
json.dump(R,open(OUT+'pixel_checks.json','w'),indent=1,default=float)
for k,v in R['seasons'].items():
    print(k,{kk:v[kk] for kk in ['n_events','target_inside_stamp','centroid_dist_to_target_px','centroid_dist_to_host_px','snr_at_target_pixel','snr_at_host_pixel','max_snr','max_snr_pixel_yx','chi2_null','fit_at_target','fit_at_host','best_grid']})
