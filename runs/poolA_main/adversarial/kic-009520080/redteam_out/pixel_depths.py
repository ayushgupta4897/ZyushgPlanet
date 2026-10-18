# Independent per-pixel depth test on the originating TPF (KIC 9520088), no pipeline imports.
# For each quarter: per-pixel light curve -> per-quarter normalise -> transit-masked 1-d running median -> box depth at ephemeris.
# Report depth in e/s (missing light) and in ppm of that pixel's flux. Compare target-edge pixel vs bright-star core pixels.
import numpy as np, glob, json
from astropy.io import fits
from astropy.wcs import WCS
P=2.377014357204058; t0=132.15453895516694; dur=1.125/24
tgt=(289.01257773399,46.16274216212); host=(289.01653098528,46.1636484614); g139=(289.01697392124,46.16711309667); g209=(289.0123224849,46.16369100903)
files=sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr009520088_lc_Q011111111111111111/*lpd-targ.fits.gz'))
res=[]
for f in files:
    hd=fits.open(f); h0=hd[0].header; d=hd[1].data; ap=hd[2].data
    qn=h0['QUARTER']; t=d['TIME']; fl=d['FLUX']; qual=d['QUALITY']
    ok=np.isfinite(t)&(qual==0)&np.isfinite(fl).all(axis=(1,2))
    t=t[ok]; fl=fl[ok]
    w=WCS(hd[2].header)
    def px(radec):
        x,y=w.all_world2pix(radec[0],radec[1],0); return float(x),float(y)
    tp,hp,g9,gf=px(tgt),px(host),px(g139),px(g209)
    ph=((t-t0)/P+0.5)%1-0.5; intr=np.abs(ph*P)<dur/2; near=np.abs(ph*P)<0.75*dur
    ny,nx=fl.shape[1:]
    dep=np.full((ny,nx),np.nan); dep_ppm=np.full((ny,nx),np.nan); err=np.full((ny,nx),np.nan); med=np.full((ny,nx),np.nan)
    # running median trend per pixel (1-d window) excluding near-transit points
    order=np.argsort(t)
    for iy in range(ny):
        for ix in range(nx):
            y=fl[:,iy,ix].astype(float)
            if not np.all(np.isfinite(y)) or np.nanmedian(y)<=0: continue
            tr=np.empty_like(y)
            for i in range(len(t)):
                s=(np.abs(t-t[i])<0.5)&(~near); tr[i]=np.median(y[s]) if s.sum()>5 else np.nan
            r=y-tr
            good=np.isfinite(r); 
            sig=1.4826*np.median(np.abs(r[good]-np.median(r[good])))
            cl=good&(np.abs(r)<5*sig)
            din=r[cl&intr]; 
            if len(din)<4: continue
            dep[iy,ix]=-din.mean(); err[iy,ix]=np.std(r[cl&~near])/np.sqrt(len(din)); med[iy,ix]=np.median(y); dep_ppm[iy,ix]=dep[iy,ix]/med[iy,ix]*1e6
    snr=dep/err
    tot=np.nansum(dep[np.isfinite(dep)])
    # pixel nearest target, nearest host
    def nearest(p):
        ix,iy=int(round(p[0])),int(round(p[1])); ix=min(max(ix,0),nx-1); iy=min(max(iy,0),ny-1); return iy,ix
    ty,tx=nearest(tp); hy,hx=nearest(hp)
    r=dict(quarter=int(qn),channel=int(h0['CHANNEL']),shape=[ny,nx],target_px=[round(tp[0],2),round(tp[1],2)],host_px=[round(hp[0],2),round(hp[1],2)],g139_px=[round(g9[0],2),round(g9[1],2)],g209_px=[round(gf[0],2),round(gf[1],2)],
           n_cad=int(len(t)),n_in=int(intr.sum()),
           target_pixel=dict(iy=ty,ix=tx,flux_e_s=float(med[ty,tx]),deficit_e_s=float(dep[ty,tx]),err_e_s=float(err[ty,tx]),ppm=float(dep_ppm[ty,tx]),snr=float(snr[ty,tx])),
           host_pixel=dict(iy=hy,ix=hx,flux_e_s=float(med[hy,hx]),deficit_e_s=float(dep[hy,hx]),err_e_s=float(err[hy,hx]),ppm=float(dep_ppm[hy,hx]),snr=float(snr[hy,hx])),
           deficit_map_e_s=np.round(dep,2).tolist(),snr_map=np.round(snr,1).tolist(),ppm_map=np.round(dep_ppm,0).tolist(),flux_map=np.round(med,0).tolist(),
           total_deficit_e_s=float(tot))
    # brightest-5 pixels (host-dominated) weighted ppm
    idx=np.argsort(np.nan_to_num(med).ravel())[::-1][:5]
    r['bright5_ppm_weighted']=float(np.nansum(dep.ravel()[idx])/np.nansum(med.ravel()[idx])*1e6)
    r['bright5_deficit_e_s']=float(np.nansum(dep.ravel()[idx])); r['bright5_err_e_s']=float(np.sqrt(np.nansum(err.ravel()[idx]**2)))
    res.append(r)
    print(f"Q{qn:2d} ch{r['channel']} tgt_px={r['target_px']} host_px={r['host_px']} | target px flux {med[ty,tx]:.0f} e/s deficit {dep[ty,tx]:+.2f}+-{err[ty,tx]:.2f} ({dep_ppm[ty,tx]:+.0f} ppm, snr {snr[ty,tx]:.1f}) | host px flux {med[hy,hx]:.0f} deficit {dep[hy,hx]:+.2f}+-{err[hy,hx]:.2f} ({dep_ppm[hy,hx]:+.0f} ppm, snr {snr[hy,hx]:.1f}) | bright5 {r['bright5_deficit_e_s']:+.2f}+-{r['bright5_err_e_s']:.2f} e/s = {r['bright5_ppm_weighted']:+.0f} ppm")
    print('   SNR map (rows=y):'); 
    for row in r['snr_map']: print('   ',row)
json.dump(res,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic009520080/out/pixel_depths.json','w'),indent=1)
# summary
tsn=np.array([r['target_pixel']['snr'] for r in res]); hsn=np.array([r['host_pixel']['snr'] for r in res])
td=np.array([r['target_pixel']['deficit_e_s'] for r in res]); te=np.array([r['target_pixel']['err_e_s'] for r in res])
hd_=np.array([r['bright5_deficit_e_s'] for r in res]); he=np.array([r['bright5_err_e_s'] for r in res])
print('SUMMARY: target pixel combined deficit %.2f +- %.2f e/s (snr %.1f); bright-5 (host) pixels combined %.2f +- %.2f e/s (snr %.1f); bright-5 weighted ppm median %.0f'%(
  (td/te**2).sum()/(1/te**2).sum(), (1/te**2).sum()**-0.5, ((td/te**2).sum()/(1/te**2).sum())/((1/te**2).sum()**-0.5),
  (hd_/he**2).sum()/(1/he**2).sum(), (1/he**2).sum()**-0.5, ((hd_/he**2).sum()/(1/he**2).sum())/((1/he**2).sum()**-0.5), np.median([r['bright5_ppm_weighted'] for r in res])))
