"""Own difference-image localisation for kic-006367225 on the KIC 6367220 TPFs (read as data only)."""
import numpy as np, json, glob, sys, warnings
warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
P=12.971527521795288; T0=142.75300012519583; HW=1.125/24  # in-transit half-width (h->d): box 2.25 h
OOT=(2.5/24, 8.0/24)
SRC={'target':(292.86005384301,41.76715441206),'G20.5':(292.85968261898,41.76810611087),'bright':(292.85822280656,41.76535037522),'G20.9':(292.85724316874,41.76548057458)}
files=sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr006367220_lc_Q011111111111111111/*targ.fits.gz'))
def gauss(p,X,Y):
    a,x0,y0,s,c=p; return a*np.exp(-((X-x0)**2+(Y-y0)**2)/(2*s*s))+c
def fitg(img,err,X,Y,p0,fix=None):
    m=np.isfinite(img)&np.isfinite(err)&(err>0)
    def res(q):
        p=list(q) if fix is None else [q[0],fix[0],fix[1],q[1],q[2]]
        return ((gauss(p,X,Y)-img)/err)[m]
    q0=p0 if fix is None else [p0[0],p0[3],p0[4]]
    r=least_squares(res,q0,method='lm' if m.sum()>len(q0) else 'trf')
    p=list(r.x) if fix is None else [r.x[0],fix[0],fix[1],r.x[1],r.x[2]]
    return p,float((r.fun**2).sum()),int(m.sum())
out={'quarters':[]}; figs=[]
for fn in files:
    h=fits.open(fn); hd=h[0].header; q=hd['QUARTER']; d=h[1].data; ap=h[2]
    t=d['TIME']; fl=d['FLUX']; fe=d['FLUX_ERR']; qual=d['QUALITY']
    ok=np.isfinite(t)&(qual==0)&np.isfinite(fl).any(axis=(1,2)); t=t[ok]; fl=fl[ok]; fe=fe[ok]
    w=WCS(ap.header); ny,nx=fl.shape[1:]
    Y,X=np.mgrid[0:ny,0:nx].astype(float)
    pos={k:tuple(map(float,w.all_world2pix([[ra,dec]],0)[0])) for k,(ra,dec) in SRC.items()}  # (col,row) 0-based
    valid=np.isfinite(fl).all(axis=0)&(np.nanmedian(fl,axis=0)>0)
    # step detector on summed flux
    S=np.nansum(np.where(valid,fl,0),axis=(1,2)); steps=[]
    for i in range(6,len(S)-6):
        a,b=np.median(S[i-6:i]),np.median(S[i:i+6]); sig=1.4826*np.median(np.abs(np.diff(S)))/np.sqrt(2)
        if abs(b-a)>10*sig: steps.append(t[i])
    steps=np.array(steps)
    # events
    n0=int(np.ceil((t.min()-T0)/P)); n1=int(np.floor((t.max()-T0)/P))
    D=[];V=[];used=[];skipped=[]
    for n in range(n0,n1+1):
        tc=T0+n*P; dt=np.abs(t-tc); inn=dt<HW; oo=(dt>OOT[0])&(dt<OOT[1])
        if inn.sum()<3 or oo.sum()<12: skipped.append((n,'few')); continue
        if steps.size and np.min(np.abs(steps-tc))<1.0: skipped.append((n,'step')); continue
        mo=np.nanmean(fl[oo],axis=0); so=np.nanstd(fl[oo],axis=0,ddof=1); mi=np.nanmean(fl[inn],axis=0)
        D.append(mo-mi); V.append(so**2*(1/inn.sum()+1/oo.sum())); used.append(n)
    if len(D)<2: out['quarters'].append({'quarter':int(q),'n_events':len(D),'skipped':skipped}); continue
    D=np.array(D); V=np.array(V); W=1/V
    Dm=(D*W).sum(0)/W.sum(0); Em=W.sum(0)**-0.5; Dm[~valid]=np.nan; Em[~valid]=np.nan
    M=np.nanmean(fl,axis=0); M[~valid]=np.nan; Me=np.nanstd(fl,axis=0)/np.sqrt(len(fl)); Me[~valid]=np.nan
    # WCS calibration on bright star: fit gaussian to mean image
    bx,by=pos['bright']; pM,c2M,_=fitg(M,np.sqrt(M)+1,X,Y,[np.nanmax(M),bx,by,0.7,0.0])
    off=(pM[1]-bx,pM[2]-by); sig=abs(pM[3])
    cpos={k:(v[0]+off[0],v[1]+off[1]) for k,v in pos.items()}
    # deficit: free fit
    i0=np.nanargmax(Dm/Em); y0,x0=np.unravel_index(i0,Dm.shape)
    pD,c2D,npx=fitg(Dm,Em,X,Y,[np.nanmax(Dm),x0,y0,sig,0.0])
    # refit with sigma fixed to stellar sigma (free pos)
    def fit_fixed_sigma(img,err,x0,y0):
        m=np.isfinite(img)
        def res(qq): return ((qq[0]*np.exp(-((X-qq[1])**2+(Y-qq[2])**2)/(2*sig*sig))+qq[3]-img)/err)[m]
        r=least_squares(res,[np.nanmax(img),x0,y0,0.0]); return list(r.x),float((r.fun**2).sum())
    pF,c2F=fit_fixed_sigma(Dm,Em,x0,y0)
    # fixed-position fits (amp, const only)
    c2fix={}
    for k,(cx,cy) in cpos.items():
        m=np.isfinite(Dm); g=np.exp(-((X-cx)**2+(Y-cy)**2)/(2*sig*sig))
        A=np.vstack([g[m],np.ones(m.sum())]).T/Em[m][:,None]; b=Dm[m]/Em[m]
        sol,resid,_,_=np.linalg.lstsq(A,b,rcond=None); c2fix[k]=(float(((A@sol-b)**2).sum()),float(sol[0]))
    c2null=float(np.nansum((Dm/Em)**2))
    # bootstrap centroid over events
    rng=np.random.default_rng(1); bs=[]
    for _ in range(150):
        idx=rng.integers(0,len(D),len(D)); Db=(D[idx]*W[idx]).sum(0)/W[idx].sum(0); Eb=W[idx].sum(0)**-0.5; Db[~valid]=np.nan
        pb,_=fit_fixed_sigma(Db,Eb,pF[1],pF[2]); bs.append(pb[1:3])
    bs=np.array(bs); bsm=np.median(bs,0); bss=np.std(bs,0)
    # sky position of deficit (undo WCS offset)
    sky=w.all_pix2world([[pF[1]-off[0],pF[2]-off[1]]],0)[0]
    def dist(a,b): return float(np.hypot(a[0]-b[0],a[1]-b[1]))
    rec={'quarter':int(q),'module':int(hd['MODULE']),'channel':int(hd['CHANNEL']),'n_events':len(D),'events':used,'skipped':skipped,'n_steps':int(steps.size),
         'stamp':[int(ny),int(nx)],'wcs_offset_px':[float(off[0]),float(off[1])],'star_sigma_px':float(sig),
         'pos_corrected_colrow':{k:[float(a),float(b)] for k,(a,b) in cpos.items()},
         'deficit_fit_colrow':[float(pF[1]),float(pF[2])],'deficit_boot_std_px':[float(bss[0]),float(bss[1])],'deficit_amp_e_s':float(pF[0]),
         'peak_snr':float(np.nanmax(Dm/Em)),'chi2_null':c2null,'chi2_free':c2F,'chi2_fixed':{k:v[0] for k,v in c2fix.items()},'amp_fixed':{k:v[1] for k,v in c2fix.items()},
         'dist_px':{k:dist((pF[1],pF[2]),v) for k,v in cpos.items()},'deficit_sky':[float(sky[0]),float(sky[1])],'npix':npx}
    out['quarters'].append(rec); figs.append((q,M,Dm/Em,cpos,(pF[1],pF[2]),bss))
    print(f"Q{q:2d} mod{rec['module']:2d} nev={len(D)} skip={len(skipped)} steps={steps.size} off=({off[0]:+.2f},{off[1]:+.2f}) sig={sig:.2f} pkSNR={rec['peak_snr']:.1f} chi2 null/free {c2null:.0f}/{c2F:.0f} fixed: "+' '.join(f"{k}={v[0]:.0f}" for k,v in c2fix.items())+f"  dist(px): tgt={rec['dist_px']['target']:.2f} G20.5={rec['dist_px']['G20.5']:.2f} bright={rec['dist_px']['bright']:.2f}  boot std=({bss[0]:.2f},{bss[1]:.2f})")
# combine in sky coords: offset of deficit from target in arcsec, per quarter, weighted by bootstrap errors (px->arcsec 3.98)
rows=[r for r in out['quarters'] if 'deficit_sky' in r and r['peak_snr']>3]
tr,td=SRC['target']; gr,gd=SRC['G20.5']
def dsky(r): return ((r['deficit_sky'][0]-tr)*np.cos(np.radians(td))*3600,(r['deficit_sky'][1]-td)*3600)
dx=np.array([dsky(r)[0] for r in rows]); dy=np.array([dsky(r)[1] for r in rows]); e=np.array([np.hypot(*r['deficit_boot_std_px'])*3.98/np.sqrt(2) for r in rows])
wgt=1/e**2; mx=(dx*wgt).sum()/wgt.sum(); my=(dy*wgt).sum()/wgt.sum(); em=wgt.sum()**-0.5
gx=(gr-tr)*np.cos(np.radians(td))*3600; gy=(gd-td)*3600
print(f"\nCombined (quarters with peak SNR>3, n={len(rows)}): deficit offset from target = ({mx:+.2f}, {my:+.2f}) +- {em:.2f} arcsec; G20.5 star is at ({gx:+.2f},{gy:+.2f}) arcsec from target")
print(f"  distance to target {np.hypot(mx,my):.2f}\" ({np.hypot(mx,my)/em:.1f} sigma); distance to G20.5 {np.hypot(mx-gx,my-gy):.2f}\" ({np.hypot(mx-gx,my-gy)/em:.1f} sigma)")
chi_t=((dx/e)**2+(dy/e)**2).sum(); chi_g=(((dx-gx)/e)**2+((dy-gy)/e)**2).sum()
print(f"  chi2 (deficit at target) = {chi_t:.1f}, chi2 (deficit at G20.5) = {chi_g:.1f}, for {2*len(rows)} dof")
print("  per-quarter offsets (arcsec):",[(r['quarter'],round(dsky(r)[0],2),round(dsky(r)[1],2),round(ee,2)) for r,ee in zip(rows,e)])
S={k:sum(r['chi2_fixed'][k] for r in rows) for k in SRC}; print('  summed fixed-position chi2 over those quarters:',{k:round(v,1) for k,v in S.items()},' null:',round(sum(r['chi2_null'] for r in rows),1))
out['combined']={'n_quarters':len(rows),'offset_from_target_arcsec':[mx,my],'err_arcsec':em,'G20.5_offset_arcsec':[gx,gy],'dist_target':float(np.hypot(mx,my)),'dist_G20.5':float(np.hypot(mx-gx,my-gy)),'chi2_at_target':float(chi_t),'chi2_at_G20.5':float(chi_g),'sum_chi2_fixed':S}
json.dump(out,open('pixel_own.json','w'),indent=1,default=float)
# figure
n=len(figs); fig,axs=plt.subplots(2,n,figsize=(2.6*n,5.6))
for i,(q,M,S_,cpos,pf,bss) in enumerate(figs):
    axs[0,i].imshow(M,origin='lower',cmap='gray'); axs[0,i].set_title(f'Q{q} mean',fontsize=8)
    axs[1,i].imshow(S_,origin='lower',cmap='RdBu_r',vmin=-6,vmax=6); axs[1,i].set_title('deficit SNR',fontsize=8)
    for ax in axs[:,i]:
        ax.plot(*cpos['target'],'o',mfc='none',mec='lime',ms=10); ax.plot(*cpos['G20.5'],'s',mfc='none',mec='orange',ms=9); ax.plot(*cpos['bright'],'x',color='cyan',ms=6)
        ax.errorbar(pf[0],pf[1],xerr=bss[0],yerr=bss[1],fmt='+',color='magenta',ms=10)
        ax.set_xticks([]); ax.set_yticks([])
fig.suptitle('kic-006367225 own difference images: lime circle=target, orange square=G20.5 star, cyan x=KIC 6367220, magenta=fitted deficit (WCS calibrated on bright star)',fontsize=8)
plt.tight_layout(); plt.savefig('pixel_own.png',dpi=110)
