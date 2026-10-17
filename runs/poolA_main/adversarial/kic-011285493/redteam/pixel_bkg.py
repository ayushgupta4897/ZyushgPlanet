import warnings; warnings.filterwarnings('ignore')
import glob, numpy as np, json
from astropy.io import fits
from astropy.wcs import WCS
ra_t,dec_t=285.1259332116789,49.016594587983434; ra_n,dec_n=285.12944145127454,49.01770673400924
P=10.790444; t0=142.33451; half_in=0.12; oot=(0.35,1.2)
rows=[]; allpix=[]
for f in sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr011285496_lc_*/*targ.fits.gz')):
    h=fits.open(f); d=h[1].data; h0=h[0].header; w=WCS(h[2].header)
    t=d['TIME']; fl=d['FLUX']; bk=d['FLUX_BKG']; q=d['QUALITY']
    m=np.isfinite(t)&(q==0)&np.isfinite(np.nansum(fl,axis=(1,2)))
    t,fl,bk=t[m],fl[m],bk[m]
    xt,yt=w.all_world2pix(ra_t,dec_t,0); xn,yn=w.all_world2pix(ra_n,dec_n,0)
    ny,nx=fl.shape[1:]; yy,xx=np.mgrid[0:ny,0:nx]
    rt=np.hypot(xx-xt,yy-yt); rn=np.hypot(xx-xn,yy-yn)
    valid=np.isfinite(fl).all(axis=0)&(np.nanmedian(fl,axis=0)>0)
    ph=((t-t0+0.5*P)%P)-0.5*P; ep=np.round((t-t0)/P).astype(int)
    defs=[]; bdefs=[]
    for e in np.unique(ep):
        s=(ep==e)&(np.abs(ph)<half_in); o=(ep==e)&(np.abs(ph)>oot[0])&(np.abs(ph)<oot[1])
        if s.sum()<4 or o.sum()<20: continue
        # per-pixel local linear baseline from out-of-transit
        tt=t[o]; A=np.c_[np.ones(o.sum()),tt-tt.mean()]
        coef=np.linalg.lstsq(A,fl[o].reshape(o.sum(),-1),rcond=None)[0]
        base=(np.c_[np.ones(s.sum()),t[s]-tt.mean()]@coef).reshape(s.sum(),ny,nx)
        defs.append(np.nanmean(base-fl[s],axis=0))
        cb=np.linalg.lstsq(A,bk[o].reshape(o.sum(),-1),rcond=None)[0]
        bb=(np.c_[np.ones(s.sum()),t[s]-tt.mean()]@cb).reshape(s.sum(),ny,nx)
        bdefs.append(np.nanmean(bb-bk[s],axis=0))
    defs=np.array(defs); bdefs=np.array(bdefs); n=len(defs)
    D=np.nanmean(defs,axis=0); E=np.nanstd(defs,axis=0)/np.sqrt(n)
    Bd=np.nanmean(bdefs,axis=0)
    img=np.nanmedian(fl,axis=0)
    tgt=valid&(rt<1.5); nbr=valid&(rn<1.5)&~tgt; far=valid&(rt>2.5)&(rn>2.5)
    def agg(mask): 
        return float(np.sum(D[mask])), float(np.sqrt(np.sum(E[mask]**2))), int(mask.sum()), float(np.mean(D[mask]))
    tot=agg(valid); tg=agg(tgt); nb=agg(nbr); fr=agg(far)
    # linear fit deficit = a + b*img over valid pixels
    A=np.c_[np.ones(valid.sum()),img[valid]]; W=1/E[valid]**2
    cov=np.linalg.inv(A.T@(A*W[:,None])); ab=cov@(A.T@(W*D[valid]))
    print(f"Q{h0['QUARTER']:2d} ch{h0['CHANNEL']:3d} nev={n:2d} valid={valid.sum()} | total deficit {tot[0]:6.2f}+-{tot[1]:.2f} e/s | target r<1.5px ({tg[2]}px): {tg[0]:5.2f}+-{tg[1]:.2f} | nbr core ({nb[2]}px): {nb[0]:5.2f}+-{nb[1]:.2f} | far ({fr[2]}px): {fr[0]:5.2f}+-{fr[1]:.2f}, per px {fr[3]:.3f} | fit a={ab[0]:.3f}+-{np.sqrt(cov[0,0]):.3f} e/s/px, b={ab[1]*1e6:.1f}+-{np.sqrt(cov[1,1])*1e6:.1f} ppm | FLUX_BKG deficit/px {np.nanmean(Bd[valid]):.3f}")
    rows.append(dict(Q=int(h0['QUARTER']),ch=int(h0['CHANNEL']),nev=n,total=tot,target=tg,nbr=nb,far=fr,a=float(ab[0]),a_err=float(np.sqrt(cov[0,0])),b_ppm=float(ab[1]*1e6),b_err=float(np.sqrt(cov[1,1])*1e6),bkg_def_px=float(np.nanmean(Bd[valid]))))
    allpix.append((D[valid],E[valid],img[valid],rt[valid],rn[valid]))
json.dump(rows,open('pixel_bkg.json','w'),indent=1)
# combined
D=np.concatenate([a[0] for a in allpix]); E=np.concatenate([a[1] for a in allpix]); I=np.concatenate([a[2] for a in allpix]); RT=np.concatenate([a[3] for a in allpix]); RN=np.concatenate([a[4] for a in allpix])
W=1/E**2
for lab,m in [('all',np.ones(len(D),bool)),('target r<1.5',RT<1.5),('nbr r<1.5',RN<1.5),('far >2.5 both',(RT>2.5)&(RN>2.5)),('far >3.5 both',(RT>3.5)&(RN>3.5))]:
    wm=np.sum(W[m]*D[m])/np.sum(W[m]); we=1/np.sqrt(np.sum(W[m]))
    print(f'ALL QUARTERS {lab:15s}: weighted mean deficit per pixel {wm:.4f} +- {we:.4f} e/s  (n px-quarters {m.sum()})')
A=np.c_[np.ones(len(D)),I]; cov=np.linalg.inv(A.T@(A*W[:,None])); ab=cov@(A.T@(W*D))
print(f'ALL: deficit = a + b*flux: a = {ab[0]:.4f} +- {np.sqrt(cov[0,0]):.4f} e/s/px ; b = {ab[1]*1e6:.1f} +- {np.sqrt(cov[1,1])*1e6:.1f} ppm')
# model comparison: uniform offset vs gaussian at target (sigma 0.9 px) vs scaled image
G=np.exp(-0.5*(RT/0.9)**2)
def fit(cols):
    A=np.array(cols).T; cov=np.linalg.inv(A.T@(A*W[:,None])); p=cov@(A.T@(W*D)); chi=np.sum(W*(D-A@p)**2); return p,np.sqrt(np.diag(cov)),chi
chi0=np.sum(W*D**2)
for lab,cols in [('uniform',[np.ones(len(D))]),('gauss@target',[G]),('scaled image',[I]),('uniform+gauss',[np.ones(len(D)),G]),('uniform+image',[np.ones(len(D)),I])]:
    p,e,chi=fit(cols); print(f'{lab:15s} chi2={chi:8.1f} (null {chi0:.1f}) params',np.round(p,4),'+-',np.round(e,4))
