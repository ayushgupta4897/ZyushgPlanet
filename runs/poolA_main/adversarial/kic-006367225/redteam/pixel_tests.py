"""Method controls: (A) phase-shifted null at +0.5 P; (B) injected 2.8 e/s Gaussian deficit at the target's WCS-corrected position at a different epoch; (C) same injection at the G20.5 position."""
import numpy as np, json, glob, warnings, sys
warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
P=12.971527521795288; T0=142.75300012519583; HW=1.125/24; OOT=(2.5/24,8.0/24)
SRC={'target':(292.86005384301,41.76715441206),'G20.5':(292.85968261898,41.76810611087),'bright':(292.85822280656,41.76535037522)}
files=sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr006367220_lc_Q011111111111111111/*targ.fits.gz'))
def run(t0, inject_at=None, inject_total_e_s=2.8, tag=''):
    rows=[]
    for fn in files:
        h=fits.open(fn); hd=h[0].header; q=hd['QUARTER']; d=h[1].data
        t=d['TIME']; fl=d['FLUX'].astype(float); qual=d['QUALITY']
        ok=np.isfinite(t)&(qual==0)&np.isfinite(fl).any(axis=(1,2)); t=t[ok]; fl=fl[ok].copy()
        w=WCS(h[2].header); ny,nx=fl.shape[1:]; Y,X=np.mgrid[0:ny,0:nx].astype(float)
        pos={k:tuple(map(float,w.all_world2pix([[ra,dec]],0)[0])) for k,(ra,dec) in SRC.items()}
        valid=np.isfinite(fl).all(axis=0)&(np.nanmedian(fl,axis=0)>0)
        M=np.nanmean(fl,axis=0); M[~valid]=np.nan
        bx,by=pos['bright']; m=np.isfinite(M)
        r=least_squares(lambda p:((p[0]*np.exp(-((X-p[1])**2+(Y-p[2])**2)/(2*p[3]**2))+p[4]-M)/(np.sqrt(M)+1))[m],[np.nanmax(M),bx,by,0.7,0])
        off=(r.x[1]-bx,r.x[2]-by); sig=abs(r.x[3]); cpos={k:(v[0]+off[0],v[1]+off[1]) for k,v in pos.items()}
        if inject_at:
            cx,cy=cpos[inject_at]; g=np.exp(-((X-cx)**2+(Y-cy)**2)/(2*sig*sig)); g*=inject_total_e_s/g.sum()
            ph=((t-t0)/P+0.5)%1-0.5; intr=np.abs(ph)*P<HW
            fl[intr]-=g[None,:,:]
        n0=int(np.ceil((t.min()-t0)/P)); n1=int(np.floor((t.max()-t0)/P)); D=[];V=[]
        for n in range(n0,n1+1):
            tc=t0+n*P; dt=np.abs(t-tc); inn=dt<HW; oo=(dt>OOT[0])&(dt<OOT[1])
            if inn.sum()<3 or oo.sum()<12: continue
            mo=np.nanmean(fl[oo],axis=0); so=np.nanstd(fl[oo],axis=0,ddof=1); mi=np.nanmean(fl[inn],axis=0)
            D.append(mo-mi); V.append(so**2*(1/inn.sum()+1/oo.sum()))
        if len(D)<2: continue
        D=np.array(D); V=np.array(V); W=1/V; Dm=(D*W).sum(0)/W.sum(0); Em=W.sum(0)**-0.5; Dm[~valid]=np.nan; Em[~valid]=np.nan
        m=np.isfinite(Dm); i0=np.nanargmax(Dm/Em); y0,x0=np.unravel_index(i0,Dm.shape)
        res=lambda qq:((qq[0]*np.exp(-((X-qq[1])**2+(Y-qq[2])**2)/(2*sig*sig))+qq[3]-Dm)/Em)[m]
        rf=least_squares(res,[np.nanmax(Dm),x0,y0,0]); c2f=float((rf.fun**2).sum())
        c2fix={}
        for k,(cx,cy) in cpos.items():
            g=np.exp(-((X-cx)**2+(Y-cy)**2)/(2*sig*sig)); A=np.vstack([g[m],np.ones(m.sum())]).T/Em[m][:,None]; b=Dm[m]/Em[m]
            sol=np.linalg.lstsq(A,b,rcond=None)[0]; c2fix[k]=float(((A@sol-b)**2).sum())
        rng=np.random.default_rng(2); bs=[]
        for _ in range(100):
            idx=rng.integers(0,len(D),len(D)); Db=(D[idx]*W[idx]).sum(0)/W[idx].sum(0); Eb=W[idx].sum(0)**-0.5; Db[~valid]=np.nan
            rb=least_squares(lambda qq:((qq[0]*np.exp(-((X-qq[1])**2+(Y-qq[2])**2)/(2*sig*sig))+qq[3]-Db)/Eb)[m],[rf.x[0],rf.x[1],rf.x[2],0]); bs.append(rb.x[1:3])
        bss=np.std(np.array(bs),0)
        sky=w.all_pix2world([[rf.x[1]-off[0],rf.x[2]-off[1]]],0)[0]
        rows.append({'q':int(q),'pkSNR':float(np.nanmax(Dm/Em)),'chi2_null':float(np.nansum((Dm/Em)**2)),'chi2_free':c2f,'chi2_fixed':c2fix,'sky':[float(sky[0]),float(sky[1])],'err_px':[float(bss[0]),float(bss[1])],'dist_px':{k:float(np.hypot(rf.x[1]-v[0],rf.x[2]-v[1])) for k,v in cpos.items()}})
    tr,td=SRC['target']; gr,gd=SRC['G20.5']; gx=(gr-tr)*np.cos(np.radians(td))*3600; gy=(gd-td)*3600
    use=[r for r in rows if r['pkSNR']>3]
    if use:
        dx=np.array([(r['sky'][0]-tr)*np.cos(np.radians(td))*3600 for r in use]); dy=np.array([(r['sky'][1]-td)*3600 for r in use]); e=np.array([np.hypot(*r['err_px'])*3.98/np.sqrt(2) for r in use]); wg=1/e**2
        mx=(dx*wg).sum()/wg.sum(); my=(dy*wg).sum()/wg.sum(); em=wg.sum()**-0.5
        summ={'n_quarters_snr3':len(use),'offset_from_target_arcsec':[float(mx),float(my)],'err':float(em),'dist_target':float(np.hypot(mx,my)),'dist_G20.5':float(np.hypot(mx-gx,my-gy)),
              'sum_chi2_fixed':{k:sum(r['chi2_fixed'][k] for r in use) for k in SRC},'sum_chi2_null':sum(r['chi2_null'] for r in use),'sum_chi2_free':sum(r['chi2_free'] for r in use),'median_pkSNR':float(np.median([r['pkSNR'] for r in rows]))}
    else: summ={'n_quarters_snr3':0,'median_pkSNR':float(np.median([r['pkSNR'] for r in rows])),'max_pkSNR':float(max(r['pkSNR'] for r in rows))}
    print(f'[{tag}]',json.dumps(summ))
    return {'summary':summ,'quarters':rows}
out={}
out['null_halfP']=run(T0+0.5*P,tag='NULL phase+0.5')
out['null_quarterP']=run(T0+0.25*P,tag='NULL phase+0.25')
out['inject_target']=run(T0+0.37*P,inject_at='target',tag='INJECT 2.8 e/s at target, epoch+0.37P')
out['inject_G20.5']=run(T0+0.37*P,inject_at='G20.5',tag='INJECT 2.8 e/s at G20.5, epoch+0.37P')
json.dump(out,open('pixel_tests.json','w'),indent=1)
