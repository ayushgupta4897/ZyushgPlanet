import warnings; warnings.filterwarnings('ignore')
import glob, numpy as np, json
from astropy.io import fits
from astropy.wcs import WCS
ra_t,dec_t=285.1259332116789,49.016594587983434; ra_n,dec_n=285.12944145127454,49.01770673400924
P=10.790444; T0=142.33451; half_in=0.12
res={}
for lab,shift in [('primary',0.0),('secondary',P/2)]:
    per_q=[]
    for f in sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr011285496_lc_*/*targ.fits.gz')):
        h=fits.open(f); d=h[1].data; h0=h[0].header; w=WCS(h[2].header)
        t=d['TIME']; fl=d['FLUX']; bk=d['FLUX_BKG']; q=d['QUALITY']
        m=np.isfinite(t)&(q==0)&np.isfinite(np.nansum(fl,axis=(1,2))); t,fl,bk=t[m],fl[m],bk[m]
        xt,yt=w.all_world2pix(ra_t,dec_t,0); xn,yn=w.all_world2pix(ra_n,dec_n,0)
        ny,nx=fl.shape[1:]; yy,xx=np.mgrid[0:ny,0:nx]; rt=np.hypot(xx-xt,yy-yt); rn=np.hypot(xx-xn,yy-yn)
        valid=np.isfinite(fl).all(axis=0)&(np.nanmedian(fl,axis=0)>0)
        aps={'target':valid&(rt<1.5),'far':valid&(rt>2.5)&(rn>2.5),'nbr':valid&(rn<1.5)&~(rt<1.5),'all':valid}
        ph=((t-T0-shift+0.5*P)%P)-0.5*P; ep=np.round((t-T0-shift)/P).astype(int)
        ev={k:[] for k in aps}; ev['bkg']=[]
        for e in np.unique(ep):
            s=(ep==e)&(np.abs(ph)<half_in); o=(ep==e)&(np.abs(ph)>0.35)&(np.abs(ph)<1.2)
            if s.sum()<4 or o.sum()<20: continue
            tt=t[o]; A=np.c_[np.ones(o.sum()),tt-tt.mean()]; As=np.c_[np.ones(s.sum()),t[s]-tt.mean()]
            for k,ap in aps.items():
                y=np.nansum(fl[:,ap],axis=1); c=np.linalg.lstsq(A,y[o],rcond=None)[0]; ev[k].append(np.mean(As@c-y[s]))
            y=np.nansum(bk[:,valid],axis=1); c=np.linalg.lstsq(A,y[o],rcond=None)[0]; ev['bkg'].append(np.mean(As@c-y[s]))
        n=len(ev['target'])
        if n<2: continue
        row={'Q':int(h0['QUARTER']),'n':n}
        for k in ev: a=np.array(ev[k]); row[k]=(float(a.mean()),float(a.std(ddof=1)/np.sqrt(n)))
        per_q.append(row)
        print(f"{lab} Q{row['Q']:2d} n={n} target {row['target'][0]:5.2f}+-{row['target'][1]:.2f} | far {row['far'][0]:5.2f}+-{row['far'][1]:.2f} | nbr {row['nbr'][0]:6.2f}+-{row['nbr'][1]:.2f} | all {row['all'][0]:6.2f}+-{row['all'][1]:.2f} | bkg {row['bkg'][0]:5.2f}+-{row['bkg'][1]:.2f}")
    comb={}
    for k in ['target','far','nbr','all','bkg']:
        d_=np.array([r[k][0] for r in per_q]); e_=np.array([r[k][1] for r in per_q]); wgt=1/e_**2
        comb[k]=(float(np.sum(wgt*d_)/np.sum(wgt)),float(1/np.sqrt(np.sum(wgt))),float(np.sum((d_-np.sum(wgt*d_)/np.sum(wgt))**2*wgt)/(len(d_)-1)))
        print(f'{lab} COMBINED {k:7s}: {comb[k][0]:6.3f} +- {comb[k][1]:.3f} e/s  (chi2/dof across quarters {comb[k][2]:.2f})')
    res[lab]={'per_q':per_q,'combined':comb}
json.dump(res,open('pixel_events.json','w'),indent=1)
