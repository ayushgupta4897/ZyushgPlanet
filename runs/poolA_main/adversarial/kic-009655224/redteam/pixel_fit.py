# 3-source Gaussian-PRF amplitude fit to per-quarter local-difference images (same differencing as pixel_check.py),
# summed over quarters; plus the pipeline's own chi2_gauss_at_source sums from pixel.json.
import numpy as np, glob, json
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import nnls
P=6.703464464595031; t0=265.183161835447; dur=2.4375/24
srcs={'target':(294.28267299472,46.34150351561),'nb_G15.9':(294.28296961419,46.340060103966586),'host_G16.4':(294.2852862259,46.33963206258)}
files=sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr009655236_lc_Q000111011101111000/*.fits.gz'))
def gauss(xx,yy,x0,y0,s): 
    g=np.exp(-((xx-x0)**2+(yy-y0)**2)/(2*s*s)); return g/ g.sum() if False else g/(2*np.pi*s*s)  # normalised to unit total over infinite plane
out={}
for sig in [0.6,0.8,1.0]:
    AtA=np.zeros((3,3)); Atb=np.zeros(3); chi2={'null':0,'target_only':0,'nb_only':0,'host_only':0,'all3':0}; n=0
    per_q={}
    for fn in files:
        h=fits.open(fn); Q=h[0].header['QUARTER']; d=h[1].data
        t=d['TIME']; fl=d['FLUX']; qu=d['QUALITY']
        goodpx=np.isfinite(fl).mean(axis=0)>0.9
        ok=np.isfinite(t)&(qu==0)&np.isfinite(fl[:,goodpx]).all(axis=1); t,fl=t[ok],fl[ok]; fl=np.where(goodpx[None],fl,np.nan)
        w=WCS(h[2].header); ny,nx=fl.shape[1:]; yy,xx=np.mgrid[0:ny,0:nx]
        pos={k:tuple(map(float,w.all_world2pix(ra,dec,0))) for k,(ra,dec) in srcs.items()}
        n0=int(np.floor((t.min()-t0)/P)); n1=int(np.ceil((t.max()-t0)/P)); diffs=[]
        for k in range(n0,n1+1):
            tc=t0+k*P; dt=t-tc; ib=np.abs(dt)<dur/2*0.9; ob=(np.abs(dt)>dur/2*1.5)&(np.abs(dt)<dur/2*1.5+0.25)
            if ib.sum()<3 or ob.sum()<8: continue
            A=np.vstack([np.ones(ob.sum()),dt[ob]]).T; coef=np.linalg.lstsq(A,fl[ob].reshape(ob.sum(),-1),rcond=None)[0]
            model=(np.vstack([np.ones(ib.sum()),dt[ib]]).T@coef).reshape(ib.sum(),ny,nx); diffs.append((model-fl[ib]).mean(axis=0))
        diffs=np.array(diffs); nev=len(diffs)
        if nev<3: continue
        D=diffs.mean(axis=0); E=diffs.std(axis=0,ddof=1)/np.sqrt(nev); v=np.isfinite(D)&goodpx&(E>0)
        M=np.array([gauss(xx,yy,*pos[k],sig)[v] for k in srcs]).T  # npix x 3
        b=D[v]/E[v]; Mw=M/E[v][:,None]
        AtA+=Mw.T@Mw; Atb+=Mw.T@b; n+=v.sum()
        chi2['null']+=np.sum(b**2)
        for name,cols in [('target_only',[0]),('nb_only',[1]),('host_only',[2]),('all3',[0,1,2])]:
            a,_,_,_=np.linalg.lstsq(Mw[:,cols],b,rcond=None); chi2[name]+=np.sum((b-Mw[:,cols]@a)**2)
        a3=np.linalg.lstsq(Mw,b,rcond=None)[0]; per_q[int(Q)]=[float(x) for x in a3]
    cov=np.linalg.inv(AtA); a=cov@Atb; err=np.sqrt(np.diag(cov))
    print(f'sigma={sig} px: joint 3-source amplitudes (e/s, total light lost): target {a[0]:.2f}+/-{err[0]:.2f}, nb_G15.9 {a[1]:.2f}+/-{err[1]:.2f}, host {a[2]:.2f}+/-{err[2]:.2f}; corr(t,nb)={cov[0,1]/err[0]/err[1]:.2f}; npix={n}')
    print('   summed chi2: null %.1f | target-only %.1f | nb-only %.1f | host-only %.1f | all3 %.1f'%tuple(chi2[k] for k in ['null','target_only','nb_only','host_only','all3']))
    out[sig]={'amp':list(map(float,a)),'err':list(map(float,err)),'chi2':chi2,'per_quarter_amp':per_q}
json.dump(out,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic009655224/work/pixel_fit.json','w'),indent=1)
# pipeline pixel.json chi2 sums per source
d=json.load(open('/home/user/ZyushgPlanet/runs/poolA_main/pixel/kic-009655224/pixel.json'))
sums={}
for q in d['quarters']:
    if 'chi2_gauss_at_target' not in q: continue
    sums['null']=sums.get('null',0)+q['chi2_null']; sums['scaled_image']=sums.get('scaled_image',0)+q['chi2_scaled_image']
    for nb in q['neighbours']:
        if nb['chi2_gauss_at_source'] is not None: sums[str(nb['source_id'])[-6:]+'_G%.1f'%nb['G']]=sums.get(str(nb['source_id'])[-6:]+'_G%.1f'%nb['G'],0)+nb['chi2_gauss_at_source']
print('\npipeline pixel.json summed chi2 over all 10 quarters:'); [print('  %-20s %.1f'%(k,v)) for k,v in sorted(sums.items(), key=lambda kv: kv[1])]
