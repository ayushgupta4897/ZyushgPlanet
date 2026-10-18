import numpy as np, glob, json
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
P=15.31335; T0=144.28486; DUR=3.75/24
OUT='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic003831969'
tgt=(285.0768009047113,38.93627480915959)
nb={'host_G15.0':(285.0794030780475,38.93540094183849),'N11_G19.9':(285.07595958231,38.9392656287),'S11_G19.2':(285.07586231462,38.93317780568),'SSW13_G20.5':(285.07461272159,38.93297109858),'NNW16_G19.2':(285.07485634722,38.94034566012),'NNW22_G17.3':(285.07361496664,38.94183189401),'E_G20.2':(285.08305412963,38.93176463934)}
files=sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr003831977_lc_Q011111111111111111/*.fits.gz'))
seasons={}
for f in files:
    h=fits.open(f); q=h[0].header['QUARTER']; s=h[0].header['SEASON']
    d=h[1].data; t=d['TIME']; fl=d['FLUX']; qual=d['QUALITY']
    ok=(qual==0)&np.isfinite(t)&np.isfinite(fl.sum(axis=(1,2)))
    t,fl=t[ok],fl[ok]
    ph=((t-T0)/P+0.5)%1-0.5; w=DUR/P
    # per-event differenced: for each event, in = |ph|<0.4w (core), out = 1.2w<|ph|<3w within +-1.5d local window
    ev=np.round((t-T0)/P).astype(int)
    diffs=[]; wts=[]
    for e in np.unique(ev):
        m=ev==e; tt=t[m]; ff=fl[m]; pp=(tt-(T0+e*P))/P
        i=np.abs(pp)<0.4*w; o=(np.abs(pp)>1.2*w)&(np.abs(pp)<4*w)
        if i.sum()<3 or o.sum()<8: continue
        # linear baseline per pixel from out points
        x=tt-(T0+e*P)
        A=np.vstack([np.ones(o.sum()),x[o]]).T
        coef,_,_,_=np.linalg.lstsq(A,ff[o].reshape(o.sum(),-1),rcond=None)
        base_in=(np.vstack([np.ones(i.sum()),x[i]]).T@coef).reshape(i.sum(),*ff.shape[1:])
        dimg=(base_in-ff[i]).mean(axis=0)   # missing light (e/s) per pixel
        resid=ff[o]-(A@coef).reshape(o.sum(),*ff.shape[1:])
        var=resid.var(axis=0)*(1/i.sum()+1/o.sum())
        diffs.append(dimg); wts.append(1/np.maximum(var,1e-6))
    if len(diffs)==0:
        print('skip Q',q,'no usable events'); continue
    diffs=np.array(diffs); wts=np.array(wts)
    D=(diffs*wts).sum(0)/wts.sum(0); Derr=1/np.sqrt(wts.sum(0))
    wcs=WCS(h[2].header)
    tx,ty=wcs.all_world2pix(tgt[0],tgt[1],0)
    npx={k:[float(x) for x in wcs.all_world2pix(v[0],v[1],0)] for k,v in nb.items()}
    med=np.nanmedian(fl,axis=0)
    print('Q',q,'season',s,'events',len(diffs)); seasons.setdefault(s,[]).append(dict(q=q,D=D,Derr=Derr,tx=float(tx),ty=float(ty),npx=npx,med=med,nev=len(diffs)))
# fit per season: stack quarters of equal shape by co-adding with inverse-variance; fit Gaussian with free centre
def model(p,yy,xx):
    A,x0,y0,sig=p; return A*np.exp(-((xx-x0)**2+(yy-y0)**2)/(2*sig**2))
summary={}
fig,axes=plt.subplots(2,4,figsize=(18,8))
for si,(s,lst) in enumerate(sorted(seasons.items())):
    shapes=set(x['D'].shape for x in lst)
    lst2=[x for x in lst if x['D'].shape==max(shapes,key=lambda sh:sum(1 for y in lst if y['D'].shape==sh))]
    W=np.array([1/x['Derr']**2 for x in lst2]); Ds=np.array([x['D'] for x in lst2])
    D=(Ds*W).sum(0)/W.sum(0); E=1/np.sqrt(W.sum(0))
    ok=np.isfinite(D)&np.isfinite(E)
    yy,xx=np.indices(D.shape)
    tx=np.mean([x['tx'] for x in lst2]); ty=np.mean([x['ty'] for x in lst2])
    def resid(p): return ((model(p,yy,xx)-D)/E)[ok]
    best=None
    for x0 in np.arange(-2,D.shape[1]+2,0.5):
        for y0 in np.arange(-2,D.shape[0]+2,0.5):
            r=least_squares(resid,[D[ok].max()*2,x0,y0,0.7],bounds=([-50,-4,-4,0.4],[200,D.shape[1]+4,D.shape[0]+4,1.5]))
            if best is None or r.cost<best.cost: best=r
    chi2=2*best.cost; p=best.x
    # error from jacobian
    J=best.jac; cov=np.linalg.pinv(J.T@J); perr=np.sqrt(np.diag(cov))
    # chi2 with centre fixed at target
    def resid_t(pp): return resid([pp[0],tx,ty,pp[1]])
    rt=least_squares(resid_t,[p[0],0.7],bounds=([-50,0.4],[200,1.5])); chi2_t=2*rt.cost
    # chi2 with centre at each neighbour
    chi2_nb={}
    for k,v in lst2[0]['npx'].items():
        def rr(pp): return resid([pp[0],v[0],v[1],pp[1]])
        rn=least_squares(rr,[p[0],0.7],bounds=([-50,0.4],[200,1.5])); chi2_nb[k]=round(2*rn.cost,1)
    chi2_null=float((((D)/E)[ok]**2).sum())
    # total deficit within stamp & in target-centred gaussian
    pxscale=3.98
    dist_t=np.hypot(p[1]-tx,p[2]-ty)*pxscale
    summary[int(s)]=dict(quarters=[x['q'] for x in lst2],shape=list(D.shape),target_px=[round(tx,2),round(ty,2)],fit_centre_px=[round(p[1],2),round(p[2],2)],fit_centre_err_px=[round(perr[1],2),round(perr[2],2)],fit_sigma_px=round(p[3],2),fit_amp=round(p[0],2),
        chi2_free=round(chi2,1),chi2_at_target=round(chi2_t,1),chi2_null=round(chi2_null,1),npix=int(ok.sum()),chi2_at_neighbours=chi2_nb,
        fit_to_target_arcsec=round(dist_t,2),neighbours_px={k:[round(v[0],2),round(v[1],2)] for k,v in lst2[0]['npx'].items()},
        deficit_sum_e_s=round(float(np.nansum(D)),2),deficit_sum_err=round(float(np.sqrt(np.nansum(E**2))),2),
        stamp_flux_at_target_px=round(float(np.nanmean([x['med'][int(round(ty)),int(round(tx))] for x in lst2 if 0<=int(round(ty))<x['med'].shape[0] and 0<=int(round(tx))<x['med'].shape[1]])),1))
    ax=axes[0,si]; im=ax.imshow(D/E,origin='lower',cmap='RdBu_r',vmin=-5,vmax=5); ax.set_title(f'season {s} Q{[x["q"] for x in lst2]} deficit SNR')
    ax.plot(tx,ty,'go',mfc='none',ms=14,mew=2); ax.plot(p[1],p[2],'k+',ms=14,mew=2)
    for k,v in lst2[0]['npx'].items(): ax.plot(v[0],v[1],'x',color='orange'); ax.text(v[0],v[1],k,fontsize=6,color='orange')
    ax.set_xlim(-3,D.shape[1]+2); ax.set_ylim(-3,D.shape[0]+2)
    ax=axes[1,si]; ax.imshow(np.log10(np.clip(lst2[0]['med'],1,None)),origin='lower',cmap='gray'); ax.plot(tx,ty,'go',mfc='none',ms=14,mew=2); ax.set_title('log out-of-transit image')
    for k,v in lst2[0]['npx'].items(): ax.plot(v[0],v[1],'x',color='orange')
    ax.set_xlim(-3,D.shape[1]+2); ax.set_ylim(-3,D.shape[0]+2)
plt.tight_layout(); plt.savefig(OUT+'/pix_free.png',dpi=110)
json.dump(summary,open(OUT+'/pix_free.json','w'),indent=1)
for s,v in summary.items(): print(s,json.dumps(v))
