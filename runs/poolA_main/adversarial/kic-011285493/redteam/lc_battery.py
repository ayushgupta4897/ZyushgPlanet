import numpy as np, json, sys
from astropy.io import fits
from astropy.timeseries import LombScargle
np.random.seed(1)
F='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-011285493_kepler_v1.0_lc.fits'
h=fits.open(F); d=h['LIGHTCURVE_STITCHED'].data
t=d['TIME']; f=d['FLUX']; fe=d['FLUX_ERR']; sap=d['SAP_FLUX']; q=d['SAP_QUALITY']; Q=d['QUARTER']
m=np.isfinite(t)&np.isfinite(f)&(q==0)&np.isfinite(sap)
t,f,fe,sap,Q=t[m],f[m],fe[m],sap[m],Q[m]
P=10.790158; t0=843.73888; T14=5.91/24
# per-quarter normalise
fn=f.copy(); sn=sap.copy()
for qq in np.unique(Q):
    s=Q==qq; fn[s]/=np.nanmedian(f[s]); sn[s]/=np.nanmedian(sap[s])
def phase(t): return ((t-t0+0.5*P)%P)-0.5*P
ph=phase(t)
intr=np.abs(ph)<0.75*T14
def detrend(y,win=1.0,mask=None):
    out=np.full_like(y,np.nan)
    gaps=np.where(np.diff(t)>0.5)[0]; edges=np.concatenate([[0],gaps+1,[len(t)]])
    for a,b in zip(edges[:-1],edges[1:]):
        tt=t[a:b]; yy=y[a:b]; mm=None if mask is None else mask[a:b]
        for i in range(len(tt)):
            s=np.abs(tt-tt[i])<win/2
            if mm is not None: s&=~mm
            out[a+i]=yy[i]/np.nanmedian(yy[s]) if s.sum()>5 else np.nan
    return out
fd=detrend(fn,1.0,intr); sd=detrend(sn,1.0,intr)
ok=np.isfinite(fd)&np.isfinite(sd)
t,fd,sd,Q,ph=t[ok],fd[ok],sd[ok],Q[ok],ph[ok]
r=(fd-1)*1e6; rs=(sd-1)*1e6
mad=lambda x:1.4826*np.median(np.abs(x-np.median(x)))
res={}
# ---- box depth with random-phase noise
def boxdepth(y,ph,w):
    i=np.abs(ph)<w/2; o=(np.abs(ph)>w)&(np.abs(ph)<3*w)
    return np.mean(y[o])-np.mean(y[i]), i.sum()
def rand_noise(y,w,n=1500):
    ds=[]
    for k in range(n):
        tt0=np.random.uniform(0,P); p2=((t-tt0+0.5*P)%P)-0.5*P
        dd,_=boxdepth(y,p2,w); ds.append(dd)
    return mad(np.array(ds))
dep,nin=boxdepth(r,ph,T14); noise=rand_noise(r,T14)
res['box_T14']=dict(depth_ppm=float(dep),n_in=int(nin),rand_noise_ppm=float(noise),snr=float(dep/noise))
print('T14 box depth',dep,'+-',noise, 'SNR',dep/noise, 'n_in',nin)
# ---- per season / quarter
seas={}
for qq in np.unique(Q):
    s=Q==qq; dd,n=boxdepth(r[s],ph[s],T14); 
    o=(np.abs(ph[s])>T14)&(np.abs(ph[s])<3*T14)
    err=np.std(r[s][o])*np.sqrt(1/max(n,1)+1/max(o.sum(),1))
    seas[int(qq)]=(float(dd),float(err),int(n))
print('per quarter (ppm, white err, n_in):'); 
for k,v in seas.items(): print(k, [round(x,0) for x in v], 'season',k%4)
res['per_quarter']=seas
for sidx in range(4):
    pass
for sidx in range(4):
    s=((Q-2)%4)==sidx
    if s.sum()<100: continue
    dd,n=boxdepth(r[s],ph[s],T14); o=(np.abs(ph[s])>T14)&(np.abs(ph[s])<3*T14)
    err=np.std(r[s][o])*np.sqrt(1/n+1/o.sum())
    print('season',sidx,'quarters',sorted(set(Q[s].astype(int))),'depth',round(dd),'+-',round(err),'n_in',n)
    res[f'season{sidx}']=dict(depth=float(dd),err=float(err),n_in=int(n))
# ---- shape: trapezoid chi2 grid with cadence integration (T23/T14)
def trap(ph,dep,T14,ratio):
    T23=ratio*T14; tau=(T14-T23)/2
    x=np.abs(ph); y=np.zeros_like(x)
    y[x<=T23/2]=-dep
    s=(x>T23/2)&(x<T14/2); y[s]=-dep*(T14/2-x[s])/max(tau,1e-9)
    return y
def integ(ph,dep,T14,ratio,n=7):
    dt=29.4/60/24; out=np.zeros_like(ph)
    for k in np.linspace(-dt/2,dt/2,n): out+=trap(ph+k,dep,T14,ratio)
    return out/n
w=np.abs(ph)<0.6
sig=mad(r[np.abs(ph)>0.3])
best=None; grid={}
for ratio in [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.0]:
    for TT in np.arange(3.0,8.6,0.25)/24:
        A=integ(ph[w],1.0,TT,ratio)
        dep_=np.sum(A*r[w])/np.sum(A*A)
        chi=np.sum((r[w]-dep_*A)**2)/sig**2
        if ratio not in grid or chi<grid[ratio][0]: grid[ratio]=(chi,TT*24,-dep_)
        if best is None or chi<best[0]: best=(chi,ratio,TT*24,-dep_)
print('shape grid best', best)
for k,v in grid.items(): print(' T23/T14=%.1f  dchi2=%.2f  T14=%.2f h depth=%.0f'%(k,v[0]-best[0],v[1],v[2]))
res['shape']=dict(best=best,grid={str(k):v for k,v in grid.items()},sigma_ppm=float(sig))
# ---- secondary scan over all phases (T14 box), and at 0.5
phs=np.linspace(-0.5,0.5,201)*P
sec=[]
for pc in phs:
    p2=((t-t0-pc+0.5*P)%P)-0.5*P
    dd,_=boxdepth(r,p2,T14); sec.append(dd)
sec=np.array(sec)
off=np.abs(phs)>1.5*T14
p5=((t-t0-P/2+0.5*P)%P)-0.5*P; d5,_=boxdepth(r,p5,T14)
print('secondary at 0.5: %.0f ppm (%.1f sigma); max off-transit dip %.0f ppm at phase %.3f; min (bump) %.0f'%(d5,d5/noise,sec[off].max(),phs[off][np.argmax(sec[off])]/P,sec[off].min()))
res['secondary']=dict(ph05_ppm=float(d5),sigma=float(d5/noise),max_offtransit_ppm=float(sec[off].max()),at_phase=float(phs[off][np.argmax(sec[off])]/P))
# ---- odd/even
ep=np.round((t-t0)/P).astype(int)
for lab,s in [('odd',ep%2==1),('even',ep%2==0)]:
    dd,n=boxdepth(r[s],ph[s],T14); print(lab,round(dd),n)
# ---- LS both polarities on out-of-transit, detrended flux is high-passed; also on per-quarter-normalised raw
oo=np.abs(ph)>T14
for lab,y in [('detrended',r[oo]),('raw_qnorm',(fn[ok][oo]-1)*1e6)]:
    ls=LombScargle(t[oo],y); fr,pw=ls.autopower(minimum_frequency=1/50,maximum_frequency=1/0.2,samples_per_peak=5)
    i=np.argsort(pw)[::-1][:5]
    print(lab,'LS top:',[(round(1/fr[j],4),round(pw[j],5)) for j in i])
    # amplitude of best
    res['ls_'+lab]=[(float(1/fr[j]),float(pw[j])) for j in i]
    for Ptest in [P,P/2,2*P,P/3]:
        print('   power at',round(Ptest,3),':',round(float(ls.power(1/Ptest)),6))
# inverted (2-f) LS is identical for LS (linear), so skip; instead check for dips vs bumps asymmetry: skewness of residuals
print('skew of residuals', float(((r-np.mean(r))**3).mean()/np.std(r)**3))
# ---- event times vs quarter edges / gaps
gaps=np.where(np.diff(t)>0.5)[0]; gap_t=np.concatenate([[t[0]],(t[gaps]+t[gaps+1])/2,[t[-1]]])
ev=[]
for e in np.unique(ep):
    tc=t0+e*P; s=np.abs(t-tc)<T14/2
    if s.sum()<4: continue
    o=(np.abs(t-tc)>T14)&(np.abs(t-tc)<3*T14)
    if o.sum()<10: continue
    dd=np.mean(r[o])-np.mean(r[s]); err=np.std(r[o])*np.sqrt(1/s.sum()+1/o.sum())
    gd=np.min(np.abs(gap_t-tc)); qq=int(Q[s][0])
    ev.append((int(e),float(tc),float(dd),float(err),int(s.sum()),float(gd),qq))
ev=np.array(ev,dtype=object)
dd=np.array([x[2] for x in ev]); ee=np.array([x[3] for x in ev]); gd=np.array([x[5] for x in ev])
wm=np.sum(dd/ee**2)/np.sum(1/ee**2)
print('n events',len(ev),'weighted mean',round(wm),'chi2/dof',round(np.sum((dd-wm)**2/ee**2)/(len(dd)-1),2))
for cut in [0.5,1.0,2.0]:
    s=gd<cut
    print(f' events within {cut} d of gap/edge: {s.sum()}  depth near-gap {np.sum(dd[s]/ee[s]**2)/np.sum(1/ee[s]**2):.0f}  far {np.sum(dd[~s]/ee[~s]**2)/np.sum(1/ee[~s]**2):.0f}')
# top 5 events
i=np.argsort(dd/ee)[::-1][:6]
print('strongest events (epoch, tc, depth, err, nin, gapdist, Q):'); [print(ev[j]) for j in i]
# drop top 3 strongest -> depth
keep=np.ones(len(dd),bool); keep[i[:3]]=False
print('weighted mean without 3 strongest', round(np.sum(dd[keep]/ee[keep]**2)/np.sum(1/ee[keep]**2)))
# ---- SAP absolute deficit vs PSF
sapmed=np.nanmedian(sap); psfmed=np.nanmedian(f)
dS,_=boxdepth(rs,ph,T14); nS=rand_noise(rs,T14,800)
print('PSF median %.0f e/s, SAP median %.0f e/s; PSF depth %.0f ppm -> %.2f e/s; SAP depth %.0f +- %.0f ppm -> %.2f e/s'%(psfmed,sapmed,dep,dep*1e-6*psfmed,dS,nS,dS*1e-6*sapmed))
res['sap']=dict(psf_med=float(psfmed),sap_med=float(sapmed),psf_depth_ppm=float(dep),sap_depth_ppm=float(dS),sap_err_ppm=float(nS))
json.dump(res,open('lc_battery.json','w'),indent=1,default=str)
