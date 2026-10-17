"""(a) Null calibration on the real LC: inverted + quarter-scrambled copies through the identical clean+BLS.
   (b) Throughput: 4-worker ProcessPool over synthetic Kepler-like and TESS-like stars (wotan + BLS + stats).
   (c) TLS use_threads 4 vs 1 on the real cleaned LC."""
import os; os.environ['OMP_NUM_THREADS']='1'; os.environ['OPENBLAS_NUM_THREADS']='1'
import time, json, sys, numpy as np, warnings; warnings.filterwarnings('ignore')
from concurrent.futures import ProcessPoolExecutor
from wotan import flatten
from astropy.timeseries import BoxLeastSquares
D='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling'
def dur_days(P,rho): return 13/24*(P/365.25)**(1/3)*rho**(-1/3)
def grid(base,rho,pmin=0.5,pmax=30,os_=3):
    f=[];fq=1/pmin
    while fq>1/pmax:
        f.append(fq);P=1/fq;fq-=0.5*dur_days(P,rho)/(os_*base*P)
    return np.sort(1/np.array(f))
def clean(t,fl,seg,window=1.0):
    flat=np.full_like(fl,np.nan)
    for s in np.unique(seg):
        m=seg==s
        if m.sum()>50: flat[m]=flatten(t[m],fl[m],method='biweight',window_length=window,break_tolerance=0.5)
    ok=np.isfinite(flat); r=flat-1; sd=1.4826*np.median(np.abs(r[ok]-np.median(r[ok]))); ok&=(r<3*sd)&(r>-10*sd)
    return t[ok],flat[ok]-1,np.full(ok.sum(),sd),seg[ok]
def sde_of(power):
    p=np.asarray(power); return (p-np.median(p))/np.std(p)
def bls_search(t,y,dy,rho,pmax=30,os_=3,ndur=6,oversample=3):
    base=t.max()-t.min(); periods=grid(base,rho,0.5,min(pmax,base/3),os_)
    durs=np.geomspace(max(0.03,0.35*dur_days(0.5,rho)),min(0.3,1.3*dur_days(pmax,rho)),ndur)
    bls=BoxLeastSquares(t,y,dy); r=bls.power(periods,durs,objective='snr',method='fast',oversample=oversample)
    s=sde_of(r.power); k=int(np.argmax(s)); st=bls.compute_stats(r.period[k],r.duration[k],r.transit_time[k])
    return dict(P=float(r.period[k]),t0=float(r.transit_time[k]),dur_h=float(r.duration[k]*24),depth_ppm=float(r.depth[k]*1e6),snr=float(r.depth[k]/r.depth_err[k]),power=float(r.power[k]),sde=float(s[k]),n_periods=len(periods),n_transits=int(np.sum(st['per_transit_count']>0)))
# ---------- (a) null runs on the real LC
RAW=None
def load_raw():
    import glob; from astropy.io import fits
    f=glob.glob(f'{D}/data/**/*kic-005689121*lc.fits',recursive=True)[0]
    T=[];F=[];Q=[]
    with fits.open(f) as h:
        for hd in h:
            if not hd.name.startswith('LIGHTCURVE_Q'): continue
            d=hd.data; t=np.asarray(d['TIME'],float); fl=np.asarray(d['FLUX'],float); e=np.asarray(d['FLUX_ERR'],float)
            good=(d['SAP_QUALITY']==0)&np.isfinite(t)&np.isfinite(fl)&np.isfinite(e)&(e>0); med=np.nanmedian(fl[good])
            T.append(t[good]);F.append(fl[good]/med);Q.append(np.full(good.sum(),int(hd.header['QUARTER'])))
    t=np.concatenate(T);fl=np.concatenate(F);q=np.concatenate(Q);o=np.argsort(t); return t[o],fl[o],q[o]
def null_job(job):
    kind,seed=job; t,fl,q=load_raw(); rng=np.random.default_rng(seed)
    if kind=='inverted': fl=2-fl
    elif kind=='scrambled':
        # keep each quarter's time stamps, permute which quarter's flux series goes where (resampled to length) and circularly shift within quarter
        qs=np.unique(q); perm=rng.permutation(qs); newf=fl.copy()
        for a,b in zip(qs,perm):
            ma=q==a; mb=q==b; src=fl[mb]; n=ma.sum()
            src=np.interp(np.linspace(0,1,n),np.linspace(0,1,len(src)),src); newf[ma]=np.roll(src,rng.integers(n))
        fl=newf
    t1=time.time(); tt,y,dy,seg=clean(t,fl,q); r=bls_search(tt,y,dy,0.505/0.507**3); r.update(kind=kind,seed=seed,seconds=round(time.time()-t1,1)); return r
# ---------- (b) throughput synthetic
def synth_kepler(seed):
    rng=np.random.default_rng(seed); t=np.arange(0,1460,29.4244/1440.); keep=np.ones(len(t),bool)
    for g0 in np.arange(0,1460,31.): keep&=~((t>g0)&(t<g0+1.0))
    for q0 in np.arange(0,1460,93.): keep&=~((t>q0)&(t<q0+3.0))
    t=t[keep]; q=(t//93).astype(int); fl=1+rng.normal(0,0.0012,len(t))+0.002*np.sin(2*np.pi*t/rng.uniform(5,40))
    return t,fl,q
def synth_tess(seed):
    rng=np.random.default_rng(seed); parts=[]
    for k,start in enumerate([0,27,365,392,730]):  # 5 sectors over 2 years, 10-min FFI cadence
        tt=np.arange(start,start+27,10/1440.); tt=tt[(tt-start<13)|(tt-start>14)]; parts.append((tt,np.full(len(tt),k)))
    t=np.concatenate([p[0] for p in parts]); q=np.concatenate([p[1] for p in parts]); fl=1+rng.normal(0,0.0015,len(t))+0.003*np.sin(2*np.pi*t/rng.uniform(2,20))
    return t,fl,q
def tp_job(job):
    kind,seed=job; t,fl,q=(synth_kepler if kind=='kepler' else synth_tess)(seed); t1=time.time()
    tt,y,dy,seg=clean(t,fl,q,window=1.0 if kind=='kepler' else 0.75)
    if kind=='tess':  # bin to 20 min for the coarse search as eumemic/Chopra do
        w=20/1440.; key=np.floor((tt-tt[0])/w).astype(int); u,inv,cnt=np.unique(key,return_inverse=True,return_counts=True)
        tt=np.bincount(inv,tt)/cnt; y=np.bincount(inv,y)/cnt; dy=np.full(len(u),dy[0]/np.sqrt(2))
    r=bls_search(tt,y,dy,3.9); r.update(kind=kind,seed=seed,n_points=len(tt),seconds=round(time.time()-t1,1)); return r
if __name__=='__main__':
    out={}
    jobs=[('inverted',0)]+[('scrambled',s) for s in range(1,8)]+[('real',99)]
    t1=time.time()
    with ProcessPoolExecutor(4) as ex: res=list(ex.map(null_job,jobs))
    out['null']=res; print('NULL RUNS (%.0fs wall, 4 workers):'%(time.time()-t1)); 
    for r in res: print(json.dumps(r))
    for kind,n in [('kepler',8),('tess',8)]:
        t1=time.time()
        with ProcessPoolExecutor(4) as ex: res=list(ex.map(tp_job,[(kind,s) for s in range(n)]))
        wall=time.time()-t1; out[kind]=dict(wall=wall,per_star_cpu=np.mean([r['seconds'] for r in res]),n_points=res[0]['n_points'],n_periods=res[0]['n_periods'])
        print(f'THROUGHPUT {kind}: {n} stars, 4 workers: wall {wall:.0f}s -> {n/wall*60:.1f} stars/min; per-star CPU {out[kind]["per_star_cpu"]:.1f}s; {res[0]["n_points"]} pts, {res[0]["n_periods"]} periods')
    # (c) TLS threads
    from transitleastsquares import transitleastsquares
    d=np.load(f'{D}/data/kic5689121_clean.npz'); t=d['time']; f=d['flux']; e=d['flux_err']
    for th in (4,1):
        t1=time.time(); m=transitleastsquares(t,f,e); r=m.power(period_min=1,period_max=30,R_star=0.507,M_star=0.505,R_star_min=0.3,R_star_max=0.8,M_star_min=0.3,M_star_max=0.8,u=[0.2125,0.3984],use_threads=th,show_progress_bar=False,verbose=False)
        print(f'TLS real LC use_threads={th}: {time.time()-t1:.0f}s, SDE {r.SDE:.1f}, P {r.period:.5f}'); out[f'tls_threads_{th}']=time.time()-t1
    json.dump(out,open(f'{D}/null_throughput.json','w'),indent=1,default=float)
