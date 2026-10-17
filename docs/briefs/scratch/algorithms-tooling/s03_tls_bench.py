"""TLS 2.0.1 timing benchmark on a synthetic 4-yr Kepler-like LC (30-min cadence, ~65k points)."""
import os, sys, time, json, numpy as np
os.environ['OMP_NUM_THREADS']='1'
from transitleastsquares import transitleastsquares
rng=np.random.default_rng(1)
# 4 years at 29.4-min cadence with ~17 quarterly gaps -> ~65k points (Kepler Q1-Q17 ~ 65k LC cadences)
t=np.arange(0,1460,29.4244/1440.)
# remove ~8% as gaps (monthly downlinks + quarter rolls)
keep=np.ones(len(t),bool)
for g0 in np.arange(0,1460,31.): keep&=~((t>g0)&(t<g0+1.0))
for q0 in np.arange(0,1460,93.): keep&=~((t>q0)&(t<q0+3.0))
t=t[keep]
sig=0.0012 # 1200 ppm per 30-min for a Kp~16 M dwarf (Kuber's star rms)
flux=1+rng.normal(0,sig,len(t))
# inject 13.77-d box, 2000 ppm, 2.3 h
P=13.769448; T0=t[0]+5.3; dur=2.3/24
phase=((t-T0+P/2)%P)-P/2
flux[np.abs(phase)<dur/2]-=0.0020
print('N points',len(t),'baseline',t.max()-t.min())
cfg=sys.argv[1]
res={}
def run(label,t,f,**kw):
    model=transitleastsquares(t,f)
    t1=time.time()
    r=model.power(R_star=0.507,M_star=0.505,R_star_min=0.3,R_star_max=0.8,M_star_min=0.3,M_star_max=0.8,
                  show_progress_bar=False,verbose=False,**kw)
    dt=time.time()-t1
    out=dict(seconds=round(dt,1),n_periods=len(r.periods),period=float(r.period),SDE=float(r.SDE),snr=float(r.snr),depth=float(r.depth),duration_hr=float(r.duration*24),transit_count=int(r.transit_count))
    print(label,json.dumps(out),flush=True); res[label]=out
if cfg=='A': run('A_full_1-30d_2threads',t,flux,period_min=1,period_max=30,use_threads=2)
if cfg=='B': run('B_full_1-10d_2threads',t,flux,period_min=1,period_max=10,use_threads=2)
if cfg=='C':
    # bin to 1 hr
    w=1/24.; key=np.floor((t-t[0])/w).astype(int); u,inv,cnt=np.unique(key,return_inverse=True,return_counts=True)
    tb=np.bincount(inv,t)/cnt; fb=np.bincount(inv,flux)/cnt
    print('binned N',len(tb)); run('C_binned1h_1-30d_2threads',tb,fb,period_min=1,period_max=30,use_threads=2)
if cfg=='D': run('D_full_1-30d_1thread_oversamp2',t,flux,period_min=1,period_max=30,use_threads=1,oversampling_factor=2,duration_grid_step=1.15)
if cfg=='E': run('E_full_1-30d_2threads_backend_fused',t,flux,period_min=1,period_max=30,use_threads=2,backend='fused')
json.dump(res,open(f'/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/tls_bench_{cfg}.json','w'),indent=1)
