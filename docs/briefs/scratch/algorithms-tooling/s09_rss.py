"""Peak RSS of one TLS worker (use_threads=1, 1-30 d) and one astropy-BLS worker on the real cleaned LC."""
import os, sys, resource, time, numpy as np, warnings; warnings.filterwarnings('ignore')
os.environ['OMP_NUM_THREADS']='1'
D='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/data'
d=np.load(f'{D}/kic5689121_clean.npz'); t=d['time']; f=d['flux']; e=d['flux_err']
which=sys.argv[1]
t1=time.time()
if which=='tls':
    from transitleastsquares import transitleastsquares
    r=transitleastsquares(t,f,e).power(period_min=1,period_max=30,R_star=0.507,M_star=0.505,R_star_min=0.3,R_star_max=0.8,M_star_min=0.3,M_star_max=0.8,u=[0.2125,0.3984],use_threads=1,show_progress_bar=False,verbose=False)
else:
    from astropy.timeseries import BoxLeastSquares
    rho=0.505/0.507**3
    def dur_days(P,rho): return 13/24*(P/365.25)**(1/3)*rho**(-1/3)
    base=t.max()-t.min(); fr=[];fq=2.0
    while fq>1/30: fr.append(fq);P=1/fq;fq-=0.5*dur_days(P,rho)/(3*base*P)
    periods=np.sort(1/np.array(fr)); r=BoxLeastSquares(t,f-1,e).power(periods,np.geomspace(0.03,0.25,6),objective='snr',oversample=3)
ru=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
print(f'{which}: wall {time.time()-t1:.0f}s, peak RSS (main process) {ru:.0f} MB')
if which=='tls':
    ch=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024; print(f'tls: peak RSS of largest child (pool worker) {ch:.0f} MB')
