"""batman + emcee fit of the cleaned KBONUS LC at the candidate ephemeris (30-min exposure integration), then LEO-vetter flux metrics."""
import os; os.environ['OMP_NUM_THREADS']='1'
import time, json, numpy as np, warnings; warnings.filterwarnings('ignore')
import batman, emcee
D='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/data'
d=np.load(f'{D}/kic5689121_clean.npz'); t=d['time']; f=d['flux']; e=d['flux_err']
P0=13.769448; T0=2455630.245431-2454833  # BKJD
ph=((t-T0+P0/2)%P0)-P0/2; sel=np.abs(ph)<0.4  # +-0.4 d around transits
t,f,e=t[sel],f[sel],e[sel]; print('points within 0.4 d of transit', len(t))
params=batman.TransitParams(); params.t0=T0; params.per=P0; params.rp=0.046; params.a=39.7; params.inc=89.3; params.ecc=0; params.w=90; params.limb_dark='quadratic'; params.u=[0.35,0.25]  # ~Claret 2011 Kepler band for Teff~3600 logg~4.8 (approx)
m=batman.TransitModel(params,t,supersample_factor=7,exp_time=1765.464/86400)  # 29.42-min Kepler LC exposure
G=6.674e-11
def rho_from_aRs(aRs,P): return 3*np.pi*aRs**3/(G*(P*86400)**2)/1408.0  # solar units (rho_sun=1408 kg/m3)
def lnprob(th):
    t0,per,rp,aRs,b,f0,lnjit=th
    if not (0<rp<0.3 and 1<aRs<200 and 0<=b<1+rp and -12<lnjit<-4 and abs(per-P0)<0.01): return -np.inf
    inc=np.degrees(np.arccos(b/aRs)); params.t0,params.per,params.rp,params.a,params.inc=t0,per,rp,aRs,inc
    mod=m.light_curve(params)*f0; s2=e**2+np.exp(2*lnjit)
    return -0.5*np.sum((f-mod)**2/s2+np.log(s2))
ndim=7; nw=32; p0=np.array([T0,P0,0.046,39.7,0.5,1.0,np.log(3e-4)])
pos=p0+np.array([1e-3,1e-5,2e-3,3,0.1,1e-4,0.1])*np.random.default_rng(1).normal(size=(nw,ndim))
pos[:,4]=np.abs(pos[:,4])
t1=time.time(); sampler=emcee.EnsembleSampler(nw,ndim,lnprob); sampler.run_mcmc(pos,2500,progress=False)
chain=sampler.get_chain(discard=1000,flat=True); dt=time.time()-t1
print('emcee 32 walkers x 2500 steps on %d points with supersample 7: %.0fs; acceptance %.2f'%(len(t),dt,np.mean(sampler.acceptance_fraction)))
names=['t0','per','rp_rs','a_rs','b','f0','lnjit']
q=lambda x:np.percentile(x,[16,50,84])
summ={}
for i,n in enumerate(names):
    lo,md,hi=q(chain[:,i]); summ[n]=(md,hi-md,md-lo); print(f'{n}: {md:.6f} +{hi-md:.6f} -{md-lo:.6f}')
rho=rho_from_aRs(chain[:,3],chain[:,1]); lo,md,hi=q(rho); print('stellar density from transit: %.2f +%.2f -%.2f rho_sun (M dwarf 0.505 Msun/0.507 Rsun -> %.2f)'%(md,hi-md,md-lo,0.505/0.507**3))
rp=chain[:,2]*0.507*109.1; lo,md,hi=q(rp); print('Rp = %.2f +%.2f -%.2f Re (for R*=0.507 Rsun, no dilution correction)'%(md,hi-md,md-lo))
dep=chain[:,2]**2; lo,md,hi=q(dep); print('(Rp/R*)^2 = %.0f +%.0f -%.0f ppm'%(md*1e6,(hi-md)*1e6,(md-lo)*1e6))
json.dump({k:list(map(float,v)) for k,v in summ.items()}|{'seconds':dt},open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/fit_result.json','w'),indent=1)
# ---- LEO-vetter flux-level metrics
print('\n=== LEO-vetter 1.2.0 flux-level vetting')
from leo_vetter.main import TCELightCurve
from leo_vetter.thresholds import check_thresholds
d=np.load(f'{D}/kic5689121_clean.npz'); t=d['time']; flux=d['flux']; err=d['flux_err']
raw=flux.copy()  # we only have the detrended series handy; 'raw' is used for plots/robustness checks
star=dict(mass=0.505,rad=0.507,Teff=3633,logg=4.73,u1=0.35,u2=0.25,Tmag=14.9,Vmag=17.5,tic=120686465,rho=0.505/0.507**3)
t1=time.time()
tlc=TCELightCurve(120686465,t,raw,flux,err,P0,T0,2.3/24)
try:
    tlc.compute_flux_metrics(star)
    M=tlc.metrics; print('LEO-vetter compute_flux_metrics: %.1fs; %d metrics'%(time.time()-t1,len(M)))
    keys=['MES','N_transit','SHP','CHI','DMM','Fred','sig_pri','sig_sec','sig_ter','sig_pos','max_SES','mean_chases','sine_sig','transit_b','transit_aRs','transit_RpRs','trap_qtran','dep','per','epo','dur']
    print({k:(round(float(M[k]),4) if isinstance(M.get(k),(int,float,np.floating,np.integer)) else M.get(k)) for k in keys if k in M})
    FA=check_thresholds(M,'FA'); FP=check_thresholds(M,'FP'); print('check_thresholds FA ->',FA); print('check_thresholds FP ->',FP)
except Exception as ex:
    import traceback; traceback.print_exc(); print('metrics so far', list(getattr(tlc,'metrics',{}).keys())[:40])
