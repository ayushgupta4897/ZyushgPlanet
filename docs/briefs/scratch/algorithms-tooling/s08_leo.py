import os; os.environ['OMP_NUM_THREADS']='1'
import time, numpy as np, warnings; warnings.filterwarnings('ignore')
from leo_vetter.main import TCELightCurve
from leo_vetter.thresholds import check_thresholds
D='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/data'
d=np.load(f'{D}/kic5689121_clean.npz'); t=d['time']; flux=d['flux']; err=d['flux_err']
P0=13.769448; T0=2455630.245431-2454833
star=dict(id=120686465,tic=120686465,rad=0.507,e_rad=0.016,mass=0.505,e_mass=0.021,Teff=3633,e_Teff=160,logg=4.73,rho=0.505/0.507**3,
          u1=0.2125,u2=0.3984,ra=285.907827,dec=40.941988,raDegrees=285.907827,decDegrees=40.941988,Tmag=14.9,planetData=None)
t1=time.time(); tlc=TCELightCurve(120686465,t,flux.copy(),flux,err,P0,T0,2.3/24)
tlc.compute_flux_metrics(star); M=tlc.metrics
print('compute_flux_metrics: %.1fs, %d metrics'%(time.time()-t1,len(M)))
keys=['MES','new_MES','N_transit','new_N_transit','SHP','CHI','DMM','Fred','sig_pri','sig_sec','sig_ter','sig_pos','FA1','FA2','max_SES','mean_chases','sine_sig','transit_b','transit_aRs','transit_RpRs','transit_chisqr','line_chisqr','transit_aic','line_aic','trap_qtran','q','aRs','dep','sig_w','sig_r','Rp','Rp_err','Teq']
print({k:(round(float(M[k]),4) if isinstance(M.get(k),(int,float,np.floating,np.integer)) else str(M.get(k))) for k in keys if k in M})
print('FA ->',check_thresholds(M,'FA')); print('FP ->',check_thresholds(M,'FP'))
