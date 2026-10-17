"""End-to-end on the real KBONUS-BKG LC of KIC 5689121: load -> clean -> wotan -> BLS -> TLS. 1 thread."""
import os; os.environ['OMP_NUM_THREADS']='1'
import time, json, glob, warnings, numpy as np
warnings.filterwarnings('ignore')
from astropy.io import fits
from wotan import flatten
from astropy.timeseries import BoxLeastSquares
f=glob.glob('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/data/**/*kic-005689121*lc.fits',recursive=True)[0]
T=[];F=[];E=[];Q=[];meta=[]
with fits.open(f) as h:
    for hd in h:
        if not hd.name.startswith('LIGHTCURVE_Q'): continue
        q=int(hd.header['QUARTER']); d=hd.data
        psffrac=hd.header.get('PSFFRAC'); crowd=hd.header.get('CROWDSAP'); flfrc=hd.header.get('FLFRCSAP')
        t=np.asarray(d['TIME'],float); fl=np.asarray(d['FLUX'],float); e=np.asarray(d['FLUX_ERR'],float)
        good=(d['SAP_QUALITY']==0)&np.isfinite(t)&np.isfinite(fl)&np.isfinite(e)&(e>0)
        med=np.nanmedian(fl[good])
        meta.append(dict(quarter=q,n=int(good.sum()),psffrac=psffrac,crowdsap=crowd,flfrcsap=flfrc,median_flux=float(med),
                         ptp_ppm=float(1.4826*np.median(np.abs(np.diff(fl[good])))/np.sqrt(2)/med*1e6)))
        T.append(t[good]);F.append(fl[good]/med);E.append(e[good]/med);Q.append(np.full(good.sum(),q))
print(json.dumps(meta,indent=0))
t=np.concatenate(T);fl=np.concatenate(F);e=np.concatenate(E);q=np.concatenate(Q);o=np.argsort(t);t,fl,e,q=t[o],fl[o],e[o],q[o]
print('N points',len(t),'baseline d',t.max()-t.min(),'quarters',sorted(set(q.tolist())))
# detrend per quarter, biweight window 1.0 d (= ~10x the 2.3 h duration; Kuber used 1.0 d), break_tolerance 0.5 d
t1=time.time(); flat=np.full_like(fl,np.nan); trend=np.full_like(fl,np.nan)
for qq in np.unique(q):
    m=q==qq
    flat[m],trend[m]=flatten(t[m],fl[m],method='biweight',window_length=1.0,break_tolerance=0.5,edge_cutoff=0.0,return_trend=True)
print('wotan biweight on %d pts: %.2fs'%(len(t),time.time()-t1))
ok=np.isfinite(flat)
# asymmetric clip: upward > 3 sigma (flares), downward only > 10 sigma (dropouts); never clip modest dips
res=flat-1; sd=1.4826*np.median(np.abs(res[ok]-np.median(res[ok])))
ok&=(res<3*sd)&(res>-10*sd)
t,y,q=t[ok],flat[ok]-1,q[ok]; dy=np.full_like(y,sd)
print('after clip N',len(t),'robust sigma ppm',sd*1e6)
np.savez('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/data/kic5689121_clean.npz',time=t,flux=y+1,flux_err=dy,quarter=q)
# BLS with stellar-density-aware grid (rho ~ 3.9 solar for 0.505 Msun/0.507 Rsun)
rho=0.505/0.507**3
def dur_days(P,rho): return 13/24*(P/365.25)**(1/3)*rho**(-1/3)
base=t.max()-t.min(); pmin,pmax=0.5,30.0
freqs=[];fq=1/pmin
while fq>1/pmax:
    freqs.append(fq); P=1/fq; fq-=0.5*dur_days(P,rho)/(3*base*P)
periods=np.sort(1/np.array(freqs)); durs=np.geomspace(0.03,0.25,8)
print('BLS grid: %d periods, %d durations'%(len(periods),len(durs)))
bls=BoxLeastSquares(t,y,dy)
t1=time.time(); r=bls.power(periods,durs,objective='snr',method='fast',oversample=5); dt_snr=time.time()-t1
k=int(np.argmax(r.power)); print('BLS objective=snr: %.1fs best P=%.5f d, t0=%.4f, dur=%.2f h, depth=%.0f ppm, power=%.2f'%(dt_snr,r.period[k],r.transit_time[k],r.duration[k]*24,r.depth[k]*1e6,r.power[k]))
pw=np.asarray(r.power); sde=(pw-np.median(pw))/np.std(pw); print('BLS SDE (global) of top peak %.1f; 2nd-best period %.4f'%(sde[k], r.period[np.argsort(pw)[-2]]))
t1=time.time(); r2=bls.power(periods,durs,objective='likelihood',method='fast',oversample=5); dt_ll=time.time()-t1
k2=int(np.argmax(r2.power)); print('BLS objective=likelihood: %.1fs best P=%.5f, depth %.0f ppm, dloglike %.1f -> sqrt(2*dll)=%.1f'%(dt_ll,r2.period[k2],r2.depth[k2]*1e6,r2.power[k2],np.sqrt(2*r2.power[k2])))
st=bls.compute_stats(r.period[k],r.duration[k],r.transit_time[k])
print('compute_stats keys',sorted(st.keys()))
print('depth',st['depth'],'depth_odd',st['depth_odd'],'depth_even',st['depth_even'],'depth_half',st['depth_half'],'depth_phased',st['depth_phased'],'harmonic_dll',st['harmonic_delta_log_likelihood'],'n transits with data',int(np.sum(st['per_transit_count']>0)))
# TLS on the same cleaned data, 1 thread, 1-30 d
from transitleastsquares import transitleastsquares, catalog_info
try:
    ab,mass,mass_min,mass_max,radius,radius_min,radius_max=catalog_info(KIC_ID=5689121)
    print('catalog_info(KIC_ID=5689121): u=',ab,'M=',mass,'(',mass_min,mass_max,') R=',radius,'(',radius_min,radius_max,')')
except Exception as ex: print('catalog_info failed:',repr(ex)); ab=(0.4804,0.1867)
model=transitleastsquares(t,y+1,dy)
t1=time.time()
res=model.power(period_min=1,period_max=30,R_star=0.507,M_star=0.505,R_star_min=0.3,R_star_max=0.8,M_star_min=0.3,M_star_max=0.8,u=list(ab),use_threads=1,show_progress_bar=False,verbose=False)
print('TLS 1 thread: %.1fs, %d periods; P=%.6f +- %.6f d, T0=%.5f, SDE=%.2f, SNR=%.2f, depth=%.5f (1-depth=%.0f ppm), dur=%.2f h, rp_rs=%.4f, odd_even_mismatch=%.2f sigma, transit_count=%d, distinct=%d, FAP=%s'%(time.time()-t1,len(res.periods),res.period,res.period_uncertainty,res.T0,res.SDE,res.snr,res.depth,(1-res.depth)*1e6,res.duration*24,res.rp_rs,res.odd_even_mismatch,res.transit_count,res.distinct_transit_count,res.FAP))
print('TLS top 5 peaks:',[(round(p,4),round(s,1)) for p,s in zip(res.periods[np.argsort(res.power)[::-1][:5]],np.sort(res.power)[::-1][:5])])
json.dump(dict(bls_snr_seconds=dt_snr,bls_ll_seconds=dt_ll,n_periods=len(periods),bls_P=float(r.period[k]),bls_depth_ppm=float(r.depth[k]*1e6),tls_P=float(res.period),tls_SDE=float(res.SDE),tls_snr=float(res.snr),tls_depth_ppm=float((1-res.depth)*1e6),tls_T0=float(res.T0),tls_dur_h=float(res.duration*24)),open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/e2e_result.json','w'),indent=1)
