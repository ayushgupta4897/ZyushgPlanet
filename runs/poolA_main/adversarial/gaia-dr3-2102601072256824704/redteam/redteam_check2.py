import numpy as np, json
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from astropy.time import Time
from wotan import flatten
F='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2102601072256824704_kepler_v1.0_lc.fits'
P=32.55706418774713; t0=148.82345687868073
d=fits.open(F)['LIGHTCURVE_STITCHED'].data
m=(d['SAP_QUALITY']==0)&np.isfinite(d['TIME'])&np.isfinite(d['FLUX'])
t=d['TIME'][m].astype(float); f=d['FLUX'][m].astype(float); q=d['QUARTER'][m]
for qq in np.unique(q): k=q==qq; f[k]/=np.nanmedian(f[k])
gap=np.where(np.diff(t)>0.5)[0]
print('gaps bracketing the three events (BKJD, UTC date, gap length d):')
for tc in [181.42,539.53,1418.53]:
    i=np.searchsorted(t,tc)
    # nearest gap before and after
    gb=[g for g in gap if t[g]<=tc]; ga=[g for g in gap if t[g+1]>=tc]
    s=''
    if gb: g=gb[-1]; s+=f' gap ending {t[g+1]:.3f} ({Time(t[g+1]+2454833,format="jd").iso[:10]}), len {t[g+1]-t[g]:.2f} d, {tc-t[g+1]:.3f} d before event |'
    if ga: g=ga[0]; s+=f' gap starting {t[g]:.3f} ({Time(t[g]+2454833,format="jd").iso[:10]}), len {t[g+1]-t[g]:.2f} d, {t[g]-tc:.3f} d after event'
    print(f' tc={tc} ({Time(tc+2454833,format="jd").iso[:10]}) Q{q[i]}:'+s)
res={}
for lab,ec in [('edge0',0.0),('edge0.3',0.3)]:
    fl=flatten(t,f,method='biweight',window_length=1.0,break_tolerance=0.5,edge_cutoff=ec); ok=np.isfinite(fl); tt=t[ok]; ff=fl[ok]
    sd=1.4826*np.median(np.abs(ff-1)); k2=(ff-1<4*sd)&(ff-1>-10*sd); tt=tt[k2]; ff=ff[k2]; dy=np.full_like(ff,sd)
    periods=np.exp(np.linspace(np.log(28),np.log(37),15000))
    for pol,yy in [('real',ff),('inv',2-ff)]:
        r=BoxLeastSquares(tt,yy,dy).power(periods,[0.10,0.14,0.18,0.25],objective='likelihood')
        pk=r.power; i=np.argmax(pk); kP=np.argmin(np.abs(periods-P))
        sde=lambda j:(pk[j]-np.median(pk))/np.std(pk)
        print(f'BLS {lab} {pol}: top P {periods[i]:.4f} SDE {sde(i):.2f} depth {r.depth[i]*1e6:.0f} ppm | at P_pipe SDE {sde(kP):.2f} depth {r.depth[kP]*1e6:.0f}+-{r.depth_err[kP]*1e6:.0f} ppm (formal)')
        res[f'{lab}_{pol}']=dict(topP=float(periods[i]),topSDE=float(sde(i)),SDE_at_P=float(sde(kP)),depth_at_P=float(r.depth[kP]*1e6))
# host consistency
G=6.674e-11; Msun=1.989e30; Rsun=6.957e8; AU=1.496e11
for R,M,lab in [(0.87,0.95,'Gaia-based R*=0.87 (primary only)'),(0.99,1.0,'KIC R*=0.99 (blend)')]:
    a=(G*M*Msun*(P*86400)**2/(4*np.pi**2))**(1/3)
    T14=P*24/np.pi*np.arcsin(R*Rsun/a)
    print(f'{lab}: a/R*={a/(R*Rsun):.1f}, central T14 = {T14:.2f} h (observed 3.4-3.6 h -> b~{np.sqrt(max(0,1-(3.5/T14)**2)):.2f})')
fr_pri=10**(-0.4*15.276)/(10**(-0.4*15.276)+10**(-0.4*15.762)); print('primary flux fraction of the 0.74" pair:',round(fr_pri,3))
for dep in [147,62]:
    print(f' observed {dep} ppm: if on primary true depth {dep/fr_pri:.0f} ppm -> Rp {np.sqrt(dep/fr_pri*1e-6)*0.87*109.1:.2f} Re; if on companion {dep/(1-fr_pri):.0f} ppm -> Rp {np.sqrt(dep/(1-fr_pri)*1e-6)*0.80*109.1:.2f} Re')
json.dump(res,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-2102601072256824704/redteam_check2.json','w'),indent=1)
