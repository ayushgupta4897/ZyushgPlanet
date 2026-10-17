import numpy as np, json, pandas as pd
from astropy.io import fits
from astropy.timeseries import LombScargle, BoxLeastSquares
np.random.seed(1)
F='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-009655224_kepler_v1.0_lc.fits'
P=6.703464464595031; t0=265.183161835447; dur=2.4375/24
d=fits.getdata(F,'LIGHTCURVE_STITCHED')
t=d['TIME']; f=d['FLUX']; q=d['QUARTER']; fs=d['SAP_FLUX']; qual=d['SAP_QUALITY']
ok=np.isfinite(t)&np.isfinite(f)&(qual==0)&(f>0)
t,f,q,fs=t[ok],f[ok],q[ok],fs[ok]
ph=((t-t0)/P+0.5)%1-0.5; intr=np.abs(ph*P)<dur/2*1.3
def runmed(t,f,mask,win=1.0):
    out=np.empty_like(f)
    for i in range(len(t)):
        s=(np.abs(t-t[i])<win/2)&~mask
        out[i]=np.median(f[s]) if s.sum()>5 else np.nan
    return out
# per-quarter normalise & detrend (transit-masked)
fn=np.full_like(f,np.nan); fsn=np.full_like(f,np.nan)
for qq in np.unique(q):
    s=q==qq
    tr=runmed(t[s],f[s],intr[s]); fn[s]=f[s]/tr-1
    tr2=runmed(t[s],fs[s],intr[s]); fsn[s]=fs[s]/tr2-1
g=np.isfinite(fn); t,fn,fsn,q,ph,intr,f=t[g],fn[g],fsn[g],q[g],ph[g],intr[g],f[g]
# clip
mad=1.4826*np.median(np.abs(fn-np.median(fn))); g=(fn<4*mad)&(fn>-8*mad); t,fn,fsn,q,ph,intr,f=t[g],fn[g],fsn[g],q[g],ph[g],intr[g],f[g]
res={}
inbox=np.abs(ph*P)<dur/2
def boxdepth(fl,inb,nrand=800):
    dep=np.median(fl[~inb])-np.mean(fl[inb])
    rs=[]; 
    for k in range(nrand):
        p2=((t-t0-np.random.uniform(0,P))/P+0.5)%1-0.5; ib=np.abs(p2*P)<dur/2
        if ib.sum()>10: rs.append(np.median(fl[~ib])-np.mean(fl[ib]))
    return dep*1e6, np.std(rs)*1e6
dep,err=boxdepth(fn,inbox); res['depth_ppm']=[dep,err]; print('global depth %.0f +/- %.0f ppm'%(dep,err))
# per quarter / per season
season={3:1,7:1,11:1,15:1,4:2,8:2,12:2,5:3,9:3,13:3,6:0,10:0,14:0}
pq={}
for qq in np.unique(q):
    s=q==qq; dq=np.median(fn[s&~inbox])-np.mean(fn[s&inbox]); eq=np.std(fn[s&~inbox])/np.sqrt((s&inbox).sum())
    pq[int(qq)]=[dq*1e6,eq*1e6,int((s&inbox).sum()), season[int(qq)]]
    print('Q%2d season %d: depth %6.0f +/- %4.0f ppm  n_in=%d'%(qq,season[int(qq)],dq*1e6,eq*1e6,(s&inbox).sum()))
res['per_quarter']=pq
ps={}
for se in [0,1,2,3]:
    s=np.isin(q,[k for k,v in season.items() if v==se]); 
    if s.sum()==0: continue
    dq=np.median(fn[s&~inbox])-np.mean(fn[s&inbox]); eq=np.std(fn[s&~inbox])/np.sqrt((s&inbox).sum()); ps[se]=[dq*1e6,eq*1e6]
    print('season %d (ch %s): depth %6.0f +/- %4.0f ppm'%(se,{1:23,2:39,3:63,0:47}[se],dq*1e6,eq*1e6))
res['per_season']=ps
vals=np.array([v[0] for v in ps.values()]); errs=np.array([v[1] for v in ps.values()]); w=1/errs**2; m=np.sum(w*vals)/np.sum(w)
print('season chi2 about weighted mean %.0f: %.2f for %d dof'%(m,np.sum(w*(vals-m)**2),len(vals)-1)); res['season_chi2']=[float(np.sum(w*(vals-m)**2)),len(vals)-1]
# time trend of per-event depths
ev=pd.read_csv('/home/user/ZyushgPlanet/runs/poolA_main/vet/kic-009655224/events.csv')
ev=ev[ev.n_in>=3]
x=ev.tc_bkjd.values; y=ev.depth.values*1e6; e=ev.err.values*1e6
A=np.vstack([np.ones_like(x),(x-x.mean())/1000]).T; W=np.diag(1/e**2)
cov=np.linalg.inv(A.T@W@A); b=cov@A.T@W@y
print('per-event depth trend: %.0f +/- %.0f ppm per 1000 d (%.1f sigma)'%(b[1],np.sqrt(cov[1,1]),b[1]/np.sqrt(cov[1,1]))); res['trend_ppm_per_1000d']=[b[1],np.sqrt(cov[1,1])]
# gap proximity
for cut in [0.5,1,2]:
    s=ev.gap_dist_d.values>cut; w=1/e[s]**2; print('events with gap_dist>%.1f d: n=%d depth %.0f +/- %.0f'%(cut,s.sum(),np.sum(w*y[s])/np.sum(w),1/np.sqrt(np.sum(w))))
    s=~s; w=1/e[s]**2; 
    if s.sum()>1: print('   near gap (<%.1f d): n=%d depth %.0f +/- %.0f'%(cut,s.sum(),np.sum(w*y[s])/np.sum(w),1/np.sqrt(np.sum(w))))
from scipy.stats import spearmanr
print('spearman depth vs gap_dist:',spearmanr(ev.gap_dist_d,ev.depth))
z=(y-np.average(y,weights=1/e**2))/e; print('per-event outliers |z|>3:',int((np.abs(z)>3).sum()),'of',len(z),'; max z %.1f at tc %.2f'%(z[np.argmax(np.abs(z))],x[np.argmax(np.abs(z))]))
# depth excluding top-5 events
o=np.argsort(-y/e)[:5]; s=np.ones(len(y),bool); s[o]=False; w=1/e[s]**2; print('depth excluding 5 highest-SNR events: %.0f +/- %.0f'%(np.sum(w*y[s])/np.sum(w),1/np.sqrt(np.sum(w))))
# cadence-level: fraction of in-transit signal from deepest 5% cadences
ib=fn[inbox]; k=int(0.05*len(ib)); srt=np.sort(ib); print('in-transit cadences: n=%d; mean %.0f ppm; mean excluding deepest 5%%: %.0f ppm; excluding deepest 10%%: %.0f'%(len(ib),ib.mean()*1e6,srt[k:].mean()*1e6,srt[2*k:].mean()*1e6))
# in-transit vs out-of-transit distribution: median instead of mean
print('median in-box - median out: %.0f ppm'%((np.median(fn[~inbox])-np.median(ib))*1e6))
# Lomb-Scargle on normalised (undetrended) flux per-quarter median-normalised
fraw=np.full_like(f,np.nan)
for qq in np.unique(q):
    s=q==qq; fraw[s]=f[s]/np.median(f[s])-1
freq=np.linspace(1/100,1/0.2,200000); ls=LombScargle(t,fraw); pw=ls.power(freq)
i=np.argmax(pw); print('LS peak: P=%.4f d power %.4g amp %.0f ppm'%(1/freq[i],pw[i],np.sqrt(pw[i]*4*np.var(fraw))*1e6*0.5))
for Pt in [P,P/2,2*P,P/3]:
    j=np.argmin(np.abs(1/freq-Pt)); print('  LS power at %.3f d: %.4g (peak %.4g); FAP-ish rank %d'%(Pt,pw[j],pw[i],int((pw>pw[j]).sum())))
res['ls_peak']=[1/freq[i],float(pw[i])]
# LS on detrended flux at P window  
# V-shape: fit trapezoid on binned fold
phh=ph*P*24; s=np.abs(phh)<6; xb=phh[s]; yb=fn[s]
def trap(x,dep,T14,T23):
    y=np.zeros_like(x); a=np.abs(x)
    y[a<=T23/2]=-dep; r=(a>T23/2)&(a<T14/2); y[r]=-dep*(T14/2-a[r])/((T14-T23)/2); return y
best=[]
for tin_frac in [0.0,0.1,0.2,0.3,0.4,0.5]:
    for T14 in np.arange(1.6,3.6,0.1):
        T23=T14*(1-2*tin_frac)
        for dep in np.arange(200,1200,25)*1e-6:
            c=np.sum((yb-trap(xb,dep,T14,T23))**2); best.append((c,tin_frac,T14,dep))
best=np.array(best); c0=best[:,0].min(); sig2=np.var(yb)
print('trapezoid grid: best', best[np.argmin(best[:,0])], )
for tf in [0.0,0.1,0.2,0.3,0.4,0.5]:
    m=best[best[:,1]==tf]; print('  ingress frac %.1f: min dchi2 vs global %.2f'%(tf,(m[:,0].min()-c0)/sig2))
cf=np.sum((yb-0)**2); print('  flat model dchi2 vs best: %.1f'%((cf-c0)/sig2))
# secondary scan all phases (box of same dur), max |depth|
scan=[]
for pc in np.arange(0.1,0.9,0.01):
    p2=((t-t0)/P-pc+0.5)%1-0.5; ib=np.abs(p2*P)<dur/2; scan.append((pc,(np.median(fn[~ib])-np.mean(fn[ib]))*1e6))
scan=np.array(scan); j=np.argmax(scan[:,1]); print('deepest box at other phases: %.0f ppm at phase %.2f (err ~%.0f)'%(scan[j,1],scan[j,0],err)); res['secondary_scan_max']=[scan[j,1],scan[j,0]]
# inverted-flux BLS narrow window
bls=BoxLeastSquares(t,fn); per=np.linspace(6.69,6.72,3000)
r=bls.power(per,dur,objective='likelihood'); ri=BoxLeastSquares(t,-fn).power(per,dur,objective='likelihood')
print('narrow BLS dlogL real %.1f vs inverted %.1f'%(r.power.max(),ri.power.max()))
# depth in PSF vs SAP in e/s
sapd,sape=boxdepth(fsn,inbox,300); print('SAP depth %.0f +/- %.0f ppm'%(sapd,sape))
res['sap_depth_ppm']=[sapd,sape]
json.dump(res,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic009655224/work/lc_checks.json','w'),indent=1,default=float)
