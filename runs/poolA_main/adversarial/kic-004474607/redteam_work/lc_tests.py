"""Red-team light-curve tests for kic-004474607 (no pipeline/ code imported).
Outputs lc_tests.json + figures in the same directory."""
import json, sys, os, numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT=os.path.dirname(os.path.abspath(__file__))
F='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-004474607_kepler_v1.0_lc.fits'
P_PIPE,T0_PIPE,DUR_PIPE=32.94271710658947,754.5746921888925,1.95/24
P_RED,T0_RED=32.98173,754.1131
res={}
h=fits.open(F); d=h['LIGHTCURVE_STITCHED'].data
m=(d['SAP_QUALITY']==0)&np.isfinite(d['FLUX'])&np.isfinite(d['TIME'])
t=d['TIME'][m].astype(float); f=d['FLUX'][m].astype(float); fs=d['SAP_FLUX'][m].astype(float); q=d['QUARTER'][m]
# per-quarter normalise
fn=np.empty_like(f); fsn=np.empty_like(f)
for qq in np.unique(q):
    s=q==qq; fn[s]=f[s]/np.nanmedian(f[s]); fsn[s]=fs[s]/np.nanmedian(fs[s])
res['n_cad']=int(len(t)); res['quarters']=[int(x) for x in np.unique(q)]
# ---- gap census: starts of segments after gaps > 0.5 d, plus quarter starts
dt=np.diff(t); starts=[0]+list(np.where(dt>0.5)[0]+1)
gaps=[]
for i in starts:
    t_start=t[i]
    # ramp: median of first 0.3 d after start vs median of 1-3 d after
    s1=(t>=t_start)&(t<t_start+0.3); s2=(t>=t_start+1.0)&(t<t_start+3.0)
    if s1.sum()>5 and s2.sum()>20:
        ramp=(np.median(fn[s1])-np.median(fn[s2]))*1e6
    else: ramp=np.nan
    gaplen=float(t[i]-t[i-1]) if i>0 else np.nan
    gaps.append(dict(t_start=float(t_start),quarter=int(q[i]),gap_len_d=gaplen,ramp_first0p3d_ppm=float(ramp)))
res['segment_starts']=gaps
# ---- epochs of pipeline ephemeris and distance to nearest segment start (after it)
def epochs(P,t0):
    n0=int(np.floor((t.min()-t0)/P)); n1=int(np.ceil((t.max()-t0)/P))
    out=[]
    for n in range(n0,n1+1):
        tc=t0+n*P
        if (np.abs(t-tc)<DUR_PIPE/2).sum()>=2: out.append((n,tc))
    return out
ev=epochs(P_PIPE,T0_PIPE)
starts_t=np.array([g['t_start'] for g in gaps])
# ---- detrend: running median, masked transit, with/without post-gap mask
def running_median(t,y,win=1.0,mask=None):
    tr=np.empty_like(y)
    for i in range(len(t)):
        s=(np.abs(t-t[i])<win/2)
        if mask is not None: s&=~mask
        tr[i]=np.median(y[s]) if s.sum()>=5 else np.nan
    return tr
def seg_detrend(t,y,mask_tr,win=1.0):
    tr=np.full_like(y,np.nan)
    for i,i0 in enumerate(starts):
        i1=starts[i+1] if i+1<len(starts) else len(t)
        sl=slice(i0,i1)
        tr[sl]=running_median(t[sl],y[sl],win,mask_tr[sl])
    return y/tr-1
intr=np.zeros(len(t),bool)
for n,tc in ev: intr|=np.abs(t-tc)<1.5*DUR_PIPE
postgap=np.zeros(len(t),bool)
for ts in starts_t: postgap|=(t>=ts)&(t<ts+0.5)
res['n_postgap_cadences_0p5d']=int(postgap.sum())
r_all=seg_detrend(t,fn,intr)
r_sap=seg_detrend(t,fsn,intr)
good=np.isfinite(r_all)
# ---- per-event local depths (box +-dur/2 in, 0.3-1.0 d out, linear local baseline) on raw normalised flux
def local_depth(tt,yy,tc,dur):
    s_in=np.abs(tt-tc)<dur/2; s_out=(np.abs(tt-tc)>dur/2+0.05)&(np.abs(tt-tc)<1.0)
    if s_in.sum()<2 or s_out.sum()<10: return np.nan,np.nan,int(s_in.sum())
    c=np.polyfit(tt[s_out]-tc,yy[s_out],1); base=np.polyval(c,tt-tc)
    rr=yy/base-1; sig=1.4826*np.median(np.abs(rr[s_out]-np.median(rr[s_out])))
    return float(-np.mean(rr[s_in])*1e6), float(sig/np.sqrt(s_in.sum())*1e6), int(s_in.sum())
events=[]
for n,tc in ev:
    dd,ee,ni=local_depth(t,fn,tc,DUR_PIPE)
    dgap=float(np.min(np.where(starts_t<=tc, tc-starts_t, np.inf)))
    events.append(dict(epoch=n,tc=float(tc),depth_ppm=dd,err_ppm=ee,n_in=ni,d_after_segment_start_d=dgap,in_postgap_0p5d=bool(dgap<0.5)))
res['events_pipeline_ephemeris']=events
def wmean(ev):
    d=np.array([e['depth_ppm'] for e in ev]); e_=np.array([e['err_ppm'] for e in ev]); ok=np.isfinite(d)&np.isfinite(e_)
    w=1/e_[ok]**2; return float(np.sum(w*d[ok])/w.sum()), float(1/np.sqrt(w.sum())), int(ok.sum())
res['stack_all']=wmean(events); res['stack_excl_postgap']=wmean([e for e in events if not e['in_postgap_0p5d']])
res['stack_postgap_only']=wmean([e for e in events if e['in_postgap_0p5d']])
# ---- box depth on detrended flux at pipeline ephemeris, with random-phase error; with and without post-gap mask
def box_depth(tt,rr,P,t0,dur,mask):
    ph=((tt-t0+P/2)%P)-P/2; s=(np.abs(ph)<dur/2)&mask
    dep=-np.mean(rr[s])*1e6
    rng=np.random.default_rng(1); fake=[]
    for k in range(2000):
        t0r=t0+rng.uniform(0.1,0.9)*P
        phr=((tt-t0r+P/2)%P)-P/2; sr=(np.abs(phr)<dur/2)&mask
        if sr.sum()>=0.5*s.sum(): fake.append(-np.mean(rr[sr])*1e6)
    return float(dep), float(np.std(fake)), int(s.sum())
res['box_depth_psf_all']=box_depth(t,r_all,P_PIPE,T0_PIPE,DUR_PIPE,good)
res['box_depth_psf_excl_postgap']=box_depth(t,r_all,P_PIPE,T0_PIPE,DUR_PIPE,good&~postgap)
res['box_depth_sap_all']=box_depth(t,r_sap,P_PIPE,T0_PIPE,DUR_PIPE,np.isfinite(r_sap))
res['box_depth_sap_excl_postgap']=box_depth(t,r_sap,P_PIPE,T0_PIPE,DUR_PIPE,np.isfinite(r_sap)&~postgap)
# ---- control: box placed 0.14 d after EVERY segment start (same offset as the two Q13 events)
ctrl=[]
for g in gaps:
    tc=g['t_start']+0.14
    dd,ee,ni=local_depth(t,fn,tc,DUR_PIPE)
    ctrl.append(dict(t_start=g['t_start'],quarter=g['quarter'],gap_len_d=g['gap_len_d'],local_depth_ppm=dd,err_ppm=ee,n_in=ni))
res['control_box_0p14d_after_every_segment_start']=ctrl
# ---- BLS on gap-masked, transit-unmasked re-detrend (fair blind search), real vs inverted, 25-40 d narrow window
r_blind=seg_detrend(t,fn,np.zeros(len(t),bool))
gb=np.isfinite(r_blind)
def bls(tt,rr,pmin,pmax):
    base=tt.max()-tt.min()
    freqs=np.arange(1/pmax,1/pmin, 0.06/24/base/3)   # ~ dur/(3 baseline)
    per=np.sort(1/freqs)
    b=BoxLeastSquares(tt,rr+1,dy=None)
    durs=np.array([1.5,2.0,3.0,4.5])/24
    r=b.power(per,durs,objective='likelihood')
    pw=r.power; sde=(pw.max()-pw.mean())/pw.std()
    i=np.argmax(pw)
    return dict(period=float(per[i]),power=float(pw[i]),sde=float(sde),depth_ppm=float(r.depth[i]*1e6),duration_h=float(r.duration[i]*24),t0=float(r.transit_time[i]),n_periods=int(len(per))), per, pw
out_bls={}
for label,mask in [('unmasked',gb),('postgap_masked',gb&~postgap)]:
    for pol,sign in [('real',1),('inverted',-1)]:
        rres,per,pw=bls(t[mask],sign*r_blind[mask],25,40)
        out_bls[f'{label}_{pol}_25_40d']=rres
        # power at 32.94 +- 0.05
        s=np.abs(per-P_PIPE)<0.05; rres['power_near_32.94']=float(pw[s].max()); rres['sde_near_32.94']=float((pw[s].max()-pw.mean())/pw.std())
        print(label,pol,rres,flush=True)
res['bls_25_40d']=out_bls
# ---- duration / density
Rs=0.8721514356136322; Ms=0.95
G=6.674e-11; Msun=1.989e30; Rsun=6.957e8; Pd=P_PIPE*86400
a=(G*Ms*Msun*Pd**2/(4*np.pi**2))**(1/3)
T14_b0_h=Pd/np.pi*np.arcsin(Rs*Rsun/a)/3600
res['duration_density']=dict(Rstar=Rs,Mstar_assumed=Ms,a_AU=a/1.496e11,T14_b0_h=float(T14_b0_h),T14_obs_h=1.95,
   b_required=float(np.sqrt(max(0,1-(1.95/T14_b0_h)**2))),
   rp_earth_for_1100ppm=float(np.sqrt(1100e-6)*Rs*109.1), rp_earth_for_2141ppm=float(np.sqrt(2141e-6)*Rs*109.1),
   rp_earth_for_600ppm=float(np.sqrt(600e-6)*Rs*109.1))
json.dump(res,open(os.path.join(OUT,'lc_tests.json'),'w'),indent=1,default=float)
# ---- figure: raw Q13 flux around the two events and the two segment starts
fig,ax=plt.subplots(2,2,figsize=(14,7))
for k,(tc,lab) in enumerate([(T0_PIPE+13*P_PIPE,'epoch 13'),(T0_PIPE+14*P_PIPE,'epoch 14')]):
    s=(t>tc-2)&(t<tc+3)
    ax[0,k].plot(t[s],(fn[s]-1)*1e6,'.',ms=3); ax[0,k].axvspan(tc-DUR_PIPE/2,tc+DUR_PIPE/2,color='r',alpha=.3); ax[0,k].set_title(f'{lab} raw PSF flux (per-quarter normalised)'); ax[0,k].set_ylabel('ppm')
    ax[1,k].plot(t[s],(fsn[s]-1)*1e6,'.',ms=3,color='g'); ax[1,k].axvspan(tc-DUR_PIPE/2,tc+DUR_PIPE/2,color='r',alpha=.3); ax[1,k].set_title(f'{lab} raw SAP flux'); ax[1,k].set_xlabel('BKJD')
    for g in gaps:
        if tc-2<g['t_start']<tc+3:
            for a_ in ax[:,k]: a_.axvline(g['t_start'],ls=':',color='k')
plt.tight_layout(); plt.savefig(os.path.join(OUT,'q13_events_raw.png'),dpi=110)
print('done')
