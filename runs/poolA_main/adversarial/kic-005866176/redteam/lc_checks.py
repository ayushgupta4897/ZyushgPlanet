"""Red-team light-curve checks for kic-005866176 (own code; nothing imported from pipeline/).
Outputs lc_checks.json + lc_checks.png next to this file."""
import json, sys, numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares, LombScargle
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rng=np.random.default_rng(42)
OUT='/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005866176/redteam/'
LC='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005866176_kepler_v1.0_lc.fits'
HOST='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam_kic005866176/hostlc/hlsp_kbonus-bkg_kepler_kepler_kic-005866192_kepler_v1.0_lc.fits'
P=29.433016127730333; T0=172.57296158611774; DUR=1.95/24.
R={}
def load(path):
    f=fits.open(path); d=f[1].data
    cols=d.columns.names
    t=d['TIME']; fl=d['FLUX']; sap=d['SAP_FLUX']; q=d['SAP_QUALITY']; qu=d['QUARTER']
    m=np.isfinite(t)&np.isfinite(fl)&(q==0)
    return t[m],fl[m],sap[m],qu[m],cols,f
t,fl,sap,qu,cols,f=load(LC); R['columns']=list(cols); R['n_good']=int(len(t))
print('columns',cols)
# per-quarter raw medians (un-rescaled) + channel
raw={}
for hdu in f:
    if hdu.name.startswith('LIGHTCURVE_Q'):
        qn=int(hdu.name[12:]); dd=hdu.data; h=hdu.header
        mm=np.isfinite(dd['FLUX'])&(dd['SAP_QUALITY']==0) if 'SAP_QUALITY' in dd.columns.names else np.isfinite(dd['FLUX'])
        raw[qn]={'channel':h['CHANNEL'],'season':h['SEASON'],'psffrac':h['PSFFRAC'],'crowdsap':h['CROWDSAP'],'flfrcsap':h['FLFRCSAP'],'npixsap':h['NPIXSAP'],
                 'median_psf_raw':float(np.nanmedian(dd['FLUX'][mm])),'median_sap_raw':float(np.nanmedian(dd['SAP_FLUX'][mm])) if 'SAP_FLUX' in dd.columns.names else None,
                 'cols':dd.columns.names}
R['per_quarter_raw']=raw
# ---- normalise per quarter, masked running median detrend (1 d), gap-split
def phase(tt,P=P,T0=T0): return ((tt-T0)/P+0.5)%1-0.5
def detrend(t,y,win=1.0,mask_dur=1.5*DUR,P=P,T0=T0):
    y=y.copy()
    for qq in np.unique(qu):
        s=qu==qq; y[s]=y[s]/np.nanmedian(y[s])
    intr=np.abs(phase(t,P,T0)*P)<mask_dur
    # gap-split
    edges=np.where(np.diff(t)>0.5)[0]; segs=np.split(np.arange(len(t)),edges+1)
    trend=np.ones_like(y)
    for sg in segs:
        ts=t[sg]; ys=y[sg]; ms=intr[sg]
        for i,idx in enumerate(sg):
            w=(np.abs(ts-ts[i])<win/2)&(~ms)
            if w.sum()<5: w=np.abs(ts-ts[i])<win/2
            trend[idx]=np.nanmedian(ys[w])
    r=y/trend-1
    # clip upward outliers 4 sigma, downward 12 sigma
    s=1.4826*np.nanmedian(np.abs(r-np.nanmedian(r)))
    keep=(r<4*s)&(r>-12*s)
    return r,keep,s
r,keep,sig=detrend(t,fl); rs,keeps,sigs=detrend(t,sap)
R['robust_scatter_ppm']=float(sig*1e6); R['robust_scatter_sap_ppm']=float(sigs*1e6)
tk,rk=t[keep],r[keep]; rsk=rs[keep]; quk=qu[keep]
# ---- gaps
dt=np.diff(tk); gi=np.where(dt>0.3)[0]
gaps=[(float(tk[i]),float(tk[i+1])) for i in gi]
qstart=[float(tk[0])]+[float(tk[i+1]) for i in gi]
R['gaps']=gaps; R['n_gaps_gt0p3d']=len(gaps)
def time_since_gap(tt):
    # positive = time since the last gap END (post-gap), for each cadence; also time to next gap start
    starts=np.array([g[1] for g in gaps]+[tk[0]]); ends=np.array([g[0] for g in gaps]+[tk[-1]])
    since=np.array([ (x-starts[starts<=x]).min() if (starts<=x).any() else 99 for x in tt])
    until=np.array([ (ends[ends>=x]-x).min() if (ends>=x).any() else 99 for x in tt])
    return since,until
since,until=time_since_gap(tk)
# ---- mechanism test: mean detrended flux as a function of time since gap (all gaps, all data)
bins=[0,0.05,0.1,0.15,0.2,0.3,0.5,0.75,1.0,1.5,2.0,3.0]
prof=[]
for a,b in zip(bins[:-1],bins[1:]):
    w=(since>=a)&(since<b)
    prof.append({'t_since_gap_d':[a,b],'n':int(w.sum()),'mean_ppm':float(rk[w].mean()*1e6),'sem_ppm':float(rk[w].std()/np.sqrt(w.sum())*1e6),
                 'median_ppm':float(np.median(rk[w])*1e6),'sap_mean_ppm':float(rsk[w].mean()*1e6)})
R['post_gap_profile']=prof
# Also pre-gap profile
prof2=[]
for a,b in zip(bins[:-1],bins[1:]):
    w=(until>=a)&(until<b)
    prof2.append({'t_until_gap_d':[a,b],'n':int(w.sum()),'mean_ppm':float(rk[w].mean()*1e6),'sem_ppm':float(rk[w].std()/np.sqrt(w.sum())*1e6)})
R['pre_gap_profile']=prof2
# per-gap post-gap mean (first 0.15 d) to see which gaps ramp
pg=[]
for g in gaps:
    w=(tk>g[1])&(tk<g[1]+0.15)
    if w.sum()>0: pg.append({'gap_end':g[1],'gap_len_d':g[1]-g[0],'n':int(w.sum()),'mean_ppm_first0p15d':float(rk[w].mean()*1e6)})
R['per_gap_first0p15d']=pg
# ---- per-event depths (local)
def event_table(tt,rr,P=P,T0=T0,dur=DUR):
    n0=int(np.floor((tt.min()-T0)/P))-1; n1=int(np.ceil((tt.max()-T0)/P))+1
    ev=[]
    for n in range(n0,n1+1):
        tc=T0+n*P
        win=(np.abs(tt-tc)<1.0)
        if win.sum()<10: continue
        inn=np.abs(tt-tc)<dur/2; out=win&(np.abs(tt-tc)>1.5*dur)
        if inn.sum()<2 or out.sum()<8: continue
        d=rr[out].mean()-rr[inn].mean()
        e=rr[out].std()*np.sqrt(1/inn.sum()+1/out.sum())
        s,u=time_since_gap(np.array([tc]))
        ev.append({'n':n,'tc':float(tc),'depth_ppm':float(d*1e6),'err_ppm':float(e*1e6),'n_in':int(inn.sum()),'n_out':int(out.sum()),
                   'since_gap_d':float(s[0]),'until_gap_d':float(u[0]),'quarter':int(quk[np.argmin(np.abs(tt-tc))]) })
    return ev
ev=event_table(tk,rk); R['events']=ev
d=np.array([e['depth_ppm'] for e in ev]); e_=np.array([e['err_ppm'] for e in ev]); sg=np.array([e['since_gap_d'] for e in ev])
w=1/e_**2; dw=(d*w).sum()/w.sum(); chi2=(((d-dw)/e_)**2).sum()
R['events_summary']={'n':len(ev),'weighted_mean_ppm':float(dw),'chi2_const':float(chi2),'dof':len(ev)-1,'unweighted_mean_ppm':float(d.mean()),'median_ppm':float(np.median(d))}
near=sg<0.15
R['events_postgap_split']={'n_postgap(<0.15d)':int(near.sum()),'postgap_wmean_ppm':float((d[near]*w[near]).sum()/w[near].sum()) if near.any() else None,
   'others_wmean_ppm':float((d[~near]*w[~near]).sum()/w[~near].sum()),'others_chi2':float((((d[~near]-(d[~near]*w[~near]).sum()/w[~near].sum())/e_[~near])**2).sum()),'others_dof':int((~near).sum()-1),
   'others_unweighted_mean_ppm':float(d[~near].mean()),'others_sem_ppm':float(d[~near].std()/np.sqrt((~near).sum()))}
# ---- box depth with random-phase errors, with/without post-gap excision
def boxdepth(tt,rr,P,T0,dur,ntrial=400):
    ph=phase(tt,P,T0)*P
    inn=np.abs(ph)<dur/2
    dep=-(rr[inn].mean())  # rr is detrended to ~0 outside; use out mean too
    out=np.abs(ph)>1.5*dur
    dep=rr[out].mean()-rr[inn].mean()
    tr=[]
    for k in range(ntrial):
        t0r=T0+rng.uniform(0.05,0.95)*P
        phr=phase(tt,P,t0r)*P; i2=np.abs(phr)<dur/2
        if i2.sum()>5: tr.append(rr[out].mean()-rr[i2].mean())
    tr=np.array(tr)
    return float(dep*1e6),float(tr.std()*1e6),int(inn.sum()),tr
dep,err,nin,trials=boxdepth(tk,rk,P,T0,DUR)
R['box_depth_all']={'depth_ppm':dep,'err_random_phase_ppm':err,'n_in':nin,'snr':dep/err,'random_trials_max_ppm':float(trials.max()*1e6),'frac_random_trials_ge_depth':float((trials*1e6>=dep).mean())}
ex=since>0.3   # excise first 0.3 d after every gap (and quarter starts)
dep2,err2,nin2,trials2=boxdepth(tk[ex],rk[ex],P,T0,DUR)
R['box_depth_postgap_excised_0p3d']={'depth_ppm':dep2,'err_random_phase_ppm':err2,'n_in':nin2,'snr':dep2/err2,'n_removed':int((~ex).sum()),'frac_random_trials_ge_depth':float((trials2*1e6>=dep2).mean())}
ex2=(since>0.3)&(until>0.3)
dep3,err3,nin3,trials3=boxdepth(tk[ex2],rk[ex2],P,T0,DUR)
R['box_depth_gap_excised_0p3d_both_sides']={'depth_ppm':dep3,'err_random_phase_ppm':err3,'n_in':nin3,'snr':dep3/err3,'n_removed':int((~ex2).sum())}
depS,errS,ninS,_=boxdepth(tk,rsk,P,T0,DUR); depS2,errS2,ninS2,_=boxdepth(tk[ex],rsk[ex],P,T0,DUR)
R['sap_box_depth']={'all':{'depth_ppm':depS,'err_ppm':errS,'snr':depS/errS},'postgap_excised':{'depth_ppm':depS2,'err_ppm':errS2,'snr':depS2/errS2}}
# ---- random-ephemeris calibration: how often does a random t0 at this P produce >= 2 events within 0.15 d after a gap, and what depth does it get
cal=[]
for k in range(1500):
    t0r=T0+rng.uniform(0,1)*P
    tcs=t0r+P*np.arange(-2,60); tcs=tcs[(tcs>tk.min())&(tcs<tk.max())]
    s_,u_=time_since_gap(tcs)
    # only events with data
    hasdata=np.array([ (np.abs(tk-tc)<DUR/2).sum()>=2 for tc in tcs])
    npost=int(((s_<0.15)&hasdata).sum())
    phr=phase(tk,P,t0r)*P; i2=np.abs(phr)<DUR/2; o2=np.abs(phr)>1.5*DUR
    dd=(rk[o2].mean()-rk[i2].mean())*1e6
    cal.append((npost,dd))
cal=np.array(cal)
R['random_ephemeris_calibration']={'n_trials':1500,'frac_with_ge1_postgap_event':float((cal[:,0]>=1).mean()),'frac_with_ge2_postgap_event':float((cal[:,0]>=2).mean()),
    'mean_depth_npost0_ppm':float(cal[cal[:,0]==0,1].mean()),'sd_depth_npost0_ppm':float(cal[cal[:,0]==0,1].std()),
    'mean_depth_npost1_ppm':float(cal[cal[:,0]==1,1].mean()) if (cal[:,0]==1).any() else None,
    'mean_depth_npost2plus_ppm':float(cal[cal[:,0]>=2,1].mean()) if (cal[:,0]>=2).any() else None,
    'sd_depth_npost2plus_ppm':float(cal[cal[:,0]>=2,1].std()) if (cal[:,0]>=2).any() else None,
    'n_npost2plus':int((cal[:,0]>=2).sum()),
    'real_ephemeris_npost':int((sg<0.15).sum()),'real_depth_ppm':dep,
    'frac_trials_depth_ge_real':float((cal[:,1]>=dep).mean())}
# ---- per season / channel depths, PSF and SAP, absolute e/s
seasons={}
for qq in np.unique(quk):
    sq=quk==qq
    dq,eq,nq,_=boxdepth(tk[sq],rk[sq],P,T0,DUR,ntrial=200); dqs,eqs,nqs,_=boxdepth(tk[sq],rsk[sq],P,T0,DUR,ntrial=200)
    sq2=sq&ex
    dq2,eq2,nq2,_=boxdepth(tk[sq2],rk[sq2],P,T0,DUR,ntrial=200)
    seasons[int(qq)]={'season':raw[int(qq)]['season'],'channel':raw[int(qq)]['channel'],'psf_depth_ppm':dq,'psf_err_ppm':eq,'n_in':nq,'sap_depth_ppm':dqs,'sap_err_ppm':eqs,
        'psf_depth_postgap_excised_ppm':dq2,'psf_err_postgap_excised_ppm':eq2,'n_in_excised':nq2,
        'psf_raw_median_e_s':raw[int(qq)]['median_psf_raw'],'abs_dip_e_s':dq*1e-6*raw[int(qq)]['median_psf_raw'],'abs_dip_err_e_s':eq*1e-6*raw[int(qq)]['median_psf_raw']}
R['per_quarter']=seasons
bys={}
for s in [0,1,2,3]:
    qs=[q for q in seasons if seasons[q]['season']==s]
    if not qs: continue
    sq=np.isin(quk,qs)
    dq,eq,nq,_=boxdepth(tk[sq],rk[sq],P,T0,DUR,ntrial=300); dqs,eqs,nqs,_=boxdepth(tk[sq],rsk[sq],P,T0,DUR,ntrial=300)
    sq2=sq&ex; dq2,eq2,nq2,_=boxdepth(tk[sq2],rk[sq2],P,T0,DUR,ntrial=300)
    bys[s]={'quarters':qs,'channel':seasons[qs[0]]['channel'],'psf_depth_ppm':dq,'psf_err_ppm':eq,'n_in':nq,'sap_depth_ppm':dqs,'sap_err_ppm':eqs,
            'psf_depth_postgap_excised_ppm':dq2,'psf_err_postgap_excised_ppm':eq2,'n_in_excised':nq2}
R['per_season']=bys
# ---- odd/even, secondary scan (post-gap excised and all)
def oddeven(tt,rr):
    ph=phase(tt,P,T0)*P; n=np.round((tt-T0)/P).astype(int)
    inn=np.abs(ph)<DUR/2; out=np.abs(ph)>1.5*DUR; base=rr[out].mean()
    o=inn&(n%2==1); e=inn&(n%2==0)
    do=(base-rr[o].mean())*1e6; de=(base-rr[e].mean())*1e6
    so=rr[out].std()/np.sqrt(o.sum())*1e6; se=rr[out].std()/np.sqrt(e.sum())*1e6
    return {'odd_ppm':float(do),'even_ppm':float(de),'sigma':float(abs(do-de)/np.hypot(so,se)),'n_odd':int(o.sum()),'n_even':int(e.sum())}
R['odd_even_all']=oddeven(tk,rk); R['odd_even_excised']=oddeven(tk[ex],rk[ex])
# phase scan
scan=[]
for phc in np.linspace(-0.5,0.5,201)[:-1]:
    ph=phase(tk[ex],P,T0); i2=np.abs(((ph-phc+0.5)%1)-0.5)*P<DUR/2
    scan.append((phc,(rk[ex][~i2].mean()-rk[ex][i2].mean())*1e6))
scan=np.array(scan); sec=scan[np.argmin(np.abs(scan[:,0]-0.5))] if False else scan[np.argmin(np.abs(np.abs(scan[:,0])-0.5))]
R['phase_scan_excised']={'secondary_ppm_at_0p5':float(sec[1]),'scan_sd_ppm_excluding_primary':float(scan[np.abs(scan[:,0])>0.05,1].std()),'deepest_other_phase':float(scan[np.abs(scan[:,0])>0.05][np.argmax(scan[np.abs(scan[:,0])>0.05,1]),0]),'deepest_other_ppm':float(scan[np.abs(scan[:,0])>0.05,1].max()),'primary_ppm':float(scan[np.argmin(np.abs(scan[:,0])),1])}
# ---- Lomb-Scargle both polarities on per-quarter-normalised flux (no detrend)
y=fl.copy()
for qq in np.unique(qu): y[qu==qq]/=np.nanmedian(y[qu==qq])
y=y[keep]
fr=np.linspace(1/100,1/0.2,200000)
ls=LombScargle(tk,y-1); pw=ls.power(fr)
pk=np.argmax(pw); 
R['lomb_scargle']={'best_period_d':float(1/fr[pk]),'best_power':float(pw[pk]),'amp_ppm':float(np.sqrt(pw[pk])*np.std(y-1)*1e6*np.sqrt(2)),
   'power_near_P':float(pw[np.abs(1/fr-P)<0.3].max()),'power_near_P_half':float(pw[np.abs(1/fr-P/2)<0.15].max()),'power_near_2P':float(pw[np.abs(1/fr-2*P)<0.6].max()),
   'top5_periods':[float(1/fr[i]) for i in np.argsort(pw)[-5:][::-1]]}
# ---- BLS 20-40 d, real and inverted, with post-gap excision
def bls(tt,rr,pmin=20,pmax=40):
    m=BoxLeastSquares(tt,rr+1)
    per=np.exp(np.linspace(np.log(pmin),np.log(pmax),20000))
    res=m.power(per,[1.5/24,2/24,2.5/24,3/24,4/24],objective='snr')
    pwr=res.power; sde=(pwr-np.median(pwr))/(1.4826*np.median(np.abs(pwr-np.median(pwr))))
    i=np.argmax(pwr); j=np.argmin(np.abs(per-P))
    jwin=np.abs(per/P-1)<0.003
    return {'top_period':float(per[i]),'top_sde':float(sde[i]),'top_depth_ppm':float(res.depth[i]*1e6),'top_dur_h':float(res.duration[i]*24),'sde_at_P':float(sde[jwin].max()),'depth_at_P_ppm':float(res.depth[jwin][np.argmax(sde[jwin])]*1e6)}
R['bls_20_40_all_real']=bls(tk,rk); R['bls_20_40_all_inverted']=bls(tk,-rk)
R['bls_20_40_excised_real']=bls(tk[ex],rk[ex]); R['bls_20_40_excised_inverted']=bls(tk[ex],-rk[ex])
# ---- host light curve at the ephemeris
th,flh,saph,quh,colsh,fh=load(HOST)
def detr_simple(t,y,qu_):
    y=y.copy()
    for qq in np.unique(qu_): s=qu_==qq; y[s]/=np.nanmedian(y[s])
    intr=np.abs(phase(t)*P)<1.5*DUR
    edges=np.where(np.diff(t)>0.5)[0]; segs=np.split(np.arange(len(t)),edges+1); trend=np.ones_like(y)
    for sg_ in segs:
        ts=t[sg_]; ys=y[sg_]; ms=intr[sg_]
        for i,idx in enumerate(sg_):
            w=(np.abs(ts-ts[i])<0.5)&(~ms)
            if w.sum()<5: w=np.abs(ts-ts[i])<0.5
            trend[idx]=np.nanmedian(ys[w])
    r=y/trend-1; s=1.4826*np.nanmedian(np.abs(r-np.nanmedian(r))); k=(r<4*s)&(r>-12*s)
    return r,k
rh,kh=detr_simple(th,flh,quh); rhs,khs=detr_simple(th,saph,quh)
dh,eh,nh,_=boxdepth(th[kh],rh[kh],P,T0,DUR); dhs,ehs,nhs,_=boxdepth(th[khs],rhs[khs],P,T0,DUR)
hostraw=float(np.nanmedian(fh['LIGHTCURVE_STITCHED'].data['FLUX']))
R['host_kic5866192']={'kepmag':fh[0].header.get('KEPMAG'),'gmag':fh[0].header.get('GMAG'),'psf_depth_ppm':dh,'psf_err_ppm':eh,'n_in':nh,'sap_depth_ppm':dhs,'sap_err_ppm':ehs,'median_psf_flux_e_s':hostraw,
   'abs_dip_2sigma_upper_e_s':(dh+2*eh)*1e-6*hostraw,'target_median_psf_flux_e_s':float(np.nanmedian(f['LIGHTCURVE_STITCHED'].data['FLUX'])),'target_abs_dip_e_s':dep*1e-6*float(np.nanmedian(f['LIGHTCURVE_STITCHED'].data['FLUX']))}
# ---- quality flags in transit
qall=f[1].data['SAP_QUALITY']; tall=f[1].data['TIME']; fin=np.isfinite(tall)&np.isfinite(f[1].data['FLUX'])
inn_all=np.abs(phase(tall[fin])*P)<DUR/2
R['quality']={'frac_flagged_all':float((qall[fin]!=0).mean()),'frac_flagged_in_transit':float((qall[fin][inn_all]!=0).mean()),'n_in_transit_total':int(inn_all.sum())}
# ---- duration / density
Rs=0.9458; Ms=1.05
a_over_R=(( (6.674e-11*Ms*1.989e30*(P*86400)**2)/(4*np.pi**2) )**(1/3))/(Rs*6.957e8)
T0c=P*24/np.pi/a_over_R  # central full duration hours (small planet)
b_for=lambda T: np.sqrt(max(0,1-(T/T0c)**2))
R['duration_density']={'a_over_Rstar':float(a_over_R),'T_central_h':float(T0c),'b_for_1p95h':float(b_for(1.95)),'b_for_3h':float(b_for(3.0)),'rp_earth_481ppm':float(np.sqrt(481e-6)*Rs*109.1),'rp_earth_270ppm':float(np.sqrt(270e-6)*Rs*109.1)}
json.dump(R,open(OUT+'lc_checks.json','w'),indent=1,default=float)
# ---- plots
fig,ax=plt.subplots(2,2,figsize=(14,9))
a=ax[0,0]; x=[np.mean(p['t_since_gap_d']) for p in prof]; a.errorbar(x,[p['mean_ppm'] for p in prof],[p['sem_ppm'] for p in prof],fmt='o-'); a.axhline(0,c='k',lw=.5); a.set_xlabel('time since gap end [d]'); a.set_ylabel('mean detrended PSF flux [ppm]'); a.set_title('post-gap profile, all gaps >0.3 d (ramps)')
a=ax[0,1]; a.errorbar([e['tc'] for e in ev],d,e_,fmt='o',c='C0'); near_i=[i for i,e in enumerate(ev) if e['since_gap_d']<0.15]; a.errorbar([ev[i]['tc'] for i in near_i],d[near_i],e_[near_i],fmt='o',c='r',label='<0.15 d after gap'); a.axhline(dw,c='C1'); a.legend(); a.set_xlabel('BKJD'); a.set_ylabel('event depth [ppm]'); a.set_title(f'per-event depths; chi2={chi2:.0f}/{len(ev)-1}')
a=ax[1,0]; a.hist(cal[cal[:,0]==0,1],bins=40,alpha=.6,label='random t0, 0 post-gap events'); a.hist(cal[cal[:,0]>=2,1],bins=40,alpha=.6,label='random t0, >=2 post-gap events'); a.axvline(dep,c='r',label='real ephemeris'); a.axvline(dep2,c='m',ls='--',label='real, post-gap excised'); a.legend(); a.set_xlabel('1.95 h box depth [ppm]'); a.set_title('random-ephemeris calibration (same P)')
a=ax[1,1]
for s,v in bys.items(): a.errorbar(s,v['psf_depth_ppm'],v['psf_err_ppm'],fmt='o',c='C0'); a.errorbar(s+0.15,v['psf_depth_postgap_excised_ppm'],v['psf_err_postgap_excised_ppm'],fmt='s',c='m'); a.errorbar(s+0.3,v['sap_depth_ppm'],v['sap_err_ppm'],fmt='^',c='C2')
a.set_xlabel('season (0: ch17, 1: ch57, 2: ch65, 3: ch25)'); a.set_ylabel('depth [ppm]'); a.set_title('per-season: PSF (o), PSF post-gap excised (s), SAP (^)'); a.axhline(0,c='k',lw=.5)
plt.tight_layout(); plt.savefig(OUT+'lc_checks.png',dpi=110)
print(json.dumps({k:R[k] for k in ['events_summary','events_postgap_split','box_depth_all','box_depth_postgap_excised_0p3d','box_depth_gap_excised_0p3d_both_sides','sap_box_depth','random_ephemeris_calibration','per_season','odd_even_all','odd_even_excised','phase_scan_excised','lomb_scargle','bls_20_40_all_real','bls_20_40_all_inverted','bls_20_40_excised_real','bls_20_40_excised_inverted','host_kic5866192','quality','duration_density','post_gap_profile','per_gap_first0p15d']},indent=1,default=float))
