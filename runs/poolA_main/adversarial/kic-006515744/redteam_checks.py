import numpy as np, glob, json, sys, csv, math
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares, LombScargle
from astropy.stats import sigma_clip
np.seterr(all='ignore')
D='/home/user/ZyushgPlanet/data/kbonus_batch/'
FN=D+'hlsp_kbonus-bkg_kepler_kepler_kic-006515744_kepler_v1.0_lc.fits'
OUT={}
f=fits.open(FN)
s=f['LIGHTCURVE_STITCHED'].data
ok=(s['SAP_QUALITY']==0)&np.isfinite(s['FLUX'])&np.isfinite(s['TIME'])
t=s['TIME'][ok]; fl=s['FLUX'][ok]; fe=s['FLUX_ERR'][ok]; sap=s['SAP_FLUX'][ok]; q=s['QUARTER'][ok]; cad=s['CADENCENO'][ok]
# simple per-quarter running median detrend (1 d window)
def rmed(t,y,w=1.0):
    out=np.empty_like(y)
    for i in range(len(t)):
        m=(t>t[i]-w/2)&(t<t[i]+w/2); out[i]=np.median(y[m])
    return out
rel=np.empty_like(fl); relsap=np.empty_like(fl)
for qq in np.unique(q):
    m=q==qq
    rel[m]=fl[m]/rmed(t[m],fl[m])-1
    relsap[m]=sap[m]/rmed(t[m],sap[m])-1
# gaps
dt=np.diff(t); gaps=[(t[i],t[i+1]) for i in np.where(dt>0.5)[0]]
def gapdist(x):
    d=[]
    for a,b in gaps: d+= [abs(x-a),abs(x-b)]
    d+=[abs(x-t[0]),abs(x-t[-1])]
    return min(d)
# per-quarter extension lookup
qext={}
for e in f[1:]:
    if e.name.startswith('LIGHTCURVE_Q'):
        qext[int(e.name[12:])]=(e.data,e.header)
def cadinfo(ti):
    i=np.argmin(abs(t-ti)); qq=q[i]; d,h=qext[qq]
    j=np.argmin(abs(d['TIME']-ti))
    good=d['SAP_QUALITY']==0
    def z(col):
        v=d[col][good]; med=np.nanmedian(v); mad=1.4826*np.nanmedian(abs(v-med)); return float((d[col][j]-med)/mad) if mad>0 else float('nan')
    return dict(time=float(t[i]),cadenceno=int(cad[i]),quarter=int(qq),channel=int(h['CHANNEL']),rel_psf_ppm=float(rel[i]*1e6),rel_sap_ppm=float(relsap[i]*1e6),
        flux_err_z=z('FLUX_ERR'),red_chi2=float(d['RED_CHI2'][j]),red_chi2_z=z('RED_CHI2'),red_chi2_med=float(np.nanmedian(d['RED_CHI2'][good])),
        sap_bkg_z=z('SAP_BKG'),cen_col_z=z('CENTROID_COLUMN'),cen_row_z=z('CENTROID_ROW'),
        psf_nova_rel_ppm=float((d['PSF_FLUX_NOVA'][j]/np.nanmedian(d['PSF_FLUX_NOVA'][good])-1)*1e6),
        gap_dist_d=float(gapdist(t[i])),
        neighbours_rel_ppm=[float(x*1e6) for x in rel[max(0,i-3):i+4]])
# A. events of pipeline ephemeris and rederive ephemeris
eph={'pipeline':(36.20754825773067,127.92048643985763,0.75/24),'rederive':(33.421805994552884,133.4822,0.48/24)}
ev={}
for k,(P,t0,dur) in eph.items():
    rows=[]
    n0=int(np.floor((t[0]-t0)/P)); n1=int(np.ceil((t[-1]-t0)/P))
    for n in range(n0,n1+1):
        tc=t0+n*P; m=abs(t-tc)<max(dur,0.0205)/2*1.0+1e-9
        m=abs(t-tc)<=0.0205  # one cadence either side ~ +/-0.49 h
        if m.sum()==0: continue
        rows.append(dict(epoch=n,tc=float(tc),n_in=int(m.sum()),psf_depth_ppm=float(-np.mean(rel[m])*1e6),sap_depth_ppm=float(-np.mean(relsap[m])*1e6),min_cad_ppm=float(np.min(rel[m])*1e6),gap_dist_d=float(gapdist(tc)),quarter=int(q[m][0])))
    ev[k]=rows
OUT['events']=ev
OUT['outlier_cadences']={str(x):cadinfo(x) for x in (200.33,1370.09)}
# the most negative cadences overall
idx=np.argsort(rel)[:10]
OUT['ten_most_negative_cadences']=[dict(time=float(t[i]),rel_ppm=float(rel[i]*1e6),sap_ppm=float(relsap[i]*1e6),quarter=int(q[i]),gap_dist=float(gapdist(t[i]))) for i in idx]
# C. BLS narrow window with and without the two cadences
def bls(tt,yy,ee,pmin=20,pmax=50):
    b=BoxLeastSquares(tt,yy,ee)
    per=np.exp(np.linspace(np.log(pmin),np.log(pmax),60000))
    r=b.power(per,[0.02,0.04,0.06,0.1],objective='likelihood')
    p=r.power; sde=(p.max()-np.median(p))/np.std(p)
    i=np.argmax(p)
    def at(P):
        j=np.argmin(abs(per-P)); return float(p[j])
    return dict(best_P=float(per[i]),sde=float(sde),depth_ppm=float(r.depth[i]*1e6),dur_d=float(r.duration[i]),t0=float(r.transit_time[i]),power_at_36_2075=at(36.2075),power_at_33_4218=at(33.4218),power_median=float(np.median(p)),power_std=float(np.std(p)))
good=np.isfinite(rel)
clip=rel<np.nanmedian(rel)+4*1.4826*np.nanmedian(abs(rel-np.nanmedian(rel)))
OUT['bls_all']=bls(t[good&clip],rel[good&clip],fe[good&clip]/fl[good&clip])
drop=(abs(t-200.33)<0.02)|(abs(t-1370.09)<0.02)
OUT['bls_drop2']=bls(t[good&clip&~drop],rel[good&clip&~drop],fe[good&clip&~drop]/fl[good&clip&~drop])
# symmetric: clip negative outliers >5 sigma (single cadence) 
mad=1.4826*np.nanmedian(abs(rel-np.nanmedian(rel)))
neg5=rel<np.nanmedian(rel)-5*mad
OUT['n_cadences_below_5mad']=int(neg5.sum()); OUT['mad_ppm']=float(mad*1e6)
OUT['bls_sym5clip']=bls(t[good&clip&~neg5],rel[good&clip&~neg5],fe[good&clip&~neg5]/fl[good&clip&~neg5])
OUT['bls_inverted_all']=bls(t[good&clip],-rel[good&clip],fe[good&clip]/fl[good&clip])
# E. Lomb-Scargle
ls=LombScargle(t[good&clip],rel[good&clip]); fr,pw=ls.autopower(minimum_frequency=1/100,maximum_frequency=1/0.1,samples_per_peak=5)
i=np.argmax(pw); OUT['lombscargle']=dict(best_period_d=float(1/fr[i]),power=float(pw[i]),fap=float(ls.false_alarm_probability(pw[i])))
# F. duration/density
Rs=0.854; Teff=5222; Ms=0.9
P=36.2075; a_AU=(Ms*(P/365.25)**2)**(1/3); a_Rs=a_AU*215.03/Rs
T14_central_h=P/np.pi/a_Rs*24*(1+0.113)
OUT['duration_physics']=dict(a_over_Rs=float(a_Rs),T14_central_h=float(T14_central_h),observed_dur_h=[0.48,0.75,1.5],implied_b_for_1h=float(np.sqrt(max(0,(1+0.113)**2-(1.0/T14_central_h*(1+0.113))**2))),rp_rs_box=float(np.sqrt(0.012763)),rp_earth=float(np.sqrt(0.012763)*Rs*109.1))
# G. ephemeris match against KOI/TCE whole-field for P (and 2P, P/2) with epoch agreement
def ephmatch(P,t0,label):
    hits=[]
    for fn,pk,tk,idk in [('koi_cumulative.csv','koi_period','koi_time0bk','kepoi_name'),('dr25_tce.csv','tce_period','tce_time0bk','kepid')]:
        for r in csv.DictReader(open('/home/user/ZyushgPlanet/data/known/'+fn)):
            try: Pc=float(r[pk]); tc=float(r[tk])
            except: continue
            for mult in (0.5,1,2,3,1/3):
                if abs(Pc/(P*mult)-1)<0.002:
                    # epoch agreement: phase of t0 in catalogue ephemeris
                    ph=((t0-tc)/Pc+0.5)%1-0.5
                    hits.append(dict(cat=fn,id=r[idk] if fn.startswith('koi') else r['kepid']+'.'+r['tce_plnt_num'],P_cat=Pc,mult=mult,t0_phase_offset_d=float(ph*Pc),ra=r['ra'],dec=r['dec']))
    return hits
OUT['ephmatch_pipeline']=[h for h in ephmatch(36.20754825773067,127.92048643985763,'pipe') if abs(h['t0_phase_offset_d'])<0.1]
OUT['ephmatch_rederive']=[h for h in ephmatch(33.4218,133.4822,'red') if abs(h['t0_phase_offset_d'])<0.1]
OUT['n_period_only_matches_pipeline']=len(ephmatch(36.20754825773067,127.92048643985763,'pipe'))
# B. cross-source test at the two cadences
cad200=OUT['outlier_cadences']['200.33']['cadenceno']; cad1370=OUT['outlier_cadences']['1370.09']['cadenceno']
res={'200.33':[], '1370.09':[]}
files=sorted(glob.glob(D+'*.fits'))
rng=np.random.default_rng(1)
for fn2 in files:
    try:
        g=fits.open(fn2); d=g['LIGHTCURVE_STITCHED'].data
        for key,cn,qq in (('200.33',cad200,2),('1370.09',cad1370,14)):
            m=(d['QUARTER']==qq)&(d['SAP_QUALITY']==0)&np.isfinite(d['FLUX'])
            if m.sum()<500: continue
            tt=d['TIME'][m]; yy=d['FLUX'][m]; cc=d['CADENCENO'][m]
            j=np.where(cc==cn)[0]
            if len(j)==0: continue
            j=j[0]
            w=(tt>tt[j]-0.5)&(tt<tt[j]+0.5); w[j]=False
            base=np.median(yy[w]); sc=1.4826*np.median(abs(yy[w]-base))
            z=(yy[j]-base)/sc if sc>0 else np.nan
            # channel
            ch=None
            try: ch=g['LIGHTCURVE_Q%d'%qq].header['CHANNEL']
            except: pass
            res[key].append((fn2.split('kepler_')[2], float(z), ch))
        g.close()
    except Exception as e:
        pass
OUT['cross_source']={}
for key in res:
    zs=np.array([r[1] for r in res[key]]); chs=[r[2] for r in res[key]]
    ch0=OUT['outlier_cadences'][key]['channel']
    same=np.array([c==ch0 for c in chs])
    OUT['cross_source'][key]=dict(n_sources=int(len(zs)),median_z=float(np.nanmedian(zs)),frac_below_m3=float(np.nanmean(zs<-3)),frac_above_3=float(np.nanmean(zs>3)),n_same_channel=int(same.sum()),median_z_same_channel=float(np.nanmedian(zs[same])) if same.sum() else None,frac_below_m3_same_channel=float(np.nanmean(zs[same]<-3)) if same.sum() else None,
        most_negative=[(r[0],round(r[1],1),r[2]) for r in sorted(res[key],key=lambda r:r[1])[:8]])
json.dump(OUT,open('/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006515744/redteam_checks.json','w'),indent=1,default=float)
print(json.dumps(OUT,indent=1,default=float))
