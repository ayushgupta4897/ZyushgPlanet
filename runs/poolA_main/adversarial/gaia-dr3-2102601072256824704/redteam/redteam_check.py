import numpy as np, json, sys
from astropy.io import fits
from astropy.timeseries import LombScargle, BoxLeastSquares
from wotan import flatten
F='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2102601072256824704_kepler_v1.0_lc.fits'
P=32.55706418774713; t0=148.82345687868073; dur=3.375/24
h=fits.open(F); d=h['LIGHTCURVE_STITCHED'].data
m=(d['SAP_QUALITY']==0)&np.isfinite(d['TIME'])&np.isfinite(d['FLUX'])&np.isfinite(d['SAP_FLUX'])
t=d['TIME'][m].astype(float); f=d['FLUX'][m].astype(float); s=d['SAP_FLUX'][m].astype(float); q=d['QUARTER'][m]
for arr in (f,s):
    for qq in np.unique(q):
        k=q==qq; arr[k]/=np.nanmedian(arr[k])
print('n',len(t),'quarters',np.unique(q))
# segments
gap=np.where(np.diff(t)>0.5)[0]; starts=np.r_[0,gap+1]; ends=np.r_[gap,len(t)-1]
edge_t=np.sort(np.r_[t[starts],t[ends]])
def edge_dist(x): return np.min(np.abs(edge_t[None,:]-np.atleast_1d(x)[:,None]),axis=1)
print('segments',len(starts))
out={}
def analyse(label,flux,edge_cutoff):
    fl,tr=flatten(t,flux,method='biweight',window_length=1.0,break_tolerance=0.5,edge_cutoff=edge_cutoff,return_trend=True)
    ok=np.isfinite(fl); tt=t[ok]; ff=fl[ok]-1
    # sigma clip positive outliers only (flares) and extreme negatives > 10 sigma
    sd=1.4826*np.median(np.abs(ff-np.median(ff))); ok2=(ff<4*sd)&(ff>-10*sd); tt=tt[ok2]; ff=ff[ok2]
    ph=((tt-t0)/P+0.5)%1-0.5; intr=np.abs(ph*P)<dur/2
    # empirical depth error: box means at random phases
    rng=np.random.default_rng(1); errs=[]
    for i in range(1500):
        p0=rng.uniform(-0.5,0.5); errs.append(-np.mean(ff[np.abs(((ph-p0+0.5)%1-0.5)*P)<dur/2]))
    e=np.std(errs); dep=-np.mean(ff[intr])*1e6
    n=np.round((tt-t0)/P).astype(int)
    ev=[]
    for k in np.unique(n[intr]):
        loc=(n==k)&(np.abs(tt-(t0+k*P))<1.0); inn=loc&intr; outn=loc&~intr&(np.abs(tt-(t0+k*P))<0.5)
        if inn.sum()>=3 and outn.sum()>=5:
            ev.append((int(k),float(t0+k*P),round(float((np.mean(ff[outn])-np.mean(ff[inn]))*1e6)),int(inn.sum()),round(float(edge_dist(t0+k*P)[0]),3)))
    ev=np.array(ev,dtype=object)
    depths=np.array([x[2] for x in ev],float); medabs=np.median(np.abs(depths))
    big=[x for x in ev if abs(x[2])>5*medabs]
    odd=-np.mean(ff[intr&(n%2==1)])*1e6; even=-np.mean(ff[intr&(n%2==0)])*1e6
    # exclude edge events (edge dist<0.15 d)
    keep=np.ones_like(ff,bool)
    for x in ev:
        if x[4]<0.15: keep&=~(n==x[0])
    dep_x=-np.mean(ff[intr&keep])*1e6
    # season depths
    seas={}
    for qq in np.unique(q):
        pass
    qq=q[ok][ok2]; 
    for sea,qs in {'ch22':[1,5,9,13,17],'ch38':[2,6,10,14],'ch62':[3,7,11,15],'ch46':[4,8,12,16]}.items():
        kk=np.isin(qq,qs)&intr; seas[sea]=(round(float(-np.mean(ff[kk])*1e6)),int(kk.sum()))
    res=dict(n=int(len(ff)),rms_ppm=round(float(np.std(ff)*1e6)),depth_ppm=round(dep),depth_err_ppm=round(e*1e6),snr=round(dep/(e*1e6),1),
             odd=round(odd),even=round(even),n_events=len(ev),big_events=[list(x) for x in big],depth_excl_edge_events=round(dep_x),
             n_events_excl=int(len([x for x in ev if x[4]>=0.15])),season_depths=seas,
             events=[list(x) for x in ev])
    print(label, {k:v for k,v in res.items() if k!='events'})
    out[label]=res
    return tt,ff
tt,ff=analyse('PSF_edge0',f,0.0)
tt2,ff2=analyse('PSF_edge0.3',f,0.3)
analyse('SAP_edge0',s,0.0)
analyse('SAP_edge0.3',s,0.3)
# Lomb-Scargle both polarities on detrended edge0
for lab,y in [('real',ff),('inverted',-ff)]:
    freq=np.linspace(1/60,1/0.2,200000); pw=LombScargle(tt,y).power(freq)
    i=np.argsort(pw)[::-1][:5]; print('LS',lab,[(round(1/freq[j],4),round(float(pw[j]),5)) for j in i])
    amp=np.sqrt(pw.max())*np.std(y)*np.sqrt(2)*1e6
    print('  approx amplitude of top peak (ppm)',round(float(amp),1))
    # power near P, P/2
    for Pt in [P,P/2,P/3]:
        k=np.abs(1/freq-Pt)<0.02*Pt; print('  max LS power near',round(Pt,2),'d:',round(float(pw[k].max()),5))
# narrow BLS 30-35 d, real and inverted, edge0 and edge0.3
for lab,(x,y) in {'edge0':(tt,ff),'edge0.3':(tt2,ff2)}.items():
    periods=np.exp(np.linspace(np.log(28),np.log(37),20000))
    for pol,yy in [('real',y),('inv',-y)]:
        bls=BoxLeastSquares(x,yy+1); r=bls.power(periods,[0.1,0.14,0.18,0.25],objective='snr')
        i=np.argmax(r.power); pk=r.power
        sde=(pk[i]-np.median(pk))/np.std(pk)
        kP=np.argmin(np.abs(periods-P))
        print('BLS',lab,pol,'top P',round(periods[i],4),'snr',round(float(pk[i]),2),'SDE-like',round(float(sde),2),'depth',round(float(r.depth[i]*1e6)),'| at P_pipe snr',round(float(pk[kP]),2),'depth',round(float(r.depth[kP]*1e6)))
# raw ramp slopes at the three events
for tc in [181.42,539.53,1418.53]:
    k=np.abs(t-tc)<0.6; 
    if k.sum()>5:
        c=np.polyfit(t[k]-tc,f[k]-1,1); print('raw PSF slope near',tc,': %.0f ppm/d, local range %.0f ppm'%(c[0]*1e6,(f[k].max()-f[k].min())*1e6), 'edge dist',round(float(edge_dist(tc)[0]),3))
json.dump(out,open(sys.argv[1],'w'),indent=1)
