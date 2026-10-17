import numpy as np, json, sys
from astropy.io import fits
P=15.31335; T0=144.28486; DUR=3.75/24
OUT='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic003831969'
def load(f, col='FLUX'):
    h=fits.open(f); d=h['LIGHTCURVE_STITCHED'].data
    t=d['TIME']; fl=d[col]; q=d['QUARTER']; ok=(d['SAP_QUALITY']==0)&np.isfinite(t)&np.isfinite(fl)
    t,fl,q=t[ok],fl[ok],q[ok]
    # per quarter normalise
    fn=np.empty_like(fl)
    for qq in np.unique(q):
        m=q==qq; fn[m]=fl[m]/np.nanmedian(fl[m])
    return t,fn,q
def detrend(t,f,win=1.0,mask=None):
    # running median excluding masked (in-transit) points
    out=np.empty_like(f); 
    src=f.copy()
    if mask is not None: src[mask]=np.nan
    for i in range(len(t)):
        m=np.abs(t-t[i])<win/2
        v=src[m]; v=v[np.isfinite(v)]
        out[i]=np.nanmedian(v) if len(v)>5 else np.nan
    return f/out-1
def phase(t): return ((t-T0)/P+0.5)%1-0.5
def boxdepth(ph,r,w=DUR/P):
    i=np.abs(ph)<w/2; o=(np.abs(ph)>w)&(np.abs(ph)<3*w)
    if i.sum()<3: return np.nan,np.nan,i.sum()
    d=np.nanmedian(r[o])-np.nanmean(r[i]); 
    e=np.nanstd(r[o])/np.sqrt(i.sum())
    return d*1e6,e*1e6,i.sum()
def randbox_err(ph,r,n=400,w=DUR/P,seed=1):
    rng=np.random.default_rng(seed); ds=[]
    for k in range(n):
        c=rng.uniform(-0.5,0.5)
        if abs(c)<2*w: continue
        pp=((ph-c)+0.5)%1-0.5
        d,_,_=boxdepth(pp,r,w); ds.append(d)
    return np.nanstd(ds)
res={}
for label,f in [('target','/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-003831969_kepler_v1.0_lc.fits'),('host3831977',OUT+'/dl/hlsp_kbonus-bkg_kepler_kepler_kic-003831977_kepler_v1.0_lc.fits')]:
    res[label]={}
    for col in ['FLUX','SAP_FLUX']:
        t,fn,q=load(f,col)
        ph=phase(t); mask=np.abs(ph)<(DUR/P)*0.8
        r=detrend(t,fn,1.0,mask)
        ok=np.isfinite(r); t,r,q,ph=t[ok],r[ok],q[ok],ph[ok]
        # sigma clip upward outliers
        s=1.4826*np.nanmedian(np.abs(r-np.nanmedian(r))); keep=r<4*s; t,r,q,ph=t[keep],r[keep],q[keep],ph[keep]
        d,e,n=boxdepth(ph,r); re=randbox_err(ph,r)
        sec,_,_=boxdepth(((ph+0.5)+0.5)%1-0.5,r)
        row={'depth_ppm':d,'err_white_ppm':e,'err_randbox_ppm':re,'n_in':int(n),'sigma_ppm':s*1e6,'secondary_ppm':sec,'median_flux_e_s':float(np.nanmedian(fits.open(f)['LIGHTCURVE_STITCHED'].data[col]))}
        # per season and per quarter
        seas={}
        for sname,qs in {'S_Q1,5,9,13,17':[1,5,9,13,17],'S_Q2,6,10,14':[2,6,10,14],'S_Q3,7,11,15':[3,7,11,15],'S_Q4,8,12,16':[4,8,12,16]}.items():
            m=np.isin(q,qs); dd,ee,nn=boxdepth(ph[m],r[m]); seas[sname]=[round(dd,0),round(ee,0),int(nn)]
        row['per_season']=seas
        pq={}
        for qq in np.unique(q):
            m=q==qq; dd,ee,nn=boxdepth(ph[m],r[m]); pq[int(qq)]=[round(dd,0),round(ee,0),int(nn)]
        row['per_quarter']=pq
        res[label][col]=row
        print(label,col,json.dumps({k:v for k,v in row.items() if k!='per_quarter'},default=float))
        print('   per quarter:',pq)
json.dump(res,open(OUT+'/lc_tests.json','w'),indent=1,default=float)
# density / duration check
G=6.674e-11; Msun=1.989e30; Rsun=6.957e8
for Rs,Ms in [(0.83,0.85),(1.0,1.0),(1.1,1.05)]:
    a=(G*Ms*Msun*(P*86400)**2/(4*np.pi**2))**(1/3)
    aR=a/(Rs*Rsun); k=0.033
    T14c=P*24/np.pi*np.arcsin(np.sqrt((1+k)**2)/aR)  # central, b=0
    print(f'Rs={Rs} Ms={Ms}: a/R*={aR:.1f} T14(b=0)={T14c:.2f} h; Rp={k*Rs*109.1:.2f} Re')
