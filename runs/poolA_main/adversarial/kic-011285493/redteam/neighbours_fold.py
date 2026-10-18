import numpy as np, glob, sys, json
from astropy.io import fits
np.random.seed(3)
P=10.790444; T0=142.33451
pc,prof=np.loadtxt('eb_profile.csv',delimiter=',',unpack=True)
tp=np.where(np.abs(pc)<0.4,prof-1,0.0); ts=np.where(np.abs(np.abs(pc)-P/2)<0.4,prof-1,0.0)
prim=1-prof.min(); 
def detrend(t,y,mask,win):
    out=np.full_like(y,np.nan); gaps=np.where(np.diff(t)>0.5)[0]; edges=np.concatenate([[0],gaps+1,[len(t)]])
    for a,b in zip(edges[:-1],edges[1:]):
        tt=t[a:b]; yy=y[a:b]; mm=mask[a:b]
        for i in range(len(tt)):
            s=(np.abs(tt-tt[i])<win/2)&~mm
            out[a+i]=yy[i]/np.nanmedian(yy[s]) if s.sum()>5 else np.nan
    return out
files=sorted(glob.glob('/home/user/ZyushgPlanet/data/kbonus_batch/*kic-0112856*')+glob.glob('/home/user/ZyushgPlanet/data/kbonus_batch/*kic-011338482*')+glob.glob('/home/user/ZyushgPlanet/data/kbonus_batch/*kic-011285855*')+glob.glob('/home/user/ZyushgPlanet/data/kbonus_batch/*kic-011285117*')+glob.glob('downloads/*.fits'))
out=[]
for F in files:
    h=fits.open(F); d=h['LIGHTCURVE_STITCHED'].data; hdr=h[0].header
    name=F.split('kepler_kepler_')[1].split('_kepler')[0]
    m=np.isfinite(d['TIME'])&np.isfinite(d['FLUX'])&(d['SAP_QUALITY']==0)
    t=d['TIME'][m]; f=d['FLUX'][m].astype(float); Q=d['QUARTER'][m]; sap=d['SAP_FLUX'][m].astype(float)
    psfmed=np.nanmedian(f)
    for qq in np.unique(Q): f[Q==qq]/=np.nanmedian(f[Q==qq])
    ph=((t-T0+0.5*P)%P)-0.5*P; mask=(np.abs(ph)<0.3)|(np.abs(np.abs(ph)-P/2)<0.3)
    r=detrend(t,f,mask,1.0)-1; ok=np.isfinite(r); t,r,ph,Q=t[ok],r[ok],ph[ok],Q[ok]
    mad=1.4826*np.median(np.abs(r-np.median(r))); keep=r<5*mad; t,r,ph,Q=t[keep],r[keep],ph[keep],Q[keep]
    Tp=np.interp(ph,pc,tp); Ts=np.interp(ph,pc,ts)
    A=np.c_[Tp,Ts]; p=np.linalg.lstsq(A,r,rcond=None)[0]
    sc=[]
    for k in range(300):
        sh=np.random.uniform(0,P); p2=((t-T0-sh+0.5*P)%P)-0.5*P; B=np.interp(p2,pc,tp); sc.append(np.sum(B*r)/np.sum(B*B))
    emp=1.4826*np.median(np.abs(np.array(sc)-np.median(sc)))
    seas=[]
    for s_ in range(4):
        m_=((Q-2)%4)==s_
        if m_.sum()<100: seas.append(None); continue
        pp=np.linalg.lstsq(np.c_[Tp[m_],Ts[m_]],r[m_],rcond=None)[0]; seas.append(round(float(pp[0]*prim*1e6)))
    row=dict(name=name,kepmag=hdr.get('KEPMAG'),gmag=hdr.get('GMAG',hdr.get('PHOT_G_MEAN_MAG')),tpforg=hdr.get('TPFORG'),psf_med_es=float(psfmed),n=int(len(r)),prim_ppm=float(p[0]*prim*1e6),err_ppm=float(emp*prim*1e6),snr=float(p[0]/emp),sec_ppm=float(p[1]*(1-prof[np.abs(np.abs(pc)-P/2)<0.3].min())*1e6),prim_es=float(p[0]*prim*psfmed),season_prim_ppm=seas)
    out.append(row)
    print(f"{name:16s} Kp={row['kepmag']} tpforg={row['tpforg']} PSFmed={psfmed:8.0f} e/s | primary-template dip {row['prim_ppm']:7.0f} +- {row['err_ppm']:4.0f} ppm (SNR {row['snr']:4.1f}) = {row['prim_es']:6.2f} e/s | secondary {row['sec_ppm']:6.0f} ppm | per season {seas}")
json.dump(out,open('neighbours_fold.json','w'),indent=1)
