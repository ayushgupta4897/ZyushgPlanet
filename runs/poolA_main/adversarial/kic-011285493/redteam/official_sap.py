import warnings; warnings.filterwarnings('ignore')
import lightkurve as lk, numpy as np, json
np.random.seed(4)
P=10.790444; T0=142.33451
pc,prof=np.loadtxt('eb_profile.csv',delimiter=',',unpack=True)
tp=np.where(np.abs(pc)<0.4,prof-1,0.0); ts=np.where(np.abs(np.abs(pc)-P/2)<0.4,prof-1,0.0); prim=1-prof.min()
def detrend(t,y,mask,win):
    out=np.full_like(y,np.nan); gaps=np.where(np.diff(t)>0.5)[0]; edges=np.concatenate([[0],gaps+1,[len(t)]])
    for a,b in zip(edges[:-1],edges[1:]):
        tt=t[a:b]; yy=y[a:b]; mm=mask[a:b]
        for i in range(len(tt)):
            s=(np.abs(tt-tt[i])<win/2)&~mm
            out[a+i]=yy[i]/np.nanmedian(yy[s]) if s.sum()>5 else np.nan
    return out
sr=lk.search_lightcurve('KIC 11285496',mission='Kepler',cadence='long'); lcs=sr.download_all(quality_bitmask='default')
res={}
for col in ['sap_flux','pdcsap_flux']:
    t=[];f=[];Q=[];med=[]
    for lc in lcs:
        lc=lc.remove_nans(); x=np.asarray(getattr(lc,col).value,float); tt=np.asarray(lc.time.value,float); ok=np.isfinite(x)&np.isfinite(tt)
        med.append(np.nanmedian(x[ok])); t.append(tt[ok]); f.append(x[ok]/np.nanmedian(x[ok])); Q.append(np.full(ok.sum(),lc.meta.get('QUARTER',-1)))
    t=np.concatenate(t); f=np.concatenate(f); Q=np.concatenate(Q); i=np.argsort(t); t,f,Q=t[i],f[i],Q[i]
    ph=((t-T0+0.5*P)%P)-0.5*P; mask=(np.abs(ph)<0.3)|(np.abs(np.abs(ph)-P/2)<0.3)
    r=detrend(t,f,mask,1.0)-1; ok=np.isfinite(r); t,r,ph,Q=t[ok],r[ok],ph[ok],Q[ok]
    mad=1.4826*np.median(np.abs(r-np.median(r))); keep=np.abs(r)<5*mad; t,r,ph,Q=t[keep],r[keep],ph[keep],Q[keep]
    Tp=np.interp(ph,pc,tp); Ts=np.interp(ph,pc,ts); p=np.linalg.lstsq(np.c_[Tp,Ts],r,rcond=None)[0]
    sc=[]
    for k in range(300):
        sh=np.random.uniform(0,P); p2=((t-T0-sh+0.5*P)%P)-0.5*P; B=np.interp(p2,pc,tp); sc.append(np.sum(B*r)/np.sum(B*B))
    emp=1.4826*np.median(np.abs(np.array(sc)-np.median(sc)))
    seas=[]
    for s_ in range(4):
        m_=((Q-2)%4)==s_
        if m_.sum()<100: seas.append(None); continue
        pp=np.linalg.lstsq(np.c_[Tp[m_],Ts[m_]],r[m_],rcond=None)[0]; seas.append(round(float(pp[0]*prim*1e6)))
    M=np.median(med)
    print(f'KIC 11285496 official {col}: median {M:.0f} e/s, n={len(r)}, 30-min MAD {mad*1e6:.0f} ppm | EB-primary template dip {p[0]*prim*1e6:.1f} +- {emp*prim*1e6:.1f} ppm (SNR {p[0]/emp:.1f}) = {p[0]*prim*M:.2f} e/s | secondary {p[1]*(1-prof[np.abs(np.abs(pc)-P/2)<0.3].min())*1e6:.1f} ppm | per season ppm {seas}')
    res[col]=dict(median_es=float(M),prim_ppm=float(p[0]*prim*1e6),err_ppm=float(emp*prim*1e6),snr=float(p[0]/emp),prim_es=float(p[0]*prim*M),season_prim_ppm=seas)
json.dump(res,open('official_sap.json','w'),indent=1)
