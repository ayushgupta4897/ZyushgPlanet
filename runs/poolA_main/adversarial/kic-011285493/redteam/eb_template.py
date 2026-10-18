import numpy as np, json
from astropy.io import fits
np.random.seed(2)
P=10.790444; T0=142.33451
te,fe_,Qe=np.load('eb_lc.npy')
def detrend(t,y,mask,win):
    out=np.full_like(y,np.nan); gaps=np.where(np.diff(t)>0.5)[0]; edges=np.concatenate([[0],gaps+1,[len(t)]])
    for a,b in zip(edges[:-1],edges[1:]):
        tt=t[a:b]; yy=y[a:b]; mm=mask[a:b]
        for i in range(len(tt)):
            s=(np.abs(tt-tt[i])<win/2)&~mm
            out[a+i]=yy[i]/np.nanmedian(yy[s]) if s.sum()>5 else np.nan
    return out
phe=((te-T0+0.5*P)%P)-0.5*P
ecl=(np.abs(phe)<0.3)|(np.abs(np.abs(phe)-P/2)<0.3)
fd=detrend(te,fe_,ecl,2.0); ok=np.isfinite(fd); te,fd,phe=te[ok],fd[ok],phe[ok]
bins=np.linspace(-0.5*P,0.5*P,int(P*24*4)+1); idx=np.digitize(phe,bins)
prof=np.array([np.nanmedian(fd[idx==i]) if (idx==i).sum()>3 else np.nan for i in range(1,len(bins))]); pc=0.5*(bins[1:]+bins[:-1])
good=np.isfinite(prof); prof=np.interp(pc,pc[good],prof[good])
prim=1-np.min(prof[np.abs(pc)<0.3]); sec=1-np.min(prof[np.abs(np.abs(pc)-P/2)<0.3])
w=pc[(np.abs(pc)<0.4)&(prof<1-0.1*prim)]
isec=np.abs(np.abs(pc)-P/2)<0.3
print(f'EB primary depth {prim:.4f}, secondary {sec:.4f}, ratio sec/prim {sec/prim:.3f}; primary T(>10% depth) {(w.max()-w.min())*24:.2f} h; secondary minimum at |phase| = {np.abs(pc[isec][np.argmin(prof[isec])]):.4f} d vs P/2 = {P/2:.4f}')
np.savetxt('eb_profile.csv',np.c_[pc,prof],delimiter=',',header='phase_d,relflux')
tp=np.where(np.abs(pc)<0.4,prof-1,0.0); ts=np.where(isec|(np.abs(np.abs(pc)-P/2)<0.4),prof-1,0.0)
# ---- target
h=fits.open('/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-011285493_kepler_v1.0_lc.fits'); d=h['LIGHTCURVE_STITCHED'].data
m=np.isfinite(d['TIME'])&np.isfinite(d['FLUX'])&(d['SAP_QUALITY']==0)
t=d['TIME'][m]; f=d['FLUX'][m].astype(float); Q=d['QUARTER'][m]
for qq in np.unique(Q): f[Q==qq]/=np.nanmedian(f[Q==qq])
ph=((t-T0+0.5*P)%P)-0.5*P
mask=(np.abs(ph)<0.3)|(np.abs(np.abs(ph)-P/2)<0.3)
r=detrend(t,f,mask,1.0)-1
ok=np.isfinite(r); t,r,ph,Q=t[ok],r[ok],ph[ok],Q[ok]
mad=1.4826*np.median(np.abs(r-np.median(r))); keep=r<5*mad; t,r,ph,Q=t[keep],r[keep],ph[keep],Q[keep]
sig=mad
Tp=np.interp(ph,pc,tp); Ts=np.interp(ph,pc,ts)
def fit(cols):
    A=np.array(cols).T; cov=np.linalg.inv(A.T@A); p=cov@(A.T@r); chi=np.sum((r-A@p)**2)/sig**2; return p,np.sqrt(np.diag(cov))*sig,chi
chi0=np.sum(r**2)/sig**2
pP,eP,cP=fit([Tp]); pS,eS,cS=fit([Ts]); pB,eB,cB=fit([Tp,Ts]); pF,eF,cF=fit([Tp+Ts])
print(f'null chi2 {chi0:.1f}')
print(f'primary-template only : scale {pP[0]*1e3:.3f} +- {eP[0]*1e3:.3f} x1e-3 -> implied dip {pP[0]*prim*1e6:.0f} ppm ; dchi2 vs null {chi0-cP:.1f}')
print(f'secondary-template only: scale {pS[0]*1e3:.3f} +- {eS[0]*1e3:.3f} x1e-3 -> implied dip {pS[0]*sec*1e6:.0f} ppm ; dchi2 {chi0-cS:.1f}')
print(f'both free             : prim scale {pB[0]*1e3:.3f}+-{eB[0]*1e3:.3f}, sec scale {pB[1]*1e3:.3f}+-{eB[1]*1e3:.3f}; ratio sec/prim scale = {pB[1]/pB[0]:.2f} (expected 1.0 for a pure copy); dchi2 {chi0-cB:.1f}')
print(f'single-scale full copy: scale {pF[0]*1e3:.3f}+-{eF[0]*1e3:.3f}; chi2 {cF:.1f} vs both-free {cB:.1f}: dchi2 = {cF-cB:.1f}')
# empirical (red-noise) error on the scale: random-phase template fits
sc=[]
for k in range(400):
    sh=np.random.uniform(0,P); p2=((t-T0-sh+0.5*P)%P)-0.5*P; A=np.interp(p2,pc,tp); sc.append(np.sum(A*r)/np.sum(A*A))
sc=np.array(sc); emp=1.4826*np.median(np.abs(sc-np.median(sc)))
print(f'empirical (random-phase) error on primary scale: {emp*1e3:.3f} x1e-3 -> implied-dip error {emp*prim*1e6:.0f} ppm; prim SNR {pP[0]/emp:.1f}; sec-scale SNR with same error {pB[1]/emp:.1f}; (prim-sec)/err = {(pB[0]-pB[1])/emp/np.sqrt(2):.1f} sigma')
# per-season primary and secondary scale (season=(Q-2)%4)
for s_ in range(4):
    m_=((Q-2)%4)==s_
    if m_.sum()<100: continue
    A=np.c_[Tp[m_],Ts[m_]]; p=np.linalg.lstsq(A,r[m_],rcond=None)[0]
    print(f' season {s_} quarters {sorted(set(Q[m_].astype(int)))}: prim dip {p[0]*prim*1e6:.0f} ppm, sec dip {p[1]*sec*1e6:.0f} ppm')
json.dump(dict(eb_prim=float(prim),eb_sec=float(sec),prim_scale=float(pP[0]),prim_scale_err_white=float(eP[0]),prim_scale_err_emp=float(emp),both=[float(x) for x in pB],chi_null=float(chi0),chi_prim=float(cP),chi_both=float(cB),chi_full=float(cF)),open('eb_template.json','w'),indent=1)
