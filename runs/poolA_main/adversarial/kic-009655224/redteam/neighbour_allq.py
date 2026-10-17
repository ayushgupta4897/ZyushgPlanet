import numpy as np, sys, glob
from astropy.io import fits
P=6.703464464595031; t0=265.183161835447; dur=2.4375/24
np.random.seed(2)
for F in sys.argv[1:]:
    d=fits.getdata(F,'LIGHTCURVE_STITCHED'); t=d['TIME']; f=d['FLUX']; fs=d['SAP_FLUX']; q=d['QUARTER']; qual=d['SAP_QUALITY']
    ok=np.isfinite(t)&np.isfinite(f)&(qual==0)&(f>0); t,f,fs,q=t[ok],f[ok],fs[ok],q[ok]
    ph=((t-t0)/P+0.5)%1-0.5; intr=np.abs(ph*P)<dur/2*1.3
    def runmed(t,f,mask,win=1.0):
        out=np.empty_like(f)
        for i in range(len(t)):
            s=(np.abs(t-t[i])<win/2)&~mask; out[i]=np.median(f[s]) if s.sum()>5 else np.nan
        return out
    fn=np.full_like(f,np.nan); fsn=np.full_like(f,np.nan); med=np.full_like(f,np.nan); meds=np.full_like(f,np.nan)
    for qq in np.unique(q):
        s=q==qq; tr=runmed(t[s],f[s],intr[s]); fn[s]=f[s]/tr-1; med[s]=np.median(f[s]); tr2=runmed(t[s],fs[s],intr[s]); fsn[s]=fs[s]/tr2-1; meds[s]=np.median(fs[s])
    g=np.isfinite(fn)&np.isfinite(fsn); t,fn,fsn,q,ph,med,meds=t[g],fn[g],fsn[g],q[g],ph[g],med[g],meds[g]
    mad=1.4826*np.median(np.abs(fn-np.median(fn))); g=(np.abs(fn)<6*mad); t,fn,fsn,q,ph,med,meds=t[g],fn[g],fsn[g],q[g],ph[g],med[g],meds[g]
    inbox=np.abs(ph*P)<dur/2
    def dep(fl,inb):
        d0=np.median(fl[~inb])-np.mean(fl[inb]); rs=[]
        for k in range(500):
            p2=((t-t0-np.random.uniform(0,P))/P+0.5)%1-0.5; ib=np.abs(p2*P)<dur/2
            if ib.sum()>10: rs.append(np.median(fl[~ib])-np.mean(fl[ib]))
        return d0,np.std(rs)
    dp,ep=dep(fn,inbox); ds,es=dep(fsn,inbox)
    print(F.split('/')[-1][:60],'quarters',sorted(set(q.astype(int))),'n_in',inbox.sum())
    print('  PSF depth at target ephemeris: %.0f +/- %.0f ppm = %.2f +/- %.2f e/s (median PSF flux %.0f e/s)'%(dp*1e6,ep*1e6,dp*np.median(med),ep*np.median(med),np.median(med)))
    print('  SAP depth at target ephemeris: %.0f +/- %.0f ppm = %.2f +/- %.2f e/s (median SAP flux %.0f e/s)'%(ds*1e6,es*1e6,ds*np.median(meds),es*np.median(meds),np.median(meds)))
    for qq in sorted(set(q.astype(int))):
        s=q==qq; dq=np.median(fn[s&~inbox])-np.mean(fn[s&inbox]); eq=np.std(fn[s&~inbox])/np.sqrt((s&inbox).sum()); print('   Q%d PSF %.0f +/- %.0f ppm (n_in %d)'%(qq,dq*1e6,eq*1e6,(s&inbox).sum()))
