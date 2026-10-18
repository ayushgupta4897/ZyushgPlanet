import warnings; warnings.filterwarnings('ignore')
import lightkurve as lk, numpy as np, json
sr=lk.search_lightcurve('KIC 11285625',mission='Kepler',cadence='long'); print(sr)
lcs=sr.download_all(quality_bitmask='default')
t=[];f=[];Q=[]
for lc in lcs:
    lc=lc.remove_nans(); x=np.asarray(lc.sap_flux.value,float); t.append(np.asarray(lc.time.value,float)); f.append(x/np.nanmedian(x)); Q.append(np.full(len(x),float(lc.meta.get('QUARTER',-1))))
t=np.concatenate(t); f=np.concatenate(f); Q=np.concatenate(Q)
np.save('eb_lc.npy',np.vstack([t,f,Q]))
Pk=5.395222; t0k=142.33451-5.395222  # EB catalogue P/2 and BJD0-P/2 so that 2P fold centres on the primary
for P,lab in [(2*Pk,'2P'),(Pk,'P')]:
    ph=((t-t0k+0.5*P)%P)-0.5*P
    bins=np.linspace(-0.5*P,0.5*P,401); idx=np.digitize(ph,bins)
    prof=np.array([np.nanmedian(f[idx==i]) if (idx==i).sum()>5 else np.nan for i in range(1,len(bins))])
    pc=0.5*(bins[1:]+bins[:-1])
    print(lab,'min at phase(h)',pc[np.nanargmin(prof)]*24,'depth',1-np.nanmin(prof))
    if lab=='2P':
        s=np.abs(pc)<0.6; print(' primary window (h, rel flux):'); print(np.round(np.c_[pc[s]*24,prof[s]],4))
        s=np.abs(np.abs(pc)-P/2)<0.6; print(' phase-0.5 window:'); print(np.round(np.c_[(np.abs(pc[s])-P/2)*24,prof[s]],4))
# odd/even at 2P: eclipses at t0k+n*Pk, n even vs odd
ph=((t-t0k+0.5*Pk)%Pk)-0.5*Pk; n=np.round((t-t0k)/Pk).astype(int)
for lab,s in [('even(n%2==0)',n%2==0),('odd',n%2==1)]:
    i=s&(np.abs(ph)<1.0/24); print(lab,'core median depth',1-np.nanmedian(f[i]),'n',i.sum())
# eclipse duration estimate at 2P fold: where prof<0.995
P=2*Pk; ph=((t-t0k+0.5*P)%P)-0.5*P
for lab,c in [('primary',0.0),('secondary',P/2)]:
    pp=((ph-c+0.5*P)%P)-0.5*P; bins=np.linspace(-0.5,0.5,121); idx=np.digitize(pp,bins)
    prof=np.array([np.nanmedian(f[idx==i]) if (idx==i).sum()>5 else np.nan for i in range(1,len(bins))]); pc=0.5*(bins[1:]+bins[:-1])
    inecl=pc[prof<0.998]; print(lab,'depth',round(1-np.nanmin(prof),4),'T(flux<0.998) h',round((inecl.max()-inecl.min())*24,2) if len(inecl) else None, 'centre h',round(pc[np.nanargmin(prof)]*24,2))
