# Are in-transit cadences preferentially near Kepler quality-flagged events (desats, tweaks, etc.)?
import numpy as np, glob
from astropy.io import fits
P=6.703464464595031; t0=265.183161835447; dur=2.4375/24
files=sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr009655236_lc_Q000111011101111000/*.fits.gz'))
T=[];Q=[]
for fn in files:
    d=fits.getdata(fn,1); ok=np.isfinite(d['TIME']); T.append(d['TIME'][ok]); Q.append(d['QUALITY'][ok])
T=np.concatenate(T); Q=np.concatenate(Q)
bits={1:'attitude tweak',2:'safe mode',4:'coarse point',8:'earth point',16:'zero crossing',32:'desaturation',64:'argabrightening',128:'cosmic ray ap',1024:'sudden sens. dropout',4096:'impulsive outlier',8192:'cosmic ray col',32768:'no fine point',65536:'no data',1048576:'thruster'}
ph=((T-t0)/P+0.5)%1-0.5; intr=np.abs(ph*P)<dur/2
rng=np.random.RandomState(0)
print('cadences %d, in-transit %d'%(len(T),intr.sum()))
for b,name in bits.items():
    fl=(Q&b)>0
    if fl.sum()==0: continue
    # fraction of in-transit cadences within +/-0.25 d of a flagged cadence, vs random phases
    tf=T[fl]; near=np.array([np.min(np.abs(tf-x))<0.25 for x in T])
    obs=near[intr].mean(); rs=[]
    for k in range(200):
        p2=((T-t0-rng.uniform(0,P))/P+0.5)%1-0.5; ib=np.abs(p2*P)<dur/2; rs.append(near[ib].mean())
    print('flag %6d %-22s n=%5d | frac in-transit within 0.25 d: %.3f vs random %.3f +/- %.3f (z=%.1f)'%(b,name,fl.sum(),obs,np.mean(rs),np.std(rs),(obs-np.mean(rs))/max(np.std(rs),1e-9)))
