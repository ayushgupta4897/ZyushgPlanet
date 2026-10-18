import numpy as np, json, sys
from astropy.io import fits
from astropy.timeseries import LombScargle
P=2.377014357204058; t0=132.15453895516694; dur=1.125/24
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-009520080_kepler_v1.0_lc.fits'
h=fits.open(f); d=h['LIGHTCURVE_STITCHED'].data
t=d['TIME']; fl=d['FLUX']; fe=d['FLUX_ERR']; q=d['QUARTER']; qual=d['SAP_QUALITY']; sap=d['SAP_FLUX']
m=np.isfinite(t)&np.isfinite(fl)&(qual==0)&(fl>0)
t,fl,fe,q,sap=t[m],fl[m],fe[m],q[m],sap[m]
out={}
# per-quarter normalise, masked running median (1 d), mask in-transit
ph=((t-t0)/P+0.5)%1-0.5
intr=np.abs(ph*P)<0.75*dur
def detrend(flux,mask):
    y=flux.copy()
    for qq in np.unique(q):
        s=q==qq; y[s]=y[s]/np.nanmedian(y[s])
    tr=np.full_like(y,np.nan)
    for i in range(len(t)):
        w=(np.abs(t-t[i])<0.5)&(~mask)
        tr[i]=np.median(y[w]) if w.sum()>5 else np.nan
    return y/tr-1
r=detrend(fl,intr); ok=np.isfinite(r)
t,r,fe,q,ph,sap,fl=t[ok],r[ok],fe[ok],q[ok],ph[ok],sap[ok],fl[ok]
sig=1.4826*np.median(np.abs(r-np.median(r)))
clip=(r<5*sig)&(r>-15*sig); t,r,fe,q,ph,sap,fl=t[clip],r[clip],fe[clip],q[clip],ph[clip],sap[clip],fl[clip]
print('N',len(t),'robust sigma ppm',sig*1e6)
# event 342 : raw look
tc=t0+342*P
w=np.abs(t-tc)<1.0
print('--- event 342 (tc BKJD %.3f, Q%d) raw PSF flux and SAP around it:'%(tc,np.median(q[w])))
for ti,ri,si,fi in zip(t[w],r[w],sap[w],fl[w]): print(f'{ti:.4f} dt_h={(ti-tc)*24:+6.2f} rel_ppm={ri*1e6:+8.0f} PSF={fi:.1f} SAP={si:.1f}')
# phase-0.5 robustness: full phase scan of box mean (box=dur) using mean, median, and 3-sigma-clipped mean
def boxstat(center):
    s=np.abs(((ph-center+0.5)%1)-0.5)*P<dur/2
    x=r[s]*1e6
    xc=x[np.abs(x-np.median(x))<3*sig*1e6]
    return x.mean(),np.median(x),xc.mean(),len(x)
centers=np.arange(0.05,0.96,0.005)
scan=np.array([boxstat(c) for c in centers])
oot=np.abs(ph)>0.1
print('--- phase scan (box = transit duration), excluding |phase|<0.05 of transit; stats of means: std %.1f ppm'%scan[:,0].std())
i=np.argmin(scan[:,0]); print('most negative (brightest) box: phase %.3f mean %.0f median %.0f clipped %.0f n %d'%(centers[i],*scan[i]))
j=np.argmax(scan[:,0]); print('most positive (deepest) box: phase %.3f mean %.0f median %.0f clipped %.0f n %d'%(centers[j],*scan[j]))
rank=(scan[:,0]<=scan[np.argmin(np.abs(centers-0.5)),0]).sum()
k=np.argmin(np.abs(centers-0.5)); print('at phase 0.500: mean %.0f median %.0f clipped %.0f n %d ; rank among %d boxes: %d'%(*scan[k],len(centers),rank))
# z-score of phase 0.5 relative to scan distribution excluding +-0.03 of transit
far=np.abs(centers-0.5)>0.03
print('z of phase-0.5 box vs scan distribution (excl. near 0.5): %.2f'%((scan[k,0]-scan[far,0].mean())/scan[far,0].std()))
# bootstrap error at phase 0.5 over events
s05=np.abs(((ph-0.5+0.5)%1)-0.5)*P<dur/2
ev=np.round((t[s05]-t0)/P).astype(int); ue=np.unique(ev)
rng=np.random.default_rng(1); bs=[]
for _ in range(2000):
    pick=rng.choice(ue,len(ue)); bs.append(np.mean(np.concatenate([r[s05][ev==e] for e in pick]))*1e6)
print('phase-0.5 box: mean %.0f ppm, event-bootstrap err %.0f ppm, n events %d'%(r[s05].mean()*1e6,np.std(bs),len(ue)))
# top contributors to the phase-0.5 brightening: per-event means
em=np.array([r[s05][ev==e].mean()*1e6 for e in ue]); 
print('per-event phase-0.5 means: median %.0f, mean %.0f, 5 most negative %s'%(np.median(em),em.mean(),np.sort(em)[:5].round(0)))
print('mean excluding 5 most negative events: %.0f'%np.sort(em)[5:].mean())
# per-season (Kepler season = quarter mod 4) phase-0.5
for sn in range(4):
    ss=s05&(((q-2)%4)==sn); print('season',sn,'phase-0.5 mean %.0f +- %.0f n %d'%(r[ss].mean()*1e6, r[ss].std()*1e6/np.sqrt(ss.sum()), ss.sum()))
# in-transit depth per season for reference
sin=np.abs(ph*P)<dur/2
for sn in range(4):
    ss=sin&(((q-2)%4)==sn); print('season',sn,'transit depth %.0f +- %.0f n %d'%(-r[ss].mean()*1e6, r[ss].std()*1e6/np.sqrt(ss.sum()), ss.sum()))
# LS both polarities on out-of-transit residuals, and on 2-f
ootm=~sin
fr,pw=LombScargle(t[ootm],r[ootm]).autopower(minimum_frequency=1/50,maximum_frequency=1/0.1,samples_per_peak=5)
i=np.argmax(pw); print('LS oot: best period %.4f d power %.4g; power at P %.4g, P/2 %.4g, 2P %.4g'%(1/fr[i],pw[i],pw[np.argmin(np.abs(fr-1/P))],pw[np.argmin(np.abs(fr-2/P))],pw[np.argmin(np.abs(fr-0.5/P))]))
amp=np.sqrt(pw[i])*np.std(r[ootm])*2; print('approx LS amplitude ppm',amp*1e6)
# odd/even timing & duration: fit box centre by scanning offsets
def boxdepth(sel,off,w=dur):
    s=sel&(np.abs(ph*P-off)<w/2); return r[s].mean()*1e6,s.sum()
ev_all=np.round((t-t0)/P).astype(int)
for lab,sel in [('odd',ev_all%2==1),('even',ev_all%2==0)]:
    offs=np.arange(-0.04,0.0401,0.002); dd=np.array([boxdepth(sel,o)[0] for o in offs])
    # parabola around min
    i=np.argmax(dd); print(lab,'best offset %.3f d (%.1f min), depth %.0f ppm'%(offs[i],offs[i]*1440,dd[i]))
    for w in (0.75,1.0,1.125,1.5,2.0):
        print('   width %.3f h depth %.0f'%(w,boxdepth(sel,0,w/24)[0]))
# SAP check at ephemeris
rs=detrend(sap,intr); oks=np.isfinite(rs); 
print('SAP depth %.0f ppm vs PSF %.0f ppm'%(-rs[oks&sin].mean()*1e6,-r[sin].mean()*1e6))
