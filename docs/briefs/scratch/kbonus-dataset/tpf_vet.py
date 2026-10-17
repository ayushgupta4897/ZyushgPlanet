"""Pixel-level vetting with the ORIGINATING Kepler TPF (KIC 5689116) at Kuber's ephemeris.
Difference image (out - in) per quarter, summed over events; flux-weighted centroid shift; compare to Gaia positions via WCS."""
import numpy as np, warnings; warnings.filterwarnings("ignore")
import lightkurve as lk
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
D="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/"
P=13.769448; T0=2455630.245431-2454833; dur=2.3/24
tgt=SkyCoord(285.907827,40.941988,unit='deg'); nb=SkyCoord(285.905233,40.942079,unit='deg')
sr=lk.search_targetpixelfile("KIC 5689116", author="Kepler", cadence="long")
print(sr)
tpfs=sr.download_all(download_dir=D)
tot_diff=None; rows=[]
for tpf in tpfs:
    q=tpf.quarter; t=tpf.time.value; F=tpf.flux.value  # (N, ny, nx) e-/s
    good=(tpf.quality==0)&np.isfinite(t)
    ph=((t-T0+P/2)%P)-P/2
    intr=good&(np.abs(ph)<dur/2*0.8); oot=good&(np.abs(ph)>dur/2+0.05)&(np.abs(ph)<0.5)
    if intr.sum()<3: continue
    # per-event local baseline: for each event, median of nearby oot cadences
    n=np.round((t-T0)/P); diffs=[]
    for k in np.unique(n[intr]):
        i=intr&(n==k); o=oot&(n==k)
        if i.sum()<2 or o.sum()<10: continue
        diffs.append(np.nanmedian(F[o],axis=0)-np.nanmean(F[i],axis=0))
    if len(diffs)<2: continue
    diff=np.nanmean(diffs,axis=0); mean=np.nanmedian(F[oot],axis=0)
    w=WCS(tpf.hdu[2].header)
    xt=w.all_world2pix([[tgt.ra.deg,tgt.dec.deg]],0)[0]; xn=w.all_world2pix([[nb.ra.deg,nb.dec.deg]],0)[0]
    ny,nx=diff.shape; Y,X=np.mgrid[0:ny,0:nx]; m=np.isfinite(diff)&np.isfinite(mean)
    # flux-weighted centroid of difference image (positive = dimming) and of mean image
    dpos=np.clip(diff,0,None)
    cd=((dpos*X)[m].sum()/dpos[m].sum(), (dpos*Y)[m].sum()/dpos[m].sum()); cm=((mean*X)[m].sum()/mean[m].sum(), (mean*Y)[m].sum()/mean[m].sum())
    pk=np.unravel_index(np.nanargmax(np.where(m,diff,-np.inf)),diff.shape)
    frac=np.nansum(diff)/np.nansum(mean)
    rows.append((q,len(diffs),frac*1e6,cd,cm,xt,xn,pk))
    print(f"Q{q:2d} ch{tpf.channel} events={len(diffs)} shape={diff.shape} summed-dip={frac*1e6:6.0f} ppm of stamp flux | diff-image centroid (x,y)=({cd[0]:.2f},{cd[1]:.2f}) | peak-dip pixel (x,y)=({pk[1]},{pk[0]}) | target WCS=({xt[0]:.2f},{xt[1]:.2f}) neighbour WCS=({xn[0]:.2f},{xn[1]:.2f}) | mean-image centroid=({cm[0]:.2f},{cm[1]:.2f})")
    if tot_diff is None: tot_diff=[]
    tot_diff.append((q,diff,mean,xt,xn))
# distance summary
dt=[np.hypot(r[3][0]-r[5][0],r[3][1]-r[5][1]) for r in rows]; dn=[np.hypot(r[3][0]-r[6][0],r[3][1]-r[6][1]) for r in rows]
print(f"\nDifference-image centroid distance to TARGET: median {np.median(dt):.2f} px; to NEIGHBOUR: median {np.median(dn):.2f} px (1 px = 3.98 arcsec); quarters where closer to target: {sum(a<b for a,b in zip(dt,dn))}/{len(dt)}")
# print one difference image
q,diff,mean,xt,xn=tot_diff[0]
np.set_printoptions(precision=0,suppress=True,linewidth=150)
print(f"\nQ{q} mean image (e-/s):\n", mean); print(f"Q{q} difference image out-in (e-/s), positive=dimming:\n", diff)
print(f"target at (x,y)=({xt[0]:.2f},{xt[1]:.2f}); neighbour at ({xn[0]:.2f},{xn[1]:.2f}); array index is [y,x]")
