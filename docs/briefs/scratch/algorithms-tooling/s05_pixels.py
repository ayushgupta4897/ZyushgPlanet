"""Pixel-level tools on the real neighbour TPF (KIC 5689116, the TPF that hosts KIC 5689121) + TESScut."""
import os, time, warnings, numpy as np; warnings.filterwarnings('ignore')
import lightkurve as lk, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/data'
t0=time.time()
sr=lk.search_targetpixelfile('KIC 5689116', author='Kepler', cadence='long', quarter=10)
print(sr); tpf=sr.download(download_dir=D, quality_bitmask='default')
print('TPF download %.1fs'%(time.time()-t0), tpf, 'shape', tpf.shape, 'channel', tpf.channel, 'col,row', tpf.column, tpf.row, 'quality_bitmask default removed', (tpf.hdu[1].data['QUALITY']!=0).sum(),'of',len(tpf.hdu[1].data))
# where are the two stars?
from astropy.coordinates import SkyCoord
import astropy.units as u
targ=SkyCoord(285.907827,40.941988,unit='deg'); nb=SkyCoord(285.905233,40.942079,unit='deg')
xy_t=tuple(float(v) for v in tpf.wcs.world_to_pixel(targ)); xy_n=tuple(float(v) for v in tpf.wcs.world_to_pixel(nb))
print('pixel coords (relative to cutout): target KIC 5689121', xy_t, ' neighbour KIC 5689116', xy_n, 'sep arcsec', targ.separation(nb).arcsec)
# 1. plot_pixels
ax=tpf.plot_pixels(aperture_mask=tpf.pipeline_mask, show_flux=True); plt.savefig(f'{D}/q10_plot_pixels.png',dpi=80); plt.close('all')
# 2. SAP light curve and centroids
lc=tpf.to_lightcurve(aperture_mask='all').remove_nans().normalize()
col,row=tpf.estimate_centroids(aperture_mask='all')
print('estimate_centroids: median col,row', np.nanmedian(col.value), np.nanmedian(row.value))
# 3. difference image at Kuber's ephemeris (T0 BJD_TDB 2455630.245431 -> BKJD = -2454833)
P=13.769448; T0=2455630.245431-2454833; dur=2.3/24
t=tpf.time.value; ph=((t-T0+P/2)%P)-P/2
intr=np.abs(ph)<dur/2; oot=(np.abs(ph)>dur)&(np.abs(ph)<3*dur)
print('Q10 cadences in transit', intr.sum(), 'out-of-transit comparison', oot.sum(), 'n transits in Q10', len(np.unique(np.round((t[intr]-T0)/P))))
flux=tpf.flux.value
diff=np.nanmean(flux[oot],axis=0)-np.nanmean(flux[intr],axis=0)   # positive where light was lost
err=np.sqrt(np.nanvar(flux[oot],axis=0)/oot.sum()+np.nanvar(flux[intr],axis=0)/intr.sum())
snr=diff/err
iy,ix=np.unravel_index(np.nanargmax(snr),snr.shape)
print('difference-image: brightest-loss pixel (x,y)=',(ix,iy),'SNR %.1f'%snr[iy,ix],'| SNR at target px %.1f | at neighbour px %.1f'%(snr[int(round(xy_t[1])),int(round(xy_t[0]))],snr[int(round(xy_n[1])),int(round(xy_n[0]))]))
w=np.clip(diff,0,None); yy,xx=np.indices(diff.shape); cx=(w*xx).sum()/w.sum(); cy=(w*yy).sum()/w.sum()
print('flux-weighted centroid of positive difference image (x,y)=(%.2f,%.2f); target at (%.2f,%.2f) neighbour at (%.2f,%.2f)'%(cx,cy,xy_t[0],xy_t[1],xy_n[0],xy_n[1]))
fig,axs=plt.subplots(1,3,figsize=(12,3.5))
for a,im,tt in zip(axs,[np.log10(np.nanmean(flux,axis=0)),diff,snr],['mean (log)','out - in (e-/s)','SNR']):
    m=a.imshow(im,origin='lower'); a.plot(*xy_t,'rx'); a.plot(*xy_n,'c+'); a.set_title(tt); fig.colorbar(m,ax=a)
fig.savefig(f'{D}/q10_diffimage.png',dpi=90)
# 4. Kepler PRF model at this channel/position
from lightkurve.prf import KeplerPRF
prf=KeplerPRF(channel=tpf.channel, shape=tpf.shape[1:], column=tpf.column, row=tpf.row)
model=prf(center_col=tpf.column+xy_t[0], center_row=tpf.row+xy_t[1], flux=1.0, scale_col=1, scale_row=1, rotation_angle=0)
print('KeplerPRF evaluate -> shape', model.shape, 'sum', model.sum().round(3), 'peak pixel', np.unravel_index(model.argmax(),model.shape))
# 5. TESScut FFI cutout for the same star
t1=time.time()
st=lk.search_tesscut('TIC 120686465', sector=14); print(st)
tc=st.download(cutout_size=11, download_dir=D)
print('TESScut 11x11 sector 14: %.1fs'%(time.time()-t1), tc, tc.shape, 'cadence (min)', np.median(np.diff(tc.time.value))*1440)
