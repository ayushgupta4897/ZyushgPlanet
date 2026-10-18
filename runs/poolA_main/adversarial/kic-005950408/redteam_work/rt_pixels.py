"""Pixel-level look at the three cadences that carry the kic-005950408 'signal' (epochs 0, 10, 29 at the
pipeline ephemeris), on the host TPF (KIC 5950394, Kp 8.45, saturated). Own code; no pipeline/ imports.
For each carrying epoch: difference image (in-box cadences minus local out-of-box median), the saturated/
bleed column geometry, the target's pixel position (from KBONUS ROW/COLUMN headers), and the summed flux
in the target's KBONUS aperture vs the whole stamp vs the bleed column, cadence by cadence.
"""
import sys, os, json, glob, warnings; warnings.filterwarnings('ignore')
import numpy as np
from astropy.io import fits
OUT = sys.argv[1]
P = 38.60671217842687; T0 = 151.4674802687019
KB = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005950408_kepler_v1.0_lc.fits'
kb = fits.open(KB)
res = {}
files = {1: glob.glob(OUT + '/tpf/mastDownload/Kepler/*/kplr005950394-2009166043257_lpd-targ.fits.gz')[0],
         5: glob.glob(OUT + '/tpf/mastDownload/Kepler/*/kplr005950394-2010174085026_lpd-targ.fits.gz')[0],
         13: glob.glob(OUT + '/tpf/mastDownload/Kepler/*/kplr005950394-2012179063303_lpd-targ.fits.gz')[0]}
epochs = {1: 0, 5: 10, 13: 29}
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axes = plt.subplots(3, 4, figsize=(18, 13))
for k, (qq, fn) in enumerate(files.items()):
    h = fits.open(fn); d = h[1].data; hd = h[1].header; ap = h[2].data
    t = d['TIME'].astype(float); fl = d['FLUX'].astype(float); qual = d['QUALITY'].astype(int); cad = d['CADENCENO'].astype(int)
    col0 = int(h[2].header['CRVAL1P']); row0 = int(h[2].header['CRVAL2P'])
    kh = kb[f'LIGHTCURVE_Q{qq}'].header; kap = kb[f'APERTURE_Q{qq}'].data
    tcol = float(kh['COLUMN']) - col0; trow = float(kh['ROW']) - row0   # target position in stamp pixel coords (x=col, y=row)
    hcol = float(hd.get('1CRV4P', col0)); # not used
    tc = T0 + epochs[qq]*P
    ok = np.isfinite(t) & np.isfinite(np.nansum(fl, axis=(1, 2))) & (np.nansum(fl, axis=(1, 2)) > 0)
    inb = ok & (np.abs(t - tc) < 0.75/48); loc = ok & (np.abs(t - tc) < 1.0) & (np.abs(t - tc) > 0.1)
    med = np.nanmedian(fl[loc], axis=0); diff = np.nanmean(fl[inb], axis=0) - med
    sat = med > 0.9*np.nanmax(med)
    # bleed column(s): columns whose median flux exceeds 1e5 e/s in more than 3 pixels
    colsum = np.nansum(med, axis=0); bleedcols = np.where((med > 1e5).sum(axis=0) >= 3)[0]
    # cadence-by-cadence time series around the epoch
    win = ok & (np.abs(t - tc) < 0.6)
    stamp_tot = np.nansum(fl[win], axis=(1, 2)); kap_tot = np.nansum(fl[win][:, kap.astype(bool)], axis=1)
    bleed_tot = np.nansum(fl[win][:, :, bleedcols], axis=(1, 2)) if len(bleedcols) else np.zeros(win.sum())
    rest_tot = stamp_tot - bleed_tot
    def rel(x): return x/np.nanmedian(x) - 1
    ser = [dict(dt_h=float((tt - tc)*24), stamp_rel_ppm=float(rel(stamp_tot)[i]*1e6), kbonus_ap_rel_ppm=float(rel(kap_tot)[i]*1e6), bleedcol_rel_ppm=float(rel(bleed_tot)[i]*1e6) if len(bleedcols) else None,
                non_bleed_rel_ppm=float(rel(rest_tot)[i]*1e6), inbox=bool(np.abs(tt - tc) < 0.75/48), cadenceno=int(cc), quality=int(qq_)) for i, (tt, cc, qq_) in enumerate(zip(t[win], cad[win], qual[win]))]
    # where is the deficit? significance-weighted
    noise = np.nanstd(fl[loc], axis=0)/np.sqrt(max(inb.sum(), 1))
    snr = diff/np.where(noise > 0, noise, np.nan)
    neg = np.where(np.isfinite(snr) & (snr < -3)); pos = np.where(np.isfinite(snr) & (snr > 3))
    res[f'Q{qq}'] = dict(epoch=epochs[qq], tc=tc, n_inbox=int(inb.sum()), stamp_shape=list(fl.shape[1:]), target_xy_px=[tcol, trow], kbonus_ap_pixels=[[int(a), int(b)] for a, b in zip(*np.where(kap))],
                         bleed_columns_x=[int(c) for c in bleedcols], n_saturated_px=int(sat.sum()), target_dist_to_nearest_bleedcol_px=float(np.min(np.abs(bleedcols - tcol))) if len(bleedcols) else None,
                         kbonus_ap_median_flux=float(np.nansum(med[kap.astype(bool)])), stamp_median_flux=float(np.nansum(med)), bleed_median_flux=float(np.nansum(med[:, bleedcols])) if len(bleedcols) else None,
                         diff_image_sum=float(np.nansum(diff)), diff_in_kbonus_ap=float(np.nansum(diff[kap.astype(bool)])), diff_in_bleedcols=float(np.nansum(diff[:, bleedcols])) if len(bleedcols) else None,
                         pixels_snr_below_m3=[[int(a), int(b), float(snr[a, b])] for a, b in zip(*neg)], pixels_snr_above_3=[[int(a), int(b), float(snr[a, b])] for a, b in zip(*pos)],
                         series=ser)
    ax = axes[k]
    im = ax[0].imshow(np.log10(np.clip(med, 1, None)), origin='lower', cmap='viridis'); ax[0].set_title(f'Q{qq} log10 median flux; x=target, box=KBONUS ap'); plt.colorbar(im, ax=ax[0])
    ax[0].plot(tcol, trow, 'rx', ms=12, mew=2)
    for a, b in zip(*np.where(kap)): ax[0].add_patch(plt.Rectangle((b-0.5, a-0.5), 1, 1, fill=False, ec='r', lw=1))
    v = np.nanmax(np.abs(np.nan_to_num(diff))); im = ax[1].imshow(diff, origin='lower', cmap='RdBu_r', vmin=-v, vmax=v); ax[1].set_title(f'epoch {epochs[qq]} diff image (in-box - local median) e/s'); plt.colorbar(im, ax=ax[1])
    ax[1].plot(tcol, trow, 'kx', ms=12, mew=2)
    im = ax[2].imshow(np.clip(snr, -8, 8), origin='lower', cmap='RdBu_r', vmin=-8, vmax=8); ax[2].set_title('diff SNR per pixel'); plt.colorbar(im, ax=ax[2]); ax[2].plot(tcol, trow, 'kx', ms=12, mew=2)
    dts = [s['dt_h'] for s in ser]
    ax[3].plot(dts, [s['kbonus_ap_rel_ppm'] for s in ser], 'o-', label='KBONUS aperture (4 px)'); ax[3].plot(dts, [s['stamp_rel_ppm'] for s in ser], 's-', label='whole stamp')
    if len(bleedcols): ax[3].plot(dts, [s['bleedcol_rel_ppm'] for s in ser], '^-', label='bleed column(s)')
    ax[3].plot(dts, [s['non_bleed_rel_ppm'] for s in ser], 'v-', label='stamp minus bleed cols')
    ax[3].axvspan(-0.375, 0.375, color='C1', alpha=0.3); ax[3].set_xlabel('hours from pipeline mid-transit'); ax[3].set_ylabel('ppm'); ax[3].legend(fontsize=8); ax[3].set_title('raw pixel sums around the epoch')
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'rt_pixels.png'), dpi=100)
json.dump(res, open(os.path.join(OUT, 'rt_pixels.json'), 'w'), indent=1)
for qq in res:
    r = res[qq]; print(qq, 'epoch', r['epoch'], 'n_inbox', r['n_inbox'], 'target_xy', r['target_xy_px'], 'bleedcols', r['bleed_columns_x'], 'dist', r['target_dist_to_nearest_bleedcol_px'], 'nsat', r['n_saturated_px'])
    print('  ap flux', r['kbonus_ap_median_flux'], 'stamp', r['stamp_median_flux'], 'bleed', r['bleed_median_flux'], '| diff sum', r['diff_image_sum'], 'in ap', r['diff_in_kbonus_ap'], 'in bleed', r['diff_in_bleedcols'])
    print('  neg px', r['pixels_snr_below_m3'][:12]); print('  pos px', r['pixels_snr_above_3'][:12])
    for s in r['series']:
        if abs(s['dt_h']) < 2.5: print('   ', {k: (round(v, 1) if isinstance(v, float) else v) for k, v in s.items()})
