"""Empirical Kepler flux zero-point from other cached KBONUS-BKG files (header GMAG vs median PSF FLUX),
to test whether kic-001577049's PSF flux is consistent with its Gaia G magnitude."""
import glob, os, json, random, sys
import numpy as np
from astropy.io import fits
D = '/home/user/ZyushgPlanet/data/kbonus_batch'
TARGET = 'hlsp_kbonus-bkg_kepler_kepler_kic-001577049_kepler_v1.0_lc.fits'
files = sorted(glob.glob(D + '/*.fits'))
random.seed(1)
sample = random.sample([f for f in files if not f.endswith(TARGET)], min(250, len(files)-1))
rows = []
def one(fn):
    with fits.open(fn, memmap=True) as f:
        h = f[0].header
        g = h.get('GMAG', np.nan); kp = h.get('KEPMAG', np.nan)
        d = f['LIGHTCURVE_STITCHED'].data
        m = (d['SAP_QUALITY'] == 0) & np.isfinite(d['FLUX']) & np.isfinite(d['SAP_FLUX'])
        if m.sum() < 100: return None
        return dict(file=os.path.basename(fn), gmag=float(g) if g is not None else np.nan,
                    kepmag=float(kp) if kp is not None else np.nan,
                    med_psf=float(np.median(d['FLUX'][m])), med_sap=float(np.median(d['SAP_FLUX'][m])))
for fn in sample:
    try:
        r = one(fn)
        if r: rows.append(r)
    except Exception as e:
        print('skip', fn, e, file=sys.stderr)
tgt = one(os.path.join(D, TARGET))
G = np.array([r['gmag'] for r in rows]); Fp = np.array([r['med_psf'] for r in rows]); Fs = np.array([r['med_sap'] for r in rows])
ok = np.isfinite(G) & (Fp > 0)
zp_psf = G[ok] + 2.5*np.log10(Fp[ok])          # G = zp - 2.5 log10(F_psf)
zp_med = float(np.median(zp_psf)); zp_mad = float(1.4826*np.median(np.abs(zp_psf - zp_med)))
pred_G_from_psf = zp_med - 2.5*np.log10(tgt['med_psf'])
pred_G_from_sap = zp_med - 2.5*np.log10(tgt['med_sap'])
excess_mag = tgt['gmag'] - pred_G_from_psf
out = dict(n_calibrators=int(ok.sum()), zp_psf_median=zp_med, zp_psf_robust_sigma=zp_mad,
           zp_psf_16_84=[float(np.percentile(zp_psf,16)), float(np.percentile(zp_psf,84))],
           target=tgt, target_G_implied_by_psf_flux=float(pred_G_from_psf),
           target_G_implied_by_sap_flux=float(pred_G_from_sap),
           target_excess_mag_psf=float(excess_mag),
           target_psf_flux_over_expected=float(10**(0.4*excess_mag)),
           target_sap_over_psf=tgt['med_sap']/tgt['med_psf'],
           sample_sap_over_psf_median=float(np.median(Fs[ok]/Fp[ok])),
           sample_sap_over_psf_16_84=[float(np.percentile(Fs[ok]/Fp[ok],16)), float(np.percentile(Fs[ok]/Fp[ok],84))],
           note='zp from other pool stars in the cache: G + 2.5log10(median PSF FLUX). Positive excess_mag = target PSF flux brighter than its Gaia G predicts.')
print(json.dumps(out, indent=1))
json.dump(dict(summary=out, calibrators=rows), open(sys.argv[1], 'w'), indent=1)
