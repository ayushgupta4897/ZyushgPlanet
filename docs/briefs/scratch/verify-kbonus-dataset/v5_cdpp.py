import numpy as np, warnings; warnings.filterwarnings("ignore")
import lightkurve as lk, inspect
from astropy.io import fits
T="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/mastDownload/HLSP/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits"
with fits.open(T) as h: q=h[1].data['SAP_QUALITY'].astype(np.int64)
bits=sorted(set(int(b) for v in np.unique(q) for b in [1<<i for i in range(32)] if v & b))
print("quality==0:", (q==0).sum(), "bits present:", bits)
print("decoded:", {b: lk.utils.KeplerQualityFlags.STRINGS.get(b) for b in bits})
print("DEFAULT_BITMASK", lk.utils.KeplerQualityFlags.DEFAULT_BITMASK, "HARD", lk.utils.KeplerQualityFlags.HARD_BITMASK, "HARDEST", lk.utils.KeplerQualityFlags.HARDEST_BITMASK)
print("bits in default mask:", [b for b in bits if b & lk.utils.KeplerQualityFlags.DEFAULT_BITMASK])
l=lk.read(T, quality_bitmask='hardest').remove_nans()
print("estimate_cdpp raw:", l.estimate_cdpp(), "| normalized:", l.normalize().estimate_cdpp(), "| normalize(unit='ppm'):", l.normalize(unit='ppm').estimate_cdpp())
ln=l.normalize(); print("normalized unit:", ln.flux.unit)
flat=ln.flatten(window_length=101); cl=flat.remove_outliers(sigma=5.0)
b=cl.bin(binsize=13); print("bin binsize=13: n", len(b), "std", np.nanstd(b.flux.value), "unit", b.flux.unit, "n finite", np.isfinite(b.flux.value).sum())
print("bin(time_bin_size=13*29.4min):", end=' '); import astropy.units as u
b2=cl.bin(time_bin_size=13*29.4245*u.min); print(len(b2), np.nanstd(b2.flux.value))
src=inspect.getsource(lk.LightCurve.estimate_cdpp); print(src[src.find('def estimate_cdpp'):][:400]); print('...'); print(src[-900:])
# official kepler neighbour to check estimate_cdpp there
import glob
f=glob.glob("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/mastDownload/Kepler/kplr005689116_lc_Q*/*llc.fits")[0]
k=lk.read(f).remove_nans(); print("official Kepler file", f.split('/')[-1], "estimate_cdpp:", k.estimate_cdpp(), "header CDPP6_0:", k.meta.get('CDPP6_0'))
print("manual bin on official:", np.nanstd(k.normalize().flatten(window_length=101).remove_outliers(sigma=5).bin(binsize=13).flux.value)*1e6)
