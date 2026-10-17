import numpy as np, warnings; warnings.filterwarnings("ignore")
import lightkurve as lk
from astropy.io import fits
T="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/mastDownload/HLSP/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits"
lc=lk.read(T); print(type(lc).__name__, len(lc), "cols", lc.columns)
print("meta AUTHOR", lc.meta.get('AUTHOR'), "FLUX_ORIGIN", lc.meta.get('FLUX_ORIGIN'), "TPFORG", lc.meta.get('TPFORG'), "QUALITY_BITMASK", lc.meta.get('QUALITY_BITMASK'))
print("time format/scale", lc.time.format, lc.time.scale, lc.time[0])
print("median flux", np.nanmedian(lc.flux.value), "median sap", np.nanmedian(lc['sap_flux'].value), "median psf_flat", np.nanmedian(lc['psf_flat_flux'].value))
lc2=lk.read(T, flux_column='sap_flux'); print("flux_column='sap_flux' -> median flux", np.nanmedian(lc2.flux.value), "FLUX_ORIGIN", lc2.meta.get('FLUX_ORIGIN'))
lc3=lc.select_flux('sap_flux'); print("select_flux('sap_flux') -> median", np.nanmedian(lc3.flux.value))
for bm in ['none','default','hard','hardest']:
    l=lk.read(T, quality_bitmask=bm); print(f"bitmask {bm:8s}: n={len(l)}")
with fits.open(T) as h: q=h[1].data['SAP_QUALITY']; print("quality==0:", (q==0).sum(), "unique bits:", sorted(set(int(b) for v in np.unique(q) for b in [1<<i for i in range(32)] if v & b)))
print("lightkurve default bitmask value:", lk.utils.KeplerQualityFlags.DEFAULT_BITMASK, "hard:", lk.utils.KeplerQualityFlags.HARD_BITMASK, "hardest:", lk.utils.KeplerQualityFlags.HARDEST_BITMASK)
print("decode of bits present:", {int(v): lk.utils.KeplerQualityFlags.decode(int(v)) for v in np.unique(q)[:12]})
try:
    lq=lk.read(T, quarter=7); print("quarter=7 worked", len(lq))
except Exception as e: print("quarter=7 ->", type(e).__name__, e)
# estimate_cdpp diagnosis
l=lk.read(T, quality_bitmask='hardest').remove_nans()
print("estimate_cdpp raw:", l.estimate_cdpp(), "normalized:", l.normalize().estimate_cdpp(), "normalized ppm unit:", l.normalize(unit='ppm').estimate_cdpp())
ln=l.normalize(); print("normalize unit:", ln.flux.unit, "flux dtype", ln.flux.dtype)
# manual: what estimate_cdpp does internally
flat=ln.flatten(window_length=101, polyorder=2, savgol_polyorder=2) if False else ln.flatten(window_length=101)
print("flatten median", np.nanmedian(flat.flux.value), flat.flux.unit)
cl=flat.remove_outliers(sigma=5.0)
print("remove_outliers n", len(cl))
b=cl.bin(binsize=13)
print("bin(binsize=13) n", len(b), "std flux", np.nanstd(b.flux.value), "unit", b.flux.unit)
print("cl.estimate_cdpp(13):", cl.estimate_cdpp(transit_duration=13))
import inspect; print(inspect.getsource(lk.LightCurve.estimate_cdpp)[-1500:])
