import numpy as np, warnings, time; warnings.filterwarnings("ignore")
from astropy.io import fits
import lightkurve as lk
T="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/mastDownload/HLSP/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits"
with fits.open(T) as h:
    d=h[1].data; good=(d['SAP_QUALITY']==0)&np.isfinite(d['FLUX'])&np.isfinite(d['TIME']); t=d['TIME'][good]; f=d['FLUX'][good]
lc=lk.LightCurve(time=t,flux=f).remove_nans().normalize()
for td in (12,13):
    print("transit_duration", td, "->", lc.estimate_cdpp(transit_duration=td))
lc_b=lk.LightCurve(time=t,flux=f).remove_nans().normalize()
print("time dtype", t.dtype, "flux dtype", f.dtype, "byteorder", f.dtype.byteorder)
# VizieR Kirk+2016
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord; import astropy.units as u
v=Vizier(columns=['KIC','Per','BJD0','Morph','Kpmag','_r'], row_limit=50)
r=v.query_region(SkyCoord(285.907827,40.941988,unit='deg'), radius=5*u.arcmin, catalog='J/AJ/151/68/catalog')
print("Kirk+2016 EBs within 5 arcmin:", r[0] if len(r) else "none")
r2=v.query_region(SkyCoord(285.907827,40.941988,unit='deg'), radius=30*u.arcmin, catalog='J/AJ/151/68/catalog'); print("within 30 arcmin:", len(r2[0]) if len(r2) else 0)
# astroquery MAST
from astroquery.mast import Observations
t0=time.time(); n=Observations.query_criteria_count(provenance_name="KBONUS-BKG"); print("KBONUS-BKG obs count:", n, f"({time.time()-t0:.1f}s)")
obs=Observations.query_criteria(provenance_name="KBONUS-BKG", target_name="Gaia DR3 2103792461826460288"); print(obs[['obs_collection','project','dataproduct_type','target_name','t_exptime','obsid']])
pl=Observations.get_product_list(obs); print(pl[['productFilename','size','dataURI']])
sr=lk.search_lightcurve(SkyCoord(285.907827,40.941988,unit='deg'), author='KBONUS-BKG', radius=20*u.arcsec); print(sr)
