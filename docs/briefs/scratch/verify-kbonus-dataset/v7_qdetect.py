import numpy as np, pandas as pd, glob, warnings; warnings.filterwarnings("ignore")
from astropy.io import fits
import lightkurve as lk
src=pd.read_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/source_catalog.parquet")
S="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/sample/"
for f in sorted(glob.glob(S+"*.fits")):
    with fits.open(f) as h:
        g=h[0].header['GAIAID']; qd=h[0].header['QDETECT']
        qs=sorted(int(x.name.split('_Q')[1]) for x in h if x.name.startswith('LIGHTCURVE_Q'))
        row=src[src.gaia_designation==g].iloc[0]
        hdus=''.join('1' if i in qs else '0' for i in range(18))
        flag=f"avail={row.psf_avail} flag={row.psf_fraction_flag} QDETECT={qd} HDUs={hdus}"
        same = (row.psf_avail==row.psf_fraction_flag)
        if not same or qd!=hdus or qd!=row.psf_avail:
            print(g, flag, "| stitched quarters present:", sorted(set(h[1].data['QUARTER'])))
            # which quarters PSFFRAC
            print("   per-HDU PSFFRAC:", {x.header['QUARTER']: round(x.header['PSFFRAC'],3) for x in h if x.name.startswith('LIGHTCURVE_Q')})
print("done; stars where all equal are not printed")
# Now reproduce the previous agent's estimate_cdpp == 0 claim
T="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/mastDownload/HLSP/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits"
with fits.open(T) as h:
    d=h[1].data; g=d['SAP_QUALITY']==0; t=d['TIME'][g]; f=d['FLUX'][g]
lc=lk.LightCurve(time=t,flux=f).remove_nans().normalize()
flat=lc.flatten(window_length=101, polyorder=2).remove_outliers(sigma=5)
print("prev-agent path: flat.estimate_cdpp(13) =", flat.estimate_cdpp(transit_duration=13), "| lc.estimate_cdpp() =", lc.estimate_cdpp(), "| flux unit:", lc.flux.unit, "flat unit:", flat.flux.unit)
lc2=lk.LightCurve(time=t,flux=f).remove_nans(); print("unnormalized raw arrays: estimate_cdpp =", lc2.estimate_cdpp())
import astropy.units as u
lc3=lk.LightCurve(time=t,flux=f*u.electron/u.s).remove_nans(); print("with e-/s unit: estimate_cdpp =", lc3.estimate_cdpp())
