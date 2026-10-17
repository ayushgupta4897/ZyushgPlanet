import numpy as np, warnings; warnings.filterwarnings("ignore")
from astropy.io import fits
import lightkurve as lk, astropy
print("lightkurve", lk.__version__, "astropy", astropy.__version__, "numpy", np.__version__)
T="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/mastDownload/HLSP/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits"
with fits.open(T) as h:
    h.info()
    hd=h[0].header
    for k in ['DOI','HLSPID','HLSPVER','LICENSE','LICENURL','ORIGIN','VERSION','PROCVER','DATE','TIMESYS','TIMEUNIT','OBJECT','KEPLERID','TPFORG','GAIAID','RA_OBJ','DEC_OBJ','EQUINOX','PMRA','PMDEC','PARALLAX','GMAG','BPMAG','RPMAG','KEPMAG','QDETECT','XPOSURE','MJD-BEG','MJD-END','APERTURE','TELESCOP','INSTRUME','RADESYS','AUTHOR','CHECKSUM','DATASUM']:
        if k in hd: print(f"{k:9s} = {hd[k]!r}  / {hd.comments[k]}")
    print("all primary keys:", list(hd.keys()))
    print("stitched cols:", h[1].columns.names)
    print("stitched header keys:", [k for k in h[1].header.keys() if k not in ('XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','TFIELDS')][:60])
    for k in h[1].header: 
        if k.startswith('TUNIT') or k.startswith('TTYPE'): print(k, h[1].header[k])
    q=h['LIGHTCURVE_Q7']; print("Q7 cols:", q.columns.names); 
    print("Q7 header:"); print(repr(q.header)[:4000])
    print("APERTURE_Q7 shape", h['APERTURE_Q7'].data.shape, "sum", h['APERTURE_Q7'].data.sum(), "NPIXSAP", q.header['NPIXSAP'])
    print("APERTURE_Q7 header:"); print(repr(h['APERTURE_Q7'].header)[:2000])
    d=h[1].data
    print("N rows", len(d), "time min/max", np.nanmin(d['TIME']), np.nanmax(d['TIME']))
    print("per-quarter season/module/channel:", [(x.header['QUARTER'],x.header['SEASON'],x.header['MODULE'],x.header['OUTPUT'],x.header['CHANNEL']) for x in h if x.name.startswith('LIGHTCURVE_Q')])
