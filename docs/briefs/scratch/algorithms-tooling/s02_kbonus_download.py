import time, numpy as np, lightkurve as lk, astropy.units as u, warnings
warnings.filterwarnings('ignore')
from astropy.coordinates import SkyCoord
from astropy.io import fits
t0=time.time()
c=SkyCoord(285.907827,40.941988,unit='deg')
sr=lk.search_lightcurve(c, author='KBONUS-BKG', radius=20*u.arcsec)
sr=sr[sr.target_name=='Gaia DR3 2103792461826460288']
lcc=sr.download_all(download_dir='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/data')
print('download took %.1fs'%(time.time()-t0), type(lcc), len(lcc))
lc=lcc[0]
print('type', type(lc), 'n', len(lc), 'columns', lc.columns)
print('meta keys sample', [k for k in list(lc.meta.keys())[:60]])
print('FLUX_ORIGIN', lc.meta.get('FLUX_ORIGIN'), 'time format', lc.time.format, lc.time.scale, lc.time[0], lc.time[-1])
print('QUARTERs present', np.unique(lc['quarter'].value) if 'quarter' in lc.columns else 'no quarter col')
import glob
f=glob.glob('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/data/**/*.fits', recursive=True)[0]
print(f)
with fits.open(f) as h:
    h.info()
    for i,hdu in enumerate(h):
        print('--- HDU',i, hdu.name)
        if hasattr(hdu,'columns') and hdu.columns is not None:
            print(hdu.columns.names)
    print({k:h[0].header[k] for k in h[0].header if k in ('OBJECT','KEPLERID','RA_OBJ','DEC_OBJ','GMAG','RUWE','TPFORG','PSFFRAC','CROWDSAP','FLFRCSAP','KEPMAG','TEFF','RADIUS','DOI','HLSPLEAD','HLSPVER','DATE')})
    print(repr(h[1].header)[:3000])
