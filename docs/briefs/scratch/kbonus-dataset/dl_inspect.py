import lightkurve as lk, astropy.units as u, warnings, os, sys
warnings.filterwarnings("ignore")
from astropy.coordinates import SkyCoord
from astropy.io import fits
dl = "/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads"
sr = lk.search_lightcurve(SkyCoord(285.907827,40.941988,unit='deg'), author='KBONUS-BKG', radius=20*u.arcsec)
paths=[]
for i in range(len(sr)):
    lc = sr[i].download(download_dir=dl)
    print("downloaded", lc.meta.get('FILENAME'), type(lc).__name__)
    paths.append(lc.meta.get('FILENAME'))
# inspect target file
f = paths[0]
print("FILE", f, os.path.getsize(f)/1e6, "MB")
with fits.open(f) as h:
    h.info()
    print("\n=== PRIMARY HEADER ===")
    for k,v in h[0].header.items():
        if k not in ('COMMENT','HISTORY'):
            print(f"{k:10s}= {v!r}  / {h[0].header.comments[k]}")
    print("\n=== HDU1 header (LIGHTCURVE_STITCHED) ===")
    for k,v in h[1].header.items():
        if k not in ('COMMENT','HISTORY'):
            print(f"{k:10s}= {v!r}  / {h[1].header.comments[k]}")
    print("\n=== HDU1 columns ===")
    print(h[1].columns)
    d=h[1].data
    import numpy as np
    print("N cadences", len(d), "time range", d['TIME'].min(), d['TIME'].max(), "quarters", sorted(set(d['QUARTER'])))
    # find first per-quarter LC HDU and aperture HDU
    names=[x.name for x in h]
    print("HDU names", names)
    for n in names:
        if n.startswith('LIGHTCURVE_Q'):
            k=names.index(n)
            print(f"\n=== {n} header ===")
            for kk,v in h[k].header.items():
                if kk not in ('COMMENT','HISTORY') and not kk.startswith(('TTYPE','TFORM','TUNIT','TDISP','TDIM')):
                    print(f"{kk:10s}= {v!r}  / {h[k].header.comments[kk]}")
            print(h[k].columns)
            break
    for n in names:
        if n.startswith('APERTURE_Q'):
            k=names.index(n)
            print(f"\n=== {n} header ===")
            for kk,v in h[k].header.items():
                if kk not in ('COMMENT','HISTORY'):
                    print(f"{kk:10s}= {v!r}  / {h[k].header.comments[kk]}")
            print("aperture shape", h[k].data.shape, "sum", h[k].data.sum(), "unique", np.unique(h[k].data))
            print(h[k].data)
            break
