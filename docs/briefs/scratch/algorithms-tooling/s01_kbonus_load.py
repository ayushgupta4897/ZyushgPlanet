import time, numpy as np, lightkurve as lk, astropy.units as u
from astropy.coordinates import SkyCoord
t0=time.time()
c=SkyCoord(285.907827,40.941988,unit='deg')
sr=lk.search_lightcurve(c, author='KBONUS-BKG', radius=20*u.arcsec)
print(sr)
print('n products', len(sr), 'search took %.1fs'%(time.time()-t0))
print(sr.table.colnames)
import pandas as pd
df=sr.table.to_pandas()[['target_name','mission','author','exptime','distance','dataURL' if 'dataURL' in sr.table.colnames else 'productFilename']]
print(df.to_string())
