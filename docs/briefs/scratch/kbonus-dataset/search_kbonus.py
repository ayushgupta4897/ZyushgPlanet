import lightkurve as lk, astropy.units as u
from astropy.coordinates import SkyCoord
import warnings; warnings.filterwarnings("ignore")
sr = lk.search_lightcurve(SkyCoord(285.907827,40.941988,unit='deg'), author='KBONUS-BKG', radius=20*u.arcsec)
print(sr)
t = sr.table
cols = [c for c in ['target_name','mission','exptime','distance','dataURL','productFilename','obs_id','project','provenance_name','t_min','t_exptime','sequence_number','dataproduct_type','description','obs_collection'] if c in t.colnames]
import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',120)
print(t[cols].to_pandas().to_string())
