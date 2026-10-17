import warnings; warnings.filterwarnings('ignore')
from astroquery.mast import Observations
import numpy as np
from astropy.time import Time
for s in [105,106,107,108,109]:
    try:
        o=Observations.query_criteria(obs_collection='TESS',sequence_number=s,provenance_name='SPOC',dataproduct_type='timeseries',t_exptime=[119,121])
        n=len(o)
        if n:
            rel=Time(np.nanmin(o['t_obs_release']),format='mjd').iso; tmax=Time(np.nanmax(o['t_max']),format='mjd').iso
            print(s,'n',n,'t_max',tmax,'release',rel)
        else: print(s,'n 0')
    except Exception as e: print(s,'ERR',e)
