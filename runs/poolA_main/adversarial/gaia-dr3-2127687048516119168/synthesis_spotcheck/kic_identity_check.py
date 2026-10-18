# Synthesis-stage check: is Gaia DR3 2127687048516119168 itself a Kepler target? (data/known manifest 2026-10-17T15:46Z)
import pandas as pd, numpy as np
ra0,dec0=289.5183982339101,46.01004525346003
ks=pd.read_csv('data/known/keplerstellar.csv')
d=np.hypot((ks.ra-ra0)*np.cos(np.radians(dec0)),ks.dec-dec0)*3600; m=d<20
print('keplerstellar within 20 arcsec:'); print(ks[m].assign(sep_arcsec=d[m]).sort_values('sep_arcsec').to_string())
tk=pd.read_csv('data/known/kepler_target_kepids.csv')
print('\nin kepler_target_kepids (observed >=1 quarter):', sorted(set(ks[m].kepid)&set(tk.iloc[:,0])))
for name in ['koi_cumulative.csv','dr25_tce.csv']:
    c=pd.read_csv('data/known/'+name); dd=np.hypot((c.ra-ra0)*np.cos(np.radians(dec0)),c.dec-dec0)*3600
    print(f'{name}: rows within 180 arcsec: {(dd<180).sum()}')
