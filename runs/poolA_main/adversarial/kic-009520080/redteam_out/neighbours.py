import numpy as np, os, subprocess, json
from astropy.io import fits
from astropy.table import Table
cat=Table.read('/home/user/ZyushgPlanet/data/kbonus-catalogs/hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits')
print(cat.colnames)
ids=['2127366575242177536','2127366575242178560','2127366575242178816','2127366472162960768']
gcol=[c for c in cat.colnames if 'gaia' in c.lower()][0]
rows={}
for i in ids:
    sel=np.array([str(x).strip().endswith(i) for x in cat[gcol]])
    if sel.sum()==0: print(i,'NOT in KBONUS catalogue'); continue
    r=cat[sel][0]; rows[i]={c:(str(r[c]) if not isinstance(r[c],(np.floating,float,np.integer,int)) else float(r[c])) for c in cat.colnames}; print(i,{c:rows[i][c] for c in cat.colnames if any(k in c.lower() for k in ('kic','mag','psffrac','crowd','flfr','npix','tpf','name','id','ra','dec'))})
json.dump(rows,open('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic009520080/out/neighbour_catalog_rows.json','w'),indent=1,default=str)
