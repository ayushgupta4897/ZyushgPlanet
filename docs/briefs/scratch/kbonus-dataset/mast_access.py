import warnings; warnings.filterwarnings("ignore")
from astroquery.mast import Observations
import time
t0=time.time()
n=Observations.query_criteria_count(provenance_name="KBONUS-BKG")
print("query_criteria_count(provenance_name='KBONUS-BKG') =", n, f"({time.time()-t0:.1f}s)")
obs=Observations.query_criteria(provenance_name="KBONUS-BKG", target_name="Gaia DR3 2103792461826460288")
print(obs['obs_collection','provenance_name','project','target_name','obs_id','dataURL','s_ra','s_dec','t_min','t_max','filters','dataproduct_type'])
prods=Observations.get_product_list(obs)
print(prods['productFilename','dataURI','size','productType','description'])
# cone query example
obs2=Observations.query_criteria(provenance_name="KBONUS-BKG", coordinates="285.907827 40.941988", radius="0.01 deg")
print("cone 36 arcsec:", len(obs2), list(obs2['target_name']))
