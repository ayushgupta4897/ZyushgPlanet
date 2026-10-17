import pandas as pd, numpy as np, time, warnings; warnings.filterwarnings("ignore")
from astroquery.gaia import Gaia
from astropy.table import Table
m=pd.read_parquet("mstars_joined.parquet")
ids=m.gaia_designation.str.split().str[-1].astype(np.int64)
up=Table({'source_id':ids.values})
t0=time.time()
q="""SELECT u.source_id, g.parallax, g.parallax_error, g.ruwe, g.phot_g_mean_mag, g.bp_rp, g.pmra, g.pmdec,
            g.non_single_star, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.phot_variable_flag, g.teff_gspphot, g.radius_gspphot=0 as dummy
     FROM tap_upload.up AS u JOIN gaiadr3.gaia_source AS g ON g.source_id=u.source_id"""
q=q.replace(", g.radius_gspphot=0 as dummy","")
job=Gaia.launch_job_async(q, upload_resource=up, upload_table_name="up", dump_to_file=False)
r=job.get_results().to_pandas()
print(f"Gaia async upload-join returned {len(r)} rows in {time.time()-t0:.0f}s")
r.to_parquet("mstars_gaia.parquet")
print(r.describe().T.to_string())
r['dist_pc']=1000/r.parallax
print("\nparallax>0 & parallax_over_error>5:", ((r.parallax>0)&(r.parallax/r.parallax_error>5)).sum())
for d in (50,100,150,200,300,500):
    print(f"  within {d} pc (plx>1000/{d} mas, plx/err>5): {((r.parallax>1000/d)&(r.parallax/r.parallax_error>5)).sum()}")
print("ruwe<1.4:", (r.ruwe<1.4).sum(), " non_single_star>0:", (r.non_single_star>0).sum())
