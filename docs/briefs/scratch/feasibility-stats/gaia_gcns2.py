from astroquery.gaia import Gaia
import time, sys, numpy as np
Gaia.ROW_LIMIT=-1
t0=time.time()
sql = """SELECT source_id, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, dist_50, ra, dec
FROM external.gaiaedr3_gcns_main_1
WHERE dist_50 < 0.1 AND phot_bp_mean_mag - phot_rp_mean_mag > 1.84
 AND phot_g_mean_mag + 5*LOG10(1000*dist_50) - 10 > 8.0"""
try:
    j = Gaia.launch_job_async(sql); r=j.get_results()
    print("rows", len(r), "time", time.time()-t0)
    r.write('dl/gcns_mdwarfs_100pc.csv', format='csv', overwrite=True)
except Exception as e:
    print("GCNS query failed:", repr(e)[:2000])
