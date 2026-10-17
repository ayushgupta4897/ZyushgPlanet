from astroquery.gaia import Gaia
import time, sys
Gaia.ROW_LIMIT=-1
t0=time.time()
# GCNS: Gaia Catalogue of Nearby Stars (Smart+2021), table external.gaiaedr3_gcns_main_1
sql = """SELECT 
 SUM(CASE WHEN phot_g_mean_mag<13 THEN 1 ELSE 0 END) AS n_g13,
 SUM(CASE WHEN phot_g_mean_mag<14 THEN 1 ELSE 0 END) AS n_g14,
 SUM(CASE WHEN phot_g_mean_mag<15 THEN 1 ELSE 0 END) AS n_g15,
 SUM(CASE WHEN phot_g_mean_mag<16 THEN 1 ELSE 0 END) AS n_g16,
 COUNT(*) AS n_all
FROM external.gaiaedr3_gcns_main_1
WHERE dist_50 < 0.1 AND phot_bp_mean_mag - phot_rp_mean_mag > 1.84
 AND phot_g_mean_mag + 5*LOG10(1000*dist_50) - 10 > 8.0"""
# bp_rp>1.84 ~ M0 boundary (Pecaut & Mamajek: M0V BP-RP ~1.84); abs G > 8 removes giants
try:
    j = Gaia.launch_job(sql); r=j.get_results(); print(r)
except Exception as e:
    print("GCNS query failed:", e)
sql2 = """SELECT 
 SUM(CASE WHEN dist_50<0.03 THEN 1 ELSE 0 END) AS n_30pc,
 SUM(CASE WHEN dist_50<0.05 THEN 1 ELSE 0 END) AS n_50pc,
 COUNT(*) AS n_100pc
FROM external.gaiaedr3_gcns_main_1
WHERE dist_50 < 0.1 AND phot_bp_mean_mag - phot_rp_mean_mag > 1.84
 AND phot_g_mean_mag + 5*LOG10(1000*dist_50) - 10 > 8.0"""
try:
    j = Gaia.launch_job(sql2); r=j.get_results(); print(r)
except Exception as e:
    print("GCNS query2 failed:", e)
print("time", time.time()-t0)
