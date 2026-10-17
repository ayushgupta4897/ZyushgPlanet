from astroquery.gaia import Gaia
import time
base = """SELECT COUNT(*) AS n FROM external.gaiaedr3_gcns_main_1
WHERE dist_50 < {d} AND phot_bp_mean_mag - phot_rp_mean_mag > 1.84
 AND phot_g_mean_mag + 5*LOG10(1000*dist_50) - 10 > 8.0 AND phot_g_mean_mag < {g}"""
print("GCNS M-dwarf candidates (BP-RP>1.84, M_G>8): counts by distance and G limit")
print("dist_pc  G<13  G<14  G<15  G<16  G<17  all")
for d in [0.03, 0.05, 0.1]:
    row=[]
    for g in [13,14,15,16,17,25]:
        r=Gaia.launch_job(base.format(d=d,g=g)).get_results()
        row.append(int(r['n'][0]))
    print(int(d*1000), *row, flush=True)
# TESS mag proxy: for M dwarfs T ~ G - 0.7 (Stassun+2019 relation for BP-RP~2.5 gives T-G ~ -0.75). Also count by G<15.7 (T<15), G<14.7 (T<14), G<13.7 (T<13), G<12.7 (T<12)
print("dist_pc  G<12.7(T<~12)  G<13.7(T<~13)  G<14.7(T<~14)  G<15.7(T<~15)")
for d in [0.03,0.05,0.1]:
    row=[int(Gaia.launch_job(base.format(d=d,g=g)).get_results()['n'][0]) for g in [12.7,13.7,14.7,15.7]]
    print(int(d*1000), *row, flush=True)
