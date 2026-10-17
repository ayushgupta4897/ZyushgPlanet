import warnings; warnings.filterwarnings("ignore")
from astroquery.gaia import Gaia
j=Gaia.launch_job("SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, ruwe, phot_g_mean_mag, bp_rp, teff_gspphot, non_single_star, ipd_frac_multi_peak, radial_velocity FROM gaiadr3.gaia_source WHERE source_id IN (2103792461826460288, 2103792427466721792)")
r=j.get_results(); print(r)
for row in r:
    print(row['source_id'], "pmra", row['pmra'], "pmdec", row['pmdec'], "-> /1000:", row['pmra']/1000, row['pmdec']/1000, "dist pc", 1000/row['parallax'])
