import warnings; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
ra,dec=285.1259461685635,49.01660406817707
job=Gaia.launch_job(f"""SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, ruwe, phot_g_mean_mag, bp_rp, phot_bp_rp_excess_factor,
 teff_gspphot, logg_gspphot, mh_gspphot, distance_gspphot, ag_gspphot, ebpminrp_gspphot, radial_velocity, non_single_star, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude,
 phot_variable_flag, DISTANCE(POINT({ra},{dec}),POINT(ra,dec))*3600 AS sep
 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT({ra},{dec}),CIRCLE(ra,dec,30./3600)) ORDER BY sep""")
t=job.get_results(); t.pprint_all()
t.write('gaia_30as.csv',overwrite=True)
# astrophysical params for target and bright neighbour
for sid in [2131988887763678720,2131988990842894976]:
    j=Gaia.launch_job(f"SELECT source_id, radius_gspphot, radius_gspphot_lower, radius_gspphot_upper, lum_flame, radius_flame, mass_flame, teff_gspphot, logg_gspphot, ag_gspphot, classprob_dsc_combmod_star, spectraltype_esphs FROM gaiadr3.astrophysical_parameters WHERE source_id={sid}")
    r=j.get_results(); print(sid); r.pprint_all()
s=Simbad(); s.add_votable_fields('otype','sp','flux(V)','ids')
try:
    r=s.query_region(SkyCoord(ra,dec,unit='deg'),radius=40*u.arcsec); print('SIMBAD 40as:'); print(r)
except Exception as e: print('SIMBAD err',e)
for name in ['KIC 11285493','KIC 11285496']:
    try: r=s.query_object(name); print(name, r)
    except Exception as e: print(name,'simbad err',e)
