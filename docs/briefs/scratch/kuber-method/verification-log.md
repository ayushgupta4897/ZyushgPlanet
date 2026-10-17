# kuber-method: sandbox verification log (2026-10-17)

Everything below was run in /home/user/ZyushgPlanet/.venv (Python 3.13, astropy 8.0.2, lightkurve 2.6.0) against live services. Downloaded data live in ./data/ (untrusted inputs; never executed).

## Kepler Bonus catalogue downloads and hashes (match evidence/codex/database_expansion/kbonus/*provenance.json and analysis/editorial-audit-20261013/verification.json)

- https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits
  2,390,400 bytes, Last-Modified 2023-08-04, sha256 406563951c367bb1ffe3c20f83c4684a4083a2fe2f3d7c36ad9a663b5f6977ef
  29,800 rows; columns: gaia_designation, ra, dec, kic, fname
- https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits
  159,022,080 bytes, sha256 2e6d5a87cb90f941bdc207aef93d61847df228ae7d56ffd38d3ab9255dd30fcf
  606,900 rows; columns: gaia_designation, ra, dec, sap_mean_flux, sap_mean_flux_err, psf_mean_flux, psf_mean_flux_err,
  flfrcsap, crowdsap, npixsap, psffrac, pertrati, pertstd, psf_avail, sap_avail, psf_fraction_flag, phot_g_mean_mag,
  phot_bp_mean_mag, phot_rp_mean_mag, tpforg, kic, kic_sep, kep_mag, fname

## Selection counts reproduced independently (2.8 s)

join M-star ids -> 29,800 matched; psffrac>0.5 -> 9,204; crowdsap>0.7 -> 6,477; psf_avail has >=8 '1' characters -> 2,495; kic != tpforg -> 1,836.
Candidate Gaia DR3 2103792461826460288 is in the eligible set with G=16.2824, rank 51 (0-based) by ascending G among the 1,836 (before primary-target exclusion).
psf_avail for the candidate: 011110111011101110 (index = quarter 0..17; quarters 1-4,6-8,10-12,14-16 available = 13 quarters).

## Candidate light-curve product (MAST)

lightkurve: lk.search_lightcurve(SkyCoord(285.907827,40.941988,unit='deg'), author='KBONUS-BKG', radius=20*u.arcsec) -> 2 products
  target_name 'Gaia DR3 2103792461826460288', obs_id hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0,
  dataURL mast:HLSP/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits, t_exptime 1765.464 s
  neighbour 'Gaia DR3 2103792427466721792' (KIC 5689116) at 3.1" from the search position (7.1" from the target).
Direct download: https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits
  8,847,360 bytes, sha256 44980afb7843a4507c898eb1ce5d5733029a0f4658f6b96d6196db49df37c06f (matches Kuber's recorded fits_sha256)
FITS structure: 28 HDUs: PRIMARY, LIGHTCURVE_STITCHED, then LIGHTCURVE_Q{1,2,3,4,6,7,8,10,11,12,14,15,16} each followed by APERTURE_Q{n}.
Primary header: OBJECT 'KIC 5689121', KEPLERID 5689121, GAIAID 'Gaia DR3 2103792461826460288', RA_OBJ 285.9078361543101, DEC_OBJ 40.94198476904723,
  XPOSURE 1765.464, TIMESYS TDB, TPFORG 5689116, KEPMAG 16.332, GMAG 16.2824, DOI 10.17909/7jbr-w430, HLSPVER V1.0.
Per-quarter columns: CADENCENO, TIME (BKJD = BJD-2454833), FLUX, FLUX_ERR (PSF), SAP_FLUX, SAP_FLUX_ERR, CENTROID_COLUMN, CENTROID_ROW,
  SAP_QUALITY, SAP_BKG, RED_CHI2, PSF_FLUX_NOVA, PSF_FLUX_ERR_NOVA.
Q1 header: QUARTER 1, PSFFRAC 0.8149, CROWDSAP 0.8131, FLFRCSAP 0.9129, NPIXSAP 6, CHANNEL 6, MODULE 3, OUTPUT 2, SEASON 3; 1,624 cadences, TIME 131.51-164.98 BKJD,
  SAP_QUALITY==0 fraction 0.893; median FLUX 4126 e-/s, median SAP_FLUX 5070 e-/s, median FLUX_ERR 3.29, median SAP_FLUX_ERR 3.38.

## Other live checks

- CDS Barber et al. 2022 table: https://cdsarc.cds.unistra.fr/ftp/J/AJ/164/88/table5.dat 80,264 bytes, Last-Modified 2022-11-08 (HTTP 200).
- PyPI current versions (2026-10-17): triceratops 1.1.0 (2026-09-02), leo-vetter 1.2.0 (2026-07-29), wotan 2.0 (2026-10-11), transitleastsquares 2.0.1 (2026-10-11),
  lightkurve 2.6.0 (2026-04-16), batman-package 2.5.3 (2025-05-05), emcee 3.1.6 (2024-04-19), pytransit 2.9.3 (2026-10-15), numba 0.68.0, astropy 8.0.2 (2026-10-15).
- Upstream pins resolve: paraschopra/tess-planet-candidates 74d5624 (2026-10-09 14:15 +0530), eumemic/exoplanets 2cf0b11 (2026-10-09 23:41 -0700), cloned read-only under ../untrusted/upstream/.
- planethunt git history (after `git fetch --unshallow`): 60 commits, tags evidence-20261010 (2026-10-10 20:32 +0530), evidence-20261013 (01:06 +0530), evidence-20261013-r2 (01:17 +0530).
- x.com tweet pages return HTTP 402/451 from this sandbox; tweet text taken from the task statement.
