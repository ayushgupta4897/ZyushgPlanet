# Kepler Bonus background light curves (KBONUS-BKG) as a planet-hunting dataset — practitioner brief

Date of checks: 2026-10-17. Everything marked **[sandbox]** was run in `/home/user/ZyushgPlanet/.venv` (lightkurve 2.6.0, astropy 8.0.2, astroquery 0.4.11); scripts and outputs are under `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/`. Everything marked **[Kuber]** comes from reading (never executing) the read-only clone of `Kuberwastaken/planethunt`. Inferences are labelled as such.

## 0. Executive summary

* KBONUS-BKG is a single v1.0 release (MAST HLSP, released 2023-11-01, DOI 10.17909/7jbr-w430, CC BY 4.0) of **606,900 Kepler long-cadence light curves** extracted with PSF ("Linearized Field Deblending", psfmachine) photometry from **all 17 Kepler quarters' TPFs** using Gaia DR3 (G<19) as the input catalogue: 406,548 are *background* sources the Kepler pipeline never processed, 200,352 are the Kepler targets re-extracted (paper numbers). No v1.1/v2 or K2/TESS follow-up data release exists as of today ([MAST page](https://archive.stsci.edu/hlsp/kbonus-bkg), [paper arXiv:2310.17733](https://arxiv.org/abs/2310.17733), v1 only, AJ 166, 265).
* Precision is **comparable to, but slightly worse than, the Kepler pipeline for the same star**: for KIC 5689116 (G=15.48, a real Kepler target) I measure 6-h CDPP 184 ppm (KBONUS PSF, stitched; 133 ppm on the pre-flattened column) vs 121 ppm for official PDCSAP **[sandbox]**. For Kuber's target (G=16.28) the PSF 6-h CDPP is ~475 ppm (per quarter 337–912 ppm). A 28-star clean background sample gives median per-quarter 6-h CDPP (PSF flux) of ~90/75/80 ppm at G=12–13/13–14/14–15, ~300 ppm at G=15–16, ~330 ppm at 16–17, ~690 ppm at 17–18 and ~1600 ppm at 18–19 **[sandbox]** (small sample; the paper says "<100 ppm for G<16", "10s ppm" for G<14, "100s ppm" for G=16–18).
* The **"M-star catalogue" (29,800 entries) is faint and distant**: median G 18.47, Gaia teff_gspphot 2616–3900 K, only 182 within 100 pc (and of those only 6 are non-Kepler-targets with PSFFRAC>0.5 and ≥8 quarters) **[sandbox]**. Kuber's 1,800 eligible stars span G=14.19–19.0 (rank 64 ends at G=16.88, rank 576 at G=18.21); only 100 of the 1,800 are within 300 pc. His search covered 1,736 of them. The remaining 28,000 M-star entries fail PSFFRAC>0.5 (20,596), CROWDSAP>0.7 (2,727), ≥8 quarters (3,982) or are Kepler targets (695).
* **The bigger, better hunting ground is not the M-star catalogue.** Among the 412,033 non-Kepler-target sources, 23,171 pass PSFFRAC>0.5 & CROWDSAP>0.7 & ≥8 PSF quarters; a Gaia DR3 cross-match shows 8,463 of them are dwarfs (M_G>4, parallax/error>5), of which ~3,600 have G<17 (noise ≤ ~350 ppm per 6 h), and 140 are within 300 pc **[sandbox]**. Nobody has published a systematic transit search of the KBONUS-BKG light curves themselves (checked: the KBonus paper, its ~13 citing papers via Crossref/OpenAlex counts, Bienias & Szabó 2025 who searched raw pixels, Barber+2022 who searched only Kepler-target PDCSAP curves, Claytor & Tayar 2025 who did rotation).
* Hard-won dataset pitfalls (all verified in the sandbox): the catalogue column names differ from the documentation (`tpforg`, `kep_mag`, `fname`); the stitched `FLUX` is level-matched across quarters while per-quarter extension `FLUX` is raw e-/s (±9% jumps); `PSF_FLAT_FLUX` is a 2-day b-spline flatten that was *not* transit-masked (it shaves Kuber's 2.2 ppt dip to ~1.7 ppt); `FLATTEN_MASK` is an outlier mask, not a transit mask; `SAP_QUALITY` is the raw Kepler TPF QUALITY word and lightkurve's *default* bitmask drops none of the bits present (use `quality==0` or `'hardest'`); `CENTROID_COLUMN/ROW` are the Gaia position plus a scene-wide jitter vector identical for every star in the TPF (useless for localisation, as Kuber found); the primary-header `PMRA/PMDEC` are labelled mas/yr but are Gaia values ÷1000 (arcsec/yr); lightkurve 2.6.0 ignores `flux_column=` for these files and rejects `quarter=`; the MAST page says "~20,000" M dwarfs while the file and the paper's Appendix F say 29,800.

---

## 1. The paper and the product

### 1.1 Reference chain
* Martínez-Palomera, Hedges & Dotson 2023, "Kepler Bonus: Light Curves of Kepler Background Sources", AJ 166, 265, DOI [10.3847/1538-3881/ad0727](https://doi.org/10.3847/1538-3881/ad0727), [arXiv:2310.17733](https://arxiv.org/abs/2310.17733) (v1 2023-10-26; no later arXiv version). Published 2023-11-24 (Crossref). Crossref `is-referenced-by-count` = 13, OpenAlex `cited_by_count` = 12 (both checked 2026-10-17).
* Method paper: Hedges et al. 2021, "Linearized Field Deblending", AJ 162, 107 (LFD). Software: [psfmachine](https://github.com/SSDataLab/psfmachine) (PyPI latest 1.1.4, 2023-12-15; the FITS headers record `ORIGIN='PSFMACHINE'`, `VERSION='1.1.2'`, `PROCVER='2.0.9'`, `DATE='2022-07-25'`, checksums updated 2023-03-31 **[sandbox]**). Configuration/notebooks: [jorgemarpa/KBonus](https://github.com/jorgemarpa/KBonus) (MIT; notebooks `tpf_stack_run.ipynb`, `light_curve_example.ipynb`, `kepler_686_example.ipynb`, `light_curve_caveats.ipynb`).
* Sister product (different pixels): KBONUS-APEXBA — aperture photometry of 9,327 sources in the EXBA "super-aperture" background masks, Martínez-Palomera et al. 2022, AJ 163, 93, [arXiv:2112.06683](https://arxiv.org/abs/2112.06683), HLSP [kbonus-apexba](https://archive.stsci.edu/hlsp/kbonus-apexba), DOI 10.17909/t9-d5wy-e535 (69 EBs, one planet candidate around Gaia EDR3 2077240046296834304, 0.81 R_J). Not the same data as KBONUS-BKG.

### 1.2 What was done (from the paper text, [arXiv HTML](https://arxiv.org/html/2310.17733v1))
* Input: 204,933 long-cadence TPFs from MAST, quarters Q0–Q17 (short cadence deferred). Processing in batches of ~200 spatially sorted TPFs per quarter/channel (~5,000 pixel time series, ~400 sources per batch).
* Source list: Gaia DR3 to G=19 within each TPF cutout + 16″ (1.4 M sources queried, proper motions propagated); keep sources within 4″ of the TPF edge; remove G<10 (saturation); of pairs within 1″ keep only the brighter. Result 606,900 sources.
* PSF model: built from Kepler FFIs (Gaia EDR3 to G=20, ~12,000 sources/100,000 pixels per CCD), evaluated on a grid 10× finer than a pixel; a time-dependent "perturbed" PSF is built from the first four CBVs smoothed with a 2-day b-spline on 90 time bins (replacing an earlier cubic time polynomial that created spurious trends). Background: cubic B-splines in time (knots every 2 h) and pixel row (knots every 6 px), fit to background pixels (`kbackground` package), motivated by rolling-band noise.
* Flux priors: Gaussian with mean F_G (from Gaia G) and σ = 10√F_G; sources that fit negative get their prior width halved (and neighbours within 5″), iterated 3 times; remaining negatives dropped (<5%). Saturated pixels (>1.2×10⁵ e⁻/s) plus up to 3 neighbours masked; pixels within 800″ of G≤8 stars masked.
* Aperture (SAP) photometry is derived from the evaluated PSF ("optimal" mode; `SAP='optimal'` in the quarter headers).
* Demonstrations: (i) deblending KOI-770.01 (P=1.506 d, 2211 ppm, centroid-offset FP) — the KOI's own PSF curve shows no transit, the eclipse is on Gaia DR3 2134870879540928896, 1.36″ away and 2.6 mag fainter; KOI-909.01 (P=16.37 d) is contaminated by the EB KIC 8256044, 8″ away. (ii) PSF photometry corrects transit dilution better than aperture photometry (depth ratio median 1.003; small planet-radius shifts). The authors explicitly say the flattened stitched light curves are meant "to help the community run transit searches"; **they did not run one**.

### 1.3 Extraction metrics (header keywords / catalogue columns)
| Keyword | Meaning (paper §/App. D–E) |
|---|---|
| `PSFFRAC` | fraction of the source's PSF that falls on the TPF pixels used; ~1 = fully enclosed. **Only quarters with PSFFRAC>0.5 are stitched.** Catalogue `psffrac` = minimum over detected quarters. |
| `FLFRCSAP` | fraction of the source flux inside the SAP aperture (from the PSF model). |
| `CROWDSAP` | source flux / total flux in the SAP aperture (1 = no contamination). |
| `NPIXSAP` | pixels in the SAP aperture (verified equal to `APERTURE_Qn.sum()`). |
| `PERTRATI` | ratio of mean flux from mean-PSF vs perturbed-PSF model; far from 1 = bad perturbation fit. |
| `PERTSTD` | ratio of light-curve std from mean vs perturbed model; small = stable. |
| `ITERNEG` | 1 if the object went negative during the iteration (per-quarter header). |
| `QDETECT` (primary header) / `psf_fraction_flag` (catalogue) | 18-char string Q0..Q17: `0` no PSF detection, `1` PSFFRAC>0.5, `2` PSFFRAC<0.5 (partial, not stitched) — definition from the KBonus caveats notebook. `psf_avail` marks any detection (`1` even when the flag is `2`), so "number of quarters" from `psf_avail` over-counts usable quarters. |

### 1.4 Known systematics (paper §IV, caveats notebook, Claytor & Tayar 2025)
* Partial-PSF sources (edge of TPF, saturation masking): noisier PSF flux, underestimated FLFRCSAP/CROWDSAP; **quarter-to-quarter flux offsets come from PSFFRAC/aperture changes** as the star lands on a different channel each season.
* PSF model independent of brightness; distorted PSFs on edge channels; high-contrast blends (≥2 mag within 15″) can go negative or be over-corrected; velocity-aberration residuals for blended/partial sources; ~1% of faint neighbours removed because their light curves correlated (r>0.5, <25″) with a brighter source; red giants near the HB (~1.6%) have 10× worse CDPP; LPV amplitudes may be biased; "empty PSF columns" exist for sources that only got SAP photometry.
* Claytor & Tayar 2025 ([arXiv:2506.03248](https://arxiv.org/abs/2506.03248), ApJ accepted; poster [Zenodo 17245793](https://zenodo.org/records/17245793)), using 449,913 KBonus light curves for rotation: **up to 63% of periodic background light curves still correlate with their foreground star** (zero-lag cross-correlation vs 3,000 random pairs, p<0.15; 55% of 3,224 pairs have periods within 10%); slight over-representation of anticorrelated pairs (over-correction); a cluster of spurious ~40-d periods possibly tied to downlink systematics. Their recommendation: do not treat a background signal as independent until you have checked it against the foreground star's light curve. For transits this is exactly the "neighbour leakage" problem (§5).

### 1.5 Measured precision **[sandbox]**
Script: `cdpp_compare.py`, `sample_cdpp.py` (6-h CDPP = std of 13-cadence bins after a 101-cadence Savitzky–Golay flatten and 5σ clip; lightkurve's own `estimate_cdpp` returned 0 in this venv, so the manual estimator was used throughout).

| Star | Column | 6-h CDPP stitched | per-quarter range |
|---|---|---|---|
| KIC 5689121 (G 16.28, Kuber's target) | PSF `FLUX` | 475 ppm | 337–912 |
| | `SAP_FLUX` (diluted) | 384 ppm | 227–589 |
| | `PSF_FLAT_FLUX` | 402 ppm | – |
| KIC 5689116 (G 15.48, the TPF host) | KBONUS PSF | 184 ppm | 103–213 |
| | KBONUS `PSF_FLAT_FLUX` | 133 ppm | – |
| | official PDCSAP (lightkurve `author='Kepler'`) | 121 ppm | 96–151 (pipeline `CDPP6_0` header: 129–177) |

Sample of 28 non-target sources with PSFFRAC>0.7, CROWDSAP>0.8, ≥12 quarters (median of per-quarter 6-h CDPP, PSF flux / SAP flux): G12–13: 92/79; 13–14: 73/54; 14–15: 82/74; 15–16: 298/161; 16–17: 327/377; 17–18: 686/599; 18–19: 1598/758 ppm. Several reddened giants (BP−RP≈2) at G≈13 reach 250–600 ppm, so colour/parallax selection matters as much as magnitude. Formal `FLUX_ERR` is typically 0.6–0.8× the measured point-to-point scatter at G≈16 and larger than it at G≈18, i.e. do not trust the error column for SNR without re-estimating noise (Paras Chopra's GUIDE.md lists "underestimated noise" as a Codex failure mode — same lesson).

---

## 2. Access

### 2.1 Catalogue files (direct URLs, sizes, SHA-256 verified **[sandbox]**, identical to Kuber's `verification.json`)
```
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits   159,022,080 B  sha256 2e6d5a87cb90f941bdc207aef93d61847df228ae7d56ffd38d3ab9255dd30fcf  (606,900 rows)
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits      2,390,400 B  sha256 406563951c367bb1ffe3c20f83c4684a4083a2fe2f3d7c36ad9a663b5f6977ef  (29,800 rows)
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_kois-nns-catalog_kepler_v1.0_cat.fits    3,430,080 B  sha256 411b3ce12c6e64aa29bbd79750aacb53623620aac0b0329a1361b55f814b8264  (35,639 rows)
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_wd-catalog_kepler_v1.0_cat.fits             14,400 B  sha256 13071542df2172dd12084e4c5ca3e40e76301728ea597a801b91df6ff96b1938  (91 rows)
https://archive.stsci.edu/hlsps/kbonus-bkg/hlsp_kbonus-bkg_kepler_kepler_all_kepler_v1.0_readme.md
```
Actual source-catalogue columns (the README documents `tpf_org`, `kepmag`, `file_name`; the file has **`tpforg`, `kep_mag`, `fname`**): `gaia_designation` (string "Gaia DR3 …"), `ra`, `dec`, `sap_mean_flux[_err]`, `psf_mean_flux[_err]` (NaN when PSF photometry failed), `flfrcsap`, `crowdsap`, `npixsap`, `psffrac`, `pertrati`, `pertstd`, `psf_avail`, `sap_avail`, `psf_fraction_flag` (18-char strings), `phot_g_mean_mag`, `phot_bp_mean_mag`, `phot_rp_mean_mag` (BP/RP can be ±1e20 sentinels), `tpforg` (KIC of the TPF the source was extracted from), `kic` (0 if no KIC match; 65,501 rows), `kic_sep` (arcsec), `kep_mag`, `fname` (`kic-005689121` or `gaia-dr3-2049134501869730432`). The M-star and WD catalogues have only `gaia_designation, ra, dec, kic, fname` (+`Pwd`); the KOI-neighbour catalogue adds `nn_koi, nn_dist`.

Catalogue facts **[sandbox]**: G range 10.00–19.00; per-quarter PSF availability 212,536 (Q0) and 309–342k (Q1–Q17); number of PSF quarters per source is strongly bimodal (96,181 sources with 1 quarter; 98,401 with 17; 45,391 with 18); `psffrac>0.5` for 47%, `crowdsap>0.7` for 60%, both for 33%.

### 2.2 Light-curve files
* Naming: `hlsp_kbonus-bkg_kepler_kepler_<kic-NNNNNNNNN | gaia-dr3-NNNN…>_kepler_v1.0_lc.fits`. Path: `lcs/<first 4 digits of the 9-digit KIC or of the Gaia ID>/<full id>/<file>`. Two working URL forms (both HTTP 200, 8,847,360 B for the target):
  * `https://archive.stsci.edu/hlsps/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits`
  * `https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits`
  * Gaia-only example from the caveats notebook: `.../lcs/2051/2051766354744241536/hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2051766354744241536_kepler_v1.0_lc.fits`.
* Volume: 4.8 TB for all light curves (8–12 MB each); bundles `bundles/…kois-nns…lcf.tar.gz` 204 GB, `…mstars…` 134 GB, `…wd…` 373 MB. Bulk scripts at the HLSP root: `hlsp_kbonus-bkg_kepler_lightcurves_all_kepler_v1_bulkdl.sh` (29 MB of `curl -f --create-dirs` lines using the MAST API URL form), `…_chunkdl.zip` (50,000 files per chunk), `…_lftp.sh` (mirrors `/hlsps/kbonus-bkg/lcs`). MAST Portal caps (50,000 load / 500,000 download) are below the total.
* astroquery **[sandbox]**: `Observations.query_criteria_count(provenance_name="KBONUS-BKG")` → 606,900 (8 s). `Observations.query_criteria(provenance_name="KBONUS-BKG", target_name="Gaia DR3 2103792461826460288")` → one observation (`obs_collection='HLSP'`, `project='Kepler'`, `dataproduct_type='timeseries'`); `get_product_list` → one `_lc.fits`; cone `coordinates="285.907827 40.941988", radius="0.01 deg"` → 2 sources. lightkurve: `lk.search_lightcurve(SkyCoord(...), author='KBONUS-BKG', radius=20*u.arcsec)` returns `target_name` = Gaia designation, `exptime` 1765.464 s, one row per source (not per quarter). Kuber's pipeline used the raw `Mast.Caom.Filtered` service with filters `provenance_name=KBONUS-BKG`, `target_name=[Gaia designations]`, pagesize 50,000 **[Kuber pilot.py]**.

### 2.3 FITS layout (inspected **[sandbox]**, `dl_inspect.py`)
```
No.  Name                 Type         Dims / rows
0    PRIMARY              PrimaryHDU   48 cards: DOI, HLSPID, LICENSE, ORIGIN=PSFMACHINE, VERSION=1.1.2, PROCVER=2.0.9,
                                       TIMESYS=TDB, OBJECT='KIC 5689121', KEPLERID, RA_OBJ/DEC_OBJ (EQUINOX 2016), KEPMAG,
                                       TPFORG=5689116, GAIAID, PMRA, PMDEC (see bug below), PARALLAX [mas], GMAG, BPMAG, RPMAG,
                                       QDETECT='011110111011101110', XPOSURE=1765.464, MJD-BEG/END, APERTURE='PSF + SAP'
1    LIGHTCURVE_STITCHED  BinTable     50182 rows × 11: TIME [BJD-2454833], CADENCENO, QUARTER, FLUX, FLUX_ERR, SAP_FLUX,
                                       SAP_FLUX_ERR, PSF_FLAT_FLUX, PSF_FLAT_FLUX_ERR, SAP_QUALITY, FLATTEN_MASK   (e-/s)
2,4,… LIGHTCURVE_Qn       BinTable     13 cols: CADENCENO, TIME, FLUX (corrected PSF), FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR,
                                       CENTROID_COLUMN, CENTROID_ROW, SAP_QUALITY, SAP_BKG, RED_CHI2, PSF_FLUX_NOVA (mean-PSF), PSF_FLUX_ERR_NOVA
                                       header: SEASON, MODULE, OUTPUT, CHANNEL, QUARTER, TPF_ROW, TPF_COL, SAP='optimal', ROW, COLUMN,
                                       FLFRCSAP, CROWDSAP, NPIXSAP, PSFFRAC, PERTRATI, PERTSTD, ITERNEG
3,5,… APERTURE_Qn         ImageHDU     uint8 mask in the originating TPF's shape (e.g. (4,5) or (5,5)); sum == NPIXSAP
```
Only quarters with QDETECT=1 get `LIGHTCURVE_Q`/`APERTURE_Q` HDUs (13 for this star: Q1–4, 6–8, 10–12, 14–16).

**Time reference**: `TIME` is BKJD = BJD_TDB − 2,454,833 (lightkurve reads it as `format='bkjd', scale='tdb'`); first cadence 131.5128 BKJD = BJD 2,454,964.5128. BTJD (TESS) = BKJD − 2167 (2,457,000 − 2,454,833); e.g. Kuber's T0 = BJD 2,455,630.245431 = BKJD 797.245431 = BTJD −1369.7546. Mid-exposure, 29.424-min cadence.

**Quarter normalisation**: stitched `FLUX` per-quarter medians agree to ±0.3% (level-matched), while `LIGHTCURVE_Qn.FLUX` medians vary 3596–4309 e⁻/s (±9%) and `SAP_FLUX` 2708–7501 e⁻/s (aperture 2–7 px) across quarters. Stitched `SAP_FLUX` is also level-matched but its *contamination* changes by season, so a transit depth measured in SAP is season-dependent (Kuber's dilution test exploits this).

**Quality**: `SAP_QUALITY` is the Kepler TPF QUALITY word. For the target, 3,659/50,182 cadences are non-zero, with bits 16 (reaction-wheel zero crossing), 128 (cosmic ray in aperture), 1024 (sudden sensitivity dropout), 2048 (impulsive outlier), 8192 (cosmic ray in collateral), 131072/262144 (rolling band in aperture/full mask). lightkurve's default bitmask keeps all 50,182; `'hard'` keeps 47,971; `'hardest'` = `quality==0` keeps 46,523 (what Kuber used).

**`PSF_FLAT_FLUX` / `FLATTEN_MASK`**: the flatten is a 2-day-window b-spline applied to the stitched PSF flux with transits *not* masked. At Kuber's ephemeris the box depth is 1,696 ± 84 ppm in `PSF_FLAT_FLUX` vs ~2,165–2,348 ppm with local baselines on raw `FLUX` (Kuber) — i.e. the supplied flattening attenuates a 2.3-h, 2 ppt dip by ~20–25%. `FLATTEN_MASK==1` on 322 cadences in runs of 3–7 consecutive cadences, 21% coinciding with non-zero quality, only 2 inside transits: it is a sigma-clip/outlier mask from the flattening, not a transit mask. Recommendation: search on `FLUX` with your own detrending (wotan biweight ~0.5–1 d window, or per-quarter), and keep `PSF_FLAT_FLUX` only for quick looks.

**Header bug**: `PMRA=-0.0054`, `PMDEC=0.00285` labelled "mas/yr" are exactly Gaia DR3 pmra=−5.4019, pmdec=2.8479 mas/yr divided by 1000 (arcsec/yr). Use Gaia directly for proper motions.

### 2.4 lightkurve 2.6.0 behaviour **[sandbox]**
```python
import lightkurve as lk, numpy as np
from astropy.table import Table; from astropy.io import fits; from astropy.time import Time
lc = lk.read(path)                       # KeplerLightCurve, PSF 'flux'; columns include sap_flux, psf_flat_flux, flatten_mask, quarter
lc_sap  = lc.select_flux('sap_flux')     # WORKS.  lk.read(path, flux_column='sap_flux') silently returns PSF flux for this product
lc_flat = lc.select_flux('psf_flat_flux')
lc = lk.read(path, quality_bitmask='hardest')   # == quality==0 (46,523 of 50,182 cadences); default keeps everything
# lk.read(path, quarter=7) raises TypeError in 2.6.0 (the KBonus notebook shows it working in 2.5.0). Per-quarter HDUs:
t = Table.read(path, hdu='LIGHTCURVE_Q7'); h = fits.getheader(path, 'LIGHTCURVE_Q7')
lq = lk.LightCurve(time=Time(np.asarray(t['TIME']), format='bkjd', scale='tdb'), flux=np.asarray(t['FLUX']),
                   flux_err=np.asarray(t['FLUX_ERR']), quality=np.asarray(t['SAP_QUALITY']))
print(h['PSFFRAC'], h['CROWDSAP'], h['FLFRCSAP'], h['NPIXSAP'], h['CHANNEL'], h['SEASON'])
ap = fits.getdata(path, 'APERTURE_Q7')   # uint8 mask, TPF shape
```
`lc.meta['AUTHOR']` is reported as 'Kepler' (reader default), not 'KBONUS-BKG'; `lc.meta['TPFORG']` gives the originating TPF's KIC.

### 2.5 Working code: list all KBONUS M dwarfs with Gaia IDs, magnitudes, quarter availability, parallax
(`catstats.py`, `gaia_xmatch.py`, `hunting.py` in the scratch dir; the Gaia join of 29,800 IDs took 15 s via an anonymous async upload job.)
```python
import numpy as np, pandas as pd
from astropy.table import Table
from astroquery.gaia import Gaia
C = "catalogs/"
src = Table.read(C+"hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits").to_pandas()
ms  = Table.read(C+"hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits").to_pandas()
for df in (src, ms):
    for c in df.columns:
        if df[c].dtype == object and isinstance(df[c].iloc[0], bytes): df[c] = df[c].str.decode('utf-8')
m = src.merge(ms[['gaia_designation']], on='gaia_designation')            # 29,800 rows
m['nq_any']  = m.psf_avail.str.count('1')                                   # quarters with any PSF detection
m['nq_good'] = m.psf_fraction_flag.str.count('1')                           # quarters with PSFFRAC>0.5 (the stitched ones)
kep = set(pd.read_csv("keplerstellar_kepids.csv").kepid)                    # TAP: select distinct kepid from keplerstellar (200,038 ids)
m['is_kepler_target'] = m.kic.isin(kep)
ids = m.gaia_designation.str.split().str[-1].astype(np.int64).values
job = Gaia.launch_job_async("""SELECT u.source_id, g.parallax, g.parallax_error, g.ruwe, g.bp_rp, g.pmra, g.pmdec,
        g.non_single_star, g.ipd_frac_multi_peak, g.teff_gspphot, g.logg_gspphot
        FROM tap_upload.up AS u JOIN gaiadr3.gaia_source AS g ON g.source_id = u.source_id""",
        upload_resource=Table({'source_id': ids}), upload_table_name="up")
g = job.get_results().to_pandas()
m = m.merge(g, left_on=ids, right_on='source_id') if False else m.assign(sid=ids).merge(g, left_on='sid', right_on='source_id')
m['dist_pc'] = 1000/m.parallax
elig = m[(m.psffrac>0.5)&(m.crowdsap>0.7)&(m.nq_any>=8)&(m.kic!=m.tpforg)&(~m.is_kepler_target)] \
        .sort_values(['phot_g_mean_mag','gaia_designation'])             # reproduces Kuber's 1,800 (first 64 in his order)
```
(Chunk uploads of >~20k rows into 2–3 jobs if the Gaia server returns HTTP 500, as it did once here.)

---

## 3. Prior searches — how unsearched is it really?

### 3.1 Published work touching Kepler background sources
| Work | Data | What was searched | Relevance |
|---|---|---|---|
| Martínez-Palomera+2023 (KBonus paper) | KBONUS-BKG | No blind search; two KOI deblending demos + dilution study | Authors invite transit searches |
| Martínez-Palomera+2022 (EXBA, AJ 163, 93) | 9,327 EXBA-mask sources (different pixels) | 69 EBs, 1 planet candidate (0.81 R_J) | Shows background planet candidates are findable |
| Forró & Szabó 2020; Bienias+2021 ([arXiv:2107.04488](https://arxiv.org/abs/2107.04488), ApJS 256, 11); Forró+2022 | Raw background pixels of Kepler target apertures | 776 and 547 new (mostly short-period) EBs; 26 RR Lyrae | Pixel-level, not KBonus LCs |
| **Bienias & Szabó 2025**, A&A 694, A209, [arXiv:2501.09152](https://arxiv.org/abs/2501.09152) | Q4 pixels of ~150,000 KIC target apertures, then Q2–Q16 for candidates | Lomb–Scargle (<4 c/d, amp ≥4× local rms) plus a "≥3 dips >3.5 rms" search; Gaia host match within 12″; `Vetting` centroid tests | **7 hot-Jupiter candidates** (P 2.51–4.79 d, 0.89–1.52 R_J, hosts G≈16); stated minimum detectable depth 0.0025 (single-quarter pixels) / 0.0018 (merged), "planets below ~0.5 R_J almost entirely outside the detectable range"; they compared two hosts (KIC 4543171, KIC 9040849) against KBonus curves; psfmachine failed on KIC 9674230. **Sub-Neptunes around M dwarfs are below their floor.** |
| Barber+2022 (THYME VIII, AJ 164, 88, [arXiv:2206.08383](https://arxiv.org/abs/2206.08383)) | Kepler PDCSAP (long cadence) for 84 of 1,007 MELANGE-3 candidates + TESS | Notch pipeline; 19 targets of interest, 9 known KOIs, rest rejected/unconfirmed (KIC 6366739 open) | Table 5 flag "Kepler? = Y" means "had a Kepler light curve of sufficient quality for Notch" (CDS ReadMe). KIC 5689121 was never a Kepler target, so its Y must come from the coordinate-matched **neighbour KIC 5689116's PDCSAP** (inference). In that official PDCSAP I measure only 101 ± 27 ppm at Kuber's ephemeris **[sandbox]** (expected leakage ≤701 ppm × aperture fraction), far below anything Notch or the DR25 pipeline would flag (keplerstellar: `ntce=0`, `nkoi=0` for 5689116). So Barber's search did not test this signal. |
| Claytor & Tayar 2025 (rotation); Choi, Espinoza-Rojas & Hekker 2026, A&A 710 ([arXiv:2604.24476](https://arxiv.org/abs/2604.24476), asteroseismic binaries) | KBONUS-BKG | Rotation periods / oscillations | Not transit searches; useful systematics notes |
| Kuber Mehta 2026 **[Kuber]** | KBONUS-BKG M-star catalogue, 1,736 of 1,800 eligible | BLS 0.5–25 d, holdout quarters, inverted-flux nulls | Real first-screen passes 31/29 vs inverted 32/28; one candidate (KIC 5689121); other unresolved leads listed in his extension README (KIC 6629201 2.638 d; 12164160 19.42 d; 8743279 21.96 d; 8803342 5.010 d) |

Planet Hunters / Zooniverse: no project has used KBonus (web search 2026-10-17 finds none). No KOI/TCE/TOI/CTOI has ever been issued from KBONUS-BKG. **Conclusion: the KBONUS-BKG light curves have never been subjected to a published systematic transit search; the only overlapping work (Bienias & Szabó) is pixel-level with a ~2 ppt depth floor (giant planets only).**

### 3.2 Accounting of what remains **[sandbox]**
Source catalogue: 606,900. Kepler targets by `kic ∈ keplerstellar` (200,038 distinct kepids): 194,867; non-targets 412,033 (paper: 200,352/406,548 — the paper must define "target" slightly differently; `kic==tpforg` gives only 181,681, and 18,583 sources have `kic≠tpforg` yet are Kepler targets extracted from another target's TPF, so **do not use `kic != tpforg` as the background definition** — Kuber's "primary target" exclusion step was needed for exactly this reason).

Non-target sources passing PSFFRAC>0.5 & CROWDSAP>0.7 & ≥8 PSF quarters: **23,171** (G<17: 4,204+1,449+549+126+49+1 ≈ 6,380). By Gaia BP−RP: 881 (<0.8), 17,383 (0.8–1.4), 3,066 (1.4–1.84), 1,841 (≥1.84). Gaia DR3 join: 10,278 have parallax/error>5; of those **8,463 are dwarfs (M_G>4)** and 1,815 giants; `logg_gspphot>4` for 20,910. Dwarfs by G: 12–13: 2; 13–14: 16; 14–15: 143; 15–16: 668; 16–17: 2,806; 17–18: 3,215; 18–19: 1,613. Dwarfs by distance: 8 within 100 pc, 41 within 200 pc, 140 within 300 pc, 498 within 500 pc, 1,418 within 1 kpc.

M-star catalogue (29,800): 27,760 background by `kic≠tpforg`, 2,040 Kepler targets; median G 18.47 (25% brighter than 17.95); BP−RP 1.31–4.18; teff_gspphot 2,616–3,900 K; RUWE<1.4 for 29,105; `non_single_star>0` for 34. Kuber's cuts: psffrac>0.5 → 9,204; crowdsap>0.7 → 6,477; ≥8 quarters → 2,495; kic≠tpforg → 1,836; not a Kepler target → 1,800 (reproduced exactly, including the first-64 order). G of eligible: 14.19–19.00 (3 stars at G<15, 5 at 15–16, 67 at 16–17, 335 at 17–18, 1,390 at 18–19). Distances (plx/err>5): 2 within 100 pc, 30 within 200 pc, 100 within 300 pc. Not searched: 20,596 with psffrac≤0.5 (partial PSF; searchable but noisier), 2,727 with crowdsap≤0.7, 3,982 with <8 quarters, 695 Kepler-target M dwarfs (already searched by DR25), plus Kuber's 64 skipped/timeouts.

---

## 4. Where to hunt

1. **Nearby non-target dwarfs of any colour with clean extraction** (the real sweet spot): `is_kepler_target==False & psffrac>0.7 & crowdsap>0.8 & nq_good>=12 & parallax/parallax_error>5 & M_G>4 & G<17`. Expect 6-h CDPP 80–350 ppm; a 2 R⊕ planet on a 0.5 R☉ star (1,300 ppm) is ≥10σ per 6-h bin at G≈16 and a 1.5 R⊕ planet on a 0.3 R☉ star (2,100 ppm) is detectable to G≈17.5. Files written: `background_quality_pool.csv`, `background_quality_pool_gaia.csv` (23,171 rows with Gaia parallax/RUWE/teff/logg), `mstar_eligible_ranked.csv` (Kuber's 1,800 in his order), `mstar_nearby_nontarget_pool.csv` (68 M dwarfs <200 pc, psffrac>0.5, ≥6 quarters; nearest: Gaia DR3 2073565650248046848 at 70 pc G=14.3 but crowdsap 0.17; Gaia DR3 2105301782052654592 at 77 pc, G=14.19, psffrac 0.60, crowdsap 0.77, 18 quarters, TPF host KIC 7940546).
2. **The 20,596 partial-PSF M dwarfs** Kuber skipped: use `LIGHTCURVE_Qn` extensions with PSFFRAC in 0.3–0.5 and ITERNEG=0, and require consistency of depth across seasons (PSFFRAC varies by season, so a real transit must have a constant PSF depth while the SAP depth tracks FLFRCSAP×PSF/SAP).
3. **KOI false-positive neighbours**: the `kois-nns` bundle (35,639 LCs, 204 GB) holds every Gaia source within 30″ of every KOI. KOIs dispositioned FALSE POSITIVE with the centroid-offset flag (`koi_fpflag_co=1` in `cumulative`) have their eclipse on *some* neighbour; KBonus's own demo shows PSF photometry pins it. Re-attributing a shallow KOI FP to a fainter, smaller neighbour can turn an EB-flagged signal into a planet-sized one (depth scales with 10^(0.4Δm)). Query: `select kepid,kepoi_name,koi_period,koi_time0bk,koi_depth,ra,dec from cumulative where koi_disposition='FALSE POSITIVE' and koi_fpflag_co=1`.
4. **Re-extracted Kepler targets** (194,867 PSF curves): deblended depths change planet radii (paper §V); lower-value for *new* planets, but a cheap consistency check for any candidate.

ADQL (Gaia archive) for an independent pre-selection of Kepler-field dwarfs, then intersect by designation:
```sql
SELECT source_id, ra, dec, parallax, parallax_error, ruwe, phot_g_mean_mag, bp_rp, teff_gspphot, logg_gspphot
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), POLYGON('ICRS',279.5,36.5, 302.5,36.5, 302.5,52.5, 279.5,52.5))
  AND phot_g_mean_mag < 19 AND parallax > 3.3 AND parallax_over_error > 5
  AND phot_g_mean_mag + 5*log10(parallax/100) > 4
```
(Kepler FoV ≈ RA 280–301°, Dec 36–53°; returns ~10⁵ rows, run async.) The upload-join in §2.5 is faster when you already have the KBonus IDs.

Ranking recipe that avoids Kuber's "brightest first" bias toward distant/giant contaminated stars: score = expected transit SNR for a fiducial 2 R⊕ planet = (depth(R*) / CDPP6(G,colour)) × sqrt(N_quarters × 90 d / P × T14/6 h), with R* from Gaia M_G–BP−RP, and penalties for crowdsap<0.8, RUWE>1.4, `ipd_frac_multi_peak>5`, and for a TPF host that is a known EB/KOI.

---

## 5. Pitfalls specific to background PSF photometry

### 5.1 Aperture (SAP) vs PSF flux
SAP includes neighbour flux; `CROWDSAP` and `FLFRCSAP` change every season because the star lands on a different channel with a different aperture (2–7 px for Kuber's target). A real on-target transit therefore has **constant PSF depth** and a **season-dependent SAP depth ≈ PSF depth × (target share of aperture flux) = depth × FLFRCSAP × F_PSF / F_SAP** (Kuber's dilution test: χ² 15.4/12 for this model vs p=0.005 for a constant SAP depth **[Kuber SCIENCE.md]**). If instead the SAP depth is *larger* in the most contaminated seasons, the signal is on the contaminant. PSF and SAP are the same photons, not independent evidence.

### 5.2 The neighbour's signal leaking in — mandatory ephemeris matching
Every background source lives in a TPF centred on a Kepler target (`TPFORG`), and PSF deblending leaves residual correlation (Claytor & Tayar: up to 63% of periodic pairs). Before anything else, match the candidate period/epoch (and 1/3, 1/2, 1, 2, 3 harmonics, with timing tolerance grown over the baseline) against **(a) the TPF host, (b) every KIC within ~60″, (c) every Gaia source in the TPF footprint**. Exact queries that work today **[sandbox]**:
```sql
-- NASA Exoplanet Archive TAP, https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=...
-- (a) KOIs on the TPF host and the candidate's own KIC
select kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_time0bk,koi_depth,koi_duration,koi_fpflag_co,ra,dec
  from cumulative where kepid in (5689116,5689121)
-- DR25 threshold-crossing events (34,032 rows total) for the same kepids
select kepid,tce_plnt_num,tce_period,tce_time0bk,tce_depth,tce_duration,tce_max_mult_ev from q1_q17_dr25_tce where kepid in (5689116,5689121)
-- (b) all KOIs / TCEs in a box around the source (TCE table has no coordinates: join keplerstellar; Oracle syntax, NO 'AS' on table aliases;
--     ADQL CONTAINS/CIRCLE triggered a Cloudflare block from this sandbox, so use a RA/Dec box)
select kepid,kepoi_name,koi_disposition,koi_period,koi_time0bk,koi_depth,ra,dec from cumulative
  where ra between 285.7975 and 286.0181 and dec between 40.8587 and 41.0253
select t.kepid,t.tce_plnt_num,t.tce_period,t.tce_time0bk,t.tce_depth,t.tce_duration,s.kepmag
  from q1_q17_dr25_tce t join keplerstellar s on s.kepid=t.kepid
  where s.ra between 285.7975 and 286.0181 and s.dec between 40.8587 and 41.0253
-- note: keplerstellar has one row per delivery (q1_q12 … dr25_supp) so DISTINCT or filter st_delivname='q1_q17_dr25_stellar'
```
Kepler eclipsing binaries: the Villanova site's CSV endpoint returned HTTP 400 here; use VizieR `J/AJ/151/68/catalog` (Kirk+2016, 2,876 EBs; columns KIC, Per, BJD0, Morph, Kpmag) via `astroquery.vizier.Vizier(...).query_region(coord, radius=5*u.arcmin, catalog='J/AJ/151/68/catalog')` **[sandbox, works]**. Also check `kois-nns` (`nn_koi`, `nn_dist`) for the candidate's designation, and Gaia `non_single_star`, `ruwe`, `ipd_frac_multi_peak` for the host. For Kuber's star the 5-arcmin box returns Kepler-169 b–f (KIC 5689351) and KOI-3784.01 (KIC 5688997, FP, 23.87 d); Kepler-169 e at 13.76710058 d differs from the candidate's 13.769448 d by 32σ and the ephemerides never overlap (nearest events 67 h apart) **[Kuber]**.

### 5.3 Quarter-to-quarter offsets and roll-angle (season) effects
Kepler rolls 90° each quarter, so the same star repeats on the same channel every 4th quarter (seasons 0–3). PSFFRAC, aperture, CROWDSAP, rolling-band bits and even detection all depend on season. For Kuber's target Q1 sat on **module 3 (channel 6)**, which failed on 2010-01-09 during Q4; every later season-3 quarter (Q5, Q9, Q13, Q17) is therefore absent (QDETECT `011110111011101110`), a ~20% duty-cycle loss shared by ~1/5 of the field ([Kepler documentation; see e.g. arXiv:1202.5852, 1501.03586]). Consequences: (i) never let a BLS fold mistake season gaps for transit windows (periods near 90 d × n/4 aliases); (ii) do odd/even, depth and shape consistency **per season**, not just per quarter; (iii) use `LIGHTCURVE_Qn` headers to split by `SEASON`/`CHANNEL`.

### 5.4 Centroid columns are useless for localisation
`CENTROID_COLUMN/ROW` = Gaia position + a scene-wide jitter vector: in every quarter the target and the neighbour have **identical** centroid standard deviations (e.g. Q1 (0.0116, 0.0256) px for both) and differ only by their fixed Gaia offset **[sandbox]**; psfmachine computes them via `get_source_centroids(method="scene"/"poscor")`, i.e. the Gaia coordinates shifted by the per-cadence scene centroid or by the TPF POS_CORR columns. Kuber's `kcentroid.py` found χ²(zero shift)=25.5/26, i.e. no in-transit response at all **[Kuber]**. Do not use them for vetting.

### 5.5 Pixel-level vetting with the originating TPF (do this instead)
Download the TPF of `TPFORG` (`lk.search_targetpixelfile("KIC 5689116", author="Kepler", cadence="long").download_all()`, 13 files here), fold at the candidate ephemeris and build per-event difference images (out-of-transit median minus in-transit mean), then compare the difference-image centroid / PRF fit with the Gaia positions of all stars in the stamp via the TPF WCS (`WCS(tpf.hdu[2].header)`). Demo result **[sandbox, `tpf_vet.py`]**: in all 13 quarters the difference-image centroid lies 0.54 px (median) from KIC 5689121 vs 1.30 px from KIC 5689116 (1 px = 3.98″), and the peak-dip pixel is the target's pixel in every quarter — consistent with Kuber's PRF-fit localisation (Δχ² = 934 against the neighbour) and raw-TPF flux-weighted centroid test (χ² 43.9/24 for the target model) **[Kuber]**. Minimal code:
```python
import numpy as np, lightkurve as lk
from astropy.wcs import WCS
P, T0, dur = 13.769448, 2455630.245431-2454833, 2.3/24          # BKJD
tgt, nb = (285.907827, 40.941988), (285.905233, 40.942079)
for tpf in lk.search_targetpixelfile("KIC 5689116", author="Kepler", cadence="long").download_all():
    t, F, ok = tpf.time.value, tpf.flux.value, tpf.quality == 0
    ph = ((t - T0 + P/2) % P) - P/2; n = np.round((t - T0)/P)
    intr = ok & (abs(ph) < 0.4*dur); oot = ok & (abs(ph) > dur/2 + 0.05) & (abs(ph) < 0.5)
    diffs = [np.nanmedian(F[oot & (n==k)], 0) - np.nanmean(F[intr & (n==k)], 0)
             for k in np.unique(n[intr]) if (intr & (n==k)).sum() >= 2 and (oot & (n==k)).sum() >= 10]
    diff = np.nanmean(diffs, 0); w = WCS(tpf.hdu[2].header)
    xt, xn = w.all_world2pix([tgt], 0)[0], w.all_world2pix([nb], 0)[0]
    Y, X = np.mgrid[:diff.shape[0], :diff.shape[1]]; d = np.clip(diff, 0, None)
    cx, cy = (d*X).sum()/d.sum(), (d*Y).sum()/d.sum()
    print(tpf.quarter, len(diffs), f"centroid=({cx:.2f},{cy:.2f}) target={xt.round(2)} neighbour={xn.round(2)}")
```
Upgrade path: fit the official Kepler PRF (`lightkurve.prf.KeplerPRF(channel, shape, column, row)`) to the difference image at each Gaia position and compare χ² (Kuber's `pixel_localize.py` **[Kuber]**); also fold the *neighbour's own* KBONUS PSF curve (Kuber: 24 ± 41 ppm) and the official PDCSAP of the TPF host (here 101 ± 27 ppm, the expected small leak).

### 5.6 Other traps
* `psf_avail` counts quarters with PSFFRAC<0.5 too; use `psf_fraction_flag.count('1')` for *stitched* quarters.
* Sources with empty PSF columns but valid SAP exist (negative-flux iteration); `ITERNEG` flags them per quarter.
* Faint sources inside 15″ of a star ≥2 mag brighter are the ones most likely to be over-/under-subtracted; a high-contrast pair can show an *inverted* (brightening) copy of the neighbour's eclipse — inverted-flux null runs catch these only partially; always fold the neighbour.
* G<10 stars were removed entirely and pixels near G≤8 stars masked; bright-star fields are simply absent.
* MAST page says "~20,000" M dwarfs and the paper "about 27,500 new M-dwarf light curves"; the file has 29,800 (2,040 of them Kepler targets). Quote the file.
* The catalogue `ra/dec` are J2016 Gaia positions; Kepler observed 2009–2013, so propagate proper motion for high-PM nearby stars before pixel work (and remember the header PM unit bug).
* `FLUX_ERR` is a formal error; recompute noise empirically for SNR/MES (see §1.5).

---

## 6. Sources
* MAST HLSP page: https://archive.stsci.edu/hlsp/kbonus-bkg ; README: https://archive.stsci.edu/hlsps/kbonus-bkg/hlsp_kbonus-bkg_kepler_kepler_all_kepler_v1.0_readme.md ; directory index https://archive.stsci.edu/hlsps/kbonus-bkg/ (bundles/, catalogs/, lcs/, bulk scripts).
* Paper: https://arxiv.org/abs/2310.17733 ; HTML https://arxiv.org/html/2310.17733v1 ; DOI 10.3847/1538-3881/ad0727 ; Crossref API record (13 citations).
* psfmachine: https://github.com/SSDataLab/psfmachine (README, `src/psfmachine/tpf.py::get_source_centroids`), PyPI JSON (1.1.4, 2023-12-15), docs https://ssdatalab.github.io/psfmachine/ ; KBonus repo https://github.com/jorgemarpa/KBonus (README, `notebooks/light_curve_caveats.ipynb`).
* APEXBA: https://archive.stsci.edu/hlsp/kbonus-apexba ; https://arxiv.org/abs/2112.06683.
* Bienias & Szabó 2025: https://arxiv.org/abs/2501.09152 (HTML v1 read); Bienias+2021 https://arxiv.org/abs/2107.04488.
* Barber+2022: https://arxiv.org/abs/2206.08383 (ar5iv HTML), CDS ReadMe https://cdsarc.cds.unistra.fr/ftp/J/AJ/164/88/ReadMe, VizieR J/AJ/164/88/table5 row 546.
* Claytor & Tayar 2025: https://arxiv.org/abs/2506.03248 ; poster https://zenodo.org/records/17245793. Choi+2026: https://arxiv.org/abs/2604.24476.
* Kepler module-3 failure: Kepler Data Characteristics Handbook (Van Cleve+2016) and e.g. https://arxiv.org/abs/1202.5852, https://arxiv.org/abs/1501.03586.
* NASA Exoplanet Archive TAP: https://exoplanetarchive.ipac.caltech.edu/TAP/sync (tables `cumulative`, `q1_q17_dr25_tce`, `keplerstellar`, `TAP_SCHEMA.tables`). VizieR J/AJ/151/68 (Kirk+2016).
* Kuber's repo (read-only): `analysis/editorial-audit-20261013/{README.md,verification.json,verify.py}`, `analysis/kbonus-extension-20261011/README.md`, `analysis/validation-20261011/SCIENCE.md`, `analysis/validation-20261011/scripts/kepler/kcentroid_tpf.py`, `evidence/codex/database_expansion/kbonus/{pilot.py,pixel_localize.py,crosscheck.py,known_control.py}`, `evidence/codex/database_expansion/bienias2025_appendixA.txt`.

## 7. Scratch files written (all under `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/`)
`search_kbonus.py`, `dl_inspect.py`, `noise.py`, `cdpp_compare.py`, `flatmask.py`, `catstats.py`, `gaia_xmatch.py`, `hunting.py`, `pick_sample.py`, `sample_cdpp.py`, `mast_access.py`, `tpf_vet.py`; data: `catalogs/*.fits` (4 catalogues + readme), `downloads/mastDownload/HLSP/...` (KBONUS files for KIC 5689121 and 5689116), `downloads/mastDownload/Kepler/...` (official LCs and TPFs of KIC 5689116), `sample/*.fits` + `sample/sample_cdpp.csv` (28-star CDPP sample), `source_catalog.parquet`, `mstars_joined.parquet`, `mstars_gaia.parquet`, `keplerstellar_kepids.csv`, `mstar_eligible_ranked.csv`, `mstar_nearby_nontarget_pool.csv`, `background_quality_pool.csv`, `background_quality_pool_gaia.csv`, `kbonus_repo/light_curve_caveats.ipynb`.