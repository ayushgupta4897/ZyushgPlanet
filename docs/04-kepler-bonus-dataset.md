# 04 · Kepler Bonus as a planet-hunting dataset

> Adversarially verified research brief (2026-10-17) on the KBONUS-BKG high-level science product: what it contains, how to read the files, where its photometry breaks, who has already searched it, and how to select targets. The checker opened real FITS files, re-counted the catalogue, and re-queried the archive.

**How this repo uses it.** Pool A (`targets/poolA_kbonus_bg_dwarfs_ranked.csv`) is built from the source catalogue exactly as described here: background extractions only, PSF fraction > 0.5, crowding > 0.7, at least 8 PSF quarters, Gaia-classified dwarfs, own KIC never a Kepler target, Kuber's 1,800 removed. Light curves are read from the `LIGHTCURVE_STITCHED` extension (`FLUX` = PSF flux, BKJD time base) and re-flattened; the per-quarter extensions and their aperture masks are used in stage-2 vetting.

# Kepler Bonus background light curves (KBONUS-BKG) as a planet-hunting dataset — verified practitioner brief

Date of checks: 2026-10-17. Everything marked **[sandbox]** was run in the repository's `.venv` (Python 3.13, lightkurve 2.6.0, astropy 8.0.2, astroquery 0.4.11, numpy 2.5.4). The scripts and outputs of this verification pass (`v1_fits.py` … `v10_pdcsap.py` and their outputs, listed in section 7) and the first-pass material it builds on (catalogues, downloaded FITS, 28-star sample, Gaia joins) were kept in a scratch area outside this repository and are not published. Everything marked **[Kuber]** comes from reading (never executing) the read-only clone of `Kuberwastaken/planethunt`. Inferences are labelled.

## 0. Executive summary

* **What it is.** KBONUS-BKG is a single v1.0 MAST HLSP (released and last updated 2023-11-01; DOI 10.17909/7jbr-w430; CC BY 4.0) of **606,900 Kepler long-cadence light curves** extracted with PSF photometry (psfmachine / Linearized Field Deblending) from all Kepler quarters' TPFs with Gaia DR3 (G<19) as the input catalogue: 406,548 background sources never processed by the Kepler pipeline plus 200,352 re-extracted Kepler targets (MAST page; paper arXiv:2310.17733 v1, AJ 166, 265). No later data release exists **[sandbox: MAST page text, arXiv API]**.
* **Someone is already searching it.** A 2026 Sagan Summer Workshop poster by the KBonus team (Moquin, Wilson, Quintana, Hedges, Martínez-Palomera, Claytor, Gilbert, Rowe, Hoffman, Kostov, Ciardi) presents **TraSHDUMP**, a matched-filter TCE pipeline on the KBonus light curves, with Robovetter-style and pixel-level vetting, a first candidate **BOI-1.01** (G=17.8 M dwarf, P=3.37 d, 2.94 R⊕) and an *expected yield of several hundred new planet candidates* ([poster PDF](https://vmnexsci.ipac.caltech.edu/workshop/2026/posters/Poster_JulieMoquin_41.pdf)). No refereed paper or candidate list was found as of 2026-10-17. The only other transit work on Kepler background stars (Bienias & Szabó 2025) is pixel-level with a ~2 ppt floor (giant planets). **So: no refereed systematic search exists, but the dataset's authors have an active pipeline. Move fast, differentiate (nearby dwarfs, KOI-FP neighbours, M dwarfs below their MES threshold), and expect to be scooped on the easiest targets.**
* **Precision.** Comparable to, but slightly worse than, the Kepler pipeline for the same star: KIC 5689116 (G=15.48) gives 6-h CDPP 184 ppm (KBONUS PSF, stitched) vs 121 ppm (official PDCSAP); Kuber's G=16.28 target ~475 ppm (337–912 per quarter) **[sandbox, re-run]**. A 28-star background sample gives median per-quarter 6-h CDPP of roughly 90/75/80 ppm at G 12–15, ~300 at 15–16, ~330 at 16–17, ~690 at 17–18, ~1,600 at 18–19 (small sample, includes variables). The paper's "<100 ppm for G<16" describes the best stars, not the median.
* **The M-star catalogue is faint, distant and half-contaminated.** 29,800 entries, median G 18.47, Teff 2616–3900 K, only 182 within 100 pc (6 of them eligible non-targets). Kuber's 1,800 eligible stars span G 14.19–19.0 (rank 64 = G 16.88); only 100 are within 300 pc. **New in this pass:** 801 of the 1,800 (45%) have a catalogue PSF flux more than 2× brighter than their Gaia G allows — a signature of residual flux from the brighter TPF host, i.e. incompletely deblended light curves **[sandbox, §1.6]**.
* **The better hunting ground.** Among non-Kepler-target sources, 23,157 pass PSFFRAC>0.5 & CROWDSAP>0.7 & ≥8 usable quarters; 15,473 of those also pass the flux-consistency test (G<17: 5,466; G<16: 2,000). A Gaia DR3 cross-match of the 23k pool shows 8,463 dwarfs (M_G>4, plx/err>5), 3,635 of them with G<17 and 140 within 300 pc **[sandbox]**.
* **Dataset pitfalls (all verified):** catalogue column names differ from the docs (`tpforg`, `kep_mag`, `fname`); `psf_avail` and `psf_fraction_flag` encode different things and the usable-quarter count is their AND; flag-`2` quarters still get HDUs; the stitched `FLUX` is rescaled per quarter to the catalogue mean and is populated even where the per-quarter PSF flux is NaN; `PSF_FLAT_FLUX` is a 2-day b-spline flatten that was not transit-masked (shaves Kuber's dip from ~2.2 to ~1.7 ppt); `FLATTEN_MASK` is an outlier mask; `SAP_QUALITY` bits present are all outside lightkurve's default bitmask; `CENTROID_COLUMN/ROW` are Gaia positions plus a scene-wide jitter vector (useless for localisation); header `PMRA/PMDEC` are Gaia values ÷1000; lightkurve 2.6.0 ignores `flux_column=` and rejects `quarter=`; `estimate_cdpp` works but returns a ppm Quantity; the MAST page's "~20,000" M dwarfs and "29 MB" bulk script are wrong (29,800; 168 MB).

---

## 1. The paper and the product

### 1.1 Reference chain
* Martínez-Palomera, Hedges & Dotson 2023, "Kepler Bonus: Light Curves of Kepler Background Sources", AJ 166, 265, DOI [10.3847/1538-3881/ad0727](https://doi.org/10.3847/1538-3881/ad0727), [arXiv:2310.17733](https://arxiv.org/abs/2310.17733) (v1 2023-10-26, the only version; Crossref: published 2023-11-24, CC BY 4.0). Crossref `is-referenced-by-count` = 13, OpenAlex `cited_by_count` = 12 **[sandbox, 2026-10-17]**. Citing works enumerated before OpenAlex rate-limited: Claytor & Tayar 2025 (rotation); Wang & Esamdin 2025 AJ (KOI-1755 background EB, DOI 10.3847/1538-3881/ade718); "Exploring Long-period Architectures: Four New Planet Candidates from Kepler with Periods >342 days" (arXiv:2608.23425, a CNN single-transit search of *known* Kepler systems, not KBonus); "Why Estimating η⊕ is Difficult" (arXiv:2511.05658); THYME XIV (arXiv:2603.01313); two rotation papers (arXiv:2509.14423, 2507.08266); a Pandora simulation SPIE paper. None is a KBonus transit search.
* Method paper: Hedges et al. 2021, "Linearized Field Deblending", AJ 162, 107 (cited as [3] in the HLSP README). Software: [psfmachine](https://github.com/SSDataLab/psfmachine) (PyPI latest **1.1.4**, 2023-12-15; releases 0.1.0, 1.0.0, 1.1.0, 1.1.4; the FITS headers record `ORIGIN='PSFMACHINE'`, `VERSION='1.1.2'`, `PROCVER='2.0.9'`, `DATE='2022-07-25'`, checksums 2023-03-31 **[sandbox]**). Source lives on branch **`master`** (`main` 404s). Configuration/notebooks: [jorgemarpa/KBonus](https://github.com/jorgemarpa/KBonus) (`light_curve_caveats.ipynb` defines QDETECT).
* M-dwarf selection recipe (paper §III.2): Bentley, Tinney & Sharma 2018, [arXiv:1804.06982](https://arxiv.org/abs/1804.06982) (Gaia+WISE+2MASS photometric M-dwarf selection) plus Gaia Teff cuts → photometric candidates, not a spectroscopic sample.
* Sister product (different pixels): KBONUS-APEXBA — aperture photometry of **9,327** sources in the EXBA background masks, Martínez-Palomera et al. 2022, AJ (DOI 10.3847/1538-3881/ac4331, [arXiv:2112.06683](https://arxiv.org/abs/2112.06683)), HLSP [kbonus-apexba](https://archive.stsci.edu/hlsp/kbonus-apexba), DOI 10.17909/t9-d5wy-e535; 69 EBs and one 0.81 R_J candidate (Gaia EDR3 2077240046296834304) **[sandbox: MAST page, arXiv abstract]**.

### 1.2 What was done (paper text, `arxiv_src/x/main.tex`)
* Input: 204,933 long-cadence TPFs (l.237), Q0–Q17. Source list: Gaia DR3 to G=19 within the TPF + 4″ margin; G<10 removed (saturation); of pairs within 1″ the fainter dropped (l.255–258). Result 606,900 sources; **"200,352 corresponds to Kepler targets for which the Kepler pipeline produced light curve files"**, identified by a 2″ KIC–Gaia cross-match with proper motion (l.260–262).
* Background model: cubic b-splines in time (knots every 2 h) and pixel row (knots every 6 px) fit to background pixels (l.322). Perturbed PSF built from the first four CBVs smoothed with a 2-day b-spline (l.443). Negative-flux handling: priors for the source and neighbours within 5″ narrowed by 2×, three iterations; remaining negatives dropped (<5% of input) (l.507–518).
* Stitched product: "aperture, corrected PSF, and a flattened version of the PSF photometry … flattened with a 2-day window b-spline function, designed to better enable the community to perform transit searches" (l.546–547). Transit masking is mentioned only for the dilution comparison (l.703), consistent with the measured ~25% attenuation of Kuber's dip (§2.3).
* Only quarters with PSFFRAC ≥ 0.5 are stitched; **"We still include all the extracted quarters in the light curve FITS file, regardless of low PSFFRAC values"** (l.578–579, l.806). Verified on a file (§1.3).
* Demonstrations: KOI-770.01 and KOI-909.01 deblending; PSF-vs-PDCSAP transit-depth dilution. **They did not run a blind search** in the paper — but see §3 for their 2026 pipeline.

### 1.3 Flags, header keywords and what they really mean **[sandbox, verified on files]**
| Field | Meaning |
|---|---|
| `PSFFRAC` (per-quarter header) / `psffrac` (catalogue = minimum over detected quarters) | fraction of the source PSF on the TPF pixels. Only quarters with PSFFRAC>0.5 are stitched. Because the catalogue value is the minimum, `psffrac>0.5` implies no quarter has flag `2` (0 exceptions in 606,900 rows). |
| `FLFRCSAP`, `CROWDSAP`, `NPIXSAP` | PSF-model-derived aperture flux fraction, crowding (1 = clean) and pixel count; `APERTURE_Qn.sum() == NPIXSAP`. Both FLFRCSAP and CROWDSAP are *over*-estimated for partial sources (paper l.563–566). |
| `PERTRATI`, `PERTSTD` | mean-PSF/perturbed-PSF flux ratio and std ratio; pool quantiles 0.9983–1.0001 and 0.0002–0.0102 (1–99%), so they rarely discriminate. |
| `ITERNEG` (per-quarter header) | 1 if the source went negative during the prior-narrowing iteration. Observed effect: `PSF_FLUX_NOVA` (mean-PSF flux) is entirely NaN in ITERNEG=1 quarters while corrected `FLUX` is populated. |
| `QDETECT` (primary header) ≡ `psf_fraction_flag` (catalogue) | 18 chars Q0..Q17: `0` no PSF detection / no HDU; `1` PSFFRAC>0.5, HDU present, in the stitched FLUX; `2` PSFFRAC<0.5, **HDU present but NaN in the stitched FLUX** (KIC 3459079: Q1 PSFFRAC=0.014, flag `2`, LIGHTCURVE_Q1 has 1,624 finite FLUX values, stitched Q1 rows have 0 finite FLUX). Definition from the KBonus caveats notebook; identical to the catalogue string in every file checked. |
| `psf_avail` (catalogue) | `1` if the per-quarter `LIGHTCURVE_Qn.FLUX` column is populated. It is `1` for 1,417,775 flag-`2` quarter-slots and `0` for 358,493 flag-`1` slots (PSF fit failed; per-quarter FLUX all NaN, SAP finite, **stitched FLUX nevertheless populated**). |
| `sap_avail` | `1` if SAP exists; SAP exists without PSF for 757,125 quarter-slots. |
| **Usable stitched quarters** | `sum(psf_avail[i]=='1' and psf_fraction_flag[i]=='1')`. Counting `psf_fraction_flag` alone over-counts (25,709 vs 23,157 pool sources at ≥8); counting `psf_avail` alone is only safe when `psffrac>0.5`. |
| `psf_mean_flux` | NaN for 156,955 sources (155,632 background) that have SAP only; 0 in the quality pool. |

### 1.4 Known systematics (paper §IV, caveats notebook, Claytor & Tayar 2025)
* Partial-PSF sources (TPF edge, saturation masks) are noisier and their FLFRCSAP/CROWDSAP are biased; PSFFRAC and the aperture change with season because the star lands on a different channel each quarter.
* PSF model independent of brightness; distorted PSFs on edge channels; high-contrast blends (≥2 mag within 15″) can go negative or be over-corrected; ~1% of faint neighbours were removed because their light curves correlated (r>0.5, <25″) with a brighter source.
* Claytor & Tayar 2025 ([arXiv:2506.03248](https://arxiv.org/abs/2506.03248); poster [Zenodo 17245793](https://zenodo.org/records/17245793)), 449,913 KBonus light curves (hardest bitmask, relaxed to default for 7 stars): of 3,224 foreground–background pairs with periods for both, 55% agree within 10%; vs 3,000 random control pairs, 53% are positively and 10% negatively correlated at p<0.15, so **"as many as 63% of light curve pairs with rotation periods are likely to remain blended"**; anticorrelated pairs are slightly over-represented (over-subtraction); **"the cross-correlation statistic did not show significant trends with source separation, flux contrast, or Kbonus quality metrics such as PSFFRAC"**; a spurious ~40-d period cluster "may be related to systematics around downlink times" **[sandbox: arXiv HTML read]**. For transits this is the neighbour-leakage problem of §5.2 — and it means the KBonus quality metrics alone do not certify independence.

### 1.5 Measured precision **[sandbox]**
6-h CDPP = std of 13-cadence running means after a 101-cadence Savitzky–Golay flatten and 5σ clip (this is exactly lightkurve's `estimate_cdpp`, which **works** on these files — 472.9 ppm for the target — but returns an astropy Quantity in ppm; `float()` of it gives 4.5e-4, which the first-pass script printed as "0 ppm". Use `.value`).

| Star | Column | 6-h CDPP stitched | per quarter |
|---|---|---|---|
| KIC 5689121 (G 16.28, Kuber's target) | PSF `FLUX` | 475 ppm | 337–912 |
| | `SAP_FLUX` (diluted) | 384 | 227–589 |
| | `PSF_FLAT_FLUX` | 402 | – |
| KIC 5689116 (G 15.48, TPF host) | KBONUS PSF `FLUX` | 184 | 103–213 |
| | KBONUS `SAP_FLUX` | 159 | – |
| | KBONUS `PSF_FLAT_FLUX` | 133 | – |
| | official PDCSAP (`author='Kepler'`, 13 quarters) | 121 | 96–151 (pipeline `CDPP6_0` header 129–177; lightkurve's SG estimate lacks the Van Cleve normalisation) |

28-star sample (non-targets, PSFFRAC>0.7, CROWDSAP>0.8, ≥12 quarters; medians of per-quarter 6-h CDPP, PSF/SAP): G 12–13: 92/79; 13–14: 73/54; 14–15: 82/74; 15–16: 298/161; 16–17: 327/377; 17–18: 686/599; 18–19: 1598/758 ppm (`sample/sample_cdpp.csv`; several reddened giants at G≈13 reach 250–600 ppm). Formal `FLUX_ERR` is typically 0.6–0.8× the measured point-to-point scatter at G≈16 and larger than it at G≈18: re-estimate noise empirically for any SNR/MES.

### 1.6 NEW — flux-vs-magnitude consistency test (deblending-failure detector) **[sandbox]**
For Kepler targets with 12<G<16, `log10(psf_mean_flux) + 0.4·G` = 10.11 with 0.03 dex scatter (zero point varies 10.06–10.13 with BP−RP; F(G=12) ≈ 2.0×10⁵ e⁻/s). A PSF flux of a point source cannot legitimately exceed what its Gaia G implies by 2×, yet in the background quality pool **32.8%** of sources do (1.2% at G<14, 7.9% at 14–16, 17% at 16–17, 30% at 17–18, 47% at 18–19). The excess is leakage from the TPF host: the too-bright fraction vs contrast ΔG = G_src − G_host is 0.1% (1–2 mag), 5% (2–3), 24% (3–4), 49% (4–5), 70% (>5 mag), and the median excess equals 8% of the host's PSF flux. Example: Gaia DR3 2102612715913245696 (G=18.8, psffrac 0.77, crowdsap 0.86 — passes Kuber-style cuts) sits 7.4″ from a G=12.1 host and has psf_mean_flux 13,013 e⁻/s (1.5 dex too bright); its per-quarter PSF medians range 70–4,488 e⁻/s with ITERNEG=1 in 8 quarters, while the stitched FLUX looks like a flat 13,000 e⁻/s series. **801 of Kuber's 1,800 eligible M dwarfs (45%; 14 of his first 64) fail |Δdex|>0.3; KIC 5689121 passes (0.006 dex).** Recommended cut: `abs(log10(psf_mean_flux)+0.4*G − zp(BP−RP)) ≤ 0.3`; clean pool = 15,473 sources (G<17: 5,466; G<16: 2,000), file `background_pool_fluxconsistent.csv`. Caveat: the test also flags genuinely variable/binary sources and is calibrated on Kepler targets only.

---

## 2. Access

### 2.1 Catalogue files (direct URLs; sizes and SHA-256 verified **[sandbox]**, identical to Kuber's `verification.json`)
```
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits   159,022,080 B  sha256 2e6d5a87cb90f941bdc207aef93d61847df228ae7d56ffd38d3ab9255dd30fcf  (606,900 rows)
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits      2,390,400 B  sha256 406563951c367bb1ffe3c20f83c4684a4083a2fe2f3d7c36ad9a663b5f6977ef  (29,800 rows)
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_kois-nns-catalog_kepler_v1.0_cat.fits    3,430,080 B  sha256 411b3ce12c6e64aa29bbd79750aacb53623620aac0b0329a1361b55f814b8264  (35,639 rows; 8,073 KOIs; nn_dist ≤ 30″)
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_wd-catalog_kepler_v1.0_cat.fits             14,400 B  sha256 13071542df2172dd12084e4c5ca3e40e76301728ea597a801b91df6ff96b1938  (91 rows)
https://archive.stsci.edu/hlsps/kbonus-bkg/hlsp_kbonus-bkg_kepler_kepler_all_kepler_v1.0_readme.md                           14,178 B
```
Actual source-catalogue columns (README says `tpf_org`, `kepmag`, `file_name`): `gaia_designation` ("Gaia DR3 …"), `ra`, `dec` (J2016), `sap_mean_flux[_err]`, `psf_mean_flux[_err]`, `flfrcsap`, `crowdsap`, `npixsap`, `psffrac`, `pertrati`, `pertstd`, `psf_avail`, `sap_avail`, `psf_fraction_flag` (18-char strings), `phot_g_mean_mag`, `phot_bp_mean_mag`, `phot_rp_mean_mag` (±1e20 sentinels occur), **`tpforg`** (KIC of the TPF the source was extracted from), `kic` (0 for 65,501 rows), `kic_sep` (arcsec), **`kep_mag`**, **`fname`** (`kic-005689121` or `gaia-dr3-…`). M-star/WD catalogues: `gaia_designation, ra, dec, kic, fname` (+`Pwd`); kois-nns adds `nn_koi, nn_dist`. Catalogue facts: G 10.00–19.00; RA 279.63–301.83, Dec 36.53–52.47; per-quarter PSF availability 212,536 (Q0), 309–342k (Q1–Q17); usable-quarter count is bimodal (164,198 with 0; 76,751 with 1; 82,853 with 17; 29,670 with 18); `psffrac>0.5` 47%, `crowdsap>0.7` 60%, both 33%.

### 2.2 Light-curve files
* Naming `hlsp_kbonus-bkg_kepler_kepler_<kic-NNNNNNNNN | gaia-dr3-NNN…>_kepler_v1.0_lc.fits`; path `lcs/<first 4 digits of the 9-digit KIC or of the Gaia ID>/<full id>/<file>`. Both URL forms return HTTP 200 (8,847,360 B for the target; Last-Modified 2023-03-31):
  * `https://archive.stsci.edu/hlsps/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits`
  * `https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits`
  * Gaia-only example: `.../lcs/2051/2051766354744241536/hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2051766354744241536_kepler_v1.0_lc.fits` (11,574,720 B).
* Volume: 4.8 TB total (MAST page), 4–12 MB per file. Bulk scripts at the HLSP root **[sandbox HEAD]**: `hlsp_kbonus-bkg_kepler_lightcurves_all_kepler_v1_bulkdl.sh` **175,634,660 B** (MAST page wrongly says 29 MB; `curl -f --create-dirs` lines using the MAST API URL form), `…_chunkdl.zip` 6,431,250 B (50,000 files per script), `…_lftp.sh` 1,018 B. Bundles: `bundles/hlsp_kbonus-bkg_kepler_kepler_{kois-nns,mstars,wd}_kepler_v1.0_lcf.tar.gz` (mstars = 107,156,587,077 B ≈ 107 GB; page says 134 GB). MAST Portal caps (50,000 load / 500,000 download) are below the total.
* astroquery **[sandbox]**: `Observations.query_criteria_count(provenance_name="KBONUS-BKG")` → 606900 (4.2 s); `Observations.query_criteria(provenance_name="KBONUS-BKG", target_name="Gaia DR3 2103792461826460288")` → one observation (`obs_collection='HLSP'`, `project='Kepler'`, `dataproduct_type='timeseries'`, `t_exptime=1765.464`, obsid 186825724); `get_product_list` → the single `_lc.fits`. lightkurve: `lk.search_lightcurve(SkyCoord(285.907827,40.941988,unit='deg'), author='KBONUS-BKG', radius=20*u.arcsec)` → 2 rows (target 0.0″, neighbour listed at 3.1″ because MAST stores TPF-centre coordinates), `exptime` 1765 s, one row per source. Kuber used the raw `Mast.Caom.Filtered` service (provenance_name, target_name list, pagesize 50,000) **[Kuber pilot.py]**.

### 2.3 FITS layout (**[sandbox]**, `v1_fits.py`)
```
No.  Name                 Type         Dims / rows
0    PRIMARY              PrimaryHDU   48 cards: DOI, HLSPID, HLSPVER='V1.0', LICENSE='CC BY 4.0', ORIGIN='PSFMACHINE', VERSION='1.1.2', PROCVER='2.0.9',
                                       DATE='2022-07-25', TIMESYS='TDB', OBJECT='KIC 5689121', KEPLERID, TPFORG=5689116, GAIAID, RA_OBJ/DEC_OBJ (EQUINOX 2016.0),
                                       PMRA, PMDEC (bug, see below), PARALLAX [mas], GMAG, BPMAG, RPMAG, KEPMAG, QDETECT='011110111011101110',
                                       XPOSURE=1765.464, MJD-BEG/END, APERTURE='PSF + SAP', RADESYS='ICRS'
1    LIGHTCURVE_STITCHED  BinTable     50182 rows × 11: TIME [BJD-2454833], CADENCENO, QUARTER, FLUX, FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR,
                                       PSF_FLAT_FLUX, PSF_FLAT_FLUX_ERR, SAP_QUALITY, FLATTEN_MASK   (e-/s)
2,4,… LIGHTCURVE_Qn       BinTable     13 cols: CADENCENO, TIME, FLUX (corrected PSF), FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR, CENTROID_COLUMN, CENTROID_ROW,
                                       SAP_QUALITY, SAP_BKG, RED_CHI2, PSF_FLUX_NOVA (mean-PSF), PSF_FLUX_ERR_NOVA
                                       header: SEASON, MODULE, OUTPUT, CHANNEL, QUARTER, TPF_ROW, TPF_COL, SAP='optimal', ROW, COLUMN,
                                       FLFRCSAP, CROWDSAP, NPIXSAP, PSFFRAC, PERTRATI, PERTSTD, ITERNEG
3,5,… APERTURE_Qn         ImageHDU     uint8 mask in the originating TPF's shape ((5,4) or (5,5) here); sum == NPIXSAP
```
HDUs exist for every quarter whose flag is `1` **or `2`** (flag-`2` quarters are NaN in the stitched FLUX). For KIC 5689121: 13 quarters (Q1–4, 6–8, 10–12, 14–16).

**Time**: `TIME` is BKJD = BJD_TDB − 2,454,833 (lightkurve reads `format='bkjd', scale='tdb'`); first cadence 131.5128 BKJD = BJD 2,454,964.5128; cadence 29.424 min. BTJD = BKJD − 2167; Kuber's T0 = BJD 2,455,630.245431 = BKJD 797.245431 = BTJD −1369.754569 **[Kuber pilot.py uses BKJD−2167]**.

**Quarter normalisation**: stitched `FLUX` = per-quarter `FLUX` × (psf_mean_flux / quarter level): stitched/per-quarter ratios 0.970, 0.925, 1.086 (Q1–Q3) with std <1e-2, and the stitched median equals the catalogue `psf_mean_flux` (3991.029) exactly. Per-quarter `FLUX` medians 3596–4309 e⁻/s (±9%); `SAP_FLUX` 2708–7501 e⁻/s (aperture 2–7 px) across quarters; stitched `SAP_FLUX` likewise rescaled to `sap_mean_flux`. **Pitfall:** stitched FLUX is populated even in quarters whose per-quarter PSF FLUX is all NaN (`psf_avail` 0), and the rescaling hides quarter-to-quarter PSF-flux inconsistency of failed extractions (§1.6). Treat the stitched table as a convenience view; search and vet from `LIGHTCURVE_Qn`.

**Quality**: `SAP_QUALITY` is the Kepler TPF QUALITY word. For the target, 3,659/50,182 cadences are non-zero with bits 16 (reaction-wheel zero crossing), 128 (cosmic ray in aperture), 1024 (sudden sensitivity dropout), 2048 (impulsive outlier), 8192 (cosmic ray in collateral), 131072/262144 (rolling band in aperture/full mask). lightkurve's `DEFAULT_BITMASK` (1130799) contains none of these: default and `'none'` keep 50,182, `'hard'` 47,971, `'hardest'` = `quality==0` = 46,523 (what Kuber used; Claytor & Tayar also used `'hardest'`).

**`PSF_FLAT_FLUX` / `FLATTEN_MASK`**: at Kuber's ephemeris the box depth in `PSF_FLAT_FLUX` is 1,696 ± 84 ppm (324 in-transit cadences) vs 2,165–2,348 ppm with local baselines on raw `FLUX` **[Kuber]** — the supplied flatten attenuates a 2.3-h, 2 ppt dip by ~25%. `FLATTEN_MASK==1` on 322 cadences in runs of 3–7 consecutive cadences, 21% coinciding with non-zero quality, only 2 inside transits: an outlier mask from the flattening, not a transit mask. Search raw `FLUX` with your own detrending (wotan biweight 0.5–1 d, per quarter).

**Header bug**: `PMRA=-0.0054019`, `PMDEC=0.0028479` labelled mas/yr are Gaia DR3 pmra=−5.40189636, pmdec=2.84787494 mas/yr ÷1000 (live Gaia query **[sandbox]**). Take proper motions from Gaia.

### 2.4 lightkurve 2.6.0 behaviour **[sandbox]**
```python
import lightkurve as lk, numpy as np, astropy.units as u
from astropy.table import Table; from astropy.io import fits; from astropy.time import Time
lc = lk.read(path)                              # KeplerLightCurve; flux = PSF FLUX; columns sap_flux, psf_flat_flux, flatten_mask, quarter, sap_quality
lc_sap  = lc.select_flux('sap_flux')            # WORKS. lk.read(path, flux_column='sap_flux') silently returns PSF flux for this product
lc      = lk.read(path, quality_bitmask='hardest')   # == quality==0 (46,523/50,182); default keeps everything
cdpp    = lc.remove_nans().estimate_cdpp().value     # 472.9 (ppm); do NOT float() the Quantity
# lk.read(path, quarter=7) -> TypeError in 2.6.0 (the KBonus notebook shows it working in an older release). Per-quarter HDUs:
t = Table.read(path, hdu='LIGHTCURVE_Q7'); h = fits.getheader(path, 'LIGHTCURVE_Q7')
lq = lk.LightCurve(time=Time(np.asarray(t['TIME']), format='bkjd', scale='tdb'), flux=np.asarray(t['FLUX']),
                   flux_err=np.asarray(t['FLUX_ERR']), quality=np.asarray(t['SAP_QUALITY']))
print(h['PSFFRAC'], h['CROWDSAP'], h['FLFRCSAP'], h['NPIXSAP'], h['CHANNEL'], h['SEASON'], h['ITERNEG'])
ap = fits.getdata(path, 'APERTURE_Q7')          # uint8 mask in TPF shape
```
`lc.meta['AUTHOR']` is 'Kepler' (reader default), `lc.meta['TPFORG']` = originating KIC.

### 2.5 Working code: build the target pool correctly (**corrected**)
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
A = np.array([list(s) for s in src.psf_avail]) == '1'; F = np.array([list(s) for s in src.psf_fraction_flag]) == '1'
src['nq_usable'] = (A & F).sum(1)                       # quarters with populated PSF flux AND PSFFRAC>0.5 (the stitched ones)
kep = set(pd.read_csv("keplerstellar_kepids.csv").kepid)  # TAP: select distinct kepid from keplerstellar  (200,038 ids)
tp  = set(src.tpforg.unique())                           # 194,770 TPFs processed
src['is_kepler_target'] = src.kic.isin(kep) | src.kic.isin(tp)   # 200,318 targets / 406,582 background (paper: 200,352 / 406,548)
src['bp_rp'] = src.phot_bp_mean_mag - src.phot_rp_mean_mag; src.loc[src.bp_rp.abs() > 20, 'bp_rp'] = np.nan
# flux-vs-G consistency (zero point from Kepler targets 12<G<16, per BP-RP bin; see §1.6)
kt = src[src.is_kepler_target & src.phot_g_mean_mag.between(12,16) & src.psf_mean_flux.gt(0) & src.bp_rp.notna()]
lr0 = np.log10(kt.psf_mean_flux) + 0.4*kt.phot_g_mean_mag
bins = [-1,0.6,0.9,1.2,1.5,1.8,2.1,2.5,3.0,5]; zp = lr0.groupby(pd.cut(kt.bp_rp,bins), observed=True).median()
cent = [(a+b)/2 for a,b in zip(bins[:-1],bins[1:])]
src['flux_dex'] = np.log10(src.psf_mean_flux) + 0.4*src.phot_g_mean_mag - np.interp(src.bp_rp.fillna(1.3), cent, zp.values)
pool = src[(~src.is_kepler_target) & (src.psffrac>0.5) & (src.crowdsap>0.7) & (src.nq_usable>=8) & (src.flux_dex.abs()<=0.3)]  # 15,473
ids = pool.gaia_designation.str.split().str[-1].astype(np.int64).values
job = Gaia.launch_job_async("""SELECT u.source_id, g.parallax, g.parallax_error, g.ruwe, g.bp_rp, g.pmra, g.pmdec,
        g.non_single_star, g.ipd_frac_multi_peak, g.teff_gspphot, g.logg_gspphot
        FROM tap_upload.up AS u JOIN gaiadr3.gaia_source AS g ON g.source_id = u.source_id""",
        upload_resource=Table({'source_id': ids}), upload_table_name="up")   # chunk to <~10k rows if the server returns HTTP 500
g = job.get_results().to_pandas(); pool = pool.assign(sid=ids).merge(g, left_on='sid', right_on='source_id')
pool['MG'] = pool.phot_g_mean_mag_x + 5*np.log10(pool.parallax/100); pool['dist_pc'] = 1000/pool.parallax
# Kuber's chain (reproduces 29,800 -> 9,204 -> 6,477 -> 2,495 -> 1,836 -> 1,800; target at rank 23):
m = src.merge(ms[['gaia_designation']], on='gaia_designation')
elig = m[(m.psffrac>0.5)&(m.crowdsap>0.7)&(m.psf_avail.str.count('1')>=8)&(m.kic!=m.tpforg)&(~m.is_kepler_target)].sort_values(['phot_g_mean_mag','gaia_designation'])
```

---

## 3. Prior searches — how unsearched is it really?

| Work | Data | What was searched | Relevance |
|---|---|---|---|
| Martínez-Palomera+2023 (KBonus paper) | KBONUS-BKG | No blind search; two KOI deblending demos + dilution study | Authors invite transit searches |
| **Moquin et al. 2026, Sagan Summer Workshop poster 41** (NExScI/Caltech, 20–24 July 2026; authors incl. Hedges, Martínez-Palomera, Claytor, Quintana, Rowe, Kostov, Ciardi) | KBONUS-BKG light curves | **TraSHDUMP** matched-filter TCE pipeline ("motivated by transit-search concepts developed for the Kepler pipeline", MES threshold 7.1σ, weak-secondary test), **RECYCLEBin** (Robovetter-adapted flux vetting: odd/even, transit-vs-sinusoid), **INCInERATOr** (pixel-level depth map + PRF fit), kiauhoku/SPOTS stellar parameters, pytfit5 transit fits, TRICERATOPS-PLUS, AO imaging | First "Bonus Object of Interest" **BOI-1.01**: G=17.8 M dwarf (0.51 R☉, 0.55 M☉, <1 Gyr), P=3.37 d, ~2 h, 2.94±0.15 R⊕; "expected yield of several hundred new planet candidates"; "triple the number of planet-search stars", "M dwarfs ×7". **The dataset's own team is running the systematic search.** Not refereed; no candidate list found (2026-10-17). |
| Martínez-Palomera+2022 (APEXBA, AJ, 10.3847/1538-3881/ac4331) | 9,327 EXBA-mask sources (different pixels) | 69 EBs, 1 planet candidate (0.81 R_J) | Background planet candidates are findable |
| Forró & Szabó 2020; Bienias+2021 ([arXiv:2107.04488](https://arxiv.org/abs/2107.04488)); Forró+2022 | Raw background pixels of Kepler target apertures | 776 and 547 new (mostly short-period) EBs; RR Lyrae | Pixel-level, not KBonus LCs |
| **Bienias & Szabó 2025**, A&A, [arXiv:2501.09152](https://arxiv.org/abs/2501.09152) | Q4 pixels of the KIC target apertures ("about 150,000 pre-selected targets"), then Q2–Q16 for candidates | Lomb–Scargle periodograms per pixel + "≥3 dips > 3.5× rms" search + PDM; Gaia host search within 12″ (3 px); `Vetting` centroid tests | **7 hot-Jupiter candidates** (0.8878–1.5174 R_J, 2.5089–4.7918 d; hosts avg G 15.96, 5–23″ from the Kepler target); floors 0.0025 (Q4) / 0.0018 (Q2–Q16) fractional depth; "planets with radius <0.5 R_Jup are almost entirely outside the detectable range"; two hosts (KIC 4543171, KIC 9040849) cross-checked against KBonus curves; psfmachine produced nothing for KIC 9674230. **Sub-Neptunes around M dwarfs are below their floor.** |
| **Wang & Esamdin 2025**, AJ, DOI 10.3847/1538-3881/ade718 | Kepler TPF pixels of KOI-1755 | Pixel-level flux modelling + centroid shifts + Gaia DR3 cross-match | The 25-d KOI-1755 signal is a 6.14-d background EB — precedent for the KOI-FP-neighbour pool (§4.3) |
| Barber+2022 (THYME VIII, AJ 164, 88, [arXiv:2206.08383](https://arxiv.org/abs/2206.08383)) | Kepler PDCSAP (long cadence) of MELANGE-3 candidates + TESS | Notch pipeline | CDS table5.dat (1,006 rows) line 546: `2103792461826460288 285.90783 40.94199 16.282 0.795 M1.6 N Y 3.406 0.04`; ReadMe: identifier is Gaia **EDR3**, flags TESS?/Kepler? with Note (2) "Of sufficient quality for running the Notch transit-search pipeline". KIC 5689121 was never a Kepler target, so the Y must come from a coordinate-matched Kepler light curve, plausibly the 7″ neighbour KIC 5689116 (inference). That PDCSAP shows only 98 ± 32 ppm at Kuber's ephemeris (SAP 13 ± 53 ppm; `v10_pdcsap.py`), and keplerstellar lists `ntce=0, nkoi=0` for 5689116 — so Barber's search did not test this signal. |
| Claytor & Tayar 2025; Choi, Espinoza-Rojas & Hekker 2026 ([arXiv:2604.24476](https://arxiv.org/abs/2604.24476), red-giant asteroseismic binaries) | KBONUS-BKG | Rotation / oscillations | Not transit searches; systematics notes |
| Kuber Mehta 2026 **[Kuber]** | KBONUS M-star catalogue, 1,736 of 1,800 eligible | BLS 0.5–25 d on per-quarter `LIGHTCURVE_Q.FLUX` (quality==0, finite err>0, PSFFRAC≥0.5), holdout quarters, inverted-flux nulls | Real first-screen passes 31/29 vs inverted 32/28; one candidate (KIC 5689121); unresolved leads KIC 6629201 (2.638 d), 12164160 (19.42 d), 8743279 (21.96 d), 8803342 (5.010 d) **[extension README]**; two inverted stars (KIC 8701308, 5283912) exceed the candidate's SNR **[verification.json]** |

No KOI/TCE/TOI/CTOI has ever been issued from KBONUS-BKG; no Zooniverse project uses it (web search 2026-10-17). **Conclusion: no refereed systematic transit search of the KBONUS-BKG light curves exists, but the KBonus team's TraSHDUMP pipeline is active with a first candidate and a multi-hundred expected yield; the pixel-level search by Bienias & Szabó reaches only giant planets.**

### 3.2 Accounting of what remains **[sandbox]**
Source catalogue 606,900. Kepler targets by `kic ∈ keplerstellar ∪ {tpforg}`: 200,318 (paper 200,352); background 406,582 (paper 406,548). keplerstellar alone gives 194,867 (misses 5,451 observed targets absent from the stellar table); `kic==tpforg` gives 181,681; 18,583 Kepler targets have `kic≠tpforg` (extracted from a neighbour's TPF). 341,081 background rows have a KIC number without ever having been Kepler targets. **Never use `kic != tpforg` or "has a KIC" as the background/searched criterion.**

Background sources passing PSFFRAC>0.5 & CROWDSAP>0.7 & ≥8 usable quarters: **23,157** (G<17: 6,366); 7,091 of them have ≥1 quarter with populated stitched FLUX but NaN per-quarter PSF flux; **15,473 also pass the flux-consistency test** (G<17: 5,466; G<16: 2,000). Gaia DR3 join of the 23k pool: 10,278 with plx/err>5; **8,463 dwarfs (M_G>4)**, 1,815 giants; dwarfs by G: 12–13: 2; 13–14: 16; 14–15: 143; 15–16: 668; 16–17: 2,806; 17–18: 3,215; 18–19: 1,613; within 100/200/300/500/1000 pc: 8/41/140/498/1,916; G<17: 3,635.

M-star catalogue (29,800): 2,040 have `kic==tpforg`, 2,345 are Kepler targets by the union definition; median G 18.47 (quartiles 17.95/18.77); Teff_gspphot 2,616–3,900 K; RUWE<1.4 for 29,105; `non_single_star>0` for 34; 182 within 100 pc (6 eligible non-targets). Kuber's cuts: psffrac>0.5 → 9,204; crowdsap>0.7 → 6,477; ≥8 quarters → 2,495; kic≠tpforg → 1,836; not a Kepler target → 1,800 (reproduced exactly; same with the union definition). Eligible G 14.19–19.00 (3 at G<15, 5 at 15–16, 67 at 16–17, 335 at 17–18, 1,390 at 18–19); 2/30/100 within 100/200/300 pc; median BP−RP 2.08; **801 fail the flux-consistency test**. Not searched: 20,596 with psffrac≤0.5, 2,727 with crowdsap≤0.7, 3,982 with <8 quarters, 695 Kepler-target M dwarfs, plus 64 skipped/timeouts.

---

## 4. Where to hunt

1. **Nearby, flux-consistent non-target dwarfs of any colour** (the sweet spot): `~is_kepler_target & psffrac>0.7 & crowdsap>0.8 & nq_usable>=12 & abs(flux_dex)<=0.3 & parallax/parallax_error>5 & MG>4 & G<17`. Expect 6-h CDPP 80–350 ppm; a 2 R⊕ planet on a 0.5 R☉ star (1,300 ppm) is ≥10σ per 6-h bin at G≈16; a 1.5 R⊕ planet on a 0.3 R☉ star (2,100 ppm) is detectable to G≈17.5. Starting tables: `verify-kbonus-dataset/background_pool_fluxconsistent.csv` (15,473 rows with host G and separation) and `kbonus-dataset/background_quality_pool_gaia.csv` (23,171 rows with Gaia columns). Nearest M dwarfs in the M-star catalogue: Gaia DR3 2073565650248046848 (70 pc, G 14.3, but crowdsap 0.17) and Gaia DR3 2105301782052654592 (77 pc, G 14.19, psffrac 0.60, crowdsap 0.77, 18 quarters, TPF host KIC 7940546, Kuber's rank 1).
2. **Below TraSHDUMP's threshold.** The KBonus team's pipeline uses a Kepler-like MES ≥ 7.1σ. Multi-quarter TLS/BLS with injection-calibrated thresholds, per-season consistency and inverted-null calibration can legitimately go deeper on the flux-consistent pool — but every such candidate needs stronger pixel/statistical validation, not weaker.
3. **KOI false-positive neighbours.** `kois-nns` (35,639 LCs, 8,073 KOIs, every Gaia source within 30″). 1,889 KOIs are `FALSE POSITIVE` with `koi_fpflag_co=1` (of 4,839 FPs) **[sandbox TAP]**; their eclipse is on *some* neighbour, and re-attributing it to a fainter, smaller star rescales the depth by 10^(0.4Δm) — the KBonus KOI-770 demo and Wang & Esamdin's KOI-1755 show the method works. Query: `select kepid,kepoi_name,koi_period,koi_time0bk,koi_depth,ra,dec from cumulative where koi_disposition='FALSE POSITIVE' and koi_fpflag_co=1`.
4. **Partial-PSF M dwarfs** (20,596 with psffrac≤0.5): use `LIGHTCURVE_Qn` with PSFFRAC 0.3–0.5, ITERNEG=0, psf_avail=1, and demand season-constant PSF depth.
5. **Re-extracted Kepler targets** (200k PSF curves): deblended depths change radii (paper §V); low value for new planets, cheap consistency check.

Independent Gaia pre-selection (ADQL, async; Kepler footprint RA 279.6–301.8°, Dec 36.5–52.5° matches the catalogue):
```sql
SELECT source_id, ra, dec, parallax, parallax_error, ruwe, phot_g_mean_mag, bp_rp, teff_gspphot, logg_gspphot
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), POLYGON('ICRS',279.5,36.5, 302.5,36.5, 302.5,52.5, 279.5,52.5))
  AND phot_g_mean_mag < 19 AND parallax > 3.3 AND parallax_over_error > 5
  AND phot_g_mean_mag + 5*log10(parallax/100) > 4
```
Ranking: score = expected SNR for a fiducial 2 R⊕ planet = (depth(R*) / CDPP6(G, colour)) × sqrt(N_quarters × 90 d / P × T14/6 h), with R* from Gaia M_G–(BP−RP), penalties for crowdsap<0.8, RUWE>1.4, `ipd_frac_multi_peak>5`, |flux_dex|>0.2, host contrast >3 mag within 8″, and a TPF host that is a known EB/KOI. Brightest-first (Kuber) biases toward distant reddened giants and host-contaminated faint stars.

---

## 5. Pitfalls specific to background PSF photometry

### 5.1 Aperture (SAP) vs PSF flux
SAP includes neighbour flux; `CROWDSAP`/`FLFRCSAP` change every season (aperture 2–7 px for Kuber's target). A real on-target transit has **constant PSF depth** and a **season-dependent SAP depth ≈ depth × FLFRCSAP × F_PSF / F_SAP** (Kuber: χ² 15.4/12 for this model, p=0.005 for a constant SAP depth **[Kuber SCIENCE.md]**). If the SAP depth is *larger* in the most contaminated seasons, the signal is on the contaminant. PSF and SAP are the same photons, not independent evidence.

### 5.2 The host's signal leaking in — mandatory ephemeris matching and flux-consistency check
Every background source lives in a TPF centred on a Kepler target (`TPFORG`); residual correlation is common (Claytor & Tayar: up to 63% of periodic pairs) and a third of the quality pool carries >2× excess flux from the host (§1.6). Before anything else: (i) check `flux_dex` and the host contrast/separation; (ii) match period/epoch (and 1/3, 1/2, 1, 2, 3 harmonics, timing tolerance grown over the baseline) against the TPF host's own KBONUS PSF and official PDCSAP curves, every KIC within ~60″, and every Gaia source in the TPF footprint. Working NASA Exoplanet Archive TAP queries **[sandbox, 2026-10-17]**:
```sql
-- https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=...   (Oracle syntax: no AS on table aliases)
select kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_time0bk,koi_depth,koi_duration,koi_fpflag_co,ra,dec
  from cumulative where kepid in (5689116,5689121)                                    -- (empty here)
select kepid,tce_plnt_num,tce_period,tce_time0bk,tce_depth,tce_duration,tce_max_mult_ev
  from q1_q17_dr25_tce where kepid in (5689116,5689121)                               -- (empty; table has 34,032 rows)
select kepid,kepoi_name,koi_disposition,koi_period,koi_time0bk,koi_depth,ra,dec from cumulative
  where ra between 285.7975 and 286.0181 and dec between 40.8587 and 41.0253          -- Kepler-169 b-f + KOI-3784.01 (FP, 23.87 d)
select t.kepid,t.tce_plnt_num,t.tce_period,t.tce_time0bk,t.tce_depth,t.tce_duration,s.kepmag
  from q1_q17_dr25_tce t join keplerstellar s on s.kepid=t.kepid
  where s.ra between 285.7975 and 286.0181 and s.dec between 40.8587 and 41.0253 and s.st_delivname='q1_q17_dr25_stellar'
-- keplerstellar has 5 delivery rows per star: filter st_delivname (or DISTINCT). Spatial ADQL is unreliable here:
-- CONTAINS(POINT,CIRCLE) on cumulative -> "ORA-00904: 'HTM20': invalid identifier"; on ps -> Cloudflare block. Use RA/Dec boxes.
```
Kepler EBs: the Villanova endpoints return 400/500 (`/results/?q=&format=csv`, `/results/?kic=…&format=csv`); use VizieR `J/AJ/151/68/catalog` (Kirk+2016, 2,876 EBs) via `Vizier(columns=['KIC','Per','BJD0','Morph','Kpmag','_r']).query_region(coord, radius=5*u.arcmin, catalog='J/AJ/151/68/catalog')` **[sandbox: none within 5′ of Kuber's star, 4 within 30′]**. Also check `kois-nns` (`nn_koi`, `nn_dist`) and Gaia `non_single_star`, `ruwe`, `ipd_frac_multi_peak`. For Kuber's star the 5′ box returns Kepler-169 b–f (KIC 5689351; e: 13.76710058 d, T0 174.81089 BKJD, 684.8 ppm) and KOI-3784.01; Kepler-169 e differs from the candidate's 13.769448 d by 32σ and the ephemerides never overlap **[Kuber]**.

### 5.3 Quarter-to-quarter offsets, seasons and module 3
Kepler rolled 90° each quarter, so a star repeats on the same channel every fourth quarter (seasons 0–3). For KIC 5689121: season 0 = module 11 ch 34 (Q2, 6, 10, 14), season 1 = module 23 ch 78 (Q3, 7, 11, 15), season 2 = module 15 ch 50 (Q4, 8, 12, 16), season 3 = **module 3 ch 6 (Q1 only)**. "All 4 outputs of Module 3 failed at 17:52 UTC Jan 9, 2010 during LC CIN 12935 … 20% of the FOV suffered a one-quarter data outage every year" (Kepler Data Characteristics Handbook KSCI-19040-005 §4.6); channels 5–8, 1,022 of 4,397 Q4 cadences, ~19% of targets observed three of four quarters (Batalha et al. 2013, [arXiv:1202.5852](https://arxiv.org/abs/1202.5852)). Hence QDETECT `011110111011101110` (Q5, 9, 13, 17 absent). Consequences: never let a BLS fold mistake season gaps for transit windows (aliases near n×90/4 d); do odd/even, depth and shape consistency **per season**; split by `SEASON`/`CHANNEL` from the `LIGHTCURVE_Qn` headers; check per-quarter PSF-flux medians for consistency (a >30% spread after PSFFRAC correction flags a failed extraction).

### 5.4 Centroid columns are useless for localisation
`CENTROID_COLUMN/ROW` = Gaia position + a scene-wide jitter vector: target and neighbour have identical per-quarter centroid standard deviations (Q1 (0.0116, 0.0256) px … Q16 (0.0086, 0.0683) px for both) and differ only by their fixed Gaia offset **[sandbox]**. psfmachine `get_source_centroids`: `'poscor'` = `tpf.wcs.all_world2pix(ra,dec)` + the TPF `pos_corr1/pos_corr2` per cadence; `'scene'` = Gaia ra/dec + the common scene centroid; `'aperture'` (moments) is "ill-defined" for partial sources (docstring, `master` branch). Kuber's `kcentroid.py` found χ²(zero shift)=25.5/26 **[Kuber]**. Do not use them for vetting.

### 5.5 Pixel-level vetting with the originating TPF (do this instead)
Download the TPF of `TPFORG` (`lk.search_targetpixelfile("KIC 5689116", author="Kepler", cadence="long").download_all()`, 13 files), fold at the candidate ephemeris, build per-event difference images (out-of-transit median − in-transit mean), compare the difference-image centroid / PRF fit with the Gaia positions via the TPF WCS (`WCS(tpf.hdu[2].header)`). Result **[sandbox, re-run]**: in 13/13 quarters the difference-image centroid lies 0.54 px (median) from KIC 5689121 vs 1.30 px from KIC 5689116 (1 px = 3.98″) and the peak-dip pixel is the target's pixel every quarter — consistent with Kuber's PRF-fit (Δχ² = 934) and flux-weighted-centroid (χ² 43.9/24) tests **[Kuber]** and with what INCInERATOr does for BOI-1. Minimal code:
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
Upgrade: fit `lightkurve.prf.KeplerPRF(channel, shape, column, row)` to the difference image at each Gaia position and compare χ²; fold the TPF host's own KBONUS PSF curve (Kuber: 24 ± 41 ppm) and its official PDCSAP (98 ± 32 ppm here — the expected small leak).

### 5.6 Other traps
* `psf_avail.count('1')` includes flag-`2` quarters; `psf_fraction_flag.count('1')` includes quarters whose PSF flux is NaN. Use the AND (§1.3).
* Flag-`2` quarters have HDUs and stitched rows (NaN FLUX); don't assume HDU presence means usable.
* Stitched FLUX is rescaled to `psf_mean_flux` and filled even where the per-quarter PSF is NaN; its provenance there is undocumented — search from `LIGHTCURVE_Qn`.
* ITERNEG=1 quarters have no `PSF_FLUX_NOVA`; faint sources near a ≥2 mag brighter star within ~15″ are the ones most likely over-/under-subtracted; a high-contrast pair can show an *inverted* copy of the neighbour's eclipse — inverted-flux nulls catch this only partially; always fold the host.
* G<10 stars were removed and pixels near G≤8 stars masked; bright-star fields are absent.
* MAST page errors: "~20,000" M dwarfs (29,800), "~32,000" kois-nns files (35,639), "29 MB" bulk script (168 MB), "134 GB" M-star bundle (107 GB), per-quarter column names. Quote the files.
* Catalogue `ra/dec` are J2016; propagate proper motion (from Gaia, not the header) for high-PM stars before pixel work.
* `FLUX_ERR` is formal; recompute noise empirically.
* Gaia upload-joins of >~20k rows can return HTTP 500 — chunk to ~8–10k.

---

## 6. Sources
* MAST HLSP page https://archive.stsci.edu/hlsp/kbonus-bkg ; README https://archive.stsci.edu/hlsps/kbonus-bkg/hlsp_kbonus-bkg_kepler_kepler_all_kepler_v1.0_readme.md ; directory index https://archive.stsci.edu/hlsps/kbonus-bkg/ (bundles/, catalogs/, lcs/, bulk scripts).
* Paper https://arxiv.org/abs/2310.17733 ; DOI 10.3847/1538-3881/ad0727 ; Crossref API (13 citations); OpenAlex W4388973240 (12).
* Moquin et al. 2026 poster https://vmnexsci.ipac.caltech.edu/workshop/2026/posters/Poster_JulieMoquin_41.pdf (2026 Sagan Summer Workshop, https://vmnexsci.ipac.caltech.edu/workshop/2026/).
* psfmachine https://github.com/SSDataLab/psfmachine (branch master, `src/psfmachine/tpf.py::get_source_centroids`), PyPI JSON (1.1.4); KBonus repo https://github.com/jorgemarpa/KBonus (`notebooks/light_curve_caveats.ipynb`).
* APEXBA https://archive.stsci.edu/hlsp/kbonus-apexba ; https://arxiv.org/abs/2112.06683.
* Bienias & Szabó 2025 https://arxiv.org/abs/2501.09152 ; Bienias+2021 https://arxiv.org/abs/2107.04488 ; Wang & Esamdin 2025 DOI 10.3847/1538-3881/ade718 (press: https://english.cas.cn/newsroom/research_news/phys/202509/t20250904_1053933.shtml).
* Barber+2022 https://arxiv.org/abs/2206.08383 ; CDS https://cdsarc.cds.unistra.fr/ftp/J/AJ/164/88/ReadMe and table5.dat.
* Claytor & Tayar 2025 https://arxiv.org/abs/2506.03248 ; https://zenodo.org/records/17245793 ; Choi+2026 https://arxiv.org/abs/2604.24476 ; Bentley+2018 https://arxiv.org/abs/1804.06982.
* Kepler module 3: Data Characteristics Handbook KSCI-19040-005 https://archive.stsci.edu/kepler/manuals/Data_Characteristics.pdf §4.6 ; Batalha+2013 https://arxiv.org/abs/1202.5852.
* NASA Exoplanet Archive TAP https://exoplanetarchive.ipac.caltech.edu/TAP/sync (`cumulative`, `q1_q17_dr25_tce`, `keplerstellar`); VizieR J/AJ/151/68 (Kirk+2016); Gaia archive (astroquery.gaia).
* Kuber's repo (read-only): `analysis/editorial-audit-20261013/{README.md,verification.json,verify.py}`, `analysis/kbonus-extension-20261011/README.md`, `analysis/validation-20261011/SCIENCE.md`, `evidence/codex/database_expansion/kbonus/pilot.py`.

## 7. Verification scratch files
Kept outside this repository (not published): `v1_fits.py` (FITS layout/headers), `v2_lk.py` + `v5_cdpp.py` (lightkurve behaviour, bitmasks, estimate_cdpp), `v3_cat.py` (columns, counts, Kuber chain), `v4_flags.py` (flag cross-tabs, target definitions), `v6_gaia.py` (PM bug), `v7_qdetect.py` + `v8_flag2.py` (QDETECT/HDU semantics, `dl/…kic-003459079…lc.fits`), `v9_misc.py` (VizieR, MAST), `v10_pdcsap.py` (host PDCSAP depth), `nq_counts.parquet`, `logratio_colour.parquet`, `background_pool_fluxconsistent.csv` (15,473 rows), `pdfs/{moquin,batalha,datachar}.txt`, `claytor.txt`, `bienias.txt`, `barber_readme.txt`, `barber_table5.dat`, `psfmachine_tpf.py`, `mast_page.html`, `oa_cites.json`.