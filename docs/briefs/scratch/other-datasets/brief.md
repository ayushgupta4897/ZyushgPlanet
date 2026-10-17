# Other archival datasets: where an unreported transiting-planet candidate can still be found in October 2026, and how thoroughly each has been searched

Research brief for the ZyushgPlanet playbook. Dimension: **other-datasets**. Date of all live checks: **2026-10-17** (UTC). All counts labelled "(TAP today)", "(ExoFOP today)", "(MAST today)" or "(sandbox)" were produced by commands run in the sandbox on that date; the scripts and raw outputs are in `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/other-datasets/`.

---

## 0. Headline findings (read this first)

1. **"Searched" has three different meanings and the gaps are mostly in the second and third.** (a) *Detection* searches (TPS/BLS run on a light curve), (b) *vetting* into a public candidate list (TOI/KOI/CTOI), and (c) *novelty* catalogues you must cross-match. For TESS the detection step is near-complete to T = 13.5 (QLP, TESS-SPOC, SPOC) and the vetting step is badly bottlenecked; for faint stars (T > 13.5) even the detection step is patchy. For Kepler the primary-target search is saturated, but the Kepler Bonus background sources are a largely unsearched 400k-star pool that one Goddard team has just started mining.
2. **TESS right now (2026-10-17):** TESS is observing Sector 109 (Cycle 9, Extended Mission 3, 200-s FFIs, ~13,000 2-min and ~2,500 20-s targets). SPOC 2-min light curves and single-sector DV/TCE products are public through **Sector 107** (released 2026-10-01, 24 days after the sector ended). **No TOIs have yet been alerted from Sectors 103-107**; the newest single-sector SPOC TOIs are from Sector 102 (alerted 2026-08-06) and the newest QLP faint-star TOIs are from Sector 98 (alerted 2026-07-24). Measured lags: data public 24-39 days after sector end; SPOC single-sector TOIs 94-148 days; QLP faint-search TOIs 200-338 days. QLP Sectors 99-104 carry a known barycentric-timing error and have not been redelivered. TESS-SPOC FFI light curves exist only through Sector 85 (for southern CVZ stars, lightkurve returns TESS-SPOC only through Sector 69). TGLC on MAST stops at Sector 55; GSFC-eleanor-lite at Sector 26; T16 covers Cycle 1 only; the new **TARS** HLSP (June 2026) gives T < 16, d < 500 pc light curves for 7.5 million stars through Sector 96.
3. **The under-searched TESS regime, quantified:** the TIC holds **100,567 dwarfs with Teff < 3900 K within 100 pc and 13.5 < T < 16.5**; only **33** are TOI hosts and **43** are CTOI hosts (TIC criteria query + ExoFOP cross-match, sandbox). Only **457 of 8,151 TOIs** and **351 of 5,138 CTOIs** have T > 13.5; only 35 TOIs and 44 CTOIs have T > 15. eumemic's full-sky sweep stopped at T <= 13.5 (4,166,638 FGKM dwarfs, Sectors 1-104) plus M-dwarf slices to T <= 13.5; Chopra's covered 526 GCNS stars. The faint nearby M dwarfs are therefore the biggest TESS pool that is (i) largely unvetted, (ii) ground-followable (d < 100 pc, Ks ~ 10-12), and (iii) has deep transits for small planets.
4. **Kepler Bonus beyond Kuber's slice:** I reproduced his funnel exactly from the v1.0 catalogues (29,800 -> 9,204 -> 6,477 -> 2,495 -> 1,836). Applying the *same* quality cuts to the **non-M** background sources leaves **25,542 stars** (6,151 with G < 16; 10,417 with G < 17) that nobody has published a transit search of. Caveat: a NASA Goddard/Ames team (Moquin, Wilson, Quintana, Hedges, Martinez-Palomera, Kostov, Ciardi et al.) presented "BOI-1.01" and the TraSHDUMP matched-filter pipeline at the 2026 NExScI workshop, explicitly expecting "several hundred new planet candidates" from Kepler Bonus. You would be racing a funded team with the product's own authors on it.
5. **K2 is more searched than it looks, but has structural holes:** 361,770 unique EPIC targets across C0-C19; Kruse 2019 (C0-8, Kp 9-15, 818 candidates) and Zink 2021 (C1-8, C10-18, 747 candidates, maximum completeness 61%) are the systematic sweeps. Holes: Campaign 0 and 19 (excluded by Zink), Kp > 15, periods > ~20-40 d (80-day campaigns), and *no* Kepler-Bonus-style PSF-deblended background product exists for K2.
6. **Nothing else is competitive for a solo developer in 2026:** Gaia DR4 (2 Dec 2026) will release epoch photometry for ~2.8 billion sources but with ~70 sparse epochs per star over 66 months, and DPAC publishes its own `vari_planetary_transit` table (1,007 rows planned); ground surveys (ZTF 8-25 mmag repeatability, ASAS-SN to V ~ 17, ATLAS 2-day cadence) only reach hot Jupiters; CHEOPS is targeted (1-yr proprietary); PLATO launches March 2027; Roman launched 30 Aug 2026 but its bulge survey has not begun.
7. **Ranked recommendation (Section 10):** (1) faint (13.5 < T < 16) M dwarfs within 100 pc, all sectors, own FFI photometry + TARS/TGLC where available, with the fresh-sector TOI-lag window as a bonus; (2) the fresh-sector window itself for T <= 13.5 stars not in eumemic's dwarf cuts (giants excluded, but also Teff > 6500 K, R > 1.5 Rsun); (3) the 25,542 non-M Kepler Bonus background stars (novel, cheap, but un-validatable at G > 16 and contested); (4) K2 C0/C19/faint/long-period residue; (5) everything else only as a cross-check source.

---

## 1. Vocabulary: three kinds of "searched"

| Layer | What it means | Where the gaps are |
|---|---|---|
| **Detection** | A transit-search algorithm (TPS, BLS, TLS, QATS, ML) was run on a light curve of this star. | TESS: complete to T <= 13.5 via QLP (all sectors to 104) and SPOC/TESS-SPOC; patchy for T > 13.5 (TGLC S1-55, eleanor-lite S1-26, T16 Cycle 1, TARS S1-96 but no transit search published). Kepler: complete for 198,709 primary targets; Kepler Bonus background: only the M-dwarf slice (Kuber) and an in-progress Goddard search. K2: Kruse/Zink cover Kp <~ 15, C1-18. |
| **Vetting/publication** | The signal was examined and put on a public list (TOI, CTOI, KOI, published candidate table). | TESS: QLP TCEs below T = 10.5 were "not considered for vetting" in the prime mission (Guerrero 2021, Sec. 10); the faint-star search (Kunimoto 2022) extended this to 13.5; TESS-SPOC FFI TCEs feed TOIs only slowly (`spoc-s14-s86-b02` TOIs with T = 13.6-14.8 alerted Aug 2026). Measured lag 94-338 days. |
| **Novelty cross-match** | Catalogues you must check before claiming "unreported". | Listed in Section 8 with counts. |

The practical consequence: a signal can be present as a SPOC TCE in a DV report (detected) and still be "unreported" in the TOI sense. eumemic exploited exactly this: they "re-measured the 14,112 TCEs on 8,698 FGKM dwarfs that are not TOIs or CTOIs" (eumemic README). Kuber's star was never a target of any pipeline, so it was unreported at *all three* layers.

---

## 2. Kepler prime mission (2009-2013)

### 2.1 What exists
- **Targets and search:** The final pipeline run (DR25) searched **198,709 stellar targets** and produced **34,032 TCEs** (Twicken et al. 2016, arXiv:1604.06140). Thompson et al. 2018 (ApJS 235, 38; arXiv:1710.06758) catalogued **8,054 KOIs, 4,034 planet candidates** with periods 0.25-632 d; Robovetter completeness > 85% and reliability > 98% for P < 100 d. The cumulative KOI table today holds **9,564 rows: 2,748 CONFIRMED, 1,977 CANDIDATE, 4,839 FALSE POSITIVE** (TAP today: `select koi_disposition, count(*) from cumulative`).
- **Cadence:** 29.4-min long cadence for all targets; ~512 short-cadence (58.8 s) slots.
- **Full-frame images:** calibration FFIs "collected approximately monthly" (Montet et al. 2017, arXiv:1705.07928), i.e. a few dozen frames over the mission. Useless for transit detection; useful for FFI-level variability of all stars in the field.
- **Superstamps:** Kepler read out 200x200-pixel long-cadence superstamps of NGC 6791 and NGC 6819 (13.3 arcmin on a side) for the whole mission (Kuehn et al. 2013, arXiv:1309.6044). Published uses: variable stars (Kuehn 2013/2015) and 149 red giants via difference imaging (Espinoza-Rojas et al. 2023/2024, arXiv:2310.06757). **No published transit search of the superstamps was found** (WebSearch, 2026-10-17). These clusters are at ~2.5-4 kpc; member dwarfs are G >~ 17-19, so any candidate would be essentially un-validatable.

### 2.2 Independent searches already done (what they removed from the table)
| Search | Scope | Result |
|---|---|---|
| Huang, Bakos & Hartman 2012 (arXiv:1205.6492) | Q1-Q6, BLS + visual | 150 new periodic candidates, 7 single-transit events |
| Planet Hunters (Zooniverse); Wang et al. 2015 "PH VIII" (arXiv:1512.02559, ApJ 815, 127) | long-period, visual | 41 long-period candidates in 38 systems (17 single-transit, 14 two-transit) |
| Foreman-Mackey et al. 2016 (arXiv:1602.07848) | 39,036 bright GK dwarfs, automated long-period (1-2 transits) | 16 long-period candidates |
| Uehara et al. 2016 | 7,557 KOIs by eye | 28 single-transit events in 24 KOIs |
| Shallue & Vanderburg 2018 (AJ 155, 94) | CNN on weak TCEs around 670 known multi-planet systems | Kepler-90i, Kepler-80g |
| Herman, Zhu & Wu 2019 (arXiv:1901.01974) | re-analysis of long-period transits | corrected FM16 occurrence; candidates revisited |
| Kawahara & Masuda 2019 (arXiv:1904.04980; AJ 157, 218) | P >~ 2 yr, 1-2 transits, pixel-level vetting, Gaia | 67 candidates, 23 new |
| Kunimoto, Matthews & Ngo 2020 (arXiv:2003.04397; AJ 159, 124) | ~200,000 stars, full Q1-Q17, >= 3 transits, independent pipeline | 17 new candidates; 98.8% recovery of confirmed planets |
| Kepler EB catalogue (Kirk et al. 2016, AJ 151, 68) | all targets | 2,878 EBs (1.3% of targets) |

**What is left on Kepler primary targets:** essentially only (i) single/double-transit long-period events below the sensitivity of FM16/Kawahara-Masuda (they used bright GK samples or KOIs), and (ii) sub-threshold multi-planet additions of the Shallue-Vanderburg type around the remaining ~2,000 multi-systems. Both are low-yield and the latter has been automated repeatedly (ExoMiner on Kepler, Valizadegan 2022). I would not spend compute here.

### 2.3 Kepler Bonus (KBONUS-BKG) - the actually open Kepler pool
Facts (MAST HLSP page https://archive.stsci.edu/hlsp/kbonus-bkg, fetched 2026-10-17; Martinez-Palomera, Hedges & Dotson 2023, AJ 166, 265, arXiv:2310.17733; DOI 10.17909/7jbr-w430; CC BY 4.0):
- "606,900 light curves, from which 406,548 correspond to background sources and 200,352 are Kepler targets"; Gaia DR3 input catalogue; Linearized Field Deblending PSF photometry fitted to Kepler FFIs; CDPP "under 100 ppm for sources brighter than G = 16". In the v1.0 source catalogue I count **425,219 rows with `kic != tpforg`** (my background definition) and 181,681 with `kic == tpforg` (sandbox; the page's 406,548/200,352 split uses a different definition).
- Products: `_lc.fits` with stitched LIGHTCURVE + per-quarter LIGHTCURVE_Q + APERTURE_Q; catalogue columns `psffrac`, `crowdsap`, `flfrcsap`, `npixsap`, `pertrati`, `pertstd`, `psf_avail` (an 18-character 0/1 string per quarter, **not** an integer - parse with `.count('1')`), `sap_avail`, `psf_fraction_flag`, `tpforg`, `kic`, `kic_sep`, `kep_mag`, `phot_g/bp/rp_mean_mag`.
- Catalogues: source catalogue (159 MB FITS, 606,900 rows), M-stars catalogue (**29,800 rows** - the HLSP page text says "~20,000"; the file says 29,800, matching Kuber's number), KOIs+neighbours (~32,000), white dwarfs (~90). Download URLs: `https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_{source,mstars,kois-nns,wd}-catalog_kepler_v1.0_cat.fits` (note the `/catalogs/` path; the bare `/hlsps/kbonus-bkg/` path 404s). Light curves: `lk.search_lightcurve(coord, author='KBONUS-BKG')` or `astroquery.mast.Observations.query_criteria(provenance_name='kbonus-bkg')`; bulk via the lftp/curl scripts (4.8 TB total; the M-star bundle is 134 GB).
- **Reproduced funnel (sandbox, 2026-10-17):** M-star catalogue 29,800 -> PSFFRAC > 0.5: 9,204 -> CROWDSAP > 0.7: 6,477 -> >= 8 PSF quarters: 2,495 -> `kic != tpforg`: **1,836**. Identical to Kuber's audit numbers, so his selection is fully reproducible from public files. G-distribution of the 1,836: 2 at 13-14, 10 at 14-15, 24 at 15-16, 68 at 16-17, 341 at 17-18, 1,391 at 18-19. I.e. the pool he searched is >90% fainter than G = 17.
- **Unsearched complement (sandbox):** same cuts on the 577,100 non-M sources -> 275,659 -> 192,698 -> 132,009 -> **25,542 background sources**; G bins <13: 320, 13-14: 937, 14-15: 1,824, 15-16: 3,070, 16-17: 4,266, 17-18: 6,712, 18-19: 8,412. Median BP-RP 1.06 (K-ish). Of these, **6,151 have G < 16** where CDPP < 100 ppm applies. (107,126 primary-target rows pass the same cuts; those were searched by DR25.) Files: `kbonus/kbonus_nonM_background_survivors.csv`, `kbonus/kbonus_M_background_survivors.csv`.
- **Caveats:** `kic_sep` median 0.23 arcsec for background survivors, i.e. most "background" sources *are* KIC entries that were simply never selected as targets; Kuber's extra "exclude primary targets" step (1,836 -> 1,800) removes sources that were targets in other quarters - apply it via the KIC/KOI target lists. Blending with the originating target is the dominant false-positive mode (his candidate sits 7.1 arcsec from a 2.1x brighter star); per-quarter aperture-depth vs. flux-fraction tests are mandatory.
- **Competition:** Moquin et al., "One Man's TraSHDUMP is Another Man's Treasure: An Exoplanet Candidate in the Kepler Bonus Light Curves" (2026 NExScI workshop poster, https://vmnexsci.ipac.caltech.edu/workshop/2026/posters/Poster_JulieMoquin_41.pdf): "triple the number of planet-search stars", "an increase in the number of M dwarfs by a factor of seven", "an expected yield of several hundred new planet candidates"; BOI-1 is a G = 17.8 M dwarf with a 3.37-d, ~2-h TCE at 7.1 sigma MES, vetted with a Robovetter-like suite ("RECYCLEBin") and TRICERATOPS-PLUS. Expect a BOI catalogue from this group.
- **Sister product:** KBONUS-APEXBA (https://archive.stsci.edu/hlsp/kbonus-apexba): 9,327 sources in the Kepler EXBA background super-apertures (Martinez-Palomera et al. 2022), with 69 EBs and one 0.81 RJ candidate already reported.

---

## 3. K2 (2014-2018)

### 3.1 What exists
- **Targets:** NASA Exoplanet Archive `k2targets` holds **392,245 target-campaign rows, 361,770 unique EPIC IDs** (TAP today). Per campaign: C0 7,748; C1 21,646; C2 13,401; C3 16,375; C4 15,781; C5 25,137; C6 27,289; C7 13,261; C8 23,565; **C9 623** (microlensing superstamp campaign); C10 27,170; C11 13,607; C12 28,088; C13 21,367; C14 29,897; C15 23,278; C16 29,413; C17 32,221; C18 20,427; E (engineering) 1,951. C19 is not a separate value in this table (the archive lumps it; MAST lists C19 separately and it was cut short by fuel exhaustion).
- **Systematics:** 6-hour thruster-firing roll corrections dominate; pointing drifts; C10 lost its first days to a pointing error; C11 split into 111/112; C9 was a 3.7 deg^2 bulge superstamp at 30-min cadence for microlensing. Zink et al. 2021 attribute their 61% maximum completeness to "the heavy detrending K2 data require".
- **Light-curve products on MAST** (K2 HLSP index https://archive.stsci.edu/k2/hlsps.html and per-product pages, fetched 2026-10-17):
  - **EVEREST** v2.0 (Luger et al. 2018; DOI 10.17909/T9501J): tar bundles for c00-c08, c102, c111, c112, c12-c18 (v2.0 for C14-18 added 2019-04-19); "comparable to the original Kepler mission for targets brighter than Kp = 13, and within a factor of 2 for fainter targets". lightkurve author `EVEREST`.
  - **K2SFF** (Vanderburg & Johnson): MAST "Campaigns 0 though 14 plus the pre-campaign 0 engineering targets"; Vanderburg's CfA page lists reductions for C0-C19 (https://lweb.cfa.harvard.edu/~avanderb/k2.html). lightkurve author `K2SFF`.
  - **K2SC** (Aigrain et al. 2016): C3-8 and C102. **POLAR** (Barros): C1-6. **K2VARCAT**: C0-4. **K2GAP**: C1 giants. **KEGS**: C5,6,8 galaxies. **HALO**: ~161 saturated bright stars C4+. **PSFK2** (Libralato et al. 2016): PSF photometry with neighbour subtraction in crowded C0 cluster fields (M35/NGC 2158) to Kp ~ 24 - the only K2 analogue of KBONUS, but limited to cluster superstamps and a ground-based input list, not Gaia.
  - **K2SUPERSTAMP** (Cody et al. 2018, arXiv:1802.06354): stitched, WCS-calibrated mosaics for M35 (2,424 epochs), M67 (3,620 epochs, 400x400 px), Ruprecht 147, NGC 6530; a C0 M35/NGC 2158 pipeline found a substellar companion to a hot star (arXiv:1907.06662). No published catalogue-scale transit search of the mosaics was found.
  - **No K2 equivalent of KBONUS-BKG exists** (searched 2026-10-17): nobody has published Gaia-catalogue-driven PSF deblending of every K2 postage stamp. K2 stamps are small (typically 10-30 pixels), so the background-source yield per stamp is much smaller than Kepler's, but it is non-zero and untouched.

### 3.2 Systematic candidate searches already published
| Search | Campaigns / sample | Candidates |
|---|---|---|
| Vanderburg et al. 2016 (arXiv:1511.07820) | C0-3 | 234 candidates around 208 stars |
| Crossfield et al. 2016 (arXiv:1607.05263) | C0-4 | 197 candidates, 104 validated |
| Osborn et al. 2016 (arXiv:1512.03722) | C1-3 single transits | 7 single-transit candidates |
| Pope, Parviainen & Aigrain 2016 (arXiv:1606.01264) | C5-6 | 152 candidates |
| Mayo et al. 2018 (arXiv:1802.05277) | C0-10, Kp < 13 with spectra | 275 candidates, 149 validated |
| Exoplanet Explorers (Zooniverse; Christiansen et al. 2018, K2-138) | C12+ | citizen vetting of pipeline signals |
| Heller, Hippke & Rodenbeck 2019 (arXiv:1905.09038) | TLS on known K2 planet hosts | 17 new small planets in multi-systems |
| Dattilo et al. 2019 (arXiv:1903.10507) | AstroNet-K2 on K2 TCEs | 2 new super-Earths (K2-293b, K2-294b) |
| Kruse et al. 2019 (ApJS 244, 11) | C0-8, ~150,000 stars Kp = 9-15, QATS | 818 candidates (374 new), 100 around M dwarfs |
| Zink et al. 2021 "Scaling K2 IV" (arXiv:2109.02675; AJ 162, 259) | C1-8, C10-18, fully automated EDI-Vetter | 747 candidates (366 new); max completeness 61%; ~91% reliability |
| Christiansen et al. 2022 "Scaling K2 V" (arXiv:2203.02087) | validation of Zink candidates | 60 validated planets, 37 vetted candidates |
| NASA Exoplanet Archive `k2pandc` | all published | **4,077 rows, 1,559 distinct hosts** (TAP today) |

### 3.3 What is left in K2
- **Campaign 0 and Campaign 19** (excluded from Zink's uniform sample; C0 is covered only by Vanderburg 2016/Crossfield 2016/Kruse 2019, C19 by nobody systematically - but C19 is only ~a week of usable data).
- **Faint stars, Kp > 15:** Kruse's sample stopped at Kp = 15; Zink searched EVEREST targets (which go deeper) but completeness at Kp > 15 is low. K2 M dwarfs at Kp 15-17 within 100-200 pc are a plausible residue.
- **Long periods:** 80-day campaigns mean 2-transit detections for P > 40 d and single transits for P > 80 d; Osborn 2016 covered only C1-3 for singles.
- **Background/non-EPIC sources:** nothing published.
- **Verdict:** real but small residue; K2 fields are at ecliptic latitudes already revisited by TESS, so most bright K2 stars now have 2-10 TESS sectors - a new K2 candidate on a Kp < 13 star would likely already be a TESS TCE/TOI. Check `toi`/`ctoi` before investing.

---

## 4. TESS (2018-present)

### 4.1 Mission status on 2026-10-17 (live checks)
- **Observing:** "TESS is in Orbit 229, Sector 109. Sector 110 begins on October 31, 2026" (https://tess.mit.edu/observations/, fetched 2026-10-17). HEASARC sector table: S108 2026-09-07 to 10-05 ("Cycle 9 begins", southern hemisphere), S109 10-05 to 10-31, S110 10-31 to 11-26, S111 11-26 to 12-22, S112 12-22 to 2027-01-17, S113-114 Jan-Mar 2027 south, S115 ecliptic from 2027-03-14 (https://heasarc.gsfc.nasa.gov/docs/tess/sector.html). Cycle 8 (S97-107) was southern ecliptic pole with "extended sector length" for S97-98; S99 was interrupted by the special 3I/ATLAS pointing (sector "1751", 2026-01-15 to 01-22).
- **Extended mission:** EM3 "is expected to occur September 2025-September 2028" (Guide to TESS for EM planning, HEASARC PDF). Cadence: 200-s FFIs since EM2 (Sector 56); "TESS is currently observing ~13,000 2-minute TPFs, ~2,500 20-second TPFs, and is collecting FFIs every 200 seconds" (same guide). MAST confirms 13,000 SPOC 2-min time-series observations per sector for S100-107 (sandbox).
- **Public release lag (MAST `t_obs_release` vs sector end, sandbox):** S100 27 d (released 2026-03-27), S101 24 d (04-20), S102 31 d (05-22), S103 39 d (06-25), S105 37 d (08-17), S106 32 d (09-10), **S107 24 d (released 2026-10-01)**. No SPOC products yet for S108/109. No proprietary period.
- **Bulk-download page** (https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html): LC, TPF, fast-LC and single-sector DV scripts through **Sector 107**; FFI scripts through Sector 109 (partial); multi-sector DV through **Sectors 1-96** (`tesscurl_multisector_s0001-s0096_dv.sh`) and 14-86.
- **Data release notes page** lists DRN 142 = Sector 104 as the latest entry (page last modified 2026-09-04) - lagging the actual product releases.

### 4.2 Light-curve products: coverage, limits, access (all fetched 2026-10-17)
| Product (MAST author / provenance) | Sectors on MAST | Stars / limit | Method | Access | Notes |
|---|---|---|---|---|---|
| **SPOC** 2-min / 20-s (`SPOC`) | 1-107 | ~13,000 2-min + ~2,500 20-s GI/GTO targets per sector | aperture, PDC | `lk.search_lightcurve(name, author='SPOC')`; bulk `tesscurl_sector_N_lc.sh` | TPS runs every sector; single-sector DV to S107, multi-sector to S96. Fast LCs S27-107. |
| **TESS-SPOC** FFI HLSP (`TESS-SPOC`; Caldwell et al. 2020; DOI 10.17909/t9-wpz1-8s54) | bulk scripts 1-85 (page updated 2026-07-23); lightkurve returns max S69 for southern CVZ stars, S85 for northern | "up to 160,000 targets per sector": all 2-min targets, H <= 10, **d <= 100 pc**, T <= 13.5 | SPOC aperture + PDC on FFIs | author `TESS-SPOC`; bulk per sector | DV/TCE products for FFI targets from S36; multi-sector FFI searches since S56. **Lag ~1 year** behind observations. Faint nearby M dwarfs are *sometimes* included (TIC 142960957, T = 15.35, 49 pc: S4, S31 only) and sometimes not (TIC 130644400, T = 15.86, 58 pc: none). |
| **QLP** (`QLP`; Huang 2020a,b; Kunimoto 2021, 2022; DOI 10.17909/t9-r086-e880) | 1-104 | "TESS magnitude < 13.5" (plus some fainter 2-min targets); 24.4 M segments S1-26, 9.1 M stars S27-39, 5.7 M S40-55 | aperture, KSP detrend | author `QLP`; bulk per sector; target-list CSVs | **"A barycentric timing correction error was identified starting in Orbit 210"** (2nd half of S99) through Orbit 220 (end of S104): "spacecraft velocity was used instead of position"; "Affected sectors ... have not yet been redelivered". S74-79 earlier re-released as v02 after a +/-3 min ephemeris error. eumemic also found QLP corrected flux "dips or rises" at orbit edges in S94-104 and light curves with empty uncertainty columns. |
| **TGLC** (`TGLC`; Han & Brandt 2023, AJ 165, 71; DOI 10.17909/610m-9474) | bulk 1-55; **portal/astroquery only S1-11** ("Sectors 12-55 are not yet available via the MAST Portal, API, or astroquery") | all Gaia DR3 stars to **T = 16**; 1.3-19 M per sector; "<~2% at 16th TESS magnitude even in crowded fields" | ePSF forward model with Gaia positions | bulk curl per sector; `pip install tglc` (v0.8.1, Python 3.10-3.12, "Use a separate environment from the MIT QLP implementation because both packages install as tglc") to build new sectors from TESScut | `cal_aper_flux`/`cal_psf_flux` "may be bad for dim stars near variables". Updated 2025-12-07. |
| **GSFC-ELEANOR-LITE** (`GSFC-ELEANOR-LITE`; Powell et al. 2022; DOI 10.17909/j2yt-t417) | 1-26 (portal/astroquery S1-13 only) | all stars T < 16; "over 150 million light curves" | eleanor aperture | bulk curl | updated 2025-09-16 |
| **T16** (`T16`; Hartman et al. 2025 PASP 137, 024501; DOI 10.17909/8nxx-tw70) | Cycle 1, table shows S1-10 | T < 16; ~45 M LCs in S1-10 | image subtraction (CDIPS pipeline), raw + detrended | bulk curl; portal | **Roth et al. 2026 (arXiv:2604.18579, ApJS 284, 19) searched 83,717,159 Cycle-1 light curves: 11,554 candidates (10,091 new, 411 single-transit), P = 0.5-27 d.** Cycle 1 faint southern sky is therefore now *detected* (vetted by ML + semi-automated). |
| **TARS** (`TARS`; Boyle, Bouma & Mann 2026, arXiv:2603.05586; DOI 10.17909/2sff-fn29; released 2026-06-03) | **1-96** | TIC 8.2 stars with **T < 16 and parallax > 2 mas (d < 500 pc)**: "39 million star-sector combinations spanning 7.5 million unique target stars" | FFI extraction with the `unpopular` causal pixel model | `Observations.query_criteria(provenance_name='tars', sequence_number=N)`; author `TARS` in lightkurve | A **rotation** survey (Lomb-Scargle 0.1-20 d); no transit search published. CPM detrending can suppress or distort transits - validate with injection before trusting. This is the only multi-sector public product for T > 13.5 stars beyond S55 (Kuber's T = 15.3 star has 10 TARS sectors S14-81). |
| **DIAmante** (Montalto et al. 2020, 2022; DOI 10.17909/t9-p7k6-4b32) | 1-26 | 1.4 M FGKM dwarfs/subgiants | difference imaging | catalog CSVs; astroquery | DTARPS (Melton, Feigelson et al. 2024, AJ 167/168; VizieR J/AJ/167/203): 462 southern candidates + 310 Galactic-plane objects, ~58% new. montalto (1,076) and feigelson (208) are top CTOI submitters. |
| **CDIPS** (Bouma et al. 2019; DR6; DOI 10.17909/t9-ayd0-k727) | 1-55 | cluster/young stars Gaia RP < 16; 538,375 LCs (Cycle 3) + 680,065 (Cycle 4) | difference imaging | bulk/astroquery | young-star planet search; bouma = 87 CTOIs |
| **PATHOS** (Nardiello et al. 2019-2021; DOI 10.17909/t9-es7m-vw14) | 1-26 | open-cluster members, PSF in "the faint regime (T > 13)"; 219,256 LCs (DR2) + 150,216 (DR4) | PSF | astroquery `provenance_name='PATHOS'` | 39 + 9 candidates reported; nardiello = 72 CTOIs |
| **TASOC** (asteroseismology) | 1-26 | all stars T < 15 (per TASOC) | aperture/CBV | author `TASOC` | not transit-oriented |
| **TESScut** (FFI cutouts) | 1-109 | any position | you build photometry | https://mast.stsci.edu/tesscut/ | The only route to fresh sectors (86-109) for stars with no product; eumemic's QLP re-search and Kuber's TESS check both used cutouts. |

### 4.3 TESS candidate searches already published (beyond the TOI process)
- **Kunimoto et al. 2022 (ApJS 259, 33), "The TESS Faint-star Search":** 1,617 new TOIs from QLP in Sectors 1-26, 10.5 < T < 13.5, because manual vetting "limits vetting to planet-hosting stars brighter than T = 10.5 mag and leaves millions of potential transit signals unvetted". The faint search continues sector-by-sector (`qlp-sNN-faintsearch` TOI sources through S98 as of today).
- **Olmschenk et al. 2021 (AJ 161, 273):** CNN on FFI light curves, 181 candidates; **Salinas et al. 2025 (MNRAS 538, 2031):** transformer, 214 candidates (88 single-transit) S1-26; olmschenk (184) and salinas (207) CTOIs.
- **Planet Hunters TESS** (Eisner et al. 2021, "PHT II"): 90 candidates from S1-26, 73 single-transit; recovers 85% of TOIs > 4 Re and 51% at 3-4 Re; Malik et al. 2025 (AJ) CNN trained on volunteer scores; eisner = 182 CTOIs. Current project status not verified (Zooniverse page not checked).
- **ExoMiner++** (Valizadegan et al. 2025, arXiv:2502.09790): 147,568 2-min TCEs -> 7,330 planet-candidate classifications (1,868 TOIs, 69 CTOIs, 50 new CTOIs); **ExoMiner++ 2.0** (arXiv:2601.14877, Jan 2026) extends to FFI data. This means the 2-min TCE backlog has been machine-vetted once already.
- **TESS Ten Thousand Catalog** (Kostov et al. 2025, arXiv:2506.05631): neural network over FFIs S1-82 flagged ~1.2 M stars; 10,001 vetted EBs (7,936 new; median T = 13.8) plus ~900,000 unvetted candidates above score 0.9 - an essential false-positive cross-match (Section 8).
- **M-dwarf specific:** NEMESIS I (Feliz et al. 2021, arXiv:2103.05647; FFI M dwarfs S1-5) and NEMESIS II (poster, ~191,000 FFI light curves of M dwarfs within 100 pc, S1-26; Zenodo 12973447); Ment & Charbonneau 2023 (volume-complete mid-late M within 15 pc, 2-min, S1-42, P <~ 7 d); feliz = 16 faint CTOIs.
- **Single/duo-transits:** PHT II (73 singles), Salinas 2025 (88), Roth 2026 (411), Bass & Fabrycky 2025 (duotransits).
- **The October-2026 AI-agent wave** (for scope only): eumemic - 4,166,638 FGKM dwarfs T <= 13.5, Teff <= 6,500 K, R <= 1.5 Rsun, S1-104 (4,086,108 with TESS-SPOC/QLP/SPOC light curves); 22,448 M dwarfs 11.5 < T <= 12.5; 84,863 M dwarfs 12.5 < T <= 13.5; 9,979 nearby K/M dwarfs with S97-105; 6,430 TOI hosts; 14,112 non-TOI TCEs re-measured; 798,769 inverted stars; 645 counted candidates, 255 meeting validation thresholds bar imaging (README 2026-10-16). Chopra - 526 GCNS stars <= 60 pc never at 2-min. A third project, `tess-transit-hunter` (comdex4), posted G 249-11 on 2026-10-05 (cited in eumemic's README; not independently checked).

### 4.4 The TOI process and the measurable "TOI lag"
- Prime-mission process (Guerrero et al. 2021, arXiv:2103.12538, full text via pdftotext): "The SPOC pipeline produces, on average, about 1,200 TCEs per sector, out of about 20,000 two-minute cadence targets observed. Before triage, QLP produces 50,000 TCEs per sector (with a Tmag = 10.5 cutoff for vetting), from about half a million FFI light curves. Each sector yields about 100 total TOIs". "Astronet-Triage ingests tens of thousands of TCEs brighter than Tmag = 13.5 and selects several thousand (3000-4000) targets as potential planet candidates. Next, this number is further reduced by setting a brightness threshold of Tmag = 10.5". "TCEs from the QLP below the Tmag = 10.5 threshold were not considered for vetting". TOI sources today show three streams: `spoc-sNNN-b02` (single-sector SPOC), `spoc-s14-s86-b0A/b02` (multi-sector SPOC incl. FFI targets, T up to 14.8), `qlp-sNN-ffi` and `qlp-sNN-faintsearch`.
- **TOI table today (TAP):** 8,151 TOIs; dispositions PC 4,820, CP 815, KP 606, FP 1,315, FA 100, APC 490; latest `toi_created` 2026-10-08; 330 TOIs created in 2026 (Jan 69, Mar 23, Apr 18, Jun 133, Jul 49, Aug 35, Oct 3), 87 since 2026-07-01; median T of 2026 TOIs 12.5.
- **Lag by stream (ExoFOP TOI CSV `Date TOI Alerted` vs HEASARC sector end; sandbox):** SPOC single-sector: S95 148 d, S96 122 d, S97 116 d, S98 94 d, S99 98 d, S100 95 d, S102 107 d (alerted 2026-08-06). QLP faint-search: S95 338 d, S96 312 d, S97 257 d, S98 200 d (all alerted 2026-07-24). Nothing from S103-107 yet, even though S103 has been public since 2026-06-25. (QLP S99-104 are additionally blocked by the timing bug.)
- **Interpretation:** at any moment there are ~4-5 public sectors with SPOC TCEs that have not been turned into TOIs (currently S103-107), ~9 sectors without QLP faint-star TOIs (S99-107), and ~20 sectors (S86-107) without TESS-SPOC FFI light curves at all. Anyone who downloads the single-sector DV products (`tesscurl_sector_107_dv.sh`) gets the SPOC TCE list months before the TOI team vets it; anyone who builds FFI photometry from TESScut for S86-107 gets stars no public pipeline has touched.

### 4.5 Under-searched TESS regimes, with numbers
- **Faint-end deficit in the candidate lists (ExoFOP today):** TOIs with T > 13.5: **457** (PC 237, KP 75, FP 71, CP 45, APC 24); T > 15: 35; T > 16: 9. CTOIs: 5,138 total (814 promoted), median T 11.6, 75th percentile 12.3; T > 13.5: 351; T > 15: 44; T > 16: 9. Faint CTOIs come mostly from montalto (96), bouma (43), nardiello (35), olmschenk (30), wilson (27).
- **Nearby faint M dwarfs (TIC criteria queries, sandbox):** Teff < 3,900 K, lumclass DWARF: **100,567 within 100 pc with 13.5 < T < 16.5** (T bins 13.5-14: 21,367; 14-14.5: 21,172; 14.5-15: 19,232; 15-15.5: 16,346; 15.5-16: 12,954; 16-16.5: 9,490; median T 14.7, median R 0.23 Rsun, median Teff 3,192 K; 160 within 25 pc, 6,329 at 25-50 pc, 28,275 at 50-75 pc, 65,803 at 75-100 pc). **Only 33 are TOI hosts and 43 are CTOI hosts.** Within 50 pc and 13.5 < T < 16: 6,196; within 30 pc and T < 16: 4,861. TOIs with Teff < 4,000 K and d < 100 pc: 373, of which only 35 have T > 13.5 and none in 15-16.
- **Public light-curve coverage of this pool is incomplete:** TESS-SPOC includes "targets within 100 parsecs" up to the 160k/sector cap and only through S85; of three random faint nearby M dwarfs probed, one (T = 14.3, 55 pc, northern CVZ, a 2-min target) had SPOC/TESS-SPOC/QLP/TARS to S85-86; one (T = 15.35, 49 pc) had TESS-SPOC in only 2 of its sectors plus TARS/TGLC/eleanor-lite/T16 in S4; one (T = 15.86, 58 pc) had **no** SPOC-family product at all, only TARS (7 sectors to S88), TGLC and eleanor-lite (S6, S8). Expect to build your own photometry for a large fraction.
- **Bright stars the agent wave skipped:** eumemic's dwarf cuts (Teff <= 6,500 K, R <= 1.5 Rsun) exclude hot/early-F stars and subgiants; giants are low-value. Stars with T <= 13.5 observed only in S105-107 (first-time coverage) are inside the TOI lag but eumemic already searched S1-104/105 for dwarfs.
- **Period regime:** QLP/SPOC searches to ~ half the baseline; single-sector FFI targets (TESS-SPOC DV) use a 2-h minimum transit duration ("minimum transit duration is 2 hours for FFI data (versus 0.5 hours for two-minute data)"), so short-duration transits of M dwarfs (1-1.5 h) in FFI data are under-searched by SPOC - a specific, documented gap that favours M-dwarf FFI searches with your own BLS/TLS using 0.5-1 h durations.

---

## 5. Other datasets

| Dataset | Status 2026-10-17 | Transit-search value for a solo developer |
|---|---|---|
| **Gaia DR3 epoch photometry** | Released 2022 for ~11.7 M variable sources only; `vari_planetary_transit` has 214 sources (173 known, 41 new; Gaia-1b, Gaia-2b confirmed by RV; Panahi et al. 2022, A&A 663, A101) | Sparse (~12 epochs/yr); exhausted by DPAC. |
| **Gaia DR4** | "planned for Wednesday, 2 December 2026" (ESA Gaia news, announced 2026-05-07); 66 months of data; ~2.8 billion sources, ~400-500 TB; `epoch_photometry` 121,845,163,512 records (7.2 TB); planned `vari_planetary_transit` table with **1,007 rows** (draft data model, https://www.cosmos.esa.int/web/gaia/dr4) | A December event. ~70 G epochs per star over 5.5 yr: hot-Jupiter-depth only; DPAC will publish its own 1,007 candidates; useful mainly as (i) a novelty cross-match and (ii) a way to confirm/extend ephemerides of your TESS candidates via DataLink. Not a primary hunting ground. |
| **ZTF** (Masci et al. 2019, arXiv:1902.01872) | DR24 public via IRSA; "photometric precision (repeatability) at bright unsaturated fluxes varies between 8 and 25 millimag"; g ~ 20.8, r ~ 20.6 (5 sigma, 30 s) | Kuber's own test: 1,843 epochs gave a window depth of 1,780 +/- 6,141 ppm. Only ~1%+ transits (hot Jupiters, white-dwarf transits) are detectable; northern sky only. Cross-check/confirmation tool, not a hunting ground. |
| **ASAS-SN** (Kochanek et al. 2017, arXiv:1706.07060) | all-sky to V ~ 17, nightly; precision a few x 0.01 mag for bright stars | Same verdict; useful for EB/variability vetting (Sky Patrol light curves). |
| **ATLAS** (Tonry et al. 2018, arXiv:1802.00879) | whole sky, 2-day cadence to m ~ 19 | Same; ATLAS forced photometry server is a good vetting input. |
| **CHEOPS** | Mission archive + DACE; "By default, the proprietary period is 1 year" (max 1.5 yr); Level-2 light curves public thereafter; `dace_query`/`pycheops` | Targeted pointings of known systems; serendipitous new planets happen only around already-known hosts (e.g. extra transits in GTO targets). Not a survey. |
| **PLATO** | ESA page: launch "planned for March 2027" on Ariane 6, 26 cameras | Nothing before 2028. |
| **Roman** | Launched 2026-08-30 07:26 EDT on Falcon Heavy (Space.com, NASA); ~3-month cruise/commissioning; first images expected ~Jan 2027 (single source) | Galactic Bulge Time-Domain Survey not started; data policy: no proprietary period, but nothing to search in 2026. |
| **Kepler/K2 GO superstamps** | Kepler NGC 6791/6819 (200x200 px, full mission); K2 M35/M67/Rup 147/NGC 6530 mosaics on MAST | Open ground with no published transit catalogue, but crowded, distant (G > 17 dwarfs) and un-validatable; M67 (850 pc, solar-age) is the only scientifically attractive one and PSF-based searches already exist (Libralato et al. 2016, arXiv:1608.08539). |
| **TESS 3I/ATLAS special sector "1751"** | S1751 (2026-01-15 to 01-22) has LC/FFI/DV products | One week of ecliptic data; negligible. |

---

## 6. Access cheat-sheet (verified commands)

```python
import lightkurve as lk, astropy.units as u
from astropy.coordinates import SkyCoord
# Kepler Bonus background source (Kuber's star)
lk.search_lightcurve(SkyCoord(285.907827, 40.941988, unit='deg'), author='KBONUS-BKG', radius=20*u.arcsec)
# TESS by author; sequence_number may be masked -> guard with np.ma.is_masked
lk.search_lightcurve('TIC 150428135', author='QLP')          # to S104 for S-CVZ stars
lk.search_lightcurve('TIC 150428135', author='TESS-SPOC')    # to S69 for S-CVZ stars (bulk to S85)
lk.search_lightcurve('TIC 120686465', author='TARS')         # T<16, d<500 pc, S1-96
lk.search_lightcurve('TIC 4206066',  author='TGLC')          # indexed S1-11 only; bulk S1-55; `pip install tglc` for later
# TIC criteria query without a position (works; ~15-60 s)
from astroquery.mast import Catalogs
Catalogs.query_criteria(catalog='Tic', d=[0,100], Tmag=[13.5,16.5], Teff=[2300,3900], objType='STAR', lumclass='DWARF')
```
```bash
# NASA Exoplanet Archive TAP (sync, CSV)
curl -sS https://exoplanetarchive.ipac.caltech.edu/TAP/sync --data-urlencode "query=select count(*) from toi where st_tmag>13.5" --data-urlencode format=csv
# tables used: toi, cumulative, k2pandc, k2targets, pscomppars, keplerstellar (990,244 rows), q1_q17_dr25_tce
# ExoFOP CSVs (no auth)
curl -o toi.csv  'https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv'    # 8,151 rows, has Sectors + Date TOI Alerted
curl -o ctoi.csv 'https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv'  # 5,138 rows, has User + Promoted to TOI
# MAST bulk scripts (TESS)
https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html   # tesscurl_sector_107_{lc,tp,dv,fast-lc}.sh ; multisector s0001-s0096_dv
# HLSP bulk: QLP, TESS-SPOC, TGLC, GSFC-ELEANOR-LITE, T16, TARS, CDIPS, KBONUS-BKG pages each have per-sector curl scripts + target-list CSVs
# Kepler Bonus catalogues
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits  # 159 MB
https://archive.stsci.edu/hlsps/kbonus-bkg/catalogs/hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits  # 2.4 MB, 29,800 rows
```
Data volumes: one TESS SPOC 2-min LC ~ 2 MB/sector; QLP LC ~ 0.1-0.3 MB/sector; a full QLP sector bulk is tens of GB; KBONUS M-star bundle 134 GB; TGLC per sector 1-19 M files. Precision vs magnitude (documented): Kepler/KBONUS CDPP < 100 ppm for G < 16; TESS QLP ~ 1% at T = 13.5 per 30-min point (Huang 2020); TGLC "<~2% at 16th TESS magnitude"; eumemic measured that a 2,400 ppm dip on a T = 15.3 star becomes ~1,000 ppm in a 1.5-px aperture against 1,294 ppm scatter (Kuber's injection test) - i.e. at T ~ 15 you need either PSF photometry or many sectors.

---

## 7. Known-object catalogues for the novelty claim (cross-match all of them)

| Catalogue | Size today | Where |
|---|---|---|
| KOI cumulative | 9,564 (2,748 confirmed, 1,977 candidates, 4,839 FP) | TAP `cumulative`; also `q1_q17_dr25_tce` (34,032 TCEs) and `fpwg` false-positive table |
| Kepler EB catalogue (Villanova; Kirk et al. 2016) | 2,878 EBs | http://keplerEBs.villanova.edu ; VizieR J/AJ/151/68 ; static copy at MAST |
| Kepler Bonus KOIs+neighbours catalogue | ~32,000 files | KBONUS-BKG `kois-nns-catalog` (which background sources sit in KOI apertures) |
| K2 candidates | 4,077 rows / 1,559 hosts | TAP `k2pandc` |
| Confirmed planets | 6,451 (pscomppars) | TAP `pscomppars`/`ps` |
| TOI | 8,151 (ExoFOP has Sectors, alert dates) | TAP `toi`; ExoFOP CSV |
| CTOI | 5,138 (814 promoted) | ExoFOP CSV only |
| TESS EBs I (Prsa et al. 2022, ApJS 258, 16) | 4,584 EBs, 2-min S1-26 | https://archive.stsci.edu/hlsp/tess-ebs (CSV; CasJobs); live version at Villanova |
| TESS Ten Thousand (Kostov et al. 2025) | 10,001 vetted EBs S1-82 + ~900,000 unvetted ML candidates | arXiv:2506.05631 tables (check MAST/Zenodo for the machine-readable release) |
| Planet Hunters TESS | 90 candidates (PHT II) + CTOIs (eisner 182, conzo 108...) | ExoFOP CTOI; Eisner et al. 2021 |
| DTARPS, DIAmante, Olmschenk, Salinas, T16 Planet Hunt (11,554), Kunimoto 2022 (1,617) | published candidate tables | VizieR J/AJ/167/203; CTOI; arXiv:2604.18579 |
| Gaia DR3 | `vari_classifier_result`, `vari_summary`, `vari_planetary_transit` (214), `nss_two_body_orbit`, RUWE | Gaia archive TAP; DR4 adds 1,007 transit candidates on 2026-12-02 |
| SIMBAD | object types EB*, Pl, Candidate_Pl etc. | `astroquery.simbad` cone + otype |
| ASAS-SN/ATLAS variable catalogues; VSX | EB/variable classifications | Sky Patrol; VizieR |
| Barber et al. 2022 Table 5 and similar membership lists | "has this star's light curve been examined" questions | CDS J/AJ/164/88 (lesson from Kuber's RNAAS: check whether the star appears in *any* search's target list, not only candidate lists) |

---

## 8. Comparison table

| Dataset / pool | Approx. size | Precision | Detection search done? | Vetting done? | Validation prospects (FPP < 1.5%) | Novelty risk / competition | Verdict for a solo dev |
|---|---|---|---|---|---|---|---|
| Kepler primary targets | 198,709 | 30-100 ppm | Yes (DR25 + 8 independent) | Yes | Good (bright) but nothing left | saturated | No |
| **Kepler Bonus M-dwarf background** | 1,836 (quality-cut) | <100 ppm for G<16; most G>17 | Kuber (1,736 stars) | Kuber only | Poor: G 17-19, no AO, TESS cannot recover | Goddard TraSHDUMP team | Only as replication/extension |
| **Kepler Bonus non-M background** | 25,542 (6,151 with G<16) | <100 ppm for G<16 | **No** | No | Fair for G<15 (RV/AO possible on K dwarfs at G 13-15), poor beyond | Goddard team ("several hundred" expected) | **Yes, fast and cheap; publish quickly** |
| Kepler superstamps NGC 6791/6819 | ~10^4 stars | poor (crowded) | No transit search | No | None | low | No |
| K2 all targets | 361,770 EPIC | ~Kepler for Kp<13, x2 worse fainter | Kruse (Kp 9-15, C0-8), Zink (C1-8, 10-18; 61% complete) | Yes | Fine for bright | overlaps TESS now | Residue only (C0, C19, Kp>15, P>40 d) |
| K2 background/superstamps | unknown | poor | No | No | None | low | No |
| **TESS T<=13.5 dwarfs, S1-104** | 4.2 M | 0.1-1% | QLP/SPOC + eumemic sweep | TOI to S102/S98; eumemic 645 | Good | eumemic, Chopra, ExoMiner++ | Only via the lag window / non-dwarf cuts |
| **TESS fresh sectors (S103-107 public, S108-109 pending)** | all stars | as above | SPOC TCEs exist (DV S107); QLP S99-104 buggy | **No TOIs yet** | Good for T<13 | eumemic to S104/105; TOI team in 3-10 months | **Yes, as a standing monthly job** |
| **TESS faint nearby M dwarfs (13.5<T<16.5, d<100 pc)** | 100,567 (6,196 within 50 pc) | 0.5-2% per 200-s point; transits of 1-3 Re are 0.1-1% | Partial (TESS-SPOC subset to S85; TGLC S1-55; T16 Cycle 1; TARS LCs exist but unsearched) | 33 TOI + 43 CTOI hosts | **Good**: Ks 10-12, AO/Speckle feasible, ground photometry of 0.5-1% dips feasible, RV for M3-M5 hard but possible | NEMESIS poster; eumemic stops at 13.5 | **Yes - best combination of novelty and validatability** |
| TESS single/duo transits | all | - | PHT, Salinas, Roth, Bass | partial | Poor without period | moderate | Side-product only |
| Gaia DR4 epoch photometry | 2.8 B sources, ~70 epochs | 1-10 mmag bright | DPAC (1,007 candidates) | DPAC | Good (bright hosts) but hot Jupiters only | DPAC itself | Cross-match; opportunistic after 2026-12-02 |
| ZTF / ASAS-SN / ATLAS | 10^8-10^9 | 1-3% | Variability surveys | EB catalogues | - | - | Vetting inputs only |
| CHEOPS | targeted | 10s ppm | n/a | n/a | - | - | No |
| PLATO / Roman | future | - | - | - | - | - | No (2027+) |

---

## 9. Ranked recommendation for a solo developer with 4 CPUs and Claude Code

**Rank 1 - Faint nearby M dwarfs in TESS (13.5 < T < 16, d < 100 pc, all sectors), with FFI photometry you control.** Why: 100,567 stars, 33 TOI hosts; small planets give 0.1-1% transits (median R = 0.23 Rsun: a 1.5 Re planet is ~3,500 ppm); eumemic's and Chopra's sweeps stop at T = 13.5 and 60 pc respectively; SPOC's 2-h minimum FFI duration under-searches short M-dwarf transits; hosts are Ks ~ 10-12 so Keck/Palomar/SOAR speckle or AO, ground photometry (0.4-1 m telescopes reach 1-3 mmag on these), and reconnaissance spectroscopy are all feasible, which is what TRICERATOPS/MOLUSC needs to get below FPP 0.015 (eumemic got 255 of 645 to threshold "apart from high-resolution imaging"). How: start with the 6,196 stars within 50 pc (TIC criteria query, 14 s); for each, pull every public product (TESS-SPOC, QLP, TGLC via `tglc`, TARS, eleanor-lite) *and* build your own 200-s aperture + PRF-weighted photometry from TESScut for all sectors to S107 so coverage is uniform; detrend (wotan biweight), TLS/BLS 0.5-30 d with durations down to 0.5 h, inverted-flux and shuffled-sector nulls, LEO-Vetter-style tests, per-sector centroid/difference-image localisation (eumemic rejected 93% of LEO-Vetter passes at the pixel stage - budget for that), then TRICERATOPS with Gaia-informed MOLUSC. Expected throughput: ~5-15 s per star-sector for cutout photometry + search on one core, i.e. the 50-pc sample in ~1-2 days of 4-CPU wall time; the full 100-pc sample is a Modal/cloud job. Risk: pixel confusion in 21-arcsec pixels; mitigate with Gaia DR3 neighbour lists and the Kostov EB list.

**Rank 2 - The fresh-sector TOI-lag window as a standing job.** Each new sector becomes public 24-39 days after it ends and TOIs follow 3-5 (SPOC) to 7-11 (QLP) months later. Today S103-107 are public with zero TOIs; S108 lands ~early November, S109 ~early December. Pull `tesscurl_sector_N_dv.sh` (SPOC TCEs with DV reports) and `tesscurl_sector_N_lc.sh` the day they appear; cross-match TCEs against TOI/CTOI; run your own search on QLP (where unaffected by the timing bug) and on TESScut photometry for T > 13.5. The bright end overlaps eumemic (through S104/105) and the TOI team, so prioritise the faint and the non-dwarf (Teff > 6,500 K, R > 1.5) targets, and multi-sector combinations where a new sector pushes an old sub-threshold signal over SNR 9 (eumemic: "only 24 needed Sectors 97-104" but "47% ... already reach SNR 9" with S<=55 - new sectors mostly *confirm* rather than create). Pre-register predictions on Zenodo as Chopra did.

**Rank 3 - Kepler Bonus non-M background stars (25,542; 6,151 with G < 16).** Cheapest novel pool: public PSF light curves, 4 years of 30-min cadence, CDPP < 100 ppm for G < 16, and nobody has published a search. Run Kuber's exact pipeline (reproduced funnel above) with his inverted-flux null and his dilution/centroid tests. Trade-offs: (i) a funded Goddard team with the data's authors is doing exactly this ("several hundred new planet candidates" expected) - move fast and publish as a Zenodo/RNAAS catalogue rather than one-candidate papers; (ii) at G > 16 a candidate is un-validatable (no AO, TESS too shallow, Kuber's FPP stuck at 0.069), so restrict to G < 15-16 where a few per cent of the pool lives; (iii) the editor's objections to Kuber's RNAAS (sample construction, inverted-null parity, TESS non-detection, prior-search status) will apply to you - build the null comparison and the membership/target-list cross-check into the pipeline from day one.

**Rank 4 - K2 residue:** C0 and C19, EPIC stars with Kp > 15, and P > 40 d two-transit searches on EVEREST/K2SFF light curves; cheap because the light curves are small and the catalogues (k2pandc, 4,077) are complete. Low yield; mostly useful as a pipeline shakedown because K2 systematics are a good stress test.

**Rank 5 - Everything else as vetting inputs, not hunting grounds:** ZTF/ATLAS/ASAS-SN forced photometry for EB/variability checks; Gaia DR3 (RUWE, NSS, variability classes) and DR4 epoch photometry after 2 Dec 2026 for ephemeris extension and blend checks; CHEOPS public visits for ephemeris refinement of known hosts only.

**What not to do:** Kepler primary targets, Kepler/K2 superstamps, Gaia-only transit searches, and ground-survey-only transit searches - either saturated or physically incapable of yielding a validatable small planet.

---

## 10. Open questions and caveats
- Whether TARS light curves preserve transit depth and shape after `unpopular` CPM detrending is untested here; an injection-recovery on a few hundred TARS stars is a prerequisite before using TARS as a search product.
- TESS-SPOC inclusion of faint stars within 100 pc appears inconsistent (one of two probed faint M dwarfs had none); the per-sector target-list CSVs on the HLSP page are the authoritative way to decide which stars need home-made photometry.
- QLP S99-104 redelivery date is unknown; eumemic tracks this and marks results stale when MAST redelivers.
- The HLSP page's "~20,000" M-stars statement vs the file's 29,800 rows: the file is authoritative; the discrepancy should be noted in any paper.
- Whether the Goddard TraSHDUMP/BOI search has already covered the non-M Kepler Bonus background (the poster mentions only BOI-1); a polite email to the poster authors before investing weeks is advisable.
- Planet Hunters TESS current activity and whether its unpublished candidate list is uploaded to MAST were not verified.
- The 13,000 2-min targets/sector figure comes from the EM planning guide and the MAST observation count; 20-s target count (~2,500) is from the guide only.
- Precision numbers for ASAS-SN and ATLAS are from secondary sources (the abstracts do not state them); treat as approximate.

---

## Sources (fetched/run 2026-10-17 unless dated)
MAST HLSP pages: KBONUS-BKG https://archive.stsci.edu/hlsp/kbonus-bkg ; KBONUS-APEXBA https://archive.stsci.edu/hlsp/kbonus-apexba ; QLP https://archive.stsci.edu/hlsp/qlp ; TESS-SPOC https://archive.stsci.edu/hlsp/tess-spoc ; TGLC https://archive.stsci.edu/hlsp/tglc ; GSFC-ELEANOR-LITE https://archive.stsci.edu/hlsp/gsfc-eleanor-lite ; T16 https://archive.stsci.edu/hlsp/t16 ; TARS https://archive.stsci.edu/hlsp/tars ; DIAmante https://archive.stsci.edu/hlsp/diamante ; CDIPS https://archive.stsci.edu/hlsp/cdips ; PATHOS https://archive.stsci.edu/hlsp/pathos ; TESS-EBs https://archive.stsci.edu/hlsp/tess-ebs ; EVEREST https://archive.stsci.edu/hlsp/everest ; K2 HLSP index https://archive.stsci.edu/k2/hlsps.html ; K2SUPERSTAMP https://archive.stsci.edu/prepds/k2superstamp/ ; TESS bulk downloads https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html ; TESS DRN https://archive.stsci.edu/tess/tess_drn.html.
Mission status: https://tess.mit.edu/observations/ ; https://heasarc.gsfc.nasa.gov/docs/tess/sector.html ; Guide to TESS for EM planning https://heasarc.gsfc.nasa.gov/docs/tess/docs/Guide-to-TESS-for-EM-planning.pdf ; TOI releases https://tess.mit.edu/toi-releases/ ; ExoFOP https://exofop.ipac.caltech.edu/tess/ ; Gaia DR4 https://www.cosmos.esa.int/web/gaia/dr4 and https://www.cosmos.esa.int/web/gaia/news-2026 ; PLATO https://www.esa.int/Science_Exploration/Space_Science/Plato ; Roman launch https://www.space.com/space-exploration/launches-spacecraft/nasa-roman-space-telescope-launch-success-on-a-spacex-falcon-heavy-rocket-spectacular-launch-video ; CHEOPS proprietary data https://www.cosmos.esa.int/web/cheops-guest-observers-programme/proprietary-data.
Papers: Martinez-Palomera, Hedges & Dotson 2023 arXiv:2310.17733 ; Moquin et al. 2026 poster https://vmnexsci.ipac.caltech.edu/workshop/2026/posters/Poster_JulieMoquin_41.pdf ; Twicken et al. 2016 arXiv:1604.06140 ; Thompson et al. 2018 arXiv:1710.06758 ; Huang, Bakos & Hartman 2012 arXiv:1205.6492 ; Wang et al. 2015 arXiv:1512.02559 ; Foreman-Mackey et al. 2016 arXiv:1602.07848 ; Herman et al. 2019 arXiv:1901.01974 ; Kawahara & Masuda 2019 arXiv:1904.04980 ; Kunimoto et al. 2020 arXiv:2003.04397 ; Kirk et al. 2016 arXiv:1512.08830 ; Kuehn et al. 2013 arXiv:1309.6044 ; Espinoza-Rojas et al. 2023 arXiv:2310.06757 ; Montet et al. 2017 arXiv:1705.07928 ; Vanderburg et al. 2016 arXiv:1511.07820 ; Crossfield et al. 2016 arXiv:1607.05263 ; Osborn et al. 2016 arXiv:1512.03722 ; Pope et al. 2016 arXiv:1606.01264 ; Mayo et al. 2018 arXiv:1802.05277 ; Heller et al. 2019 arXiv:1905.09038 ; Dattilo et al. 2019 arXiv:1903.10507 ; Kruse et al. 2019 ApJS 244, 11 ; Zink et al. 2021 arXiv:2109.02675 ; Christiansen et al. 2022 arXiv:2203.02087 ; Cody et al. 2018 arXiv:1802.06354 ; Libralato et al. 2016 arXiv:1608.08539 ; Guerrero et al. 2021 arXiv:2103.12538 ; Kunimoto et al. 2022 ApJS 259, 33 ; Caldwell et al. 2020 RNAAS 4, 201 ; Huang et al. 2020a,b ; Han & Brandt 2023 AJ 165, 71 ; Powell et al. 2022 ; Hartman et al. 2025 PASP 137, 024501 ; Roth et al. 2026 arXiv:2604.18579 ; Boyle et al. 2026 arXiv:2603.05586 ; Montalto et al. 2020/2022 ; Melton et al. 2024 (DTARPS I-III; VizieR J/AJ/167/203) ; Bouma et al. 2019 ; Nardiello et al. 2019-2021 ; Prsa et al. 2022 ApJS 258, 16 ; Kostov et al. 2025 arXiv:2506.05631 ; Olmschenk et al. 2021 AJ 161, 273 ; Salinas et al. 2025 MNRAS 538, 2031 ; Eisner et al. 2021 (PHT II) ; Malik et al. 2025 AJ ; Valizadegan et al. 2025 arXiv:2502.09790 and 2026 arXiv:2601.14877 ; Feliz et al. 2021 arXiv:2103.05647 ; NEMESIS II poster Zenodo 12973447 ; Ment & Charbonneau 2023 ; Panahi et al. 2022 A&A 663, A101 ; Masci et al. 2019 arXiv:1902.01872 ; Kochanek et al. 2017 arXiv:1706.07060 ; Tonry et al. 2018 arXiv:1802.00879.
Repos: https://github.com/eumemic/exoplanets (README 2026-10-16) ; https://github.com/paraschopra/tess-planet-candidates ; https://github.com/Kuberwastaken/planethunt (read-only clone) ; https://github.com/TeHanHunter/TESS_Gaia_Light_Curve.
Sandbox outputs: `toi_exofop.csv`, `ctoi.csv`, `tic_mdwarfs_100pc_faint.csv`, `kbonus/kbonus_{M,nonM}_background_survivors.csv`, `dl/guerrero.txt`, `dl/moquin_poster.txt`, `em_guide.txt` under `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/other-datasets/`.
