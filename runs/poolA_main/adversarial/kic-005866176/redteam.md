# kic-005866176: red-team review (2026-10-17)

**Role.** Try to kill the 29.433 d first-screen signal on KIC 5866176 (Gaia DR3 2102115328638366592, G = 15.28, BP-RP = 0.81, 2.0 kpc, M_G = 3.8, R* ~ 0.95 Rsun; extracted by KBONUS-BKG from the TPF of KIC 5866192, Kp 14.50, 13.3 arcsec = 3.3 px away). Language: this is a *first-screen signal*; nothing below makes it a candidate or a planet.

**Verdict: kill.** The first-screen and stage-2 statistics (depth 481 +- 91 ppm, empirical SNR 7.54, TLS SDE 8.46; `vet.json`, run `poolA_main`) are produced by two post-gap PSF-photometry excursions on channel 17 that are anti-correlated between PSF and SAP flux and are therefore not a loss of light. With the first 0.3 d after every data gap removed, the signal is not recoverable blind (BLS 20-40 d: SDE 3.6 at 29.43 d, below the inverted-flux null top of 5.9) and the residual box depth at the pre-chosen ephemeris is 272 +- 48 ppm (empirical SNR 5.7), below the pre-registered stage-2 floor `snr_emp_pass2_min = 7.1`. Pixel-level localisation is impossible for this star (its Gaia position lies 0.2-0.4 px outside the stamp in every quarter). A ~270 ppm, ~5 sigma feature at the fixed ephemeris remains in PSF, SAP and PSF_FLAT flux; it cannot be excluded as astrophysical but it does not constitute a signal under this pipeline's registered rules and it was selected on data dominated by the artefact.

All my numbers trace to the scripts and JSON under `runs/poolA_main/adversarial/kic-005866176/redteam/` (`lc_checks.py/.json/.png`, `pixel_checks.py/.json/.png`, `event_rawflux.log/.json`, `psfflat_check.log/.json`; SHA-256 in `MANIFEST.sha256`). Nothing from `pipeline/` was imported. Inputs: the cached KBONUS file `data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005866176_kepler_v1.0_lc.fits` (39,966 good cadences, 11 quarters, `SAP_QUALITY == 0`), the 11 Kepler long-cadence TPFs of KIC 5866192 for the same quarters (MAST, downloaded 2026-10-17 to scratch), the KBONUS light curve of KIC 5866192 (MAST, scratch), `data/known/` (fetched 2026-10-17T15:46Z per `manifest.json`), the NASA Exoplanet Archive TAP (cumulative KOI, TOI; 2026-10-17), SIMBAD (2026-10-17), `runs/poolA_main/triage_all.csv`, and the independent re-derivation in `rederive/`.

## 1. Is it already known? (no)

| check | result |
|---|---|
| KIC 5866176 or the TPF host KIC 5866192 as KOI / DR25 TCE / Kepler EB (Kirk+2016) / confirmed planet | absent from all four; 5866176 is not a Kepler target at all (`kepler_target_kepids`), 5866192 is (Kp 14.50, Teff 6050 K, R 1.0 Rsun) with no TCE |
| KOI / TCE / EB / PS within 3 arcmin | KOI-506.01 = KIC 5780715 (101 arcsec, P 1.583 d, FALSE POSITIVE with centroid-offset flag), TCEs on KIC 5866029 (225 d), 5866098 (0.93, 0.93, 66.8 d), 5866138 (342 d; also the Kepler EB 5866138, 45 arcsec away), 5780782 (ten long-period TCEs, 173 arcsec). None at 29.43 d or any harmonic. PS: none. |
| ephemeris match, all KOIs/TCEs field-wide at P, 2P, P/2, 3P/2, 2P/3, 3P, P/3, 4P, P/4 (abs dP/P < 2e-3) with phase alignment | 176 period matches; only three are phase-aligned within 0.5 d: KOI/TCE 10467714-01 (58.749 d ~ 2P, 406 arcmin away, dP/P = 2e-3 which drifts 3 d over the baseline, i.e. incoherent), TCE 4843152-01 (19.598 d ~ 2P/3, 255 arcmin), TCE 6721123-01 (7.36 d ~ P/4, 474 arcmin). All on other modules; chance hits. The five cumulative KOIs with 29.38 < P < 29.49 d (TAP, 2026-10-17: K01549.01, K02123.01, K00271.02, K03265.01, K08035.01) have t0 offsets of 7-21 d from ours: none aligned. |
| TOI in a 3 arcmin box (TAP) | 0 rows |
| SIMBAD within 20 arcsec | one object: UCAC4 656-070604 = KIC 5866192, the host (plx 0.91 mas, i.e. a foreground star at ~1.1 kpc vs the target at ~2.0 kpc); nothing at the target position |
| Bienias & Szabo 2025 (arXiv:2501.09152, "Background exoplanet candidates in the original Kepler field") | title confirmed on arXiv; the candidate list as transcribed in the sibling review `adversarial/kic-011138423/redteam.md` (hosts near KIC 4459924, 8175131, 8293539, 4543171, 9674230, 9040849) does not include this field. I did not re-read the paper's table myself. |
| survey-wide (`triage_all.csv`, 6,365 stars) | six other pool stars peak within 0.3 % of 29.43 d (none pass; phase offsets 7.7-13.5 d from our t0), so there is no shared systematic ephemeris. But 349 of the 438 real first-screen passes, and 2,936 of 6,365 best periods, lie in 24-40 d: the monthly-gap systematics class named in the prereg revision log (v1.3). |

## 2. The detection is built on two post-gap PSF-fit excursions (the kill)

Per-event local depths at the pipeline ephemeris (own detrend: per-quarter median, transit-masked 1 d running median, gap-split; errors from the local out-of-transit scatter; `lc_checks.json: events`):

| event (BKJD) | quarter / channel | time after gap end | depth (ppm) |
|---|---|---|---|
| 231.44 | Q2 / ch 17 | 0.045 d (gap 230.29-231.39) | **3251 +- 607** |
| 937.83 | Q10 / ch 17 | 0.027 d (gap 936.95-937.80) | **3568 +- 229** |
| other 29 events | all | 0.65-29 d | weighted mean 275 ppm, chi2 37.0 / 28 dof; unweighted 268 +- 52 |

With the two post-gap events the per-event chi2 about a constant is 258 / 30 (the pipeline's own `event_chi2_const` is 67 / 30 and its `event_weighted_depth_ppm` is 250 vs a box depth of 481: the inconsistency was already in `vet.json`).

What the raw data do at those two epochs (`event_rawflux.log`; all in ppm of the quarter median, 3 h bins):

* **BKJD 231.48 (Q2).** Target PSF `FLUX` starts at -7830 ppm and recovers exponentially to +2700 ppm over ~2 d. Target `SAP_FLUX` does the opposite: **+18,590 ppm** decaying to +1,700. The TPF column 0 (the only column the target's PSF wing reaches; it is also the whole KBONUS SAP aperture, `APERTURE_Q2`) drifts from -3,261 to -21,954 ppm over the same 2 d, while the host's columns 2-4 drift only -450 to -1,430 ppm; `POS_CORR1` drifts monotonically from -0.023 to -0.042 px and the flux-weighted column centroid of the stamp moves from 2.348 to 2.338 px. This is a post-downlink focus/pointing settling that changes how much host light sits in column 0, and the PSF model for a star centred *off the stamp* answers it with an 8 ppt dip-and-recovery.
* **BKJD 937.85 (Q10).** Same mechanism with opposite sign: target PSF -1,523 -> +2,500 ppm step across the gap; `SAP_FLUX` +11,152 ppm decaying to 0; column 0 +22,938 -> +10,000 ppm; `POS_CORR1` 0.040 -> 0.030 px.
* A real loss of light cannot make the SAP flux go *up* by 1-2 % while the PSF flux goes down. Quarter box depths: Q2 PSF +1459 +- 135 ppm vs SAP **-1067 +- 147**; Q10 PSF +1040 +- 134 vs SAP **-992 +- 147** (`per_quarter`). The pipeline's `gap_frac_1d = 0.16` counts these as 2 of 31 events, but they carry roughly half of the depth.
* KBONUS's own `PSF_FLAT_FLUX` removes most of the Q2 excursion (local depth 337 +- 256 ppm) but not the Q10 step (2673 +- 274); its box depth at the ephemeris is 360 +- 52 ppm (`psfflat_check.log`).

Post-gap behaviour is generic in this light curve, not specific to these two gaps: the mean detrended PSF flux in the first 0.05 d after the 66 gaps > 0.3 d is -136 +- 77 ppm, and single gaps give first-0.15-d means of -3,909 (231.39), -2,438 (937.80), -1,962 (970.11), -678 (874.78), -609 (599.22) ppm (`per_gap_first0p15d`). Among 1,500 random ephemerides at the same period, 29 % land at least one event and 4.7 % at least two events inside a post-gap 0.15 d window; the two-hit subset has a box-depth scatter of 122 ppm with tails to ~430 ppm, against 52 ppm for ephemerides with no post-gap event (`random_ephemeris_calibration`). A blind search over ~1e5 trial periods therefore has many ephemerides available that convert post-gap ramps into 300-500 ppm "transits". The pipeline's own first screen reports `max_single_event_frac = 0.37` for this star.

## 3. What is left after excising 0.3 d after every gap (899 cadences, 2.2 % of the data)

| statistic | all data | post-gap 0.3 d excised |
|---|---|---|
| 1.95 h box depth at the pipeline ephemeris, PSF (error = scatter of 400 random-phase boxes) | 514 +- 50 ppm (SNR 10.2) | **272 +- 48 ppm (SNR 5.7)**; 286 +- 52 if 0.3 d is also cut before each gap |
| same, SAP_FLUX | -8 +- 56 ppm | 276 +- 55 ppm (SNR 5.0) |
| BLS 20-40 d, 1.5-4 h, own SDE: top peak (real) | 29.4318 d, SDE 6.65, 412 ppm, 3.0 h | 30.89 d, SDE 5.50; **SDE at 29.43 d = 3.59** |
| BLS 20-40 d inverted flux (2 - f): top peak | 32.99 d, SDE 5.69 | 36.78 d, SDE 5.94 |
| odd / even (1.95 h box) | 363 / 640 ppm (3.1 sigma) | 362 / 185 ppm (1.9 sigma) |
| secondary at phase 0.5 / scan over 200 phases | - | 47 ppm; scan scatter 53 ppm; deepest other phase 163 ppm at 0.40; primary 272 |
| per season, PSF (ch): S0 (17) / S1 (57) / S2 (65, Q12 only) / S3 (25) | 907 +- 76 / 298 +- 86 / 335 +- 180 / 326 +- 99 | 168 +- 74 / 299 +- 94 / 334 +- 181 / 325 +- 92 |
| per season, SAP | -665 +- 89 / 326 +- 99 / 318 +- 231 / 329 +- 124 | - |

After excision the three clean seasons agree on ~300-330 ppm and season 0 (the channel with both artefacts) gives 168 +- 74. The residual is consistent across PSF, SAP and PSF_FLAT and 0 of 1,500 random phases at the fixed period reach 272 ppm, so a small periodic dip at this ephemeris is not excluded. But (i) the ephemeris was selected by a search in which the two artefacts supplied most of the power, (ii) the 29.43 d peak is not recoverable blind in 20-40 d once the post-gap cadences are removed (SDE 3.6 < inverted-null top 5.9; the re-derivation independently found that 29.43 d is not the top blind peak even with the artefacts in and that the peak collapses to SDE 3.5-4.0 under 0.5 d / 2 d detrending), and (iii) 5.7 < 7.1, the registered stage-2 empirical-SNR floor. Under the registered rules the star fails stage 2 once the artefact cadences are removed.

## 4. Pixel level: localisation is not possible for this star

The target's Gaia position maps to x = -0.22 to -0.39 px, y = 1.2-2.1 px in the host's 6x5 (6x4 in season 3) stamp in all 11 quarters, i.e. its centre is outside the collected pixels; the host sits at x = 3.0 (`pixel_checks.json: tpfs`). The KBONUS SAP aperture is 3-4 pixels of column 0 only; `PSFFRAC` is 0.52-0.69. Season-stacked event-differenced deficit maps (own code, per-event local linear trend per pixel, inverse-variance stacking, two post-gap events excluded): deficit SNR in the pixel nearest the target (row 2, column 0) is 2.7 (S0, 8 events), 2.5 (S1, 10), 1.2 (S2, 3), 0.75 (S3, 8); the maximum |SNR| in any pixel in any season is 2.7. A shifted-PRF template at the target position improves chi2 over the null by 3.1 (S0), 10.5 (S1), 3.7 (S2), 2.6 (S3), but in S1 a template at the host fits equally well (chi2 36.8 vs 36.7, A = 3.5 +- 1.1 e/s) and the grid-search optimum lands near the host (2.75, 2.0), in S2 in the opposite corner (4.5, 3.25) and in S3 off the stamp (-1, 1). This neither localises a ~3 e/s dip (272 ppm x 11,713 e/s) nor can it: the expected per-season SNR for the target wing in column 0 is ~2-3 at best.

The host as the source is excluded directly: KIC 5866192's own KBONUS light curve folded at the ephemeris gives 19 +- 35 ppm (PSF) and 14 +- 25 ppm (SAP) over 83 in-transit cadences in Q1, 3, 5, 6, 8, 9, 11, 12, 14 (it has no Q2/Q10 data in the KBONUS product), i.e. an absolute loss < 2.1 e/s at 2 sigma against the 3.2-6.0 e/s that the target's PSF flux loses; and only a fraction of host light enters the target's PSF fit. So the two artefact events are not host eclipses either; they are the PSF model responding to the host's PRF moving in column 0.

## 5. Other checks (none kill, none rescue)

* Lomb-Scargle 0.2-100 d on per-quarter-normalised raw flux, both polarities: best 53.3 d (533 ppm, long-term trend), power near P 0.004, near P/2 0.02: no rotation alias.
* Quality flags: 8.7 % of cadences flagged in transit vs 8.7 % overall; no association. 
* Duration / density: for R* = 0.95 Rsun, M* ~ 1.05 Msun (Gaia M_G 3.8, BP-RP 0.81, RUWE 0.99, dwarf), a/R* = 43 and a central transit lasts 5.2 h; the measured 1.95 h (pipeline) or 3.0 h (re-derivation) need b = 0.93 or 0.82. Allowed but a priori unlikely, and a near-grazing geometry is exactly where a V-shaped blend is hard to exclude. Implied radius if on target: 2.3 Re at 481 ppm, 1.7 Re at 272 ppm, with a season-dependent dilution of the PSF flux that this product cannot calibrate (the star is off the stamp).
* Odd/even after excision 362 vs 185 ppm (1.9 sigma) and secondary 47 +- 53 ppm: an EB at 2P with unequal eclipses is not excluded at this SNR.
* Agreement with the re-derivation (`rederive/rederive_summary.json`): blind BLS does not return 29.433 d as the top peak; the peak is detrending-dependent; the two deep events are post-gap; depth without them 185-271 ppm at SNR 3.9-4.8. My additional findings are the SAP/PSF anti-correlation at those epochs (which makes them non-astrophysical rather than merely suspicious), the column-0 / `POS_CORR` anatomy, the random-ephemeris calibration, the off-stamp geometry, and the host exclusion.

## 6. What would settle it

1. Nothing in the Kepler data can raise a 272 +- 48 ppm box at this ephemeris above the 7.1 floor; the residual is at the noise limit of this product. Retire the signal from the survivor list and log it as a post-gap PSF-fit false positive (class: edge-of-stamp background star, channel 17).
2. For the pipeline (post-hoc, to be logged as v1.8 before any rerun): (a) excise 0.3 d after every gap before computing stage-2 depth, or require the depth with and without post-gap events to agree within 2 sigma; (b) require sign agreement of the per-quarter PSF and SAP depths; (c) flag stars whose Gaia position falls outside the stamp (`PSFFRAC` < 0.7 or off-stamp WCS position) as non-localisable. The 24-40 d class (349 of 438 first-screen passes) suggests these cuts will matter for many stars.
3. If anyone wants to pursue the residual: a season-by-season PSF re-extraction on the TPF with the target's position fixed from Gaia and a free focus term, then a fold; and a check of the Q12 and Q17 cadences individually. I do not recommend spending the effort.
