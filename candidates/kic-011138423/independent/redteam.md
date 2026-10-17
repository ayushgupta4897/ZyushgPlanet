# kic-011138423: red-team review (2026-10-17)

**Role.** Try to kill the 4.8945 d signal on KIC 11138423 (Gaia DR3 2129669944725625344, G = 17.35, BP-RP = 1.72, 894 pc, early-M dwarf, R* = 0.534 Rsun; extracted by KBONUS from the TPF of the G = 12.5 star KIC 11138418, 11.4 arcsec away). Language: this is a *signal*; nothing below makes it a candidate or a planet.

**Verdict: survives** (no kill found), with one finding that must be added to the evidence pack before any further step: the KBONUS fractional depth of this star is season-dependent because the PSF flux carries a season-dependent amount of non-target light; the light lost per transit in e/s is constant. The radius statement must be given as a range (about 2.0-2.5 Re), not as 2.1 Re.

All my numbers trace to the scripts and JSON under `runs/poolA_main/adversarial/kic-011138423/redteam/` (SHA-256 in `sha256_manifest.txt`); nothing from `pipeline/` was imported. Inputs: the cached KBONUS file `data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-011138423_kepler_v1.0_lc.fits`, the 18 cached Kepler TPFs of KIC 11138418 (`~/.cache/lightkurve/.../kplr011138418-*_lpd-targ.fits.gz`), the KBONUS light curve of KIC 11138418 (downloaded to scratch), the catalogues in `data/known/` (fetched 2026-10-17T15:46Z per `manifest.json`), the KBONUS source catalogue, SIMBAD and the NASA Exoplanet Archive TAP (TOI, PS) queried 2026-10-17.

## 1. Is it already known?

| check | result |
|---|---|
| KOI / DR25 TCE / Kepler EB (Kirk+2016) / confirmed planets within 5 arcmin | none (0 rows in each catalogue) |
| KIC 11138423 or KIC 11138418 as KOI / TCE / EB | not present (11138418 is in `keplerstellar` only: Kp 12.55, Teff 6070 K, R 1.1 Rsun, 17 quarters, no TCE) |
| ephemeris match, any KOI/TCE/EB at P, P/2, 2P, P/3, 3P, 2P/3, 3P/2, P/4, 4P (|dP/P| < 2e-3), whole field (`redteam/ephemeris_match.py`, `.log`, `.json`) | 9 KOIs and 12 TCEs within 0.2 % of P, all > 177 arcmin away; none agrees in phase (their t0 lands at phases 0.14-0.89 of our ephemeris) and none has |dP/P| < 1e-4, which a cross-talk/column contaminant would need to stay coherent over 300 epochs. Across all harmonics, 203 KOI/TCE/EB rows match within 0.2 %; exactly one object is both period-coherent and phase-aligned: KOI-70.05 = TCE 6850504-05 (Kepler-20 f, P = 19.57761 d = 4P to 1e-5, phase 0.016, depth 98 ppm, 426 arcmin away, different module). It cannot be the source: a 4P contaminant would populate only every fourth of our 253 events (88 % of events are positive and the per-event depths are consistent with a constant), and 98 ppm on a star 7 degrees away cannot reach this pixel. With ~200 matches, |dP/P| < 1e-4 and a phase window of 0.08, about one chance hit is expected. The 13 Kepler-EB harmonic matches are all > 540 arcmin away and none is coherent. |
| TOI / PS within a 3 arcmin box (Exoplanet Archive TAP, 2026-10-17) | 0 rows |
| Bienias & Szabo 2025 (arXiv:2501.09152) background hot-Jupiter candidates | 7 candidates: hosts near KIC 4459924, 8175131, 8293539, 4543171, 9674230, 9040849; none in this field |
| Kepler Bonus team (Moquin+2026 poster) | only BOI-1.01 (G = 17.8, P = 3.37 d) is public; not this star |
| SIMBAD within 30 arcsec | one entry: UCAC4 694-066714 = the host KIC 11138418, otype `*`; nothing at the target position; within 3 arcmin the only variables are ZTF J192857.62+484416.6 (EB*, 167 arcsec away) and an HB*/RG* set > 120 arcsec away |
| survey-wide systematic period | among 6,365 stars in `runs/poolA_main/triage_all.csv`, this is the only real-polarity detection within 0.2 % of 4.8944 d (also none at 2P; one inverted-null detection of another star at P/2 = 2.449 d, unrelated) |

Nothing is known at this position; no catalogued object on any channel shares the ephemeris.

## 2. The bright TPF host (the obvious blend) is excluded directly

KIC 11138418 (G = 12.5, 120,748 e/s PSF flux) folded at the signal ephemeris with the same masked detrend: **-8.1 +- 8.9 ppm (PSF), -5.0 +- 3.7 ppm (SAP)**, 939 in-transit cadences. At 2 sigma the host loses at most ~20 ppm x 1.2e5 e/s = 2.4 e/s in total; the target's PSF flux loses 1.7-1.9 e/s per transit (section 4). For the host to be the source, ~80 % of its missing light would have to be attributed to the target's PSF model, while the measured host-PRF wing at the target's distance (2.9 px) is only 0.10-0.19 % of the host flux per pixel (`pixel_checks.json: host_wing_fraction_per_pixel_at_target_distance_median`). Dead.

## 3. Other resolved neighbours (Gaia, G < 20)

Nearest sources to the target other than the host (positions from `pixel.json`, Q5 frame): G = 19.19 at 3.6 px (14 arcsec), G = 16.79 at 4.3 px (17 arcsec), G = 19.85 at 5.4 px, G = 18.69 at 5.6 px, G = 19.29 at 5.7 px, KIC 11138415 (G = 18.56, 29 arcsec; the pixel_v16 "strongest alternative" 2129669944725610624 is this star). With the measured host wing (0.1-0.2 % per pixel at 3 px) the leak of any of them into the target's PSF-fit flux is < 0.5 %; the eclipse depth each would need to produce 1.88 e/s is >= 150 % for every G > 18.5 source and 17-41 % for the G = 16.79 star. The 17-41 % case would be a deep EB 4.3 px NW of the target whose deficit would peak in the stamp's left column rows 2-4; the season stacks (section 5) peak in the target's own pixel in every season and show nothing there. None of these stars has a KBONUS light curve (404 from MAST for all six Gaia ids; the KBONUS catalogue lists only 11138415, 11138418 and 11138423 on this TPF).

## 4. The finding: season-dependent dilution of the KBONUS PSF flux

Per-quarter box depths (own detrend: per-quarter median normalisation, transit-masked 1 d running median, 4 sigma / -10 sigma clip; error = scatter of 300 random-phase boxes), grouped by Kepler season (`raw_quarter_flux.log`, `season_table.log`):

| season (quarters) | channel | PSFFRAC | raw PSF flux (e/s) | PSF depth (ppm) | absolute dip (e/s) | raw SAP flux | SAP absolute dip (e/s) |
|---|---|---|---|---|---|---|---|
| S0 (2, 6, 10, 14) | 68 | 0.61-0.64 | 1789 | **900 +- 105** | 1.61 +- 0.19 | 2740 | 1.24 +- 0.25 |
| S1 (3, 11, 15) | 28 | 0.79-0.81 | 1417 | 1243 +- 175 | 1.76 +- 0.25 | 1796 | 1.46 +- 0.23 |
| S2 (4, 8, 12, 16) | 20 | 0.91-0.92 | 1017 | 1548 +- 213 | 1.56 +- 0.21 | 4173 | 1.59 +- 0.28 |
| S3 (0, 1, 5, 9, 13, 17) | 60 | 0.87-0.92 | 1235 | **1698 +- 180** | 2.03 +- 0.21 | 2204 | 1.95 +- 0.22 |

* Fractional PSF depth: S0 vs S3 differ by 798 +- 209 ppm (3.8 sigma); per-quarter chi2 about a constant = 38.6 for 16 dof. Calibration with 60 random ephemerides (same P, random t0) on the same data gives median 13.7, 95th percentile 21.8, max 26.5: the real spread is **not** noise. The slope of depth against PSFFRAC is 2312 +- 601 ppm per unit (3.8 sigma).
* The stitched `FLUX` cannot show this: KBONUS rescales every quarter to the same mean (all 17 stitched medians are 1357 e/s). The un-rescaled per-quarter `LIGHTCURVE_Qn.FLUX` medians vary from 1017 (S2) to 1789 e/s (S0), i.e. the PSF photometry of this star contains between ~0 and ~800 e/s of light that is not the star, depending on the season's geometry (the target sits in the stamp corner; in S0 only 61 % of its PRF is inside).
* The absolute light lost per transit is constant: raw-PSF dip 1.73 e/s mean, chi2 14.4/16; raw-SAP dip 1.58 e/s, chi2 15.1/16; `PSF_FLUX_NOVA` gives the same (S0-S3 difference -0.17 +- 0.33 e/s). Pixel-level full-PRF amplitude at the target position (section 5): 2.44 +- 0.37, 1.37 +- 0.26, 1.60 +- 0.35, 2.13 +- 0.31 e/s for S0-S3 (weighted mean 1.80 e/s, chi2 7.2/3).
* Consequence: the pipeline's 1316 +- 79 ppm (vet.json) and the re-derivation's 1389 +- 95 ppm are depths relative to a season-averaged, contaminated flux. The true depth is 1.7-1.8 e/s divided by the star's own flux. The host-derived zero point (Kp 12.554 at 120,748 e/s -> 25.26) predicts ~1450-1470 e/s for Kp 17.34 / G 17.35, giving ~1200 ppm and Rp ~ 2.0 Re; the least-contaminated season (S2, 1017 e/s) gives ~1700 ppm and Rp ~ 2.45 Re. **State Rp as 2.0-2.5 Re with this caveat; do not quote 2.1 Re with a 0.1 Re error.** The rederivation's monotonic per-year rise (1112 -> 1609 ppm) is a mixture of this season effect and noise, not evidence of a changing signal.
* This also explains the rederivation's SAP/PSF depth ratio of 0.44: both photometries are diluted by the host's wings by different amounts, and `CROWDSAP` (0.90-0.998) does not describe the PSF flux. Agreed with the re-derivation: do not use CROWDSAP for anything here.

This is not a kill: a blended source would leak a *season-dependent* amount of light into the aperture and PSF fit (the geometry changes every quarter), yet the e/s dip is the same on four different CCD channels and stamp orientations, in PSF, SAP and NOVA flux alike. That is the behaviour of a source whose light is wholly inside the photometry in every season, i.e. the target or something unresolved from it.

## 5. Pixel-level localisation, independently (season-stacked, empirical PRF)

Own analysis of the 18 TPFs (`pixel_checks.py`): per event, in-transit mean image minus the median of out-of-transit cadences within +-1 d; deficits stacked per season (same pixel grid within a season: 58-71 events per stack); template = the out-of-transit image itself (dominated by the host's PRF) shifted to a trial position, amplitude fitted by weighted least squares.

| season | n events | chi2 null | chi2 PRF at target | chi2 PRF at host | best position | offset from target | sky offset (E, N) |
|---|---|---|---|---|---|---|---|
| S3 | 68 | 85.5 | 37.6 | 78.2 | (1.25, 0.75) px | 0.45 px = 1.8" | (+1.7", +1.8") |
| S0 | 71 | 156.6 | 113.8 | 147.7 | (-0.50, 0.50) | 0.74 px = 2.9" | (+1.0", +2.9") |
| S1 | 71 | 61.1 | 34.3 | 59.6 | (0.50, 0.00) | 0.42 px = 1.7" | (-0.7", +1.6") |
| S2 | 58 | 56.6 | 36.0 | 47.6 | (1.25, 1.00) | 0.53 px = 2.1" | (-0.6", -1.9") |

* Target pixel is the strongest pixel of the SNR map in every season (SNR 7.7, 7.4, 4.8, 5.7); a PRF at the target beats the null by delta chi2 = 48, 43, 27, 21 and beats a PRF at the host by 41, 34, 25, 12.
* The stamp-edge worry (the target is 0.07-0.85 px from the edge, so a source *beyond* the edge is the natural alternative): moving the PRF outward along the host-to-target line by +0.5 / +1.0 px raises chi2 by +23 / +33 (S3), +4 / +9 (S1), +6 / +11 (S2). Only S0 prefers +1 px outward (delta chi2 -14), and S0 is the season with chi2_null = 157 for 34 pixels (systematic residuals across the whole stamp; also the season with PERTSTD 0.04 in KBONUS) and the target at 0.07 px from the edge; its localisation is unreliable in the outward direction.
* The best-fit sky offsets point in different directions in different seasons (three north, one south; mean (+0.3", +1.1"), scatter 1.0" / 1.8"); a real offset source would sit in the same sky direction every season. Localisation: within ~2-3" (< 1 px) of the Gaia position of the target, consistent with the pipeline's `pixel_v16` on-target verdict (delta chi2 197 vs null, 140 vs best alternative over 13 quarters) and with `pixel.json` (median target-to-deficit 3.1").
* What this cannot do: exclude a companion or background EB within ~2-3" of the target (sub-pixel blend). That is left to statistical validation.

## 6. Light-curve kill attempts (own code, `lc_checks.log`)

| test | result |
|---|---|
| odd/even (pipeline t0 parity) | 1413 +- 134 vs 1349 +- 127 ppm: 0.35 sigma (pipeline 0.78, rederivation 0.60) |
| secondary at phase 0.5 | +51 +- 106 ppm (1.95 h box); -11 ppm with a 3x-duration box. 3 sigma limit ~370 ppm, i.e. < 27 % of the primary |
| box scan over all phases (181 trials) | largest dip elsewhere 223 ppm at phase 0.27 (2.1 sigma), largest bump -360 ppm at 0.95 (3.4 sigma, adjacent to transit egress window; not at 0.5). No second eclipse at any phase |
| Lomb-Scargle on raw normalised flux, 0.25-60 d | top peaks at 40.6 d (quarter-length systematics); power at P = 1.8e-4 (sinusoidal amplitude ~100 ppm at P, 7 % of the depth); pipeline LS period 0.607 d, amp 84 ppm. No rotation alias at P, P/2 or 2P. Inverted-flux null, identical BLS (re-derivation): SDE 8.2 vs 31.3 raw |
| per-event depths (253 events) | chi2 = 349 / 252 dof (vet.json 347); random-ephemeris calibration on the same light curve gives chi2/dof median 1.27, p84 1.44, max 1.76; the real 1.38 is at the 73rd percentile. The excess is red noise, not variable depth. Strongest event carries 3.8 % of the summed SNR^2; 88 % of events positive |
| gap proximity | 24 events within 1 d of a gap edge: 1114 +- 215 ppm; 229 others: 1239 +- 67 ppm. Gap-edge events are not driving the depth |
| per-year t0 (linear ephemeris) | best t0 shifts of +1.4, +4.3, -7.2, +1.4 min for the four years: no TTV or period drift at the cadence resolution |
| quality flags in transit | 7.6 % flagged in transit vs 8.0 % overall (bits present: 16, 128, 2048, 8192, 131072, 262144: reaction-wheel zero crossing, cosmic rays, impulsive outlier, rolling band in collateral/aperture). No excess of flagged cadences at the transit phase; flagged cadences were excluded anyway |
| PSF-fit red chi2 in vs out of transit | median ratio-1 = +0.0014 across quarters: the PSF model is not worse in transit (no in-transit blob/outlier population) |
| SAP background | `SAP_BKG` is NaN in this product; not testable here. The per-pixel analysis above uses the raw TPF flux and is unaffected |

## 7. Duration, density, radius

For R* = 0.534 Rsun and M* = 0.50-0.60 Msun: a/R* = 18.0-19.2, rho* = 3.3-3.9 rho_sun, central-transit T14 = 1.95-2.07 h, T14(b = 0.5) = 1.7-1.8 h. Measured: 1.95 h (pipeline box), 1.76 h (re-derivation box), 2.04 h (trapezoid T14). Consistent with a transit across this dwarf at b < ~0.5. T23/T14 = 0.72 (trapezoid) is cadence-limited at 30 min and is not evidence of a V shape; for Rp/R* = 0.04 the expected value is 0.92. Rp: 2.0-2.5 Re for depth 1200-1770 ppm (section 4). Gaia: M_G = 7.60 at BP-RP = 1.72 is on the main sequence; RUWE = 1.06 (no astrometric-binary signal). The signal is physically self-consistent as a sub-Neptune transit of this star.

## 8. What survives, what does not, what is left

**Killed hypotheses**
1. Signal on the bright TPF host KIC 11138418: < 20 ppm on the host; impossible by three orders of magnitude.
2. Any resolved Gaia neighbour (G < 20) within 30 arcsec: required eclipse depths >= 17 % (G = 16.8 star, 17" away) to > 100 % (all others), and the pixel deficit sits in the target's pixel, not at any of them.
3. Catalogued contaminant via ephemeris match (any harmonic, any distance): none coherent.
4. Rotation / variability alias: no LS power at P or harmonics.
5. Kepler systematics: no survey-wide period cluster at 4.894 d, no gap/quality-flag/rolling-band association, no per-event depth excess beyond red noise, stable linear ephemeris over four years, no in-transit PSF-fit degradation.
6. Already known: not in KOI, TCE, Kepler EB, PS, TOI, SIMBAD, Bienias & Szabo 2025 or the Kepler Bonus team's public candidate.

**Surviving caveats (not kills; standard items for the next stage)**
1. Season-dependent dilution of the KBONUS PSF flux (section 4): the ppm depth and the radius need the correction described; the evidence pack must carry the per-quarter/per-season table and the absolute-e/s constancy.
2. Sub-pixel blend (companion or background EB within ~2-3"): untested by any pixel data here; needs TRICERATOPS with the Gaia/pixel constraints and, ideally, high-resolution imaging (G = 17.35 makes AO hard; speckle at 4-m class is marginal).
3. Twin EB at 2P = 9.789 d with equal eclipses: odd/even agree to 0.35 sigma, which leaves room for a near-identical pair; at 9.79 d around a 0.55 Msun primary the expected T14 (2.5 h) is slightly longer than measured but not excluded. Statistical-validation territory.
4. On-target grazing stellar companion: the 3 sigma secondary limit (< 27 % of primary) and the 30-min cadence do not exclude a low-mass eclipsing companion with a cool secondary; same remedy.

**What would settle it**
* Re-extract the star with a PSF fit that includes the full host PRF and report depth per season; or simply adopt the absolute-e/s dip and an independent flux estimate for the star (Gaia G with the host-derived zero point) as the depth definition.
* TRICERATOPS run with the pixel-level 2-3" exclusion and the G < 20 neighbour list; report FPP/NFPP.
* Short-cadence or ground-based (e.g. 1-2 m class, red filter, 1.7 ppt on G = 17.35 is hard but possible for a 2 h event) transit to resolve ingress/egress shape.
* Speckle/AO imaging of the target to bound companions at 0.1-2".

## Files
* `runs/poolA_main/adversarial/kic-011138423/redteam/lc_checks.py`, `lc_checks.json`, `lc_checks.log` (per-flux-column depths, per-event chi2 and its random-ephemeris calibration, phase scan, odd/even, quality flags, LS)
* `.../redteam/season_table.json|.log` (per-quarter and per-season depths from the stitched columns; random-ephemeris calibration of the per-quarter chi2)
* `.../redteam/raw_quarter_flux.json|.log` (un-rescaled per-quarter PSF/SAP/NOVA fluxes and absolute e/s dips)
* `.../redteam/pixel_checks.py`, `pixel_checks.json|.log` (season-stacked deficit maps, empirical-PRF position fits, line scans beyond the stamp edge, host-wing fractions)
* `.../redteam/sky_offsets.json|.log` (best-fit deficit positions in sky coordinates per season)
* `.../redteam/sha256_manifest.txt`
* Scratch (not part of the pack): `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic011138423/` (includes the downloaded KBONUS light curve of KIC 11138418 and the TAP query outputs `dl/toi_near.csv`, `dl/ps_near.csv`)
