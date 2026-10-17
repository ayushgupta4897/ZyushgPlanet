# Red-team review: kic-006515744 (first-screen signal) — verdict: KILL

Date: 2026-10-17. Run id: poolA_main. Reviewer role: adversarial (independent of the stage-1/2 code; nothing under `pipeline/` imported).
Target: KIC 6515744 = Gaia DR3 2101870206266268416, G = 16.83, BP-RP 0.95, Teff 5222 K, Rstar 0.85 Rsun, d = 3.1 kpc, RUWE 1.01, TPF host KIC 6515734. Quarters with PSF data: 0,1,2,3,6,10,11,14 (channels 26 / 18 / 58 by season).

Pipeline claim under test (vet.json, rules v1.5): P = 36.2075 d, t0 = BKJD 127.9205, box depth 12 763 ppm (pass-2 masked re-detrend), duration 0.75 h, snr_emp 13.8, TLS SDE 5.4. Note that the same vet.json already records `stage2_pass: false` under the current v1.7 rules (`stage2_fails: n_events`) and `stage2_pass_v14: false` (`tls_sde`); the signal only ever passed under v1.5.

## One-line explanation

The 12 763 ppm "transit" is the average of a handful of isolated single-cadence negative excursions in the PSF flux. The two cadences that carry the depth are both Kepler/PSF-fit artefacts demonstrable in the originating target pixel file: BKJD 200.33 is the first cadence after the Q2 attitude tweak (a 27-sigma pointing excursion, dipole difference image, seen as a field-wide outlier in other KBONUS sources), and BKJD 1370.09 is a localised pixel *brightening* (+89 sigma) next to the target that the PSF fitter converted into a spurious 3.5 percent dimming while the aperture flux stayed flat. No pixel in the TPF gets fainter at that cadence. Nothing periodic remains once single-cadence outliers are treated symmetrically.

## Checks run and results

All numbers below are from the scripts and JSON in this directory unless cited to `vet.json` or `rederive/rederive_results.json`:
`redteam_checks.py` -> `redteam_checks.json` (light-curve checks, cross-source cadence test, ephemeris matching), `redteam_tpf_check.py` -> `redteam_tpf_check.json` (TPF single-cadence difference images for the two dominant cadences), `redteam_extra.py` -> `redteam_extra.json` (TPF check of the secondary contributing cadences, control for the cross-source test, per-quarter tails).

### 1. What the depth is made of (PSF flux, own 1-d running-median detrend, SAP_QUALITY == 0)
Per-epoch depths at the pipeline ephemeris, +/- 1 cadence (0.49 h) window (`redteam_checks.json: events.pipeline`):

| epoch | BKJD | Q | PSF depth ppm | SAP depth ppm | most negative cadence ppm | gap dist d |
|---|---|---|---|---|---|---|
| 0 | 127.92 | 0 | -344 | +262 | -287 | 2.3 |
| 2 | 200.34 | 2 | **35 716** | 13 850 | **-66 633** | 16.6 |
| 3 | 236.54 | 2 | 4 464 | 1 714 | -5 451 | 5.1 |
| 4 | 272.75 | 3 | 9 161 | 571 | -13 695 | 7.8 |
| 5 | 308.96 | 3 | -7 088 | 210 | +2 726 | 12.5 |
| 6 | 345.17 | 3 | 6 164 | 794 | -6 457 | 4.3 |
| 12 | 562.41 | 6 | 59 | 174 | -941 | 4.1 |
| 22 | 924.49 | 10 | -1 053 | -317 | +1 053 | 9.5 |
| 23 | 960.69 | 10 | 484 | 160 | -484 | 8.5 |
| 24 | 996.90 | 10 | -611 | -156 | +454 | 3.4 |
| 25 | 1033.11 | 11 | **26 890** | 471 | **-28 326** | **0.05** |
| 26 | 1069.32 | 11 | **20 389** | -178 | **-27 221** | 0.67 |
| 32 | 1286.56 | 14 | -428 | -316 | +50 | 2.9 |
| 33 | 1322.77 | 14 | -536 | -765 | -1 301 | 13.5 |
| 34 | 1358.98 | 14 | -880 | -154 | -305 | 12.3 |

Four epochs (200.34, 272.75, 1033.11, 1069.32) supply essentially all of the mean; ten of fifteen are consistent with zero; three are positive. Every "deep" epoch is one cadence wide. The rederivation (P = 33.42 d) is built on the same BKJD 200.33 cadence plus BKJD 1370.09 (`rederive_results.json`). The pipeline's own `events.csv` has a single usable event at 357 +/- 585 ppm and `event_dof = 0`, i.e. stage 2 had no event-consistency test to apply.

### 2. The dominant cadence BKJD 200.33 (Q2, channel 18, cadence 4473) is the first cadence after the Q2 attitude tweak
From the originating TPF `kplr006515734-2009259160929_lpd-targ.fits.gz` (`redteam_tpf_check.json: Q2_cad4473`):
- QUALITY of the preceding cadence 4472 = 98305 = NoData + NoFinePoint + **AttitudeTweak**; cadence 4473 itself is flagged 0, so neither the Kepler default bitmask nor KBONUS `SAP_QUALITY` removes it.
- POS_CORR1/2 at cadence 4473 are -6.7 sigma / **-27.1 sigma** outliers relative to the surrounding 48 cadences: a pointing excursion.
- Single-cadence difference image is a **dipole**: rows 1-2 of the stamp +15 to +21 sigma, rows 3-4 -15 to -41 sigma (16 pixels < -3 sigma, 8 pixels > +3 sigma); whole-stamp flux change -0.02 percent. That is an image shift, not a flux loss. A PSF fit at a fixed position reads a shifted image as a dimming: PSF FLUX -66 633 ppm, PSF_FLUX_NOVA -62 325 ppm, SAP -27 064 ppm, and the three following cadences sit at -5 000 ppm while the pointing recovers (`redteam_checks.json: outlier_cadences.200.33.neighbours_rel_ppm`).
- Field-wide: at the same cadence number in the 323 other cached KBONUS sources with Q2 data, 2.2 percent are below -3 sigma (vs 0.6 percent above +3 sigma) and three unrelated sources on channels 66, 13 and 84 are at -37, -29 and -12 sigma. Control: over 40 random Q2 cadences in 147 sources, the maximum number of sources below -10 sigma is **0**; at cadence 4473 it is **4** (`redteam_extra.json: cross_source_control_Q2`).

### 3. The second dominant cadence BKJD 1370.09 (Q14, channel 18, cadence 61720) is a pixel brightening, not a dimming
From `kplr006515734-2012277125453_lpd-targ.fits.gz` (`redteam_tpf_check.json: Q14_cad61720`):
- QUALITY 0 for this and all neighbouring cadences; POS_CORR normal (0.3 / 0.6 sigma).
- Difference image: **no pixel below -3 sigma**; four pixels at column 0, rows 2-4 are at +89, +72, +23, +17 sigma (adjacent to the target at column 1.5, row 4.1). Whole-stamp flux +0.35 percent; target 3x3 box +0.14 percent.
- The two KBONUS PSF estimators disagree by 15 percent at this cadence (FLUX -35 529 ppm, PSF_FLUX_NOVA **+115 501 ppm**) while SAP is +127 ppm. A localised charge/cosmic-ray event next to the target perturbed the PSF fit. There is no flux deficit anywhere in the stamp, so this cadence cannot be a transit of anything.
- 1.23 d before the end of Q14.

### 4. The remaining contributors are Q11 (noisiest quarter) edge/noise cadences with no pixel-level deficit
`redteam_extra.json: other_event_cadences`, TPF `kplr006515734-2012004120508_lpd-targ.fits.gz` (Q11, channel 58) and Q3:
- BKJD 1033.063: first good cadence after an Earth-point (downlink) gap (preceding QUALITY 98312 = NoData+NoFinePoint+EarthPoint, 16392 = DetectorAnomaly+EarthPoint). Target 3x3 box +0.01 percent; no pixel beyond +/-3 sigma.
- BKJD 1033.145, 1069.29, 1069.31: target 3x3 box within +/-0.07 percent; no pixel below -3 sigma (one at -3.8 sigma at 1033.145); 1069.29 sits next to a CosmicRayCollat-flagged cadence. PSF flux at these cadences is -27 000 to -28 000 ppm: pure PSF-fit scatter in a quarter whose detrended PSF scatter is 7 600 ppm (`rederive_results.json: per_quarter_scatter_ppm`; Q3 5 035 ppm).
- BKJD 272.73/272.75 (Q3): no pixel beyond +/-2 sigma; target box -0.06 percent.
- The single most negative cadence in the whole light curve (BKJD 349.495, -79 954 ppm PSF, -5 022 ppm SAP) is the last good cadence of Q3 with 17 pixels < -3 sigma: a quarter-end roll artefact. It does not fall on either ephemeris but shows the light curve's failure mode.

### 5. Does anything periodic survive? (BLS 20-50 d, 60 000 log-spaced periods, durations 0.02-0.1 d, likelihood objective)
`redteam_checks.json: bls_*`
- All cadences (4-sigma positive clip only): best 49.72 d, SDE 12.7, duration 0.02 d (shortest trial), depth 18 432 ppm. Power at 36.2075 d: 183; at 33.42 d: 273.
- Drop just the two cadences 200.33 and 1370.09: best period moves to 43.45 d, SDE 11.6; power at 36.2075 d falls to 124 and at 33.42 d to 59.
- Symmetric treatment (drop the 480 cadences below -5 MAD, MAD = 2 045 ppm, the same clip the pipeline applies to positive outliers): best 33.51 d, **SDE 7.4**, depth 1 963 ppm; power at 36.2075 d (15.8) is *below* the spectrum median (16.0).
- Inverted flux (2 - f), same window, all cadences: best 20.78 d, **SDE 9.7**. The real-flux SDE is therefore of the same order as its own inverted null; both are driven by single-cadence tails.
- Pipeline first screen (`survivors.csv` row for this star): SDE 11.4, snr_emp 9.0, inverted-null SDE 8.4; rederivation blind BLS: 8.2 raw / 16.7 flattened vs inverted 7.2 / 7.8.

### 6. Stellar variability
Lomb-Scargle 0.1-100 d on the detrended PSF flux: best period 1.03 d, power 0.0012, FAP 0.22 (`redteam_checks.json: lombscargle`); vet.json: 0.88 d, 238 ppm half-amplitude. No rotation signal, no alias relation to 36.2 or 33.4 d. Not a variability alias.

### 7. Duration / density / implied radius
P = 36.2 d, M ~ 0.9 Msun, R = 0.854 Rsun -> a/R* = 52, central-transit T14 = 5.9 h (`redteam_checks.json: duration_physics`; vet.json t14_max 5.7 h). Measured box durations 0.48-1.5 h are 8-12x shorter: a 1 h transit would need b > 1.0 (formally 1.10), i.e. impossible for a flat-bottomed 1.3 percent box; the 0.48 h best duration is shorter than one 29.4-min cadence. A real 0.75 h event sampled at 29.4 min cannot be single-cadence-only at every epoch. Rp/R* = 0.113 (10.5 Re) if taken at face value; irrelevant given the above.

### 8. Odd/even, secondary, V-shape
Odd/even 2.3 sigma (vet.json) / 2.0 sigma (rederive) only because the odd and even bins each contain one artefact cadence; secondary -555 ppm (-0.6 sigma). With one cadence per "event" no shape test (V vs U) is possible. These tests are uninformative, not passing.

### 9. Known objects / ephemeris contamination
- KIC 6515744 is not a Kepler target (`kepler_target_kepids.csv` lists only the host 6515734) and has no KOI, DR25 TCE or Kepler EB entry (`data/known/*`, fetched 2026-10-17T15:46 UTC). SIMBAD cone 20 arcsec returns only UCAC4 660-072021 (the G = 12.3 host, ~10 arcsec away).
- Within 3 arcmin: KIC 6515722 (40.8 arcsec), KOI-3818.01, FALSE POSITIVE, P = 1.9085 d, with nine DR25 TCEs at 1.9-72 d (incl. 35.85 d, t0 138.20, depth 74 ppm and 34.26 d). Neither 35.85 nor 34.26 d matches 36.2075 d within tolerance and the epochs disagree; these are low-SNR multiples on a noisy star, not a contaminator. KIC 6515790 (128 arcsec) TCE at 343.8 d: no relation.
- Whole-field ephemeris match (KOI + DR25 TCE, period within 0.2 percent at x1/3, x1/2, x1, x2, x3 and epoch within 0.1 d): **none** for the pipeline ephemeris (96 period-only coincidences, all fail the epoch test); for the 33.42 d rederive ephemeris only KOI-353.03 / TCE 11566064.3 at 3 x 11.1619 d (KIC 11566064, 6.5 deg away, channel differs): a chance period triple, not a same-channel contaminant.
- Bienias & Szabo 2025 (A&A 694, A209) lists seven background hot-Jupiter candidates at P = 2.5-4.8 d; no overlap with a 36 d / 33 d signal. Moquin et al. 2026 BOI-1 is P = 3.37 d. CCD channels used: 26 (Q0-1), 18 (Q2, 6, 10, 14), 58 (Q3, 11); no channel-level KOI/EB catalogue was available offline, but the question is moot once the events are shown to be non-astrophysical in the TPF.

### 10. Existing pixel runs
`runs/poolA_main/pixel_v15` and `pixel_v16` for this fname: 0 of 8 quarters conclusive, "no usable events" in Q0-3 and Q14, single-event maps in Q6/Q10/Q11 with peak S/N 0.1-1.7 and in-aperture depth -30 to -144 ppm (i.e. no deficit). Consistent with there being nothing to localise.

## Verdict: KILL

Most likely explanation: PSF-photometry single-cadence artefacts — a Q2 post-attitude-tweak pointing excursion (BKJD 200.33), a localised pixel brightening/charge event beside the target in Q14 (BKJD 1370.09), a post-downlink edge cadence and PSF-fit scatter in the noisy season-1 quarter Q11 — aligned by a 36.2 d (or 33.4 d) period search over a 1 250 d baseline with duty-cycle-limited sampling. The KBONUS `SAP_QUALITY` mask does not remove the first cadence after an attitude tweak, and the first-screen positive-only clip leaves negative single-cadence tails in place.

## What would settle it (if anyone wanted to resurrect it)
1. Re-run the first screen on this star with (a) the Kepler QUALITY bitmask (including AttitudeTweak-adjacent cadences) applied from the TPF and (b) a symmetric 5-MAD single-cadence clip; the 36.2 d peak should vanish (this review's symmetric-clip BLS: SDE 7.4, 36.2 d below median).
2. Require >= 3 epochs with >= 2 in-window cadences each and per-epoch depths consistent at < 3 sigma before any signal with duration < 2 cadences advances (the v1.7 `n_events` rule already fails it).
3. Recommend a pipeline-level hygiene check: flag any first-screen survivor whose best BLS duration equals the shortest grid value and whose depth is carried by <= 2 cadences (single-event dominance > 0.5), and cross-check the carrying cadences against the host TPF QUALITY of cadence-1 and POS_CORR z-scores. Stated as a post-hoc suggestion; not a change to prereg thresholds.

Files written: `redteam.md`, `redteam_checks.py`, `redteam_checks.json`, `redteam_tpf_check.py`, `redteam_tpf_check.json`, `redteam_extra.py`, `redteam_extra.json` in `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006515744/`.
