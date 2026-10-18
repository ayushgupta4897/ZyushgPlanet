# Red-team review: kic-005609180 (KBONUS-BKG PSF light curve, 21.106 d first-screen signal)

Date 2026-10-17. Role: adversarial review whose job is to kill the signal. Candidate language only; nothing here is a planet.
Inputs read: `runs/poolA_main/vet/kic-005609180/{vet.json,fold_events.png,events.csv}` (stage-2 pass under rules v1.7, failed v1.4 on `tls_sde`; pixel directory absent), the independent re-derivation under `runs/poolA_main/adversarial/kic-005609180/rederive/`, the cached FITS file `data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005609180_kepler_v1.0_lc.fits`, `targets/poolA_kbonus_bg_dwarfs_ranked.csv` row 4824, and `data/known/` (manifest fetched 2026-10-17T15:46 UTC).
All numbers below trace to files in `runs/poolA_main/adversarial/kic-005609180/redteam/` (produced by the inline scripts recorded in this session; numpy/astropy only, `use_threads`-free, no `pipeline/` code imported).

## Verdict: KILL

The 21.106 d signal is a coincidence of three post-gap PSF-model recovery ramps in the KBONUS stitched PSF flux (epochs 4, 36, 39), all three in season-0 quarters on channel 17 (module 7, output 1). With the first day after every data gap masked, the period is not detectable at all (BLS SDE 1.7, rank > 40 in a 15-30 d window) while the identical inverted-flux null reaches SDE 8.3. The remainder, once events within the ~2 d ramp timescale are excluded, is a 3.1-3.3 sigma PSF-only residual (201-220 +- 65 ppm) that is 1.5-1.8 sigma in SAP flux, has a per-event chi2 of 57-59 for 45-46 dof, is absent in season 2, is carried by three events, and is matched by an equal-amplitude brightening 3 h earlier in a phase scan. That is what a search-selected period looks like in red noise, not a transit.

## 1. The three dominant events are post-gap ramps, not flux deficits of the star

Gap classification (`redteam/ramp_census.txt`): epoch 4 (231.51 BKJD) is 0.113 d after a 1.10 d gap in Q2; epoch 36 (906.91) is 0.061 d after the Q9/Q10 boundary gap; epoch 39 (970.23) is 0.138 d after a 0.86 d downlink gap in Q10. In events.csv these three carry depths of 4.27, 2.42 and 4.04 ppt with `gap_dist_d` 0.11, 0.06, 0.14; the next-deepest event is 1.9 ppt.

Cadence-level comparison of the flux columns in the FITS file (`redteam/stitched_vs_quarter.txt`, all per-quarter-median normalised):

| epoch | stitched FLUX (searched) | per-quarter FLUX | PSF_FLUX_NOVA | SAP_FLUX |
|---|---|---|---|---|
| 4, first cadences after gap | -43 ppt | -32 ppt | +17 ppt | +7 ppt |
| 36 | -14 ppt | -4 ppt | -28 ppt | -6.6 ppt |
| 39 | -47 ppt | +29 to +68 ppt | +16 ppt | +8 ppt |

A flux deficit of the star would appear with the same sign and, after dilution, comparable amplitude in every column. Here the four columns disagree in sign and by up to 100 ppt. The stitched PSF flux the pipeline searched also differs from the per-quarter PSF flux by a MAD of 85-96 ppt (max 180 ppt) in Q6, Q10 and Q14, versus 0.3-5 ppt in most other quarters: the three quarters with the worst disagreement are all season 0 / channel 17, as are the quarters of all three dominant events. RED_CHI2 of the PSF fit changes by -17000 and -6400 at epochs 4 and 39 and the SAP centroid moves 0.02 px (`redteam/quarter_diag.txt`, `redteam/quarter_diag.png`): the PSF model is responding to a post-gap focus/thermal change, which is the known behaviour after Kepler monthly downlinks and quarter rolls.

Census across all 85 gaps > 0.3 d (`redteam/ramp_census.json`): 11 gaps produce a stitched-PSF ramp deeper than -5 ppt in the first 0.15 d, 7 deeper than -10 ppt. Eight of the eleven are in season-0 quarters (channel 17), three in season 1 (channel 57), none in seasons 2-3. The three largest (-38, -38, -25 ppt at 231.39, 970.09, 1337.20 BKJD) all have SAP ramps of opposite sign (+6.7, +7.2, +6.0 ppt). The ramps are therefore a generic PSF-model artefact of this light curve, strongest on channel 17, and not specific to the 21.1 d ephemeris.

## 2. The period is set by the ramps

Of the 11 big ramps, exactly three (231.39, 906.85, 970.09) fall within 0.15 d of a 21.106 d epoch (phase offsets -0.005, -0.003, -0.007). Scanning P = 1-50 d for periods that phase all three gap ends within 0.15 d gives a comb of 70 aliases, of which only 10.552 d and 21.103 d are longer than 10 d (`redteam/residual_robustness.txt`). The BLS period is therefore fixed by the gap-timing comb; 21.106 d wins over 10.55 d because the latter puts half its epochs on ordinary data.

Confirmation by masking (`redteam/masked_bls.txt`, 15-30 d BLS, 6000 periods, durations 2-6 h, likelihood objective, running-median detrend per segment, positive outliers clipped):

| data | top peak | at 21.106 d |
|---|---|---|
| PSF, all cadences | 21.1069 d, SDE 8.56, 429 ppm, 4.0 h (rank 1) | same |
| PSF, 1.0 d after every gap masked | 28.53 d, SDE 5.86 | SDE 1.66, rank > 40 |
| PSF inverted (2 - f), same mask | 23.16 d, SDE 8.31 | SDE 1.91 |
| SAP, same mask | 15.93 d, SDE 7.43 | SDE 2.34, rank > 40 |

With post-gap cadences removed the real periodogram's best peak (5.9) is below what the inverted-flux null reaches (8.3) in the same window; the 21.1 d peak is gone. This agrees with the re-derivation's 0.5 d edge-trim result (local SDE 4.9, 30 peaks above) and is stronger because 0.5 d leaves most of a ~2 d ramp in place.

## 3. What remains after removing ramp-affected events is not a transit

Own per-event depths (3.0 h box, local linear baseline over +-1 d, per-quarter-normalised flux, errors from out-of-transit MAD; `redteam/strict_threshold.txt`):

| events kept | PSF depth | PSF chi2/dof | SAP depth |
|---|---|---|---|
| all 61 | 335 +- 60 ppm | - | 199 +- 56 |
| > 1.0 d after last gap end (n=53) | 260 +- 63 (4.2 sigma) | 75/52 | 208 +- 60 (3.5 sigma) |
| > 2.0 d (n=47) | 201 +- 65 (3.1 sigma) | 59/46 | 97 +- 63 (1.5 sigma) |
| > 2.5-4.0 d (n=44-46) | 213-220 +- 65-67 (3.3 sigma) | 56-57/43-45 | 112-113 +- 63-64 (1.8 sigma) |

The step between the 1 d and 2 d thresholds is epochs 5, 20 and 67, which sit 1.3-2.0 d after gaps that have measured ramps (-6.0, -1.9, +1.1 ppt) and are three of the five highest-SNR "non-gap" events. Against the 1 d-threshold ensemble: dropping the top 1/2/3 events by SNR gives 207, 185, 161 +- 70 ppm (3.0, 2.7, 2.3 sigma); per season the depth is 386 +- 147 (S0, ch17), 505 +- 141 (S1, ch57), -110 +- 145 (S2, ch65), 233 +- 118 (S3, ch25), i.e. absent on the channel with CROWDSAP = 1.0 and present where the ramps live. A phase scan of the same ensemble in 1.5 h steps (`redteam/residual_robustness.json`, `phase_scan`) gives +254 +- 68 ppm at phase 0 and -265 +- 67 ppm (a brightening of equal significance) at -3 h, with -136 to -152 at -4.5 to -7.5 h: the "dip" is one excursion in a structured residual, not a box on a flat baseline. A random-ephemeris null of the 1 d-threshold ensemble (300 random phases at the same P, and 300 random P in 15-30 d) has sd 1.2-1.3 in SNR units (red noise) and max 3.1; the observed 3.7 is beyond it, but that is a single-trial number for a period and phase chosen by the full search and it drops to 3.1-3.3 with the ramp timescale respected and to 2.3 without the top three events. The CROWDSAP of 0.95-1.0 and FLFRCSAP 0.79-0.99 require an on-target signal to be present in SAP at >= 95 % of its PSF depth; at the 2 d threshold SAP shows 48 +- 45 % of PSF.

## 4. Other kill channels checked (none needed, all consistent with the above)

- Stellar variability / rotation aliases (`redteam/ls_physical.txt`): Lomb-Scargle 0.5-100 d on per-quarter-normalised PSF, SAP and PSF_FLAT_FLUX. No power at 21.1, 10.55, 42.2 or 7.04 d (PSF power 0.002 at 21.17 d, ~220 ppm half-amplitude ceiling; PSF_FLAT_FLUX < 10 ppm). Top raw peaks at 46, 41, 31 d are quarter-length trends. Rotation is not the origin.
- Blended/background EB: odd/even 0.09-0.25 sigma, no secondary (-0.17 to -0.55 sigma), no V-shape test possible because the folded shape is a ramp (deep ingress side, long tail; box SNR rising to the 8-10 h grid edge in the re-derivation). The non-ramp residual shows no EB signature either, but it is too weak for the test to mean anything.
- Kepler systematics: the whole signal is the post-downlink / quarter-roll thermal recovery seen through the PSF model (sections 1-2). The pipeline's own `gap_frac_1d` = 0.21, `event_weighted_depth` = 318 ppm vs `depth_pass2` = 656 ppm and `event_chi2_const` = 87/60 were already flagging this.
- Ephemeris-match contamination: channels per season are 17 (S0), 57 (S1), 65 (S2), 25 (S3). No KOI, DR25 TCE, Kepler EB (Kirk et al. 2016, J/AJ/151/68) or confirmed planet within 180 arcsec; the only TCE within 3 arcmin is KIC 5609120 (141 arcsec, P = 484.95 d, MES 8.3, unrelated). Field-wide period matches to 21.106 d or its 1/4-4x harmonics within 0.2 %: 56 KOIs, 124 TCEs, 11 EBs; none has an epoch within 0.02 in phase of ours except one 84.5 d TCE on KIC 1160891 at 4.5 deg (phase 0.012), which is a chance coincidence and has no mechanism to reach channel 17. The closest period twin, Kepler host KIC 9003401 (K02319.01, 21.1033 d), is 5.6 deg away with a 2.4 d (0.11 phase) epoch offset. KOI-6515 (KIC 5091614, EB, 21.1425 d) is 1.0 deg away, 0.17 % off in period, and would have drifted 2.5 d over the baseline.
- Already known? KIC 5609180 was never a Kepler target (not in `kepler_target_kepids.csv`); the TPF host KIC 5609183 (Kp 13.74, 12.6 arcsec) was observed in all quarters and has no KOI, TCE or EB entry; a 656 ppm leak from it would require a ~13 ppt event on a Kp 13.7 star, which the Kepler pipeline could not have missed. Live NASA Exoplanet Archive TAP queries (2026-10-17) return no KOI/TCE for either KIC and no TOI within 0.05 deg. SIMBAD within 30 arcsec returns only 2MASS J19133455+4051175 (= the TPF host). Bienias & Szabo 2025 report seven hot-Jupiter candidates at P = 2.5-4.8 d (KICs listed in `docs/02-landscape-and-open-ground.md`); neither KIC is among them and the period regime is different. The Moquin et al. 2026 poster object BOI-1.01 (P = 3.37 d) is unrelated. Novelty is not the problem; reality is.
- Duration/density: Rstar 0.93 Rsun (CSV), Teff 5650 K, BP-RP 0.83, M_G 4.19 at 3.43 kpc (plx SNR 5.6), RUWE 1.07. For M ~ 1 Msun, a = 34.6 Rs and the b = 0 maximum T14 is 4.7 h (6.5 h even for 1.3 Rsun). The pipeline's 3.0 h is allowed (Rp 2.2-2.6 Re for 483-656 ppm; 1.55 Re for 233 ppm), but the 8 h box-SNR-preferred duration of the actual signal is physically impossible for this host, as expected for ramps.

## 5. What would resurrect it (for completeness; not recommended)

Only an on-target, box-shaped signal at 21.106 d that survives with the post-gap day masked, appears in SAP at >= 95 % of its PSF depth, and has consistent per-event depths across all four channels. A pixel-level run on the originating TPF (KIC 5609183) is not warranted on the present evidence; if one is nonetheless done, it must exclude the twelve gap-edge epochs or it will localise the thermal ramp, not a star.

## Files
- `runs/poolA_main/adversarial/kic-005609180/redteam/ramp_census.{txt,json}`: post-gap ramp amplitudes in PSF and SAP at all 85 gaps, phase at 21.106 d.
- `runs/poolA_main/adversarial/kic-005609180/redteam/quarter_diag.{txt,png}`: FLUX / PSF_FLUX_NOVA / SAP / SAP_BKG / RED_CHI2 / centroids around epochs 4, 36, 39.
- `runs/poolA_main/adversarial/kic-005609180/redteam/stitched_vs_quarter.txt`: per-quarter and cadence-level comparison of the stitched PSF flux with the per-quarter PSF, NOVA and SAP columns.
- `runs/poolA_main/adversarial/kic-005609180/redteam/event_depths_null.{txt,json}`: own per-event depths and random-ephemeris nulls.
- `runs/poolA_main/adversarial/kic-005609180/redteam/residual_robustness.{txt,json}`: leave-top-k, per-season, per-year, duration scan, phase scan, alias comb.
- `runs/poolA_main/adversarial/kic-005609180/redteam/strict_threshold.{txt,json}`: residual depth vs time-since-gap threshold; big-ramp counts per season.
- `runs/poolA_main/adversarial/kic-005609180/redteam/masked_bls.{txt,json}`: 15-30 d BLS with post-gap day masked, real/inverted, PSF/SAP.
- `runs/poolA_main/adversarial/kic-005609180/redteam/ls_physical.txt`: Lomb-Scargle both columns; duration and radius bounds.
