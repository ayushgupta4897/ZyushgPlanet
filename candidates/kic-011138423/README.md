# kic-011138423 (Gaia DR3 2129669944725625344, KIC 11138423) — evidence pack

_Built 2026-10-17T21:01:00Z from run `poolA_main` at commit `6a93dbe97b24`. Status under the project's language rules: **candidate**. Nothing here is a planet; nothing is 'discovered'. Author: Aayush Gupta; analysis code written with AI assistance (Claude); every number below names the file it comes from._

## Key numbers

| quantity | value | statistic / subset | file |
|---|---|---|---|
| Period | 4.894451 d | stage-2 BLS refinement, all quarters | stage2/vet.json `P` |
| Mid-transit (BJD) | 2454958.33663 | stage-2 | stage2/vet.json `t0_bjd` |
| Box duration | 1.95 h | stage-2 box | stage2/vet.json `dur_h` |
| Depth | 1316 ppm (formal ± 79) | masked re-detrend, unweighted box mean | stage2/vet.json `depth_pass2_ppm` |
| Empirical SNR | 15.2 | depth / scatter of 80 random off-transit boxes | stage2/vet.json `snr_emp_pass2` |
| Events | 253 (87 % positive) | per-event local depths | stage2/events.csv |
| Per-event χ²/dof | 347.3 / 252 | formal errors understate red noise (controls 1.5–2.0) | stage2/vet.json |
| Odd/even | 0.78 σ | stage-2 | stage2/vet.json `odd_even_sigma` |
| Secondary at phase 0.5 | -81 ppm (-0.94 σ) | two-sided test | stage2/vet.json |
| TLS (±2 % window) | SDE 18.5, SNR 14.5, 253.0 transits | not a calibrated SDE (see prereg v1.5) | stage2/vet.json `tls` |
| R_p if on target | 2.11 R⊕ (box), 2.08 R⊕ (TLS) | R★ 0.53 R☉ colour-based | stage2/vet.json |
| Stage-2 verdict | PASS (poolA_kbonus_background_dwarfs_v1.7); v1.4: PASS | post-hoc relative to stage 1 | stage2/vet.json |
| Stage-1 BLS | P 4.89440 d, SDE 33.0, SNR_emp 12.9, 253 transits, max single-event fraction 0.02; inverted (2−f) top peak SDE 15.0 at 0.612 d | blind screen, config ccca76e882999a44 | stage1_row.json |
| Stage-1 second signal (after masking the first) | P 27.9176 d, depth 995 ppm, SDE 13.3, SNR_emp 8.1 | 0.5–40 d window | stage1_row.json `real.s2` |
| Pixel localisation (v1.6 χ²) | **on-target**: Σ(χ²null−χ²target) 197, Σ(χ²best-alt−χ²target) 140 over 13 quarters; strongest alternative 2129669944725610624 | post-hoc rule, passed control KIC 5689121 | pixels/pixel_v16.json |
| Independent re-derivation | agrees: P 4.894405 d, depth 1389 ± 95 ppm, SNR 14.7, 253 events | no shared code | independent/rederive_result.json |
| Red team | **survives**; kill reasons: none | | independent/redteam_verdict.json, independent/redteam.md |
| TRICERATOPS FPP | 0.0364 ± 0.0019 (8 seeds × 1000000 draws) | Kepler mode, TIC stellar parameters, no contrast curve | fpp/fpp_summary.json |
| TRICERATOPS NFPP | 0.00000 ± 0.00000 | nearby stars excluded only by their own light curves (fpp/fpp_summary.json `neighbour_tests`) | fpp/fpp_summary.json |
| Host (Gaia DR3) | G 17.35, BP−RP 1.72, Teff 4258 K, R★ 0.53 R☉, 894 pc, RUWE 1.06 | colour-based radius; TIC values in fpp/stars_tic.csv | target_row.json |
| Ephemeris matches (KOI/TCE/EB/confirmed within 3′, ratios 1:1,2:1,1:2,3:1,1:3, 0.2 %) | 0 (none on the same star or TPF) | catalogues dated 2026-10-17T15:46:52+00:00 | stage2/vet.json `known_matches`, known_catalogues_manifest.json |
| Pixel localisation, season-stacked (v1.8) | **on-target**: 4 seasons, 278 events: null-target 82, best-alt-target 74 | post-hoc method extension; controls: KIC 5689121 on, two background EBs offset | pixels/season_stack/pixel_season.json |
| TRICERATOPS FPP with MOLUSC companions + Gaia FLAME host | **0.0070 ± 0.0021**, NFPP 0.00000 (8 seeds) | bound companions constrained by Gaia RUWE and Gaia imaging (G ≤ 20); companion prior 14.7 % vs the default; multiplicity fraction 0.46 (conservative) | fpp_molusc/fpp_summary.json, molusc/molusc_summary.json |
| Transit fit (batman + emcee) | rp/rs 0.0381, b 0.53 (+0.28/−0.35), T14 1.91 h, Rp 2.83 R⊕ | stellar density from transit 4.24 g/cc vs Gaia FLAME 3.02 g/cc (0.4 σ) | transit_fit/fit.json |
| Field-wide ephemeris check (whole Kepler field, any channel; KOI, DR25 TCE, Kepler EB; harmonics 1:1–4:1) | 65 period coincidences within 0.1 %, 10 with a matching epoch, **0 viable** (1:1 twins, or harmonics not excluded by this signal's odd/even and secondary tests) | a 4-year baseline separates period twins whose epochs differ; harmonic twins need an odd/even or phase-0.5 signature this signal lacks | ephemeris_fieldwide.json |

## Deep single-target run (`pipeline/run.py`, `deep_run/results.json`)

| quantity | value | note |
|---|---|---|
| Blind recovery | BLS 1–60 d best P n/a d, SDE 21.9; TLS 1–50 d P n/a d, SDE n/a, SNR 15.3 | independent detrend window and grids from the batch pipeline |
| Adopted ephemeris | P 4.894397 d, t0 BJD 2454958.34715, duration 1.80 h | BLS refine on transit-masked detrended light curve |
| Agreement with stage 2 | ΔP -11.1 ppm, mid-baseline t0 offset 3.5 min; depth at the stage-2 ephemeris 1384 ± 70 ppm | |
| Box depth (pass 2) | 1453 ± 72 ppm (formal), 857 in-transit cadences | empirical SNR 18.7 (random-phase boxes: 400, std 78 ppm) |
| Odd / even | 1542 / 1365 ppm, difference 1.23 σ | |
| Secondary | phase 0.5: 23 ppm (0.33 σ); deepest other dip in a 191-phase scan 3.2 σ | scan σ scatter 1.19 |
| Inverted-flux null (masked LC) | BLS best SDE 4.9; TLS best SDE 9.8 at 5.534 d | full-range TLS SDE, not the stage-1 statistic |
| TPF host (KIC 11138418) at the ephemeris | -1 ± 4 ppm (-0.2 σ) | worst case: would need only 17 ppm if all of its light entered the target's photometry; the FPP step's aperture model gives the required depth under PSF separation (fpp/fpp_summary.json `neighbour_tests`), and the pixel verdict above is the primary host exclusion |

## What this pack does not show

- No high-resolution imaging, no radial velocities, no contrast curve: the unresolved-companion scenarios in the FPP carry their prior weight only.
- The pixel model is a Gaussian point source, not the Kepler PRF; many background targets sit at the stamp edge (quarters flagged).
- Stage-2 rules (v1.3–v1.7) and the pixel rule (v1.6) were registered after results had been seen; every such result is post-hoc relative to stage 1 (prereg_poolA.json `revision_log`).
- 'Validated' requires FPP < 0.015 and NFPP < 0.001 in a refereed paper; this pack is a pipeline product.
- Novelty holds for the catalogues listed in known_catalogues_manifest.json on their download date; the Kepler Bonus team's own candidate list is unpublished.

## Files

`SHA256SUMS` lists every file; `manifest.json` adds sizes and the git commit.
