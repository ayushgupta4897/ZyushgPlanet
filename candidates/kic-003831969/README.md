# kic-003831969 (Gaia DR3 2100154070473374976, KIC 3831969) — evidence pack

_Built 2026-10-17T18:04:56Z from run `poolA_main` at commit `0282131787ab`. Status under `CLAUDE.md` language: **candidate**. Nothing here is a planet; nothing is 'discovered'. Search, vetting, localisation, re-derivation, red team and FPP were run by Claude Code (Anthropic) agents under human direction; every number below names the file it comes from._

## Key numbers

| quantity | value | statistic / subset | file |
|---|---|---|---|
| Period | 15.313350 d | stage-2 BLS refinement, all quarters | stage2/vet.json `P` |
| Mid-transit (BJD) | 2454977.28486 | stage-2 | stage2/vet.json `t0_bjd` |
| Box duration | 3.75 h | stage-2 box | stage2/vet.json `dur_h` |
| Depth | 997 ppm (formal ± 58) | masked re-detrend, unweighted box mean | stage2/vet.json `depth_pass2_ppm` |
| Empirical SNR | 23.8 | depth / scatter of 80 random off-transit boxes | stage2/vet.json `snr_emp_pass2` |
| Events | 86 (92 % positive) | per-event local depths | stage2/events.csv |
| Per-event χ²/dof | 194.8 / 85 | formal errors understate red noise (controls 1.5–2.0) | stage2/vet.json |
| Odd/even | 0.54 σ | stage-2 | stage2/vet.json `odd_even_sigma` |
| Secondary at phase 0.5 | 1 ppm (0.02 σ) | two-sided test | stage2/vet.json |
| TLS (±2 % window) | SDE 13.3, SNR 16.7, 86.0 transits | not a calibrated SDE (see prereg v1.5) | stage2/vet.json `tls` |
| R_p if on target | 2.85 R⊕ (box), 2.67 R⊕ (TLS) | R★ 0.83 R☉ colour-based | stage2/vet.json |
| Stage-2 verdict | PASS (poolA_kbonus_background_dwarfs_v1.7); v1.4: PASS | post-hoc relative to stage 1 | stage2/vet.json |
| Stage-1 BLS | P 15.31387 d, SDE 45.8, SNR_emp 15.2, 86 transits, max single-event fraction 0.05; inverted (2−f) top peak SDE 13.2 at 0.530 d | blind screen, config ccca76e882999a44 | stage1_row.json |
| Stage-1 second signal (after masking the first) | P 0.5297 d, depth 65 ppm, SDE 13.9, SNR_emp nan | 0.5–40 d window | stage1_row.json `real.s2` |
| Pixel localisation (v1.6 χ²) | **on-target**: Σ(χ²null−χ²target) 219, Σ(χ²best-alt−χ²target) 112 over 13 quarters; strongest alternative 2100153314562872320 | post-hoc rule, passed control KIC 5689121 | pixels/pixel_v16.json |
| Independent re-derivation | agrees: P 15.313360 d, depth 1009 ± 51 ppm, SNR 19.6, 86 events | no shared code | independent/rederive_result.json |
| Red team | **survives**; kill reasons: none | | independent/redteam_verdict.json, independent/redteam.md |
| TRICERATOPS FPP | 0.0272 ± 0.0029 (8 seeds × 1000000 draws) | Kepler mode, TIC stellar parameters, no contrast curve | fpp/fpp_summary.json |
| TRICERATOPS NFPP | 0.00000 ± 0.00000 | nearby stars excluded only by their own light curves (fpp/fpp_summary.json `neighbour_tests`) | fpp/fpp_summary.json |
| Host (Gaia DR3) | G 16.53, BP−RP 0.99, Teff 5275 K, R★ 0.83 R☉, 1969 pc, RUWE 1.06 | colour-based radius; TIC values in fpp/stars_tic.csv | target_row.json |
| Ephemeris matches (KOI/TCE/EB/confirmed within 3′, ratios 1:1,2:1,1:2,3:1,1:3, 0.2 %) | 0 (none on the same star or TPF) | catalogues dated 2026-10-17T15:46:52+00:00 | stage2/vet.json `known_matches`, known_catalogues_manifest.json |
| Field-wide ephemeris check (whole Kepler field, any channel; KOI, DR25 TCE, Kepler EB; harmonics 1:1–4:1) | 59 period coincidences within 0.1 %, 0 with a matching epoch, **0 viable** (1:1 twins, or harmonics not excluded by this signal's odd/even and secondary tests) | a 4-year baseline separates period twins whose epochs differ; harmonic twins need an odd/even or phase-0.5 signature this signal lacks | ephemeris_fieldwide.json |

## Deep single-target run (`pipeline/run.py`, `deep_run/results.json`)

| quantity | value | note |
|---|---|---|
| Blind recovery | BLS 1–60 d best P n/a d, SDE 17.8; TLS 1–50 d P n/a d, SDE n/a, SNR 12.8 | independent detrend window and grids from the batch pipeline |
| Adopted ephemeris | P 15.313390 d, t0 BJD 2454977.27966, duration 3.80 h | BLS refine on transit-masked detrended light curve |
| Agreement with stage 2 | ΔP 2.6 ppm, mid-baseline t0 offset -4.8 min; depth at the stage-2 ephemeris 1037 ± 56 ppm | |
| Box depth (pass 2) | 1045 ± 56 ppm (formal), 608 in-transit cadences | empirical SNR 23.2 (random-phase boxes: 400, std 45 ppm) |
| Odd / even | 1042 / 1049 ppm, difference -0.06 σ | |
| Secondary | phase 0.5: -6 ppm (-0.12 σ); deepest other dip in a 193-phase scan 2.4 σ | scan σ scatter 0.89 |
| Inverted-flux null (masked LC) | BLS best SDE 5.1; TLS best SDE 10.0 at 1.554 d | full-range TLS SDE, not the stage-1 statistic |
| TPF host (KIC 3831977) at the ephemeris | -32 ± 19 ppm (-1.7 σ) | would need 252 ppm to host the signal |

## What this pack does not show

- No high-resolution imaging, no radial velocities, no contrast curve: the unresolved-companion scenarios in the FPP carry their prior weight only.
- The pixel model is a Gaussian point source, not the Kepler PRF; many background targets sit at the stamp edge (quarters flagged).
- Stage-2 rules (v1.3–v1.7) and the pixel rule (v1.6) were registered after results had been seen; every such result is post-hoc relative to stage 1 (prereg_poolA.json `revision_log`).
- 'Validated' requires FPP < 0.015 and NFPP < 0.001 in a refereed paper; this pack is a pipeline product.
- Novelty holds for the catalogues listed in known_catalogues_manifest.json on their download date; the Kepler Bonus team's own candidate list is unpublished.

## Files

`SHA256SUMS` lists every file; `manifest.json` adds sizes and the git commit.
