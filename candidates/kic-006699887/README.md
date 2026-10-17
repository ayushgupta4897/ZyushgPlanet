# kic-006699887 (Gaia DR3 2077554781499350144, KIC 6699887) — evidence pack

_Built 2026-10-17T18:21:46Z from run `poolA_main` at commit `87e2c804b64e`. Status under `CLAUDE.md` language: **candidate**. Nothing here is a planet; nothing is 'discovered'. Search, vetting, localisation, re-derivation, red team and FPP were run by Claude Code (Anthropic) agents under human direction; every number below names the file it comes from._

## Key numbers

| quantity | value | statistic / subset | file |
|---|---|---|---|
| Period | 15.113612 d | stage-2 BLS refinement, all quarters | stage2/vet.json `P` |
| Mid-transit (BJD) | 2454967.72786 | stage-2 | stage2/vet.json `t0_bjd` |
| Box duration | 3.75 h | stage-2 box | stage2/vet.json `dur_h` |
| Depth | 1023 ppm (formal ± 33) | masked re-detrend, unweighted box mean | stage2/vet.json `depth_pass2_ppm` |
| Empirical SNR | 20.0 | depth / scatter of 80 random off-transit boxes | stage2/vet.json `snr_emp_pass2` |
| Events | 86 (100 % positive) | per-event local depths | stage2/events.csv |
| Per-event χ²/dof | 151.6 / 85 | formal errors understate red noise (controls 1.5–2.0) | stage2/vet.json |
| Odd/even | 0.03 σ | stage-2 | stage2/vet.json `odd_even_sigma` |
| Secondary at phase 0.5 | -37 ppm (-0.73 σ) | two-sided test | stage2/vet.json |
| TLS (±2 % window) | SDE 13.7, SNR 28.1, 86.0 transits | not a calibrated SDE (see prereg v1.5) | stage2/vet.json `tls` |
| R_p if on target | 2.30 R⊕ (box), 2.23 R⊕ (TLS) | R★ 0.66 R☉ colour-based | stage2/vet.json |
| Stage-2 verdict | PASS (poolA_kbonus_background_dwarfs_v1.7); v1.4: PASS | post-hoc relative to stage 1 | stage2/vet.json |
| Stage-1 BLS | P 15.11386 d, SDE 89.2, SNR_emp 22.3, 86 transits, max single-event fraction 0.04; inverted (2−f) top peak SDE 14.1 at 0.657 d | blind screen, config ccca76e882999a44 | stage1_row.json |
| Stage-1 second signal (after masking the first) | P 23.8450 d, depth 442 ppm, SDE 33.2, SNR_emp 13.7 | 0.5–40 d window | stage1_row.json `real.s2` |
| Pixel localisation (v1.6 χ²) | **on-target**: Σ(χ²null−χ²target) 503, Σ(χ²best-alt−χ²target) 200 over 14 quarters; strongest alternative scaled_image | post-hoc rule, passed control KIC 5689121 | pixels/pixel_v16.json |
| Independent re-derivation | agrees: P 15.113718 d, depth 1086 ± 33 ppm, SNR 33.0, 86 events | no shared code | independent/rederive_result.json |
| Red team | **survives**; kill reasons: none | | independent/redteam_verdict.json, independent/redteam.md |
| TRICERATOPS FPP | 0.0240 ± 0.0208 (8 seeds × 1000000 draws) | Kepler mode, TIC stellar parameters, no contrast curve | fpp/fpp_summary.json |
| TRICERATOPS NFPP | 0.00000 ± 0.00000 | nearby stars excluded only by their own light curves (fpp/fpp_summary.json `neighbour_tests`) | fpp/fpp_summary.json |
| Host (Gaia DR3) | G 16.22, BP−RP 1.35, Teff 4542 K, R★ 0.66 R☉, 859 pc, RUWE 1.09 | colour-based radius; TIC values in fpp/stars_tic.csv | target_row.json |
| Ephemeris matches (KOI/TCE/EB/confirmed within 3′, ratios 1:1,2:1,1:2,3:1,1:3, 0.2 %) | 0 (none on the same star or TPF) | catalogues dated 2026-10-17T15:46:52+00:00 | stage2/vet.json `known_matches`, known_catalogues_manifest.json |
| Field-wide ephemeris check (whole Kepler field, any channel; KOI, DR25 TCE, Kepler EB; harmonics 1:1–4:1) | 69 period coincidences within 0.1 %, 9 with a matching epoch, **0 viable** (1:1 twins, or harmonics not excluded by this signal's odd/even and secondary tests) | a 4-year baseline separates period twins whose epochs differ; harmonic twins need an odd/even or phase-0.5 signature this signal lacks | ephemeris_fieldwide.json |

## Second signal (post-hoc, label P47)

Found by the red team outside the registered 0.5–40 d window. Stage 2 under current rules: P 47.69006 d, depth 951 ppm, SNR 21.3, 28 events (100 % positive), odd/even 2.72 σ, secondary 0.20 σ, verdict PASS; known matches 0. Pixel v1.6: inconclusive (4 quarters: null-target dchi2 70, best-alternative-target dchi2 24 (rule: both >= 25 over >= 5 quarters)). Files: second_signal_P47/.

## Deep single-target run (`pipeline/run.py`, `deep_run/results.json`)

| quantity | value | note |
|---|---|---|
| Blind recovery | BLS 1–60 d best P n/a d, SDE 33.8; TLS 1–50 d P n/a d, SDE n/a, SNR 25.5 | independent detrend window and grids from the batch pipeline |
| Adopted ephemeris | P 15.113735 d, t0 BJD 2454967.72682, duration 3.50 h | BLS refine on transit-masked detrended light curve |
| Agreement with stage 2 | ΔP 8.2 ppm, mid-baseline t0 offset 7.0 min; depth at the stage-2 ephemeris 1020 ± 31 ppm | |
| Box depth (pass 2) | 1070 ± 31 ppm (formal), 569 in-transit cadences | empirical SNR 42.0 (random-phase boxes: 400, std 25 ppm) |
| Odd / even | 1092 / 1047 ppm, difference 0.71 σ | |
| Secondary | phase 0.5: -18 ppm (-0.72 σ); deepest other dip in a 195-phase scan 2.5 σ | scan σ scatter 0.96 |
| Inverted-flux null (masked LC) | BLS best SDE 4.7; TLS best SDE 8.0 at 47.691 d | full-range TLS SDE, not the stage-1 statistic |
| TPF host (KIC 6699878) at the ephemeris | -5 ± 15 ppm (-0.3 σ) | worst case: would need only 379 ppm if all of its light entered the target's photometry; the FPP step's aperture model gives the required depth under PSF separation (fpp/fpp_summary.json `neighbour_tests`), and the pixel verdict above is the primary host exclusion |

## What this pack does not show

- No high-resolution imaging, no radial velocities, no contrast curve: the unresolved-companion scenarios in the FPP carry their prior weight only.
- The pixel model is a Gaussian point source, not the Kepler PRF; many background targets sit at the stamp edge (quarters flagged).
- Stage-2 rules (v1.3–v1.7) and the pixel rule (v1.6) were registered after results had been seen; every such result is post-hoc relative to stage 1 (prereg_poolA.json `revision_log`).
- 'Validated' requires FPP < 0.015 and NFPP < 0.001 in a refereed paper; this pack is a pipeline product.
- Novelty holds for the catalogues listed in known_catalogues_manifest.json on their download date; the Kepler Bonus team's own candidate list is unpublished.

## Files

`SHA256SUMS` lists every file; `manifest.json` adds sizes and the git commit.
