# Synthesis of the red-team pass: poolA_main survivors (2026-10-17)

Candidate language only. Nothing below is a planet; nothing is "discovered". "Candidate-grade" here is the workflow status defined for this step (independent re-derivation agrees, red team found no kill, pipeline pixel localisation on-target). Under CLAUDE.md none of these objects is yet a "candidate" in the public sense: no TRICERATOPS run exists for any of them, and the pixel rule that produced the on-target verdicts (v1.6) is itself post-hoc (see "Caveats that travel with every object").

## Pool context (run id poolA_main)

| statistic | value | source |
|---|---|---|
| stars with results | 6365 (31 skipped) | `runs/poolA_main/triage.md` |
| first-screen passes, real | 438 | `results.jsonl`, `pass_real_s1` |
| first-screen passes, inverted null (2 - f) | 95 (20 stars pass both polarities) | `results.jsonl`, `pass_inv_s1` |
| stage-2 passes, v1.4 rules | 11 | `vet_summary.csv`, `stage2_pass_v14` |
| stage-2 passes, v1.5 rules | 23 | `vet_summary.csv`, `stage2_pass` |
| red-teamed in this workflow | 4 (the verdicts synthesised here) | `runs/poolA_main/adversarial/<fname>/` |
| v1.6 pixel rule on positive controls | KIC 5689121 on-target (dchi2 663 / 388); Kepler-138 c inconclusive (target and scaled-image hypotheses coincide for a TPF host) | `runs/controls/pixel_v16/*/pixel.json` |

## Summary table (referee-sized)

Pipeline numbers from `runs/poolA_main/vet/<fname>/vet.json` (stage 2, v1.5) and `runs/poolA_main/pixel_v16/<fname>/pixel.json`; stage-1 statistics from `results.jsonl`; red-team and re-derivation numbers from `runs/poolA_main/adversarial/<fname>/`.

| rank | fname (object) | P (d) | depth (ppm, 3.75/1.125/4.5-h box) | T14 (h) | N events | SNR_emp pipeline / re-deriv. | stage-1 BLS SDE real / inverted | TLS SDE (narrow window) | stage-2 v1.5 | pixel v1.6 (dchi2 null-target / best-alt-target) | red-team localisation | host LC at ephemeris | novelty | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | kic-006699887 (KIC 6699887, G 16.2 K dwarf, 859 pc) | 15.1136 | 1023 +- 33 | 3.5-4.0 | 86 | 20.0 / 33.0 (same depth; noise estimator differs) | 89.2 / 14.1 | 13.7 | PASS | on-target (503 / 200), 14 quarters | free-position deficit 1.0 and 0.7 arcsec from Gaia position (seasons 0, 2); beyond-edge positions excluded > 3 sigma | KIC 6699878: -1.8 +- 15.8 ppm (>= 395 needed) | none in KOI/TCE/EB/confirmed/B&S 2025/Moquin 2026/TOI; never a Kepler target | candidate-grade |
| 2 | kic-003831969 (KIC 3831969, G 16.5 K/G dwarf, 1968 pc) | 15.3134 | 997 +- 58 (1009 +- 51 re-deriv.) | 4.06-4.29 (T23/T14 0.80) | 86 | 23.8 / 19.6 (quote 19.6) | 45.8 / 13.2 | 13.3 | PASS | on-target (219 / 112), 13 quarters | free-centre deficit 0.33 / 0.32 / 0.40 arcsec from target in seasons 0, 1, 3; all Gaia neighbours within 25 arcsec rejected at dchi2 > 65 | KIC 3831977: -44 +- 22 ppm PSF, +20 +- 21 SAP | none (incl. ExoFOP TOI/CTOI 2026-10-17, TIC 120425498) | candidate-grade |
| 3 | kic-009520080 (KIC 9520080, G 17.5 K dwarf) | 2.3770 | 1569 +- 62 (~2000 smearing-corrected) | 1.13-1.15 (b ~ 0.85) | 468 | 28.3 / recovered blind, SDE 108 vs inverted 7.5 | 101.5 / 19.7 | 22.2 | FAIL (secondary_2sided, -4.0 sigma; passes v1.4) | on-target (917 / 335), 17 quarters | target-pixel deficit in all 17 TPFs, 1.02 +- 0.04 e/s (24.5 sigma); host pixels -0.33 +- 0.51 e/s | KIC 9520088: -6 +- 12 ppm; KIC 9520089: -190 +- 210 ppm | none; Kepler-113 c at 2:1 is 4.8 deg away, not a contaminant | signal-needs-work |
| 4 | kic-006448637 (KIC 6448637, G 15.7 K dwarf) | 35.5873 | 846 +- 52 (897 +- 45 re-deriv.) | 4.1-4.8 | 26 | 10.0 / 20.0 (noise-scale difference) | 36.0 / 10.0 | 11.6 | PASS | inconclusive (121 / 1.8); saturated Kp 8.75 TPF host with bleed column | red-team local TPF box: deficit at target 12/12 usable quarters, centroid dx 0.00 +- 0.06, dy +0.21 +- 0.07 px; scratch result, not a pipeline product | KIC 6448650: -75 +- 121 ppm (>= 1120 needed); KIC 6448629: 971 +- 694 (1.0% needed) | none; Kepler-1715 b same channel but 20 sigma off in period and 0.84-1.20 d off in time | signal-needs-work |

All four: inverted-flux null fails stage 1 at the pre-registered thresholds (inverted SDE 14.1 / 13.2 / 19.7 / 10.0 vs real 89.2 / 45.8 / 101.5 / 36.0); all four were recovered blind by an independent re-derivation with no shared code (`adversarial/<fname>/rederive*`); none has a TRICERATOPS run.

## Why this ranking

1. **kic-006699887** has the cleanest end-to-end record: stage 2 pass, pipeline on-target pixel verdict with the largest margin over the best alternative hypothesis among the pixel-passing objects relative to its depth, an own-code free-position localisation that explicitly tested positions beyond the stamp edge (the pipeline's weakness for edge stars), a flat host light curve at > 20 sigma, SAP/PSF ratio equal to CROWDSAP, a clean 702-phase secondary scan, and a second, independently localised 47.69-d signal that disfavours a single blended EB. The residual is the generic one (unresolved EB inside ~1.5 arcsec).
2. **kic-003831969** matches it on every pipeline stage and has the tightest red-team localisation (0.3-0.4 arcsec in three independent pixel orientations) but carries two quantified caveats that lower it: a 2-sigma season dependence of the PSF depth tied to psffrac (~25% systematic on depth/radius) and per-event scatter 1.6x above white noise (chi2 195-220 / 85) so the quoted depth errors are under-estimated. Neither is evidence against the signal.
3. **kic-009520080** is physically the strongest signal (SNR_emp 28.3, 468 events, 24.5-sigma target-pixel deficit) and the red team found no false-positive mechanism, but it is a formal stage-2 FAIL under the current rules: a -4.0-sigma (pipeline) / 2.2-2.8-sigma (red team, estimator-dependent) brightening at phase 0.5. No planet produces a brightening; until an independent analysis explains it (outlier clustering, detrending, or real), it cannot be promoted. Rule (c) is post-hoc and was motivated by this very star; both facts must appear in any report.
4. **kic-006448637** survives the red team and re-derivation but the pipeline pixel product is inconclusive because the TPF host is saturated; the on-target localisation exists only as red-team scratch code. It also has the lowest conservative SNR (10.0), the fewest events (26), a 2.0-sigma odd/even under re-derivation errors, and a single Q13-edge event (epoch 28) that drives a per-year depth drift. The same-channel KOI (Kepler-1715 b) is excluded by timing, which is sufficient for crosstalk, but should be documented in the evidence pack.

## Caveats that travel with every object

- **Post-hoc rules.** Stage 2 (v1.3), the per-event chi2 rule (v1.4), the v1.5 changes and the v1.6 pixel chi2 model comparison were all registered after results had been seen (`pipeline/prereg_poolA.json`, `revision_log`). Every stage-2 and pixel result above is therefore post-hoc relative to the stage-1 registration and must be labelled as such. The v1.6 rule did pass its positive control KIC 5689121 (on-target, dchi2 663 / 388).
- **Pixel localisation of edge stars.** Background sources sit at stamp edges by construction. The pipeline's v1.5 centroid verdict is "inconclusive" for all four; the v1.6 chi2 comparison is on-target for three. The red-team Gaussian localisations are not Kepler PRF models; a PRF forward model is on every object's list.
- **Unresolved blends.** Kepler pixels cannot exclude an EB within ~1-2 arcsec for any of these (required eclipse depths 1.7-15% depending on dG). This is exactly what a documented TRICERATOPS run plus a contrast curve must price before the word "candidate" is used in public.
- **Novelty is dated.** "Not in KOI / DR25 TCE / Kepler EB / confirmed / Bienias & Szabo 2025 / Moquin et al. 2026 / ExoFOP TOI-CTOI" holds for the catalogues in `data/known/` with their 2026-10-17 download dates; the KBonus team's TraSHDUMP list is unpublished, so the correct statement is "not in any public list as of 2026-10-17".
- **Which statistic.** Where pipeline and re-derivation SNRs differ (kic-006699887 20.0 vs 33.0, kic-006448637 10.0 vs 20.0, kic-003831969 23.8 vs 19.6) the depths agree and the difference is the noise estimator (pipeline `box_noise_ppm` vs random-phase-box scatter). Quote the more conservative value in each case and say which it is.

## Per-survivor next steps

### 1. kic-006699887 — candidate-grade
- TRICERATOPS with the Kepler PRF, Gaia DR3 contrast, secondary limit (-26 +- 29 ppm) and odd/even (0.03 sigma); the only surviving scenario is an unresolved EB within ~1-1.5 arcsec at dG >= 3.
- Kepler PRF forward model of the deficit maps for all four seasons (target + host + contaminant amplitudes) to get a formal exclusion radius and to calibrate the season-dependent depth (998 / 1045 / 1278 / 1036 ppm; Rp 2.2-2.7 Re).
- Run the 47.6906-d second signal (877 +- 49 ppm, 4.5 h, 26 events, pixel-localised to the target) through stage 2 and pixel_vet as a labelled post-hoc test; it was discarded at its P/2 alias by the 0.5-40 d window.
- Audit `pixel_vet.py`'s noise model for Q2, Q13, Q14 (pipeline: no signal; red-team maps: usual target deficit); extend second-signal searches beyond 40 d before the next pool run.
- Speckle/AO imaging of the G 16.2 target inside 1 arcsec.

### 2. kic-003831969 — candidate-grade
- TRICERATOPS with the Gaia field, depth/duration, secondary limit (3-sigma < 14% of primary) and odd/even (0.3-0.5 sigma).
- Kepler PRF-model fit to the stacked difference images including the degenerate season 2 (corner pixel of a 4x5 stamp).
- Stamp-truncation-aware PSF refit for the season dependence of depth (1199 / 1089 / 836 / 772 ppm); quote depth and radius (2.7-4.0 Re) with ~25% systematics until then.
- Extinction-corrected host radius (3D dust map at ~2 kpc).
- High-resolution imaging inside ~1-2 arcsec.
- In any write-up quote SNR 19.6 and the red-noise-inflated depth error (~70-80 ppm), and state per-event chi2 195 / 85 (pipeline) and 220 / 85 (re-derivation).

### 3. kic-009520080 — signal-needs-work
- Independent re-examination of the phase-0.5 brightening (second detrending family; PSF vs SAP per quarter and season). Pipeline -222 +- 55 ppm (-4.0 sigma); red team +192 ppm mean (z 2.84), clipped +126 (z 2.24), median +51 (z 1.32). If real at ~200 ppm, downgrade; if an outlier/estimator artefact, document it and re-verdict under the recorded rules.
- Any report must carry the v1.5 `secondary_2sided` stage-2 FAIL, `stage2_pass_v14 = true`, `stage2_pass = false`, and the fact that rule (c) was added post-hoc after this star was seen.
- Ground-based transit photometry at predicted times (~2000 ppm true depth, 1.1 h, G 17.5; 1-2 m class) to separate the target from the G 20.9 star at 3.5 arcsec.
- TRICERATOPS plus AO/speckle to bound G > 21 background EBs and bound-companion EBs inside the 4-arcsec pixel.
- Document the single event n = 342 (BKJD 945.09, Q10) behind the per-event chi2 inflation (912 / 467 -> ~1.5 reduced without it).

### 4. kic-006448637 — signal-needs-work
- Rerun `pipeline/pixel_vet.py` (or a documented extension) in a saturated-host / edge-of-stamp mode (local box around the Gaia position, bleed column excluded), validated on KIC 5689121 first, so the localisation becomes a pipeline product with a run id. Until then the pixel stage is not passed.
- Quote both depth sets, with and without epoch 28 (2268 +- 349 ppm; per-year chi2 13 / 4 -> 4.5 / 4; odd/even 1.99 -> 0.41 sigma), and both SNR scales (10.0 pipeline, 20.0 re-derivation) with the reason they differ.
- Explicit 2P test at 71.175 d (EB with equal eclipses not excluded; secondary limit at P ~130 ppm, 3 sigma).
- Keep the Kepler-1715 b same-channel table (P 20 sigma off; transits 0.84-1.20 d away at every epoch) in the evidence pack.
- TRICERATOPS plus high-resolution imaging inside 1 arcsec; live cross-match against ExoFOP TOI/CTOI, VSX, Gaia NSS and TraSHDUMP when published.

## Before anything goes public (CLAUDE.md checklist, status today)
1. Stage-2 pass with human-inspected plots: objects 1, 2, 4 pass; object 3 fails v1.5. Human inspection of `fold_events.png` / `diffimages.png` is recorded by the red teams but should be repeated by the author.
2. Pixel localisation on the originating TPF: pipeline on-target for 1, 2, 3; inconclusive for 4 (red-team scratch localisation only).
3. Ephemeris match within 3 arcmin and the TPF host: clear for all four (dated 2026-10-17).
4. Independent re-derivation with no shared code: done for all four, all agree within 1 sigma on P, depth, duration.
5. Red-team review: done for all four, no kill reasons.
6. Candidate language, AI disclosure, SHA-256 evidence manifest, and a documented TRICERATOPS run: outstanding for all four. No object should be called a candidate in public until item 6 is complete.

Evidence paths: `/home/user/ZyushgPlanet/runs/poolA_main/vet/<fname>/vet.json`, `/home/user/ZyushgPlanet/runs/poolA_main/pixel_v16/<fname>/pixel.json`, `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/<fname>/{redteam.md, redteam_verdict.json, rederive_result.json, rederive/, redteam*/}`, `/home/user/ZyushgPlanet/pipeline/prereg_poolA.json` (revision_log), `/home/user/ZyushgPlanet/runs/controls/pixel_v16/`.