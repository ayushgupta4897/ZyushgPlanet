# Synthesis of the red-team pass, batch B: eight poolA_main survivors (2026-10-17)

Candidate language only. Nothing below is a planet and nothing is "discovered". "Candidate-grade" is the workflow status for this step (independent re-derivation agrees, red team found no kill, pipeline pixel localisation on target); none of the eight objects reaches it. This synthesis complements `runs/poolA_main/adversarial/synthesis_batchA.md` (four objects, two candidate-grade); seven further adversarial directories (kic-001577049, kic-003963519, kic-006367225, kic-009655224, kic-010213563, kic-011138423, kic-011285493) are outside both syntheses.

## Headline

**All eight signals are rejected.** Six were killed by their red teams; the remaining two (kic-008620445, gaia-dr3-2127687048516119168) have no red-team record, but their independent re-derivations failed to recover them blind and a synthesis-stage spot check (own code, saved with SHA-256 manifests under `adversarial/<fname>/synthesis_spotcheck/`) shows their depth is carried by the same spacecraft-wide artefact cadences that killed their siblings. Their status is **rejected, provisional** until a formal red-team record is written.

Every one of the eight is a **data-segment-edge or PSF-model artefact**, not a stellar or planetary signal:

| mechanism | objects |
|---|---|
| post-gap / post-downlink PSF ramps aligned into a period by the monthly gap comb | kic-005609180 (ch 17), kic-005866176 (ch 17, off-stamp star), kic-004474607 (Q13 starts 1182.74 / 1215.72), gaia-dr3-2127687048516119168 (same two Q13 starts), gaia-dr3-2102601072256824704 (Q2 safe mode, Q6 start, Q15 gap), kic-008620445 (Q12 ramp + Q2 step) |
| single-cadence or attitude-tweak excursions (BKJD 200.33, Q2 cadence 4473) | kic-006515744, kic-008620445 |
| saturated TPF host bleed-column slosh on one channel | kic-005950408 (channel 25, KIC 5950394, Kp 8.45) |
| pool-identity failure: target is a Kepler target, not a background dwarf | gaia-dr3-2102601072256824704 (KIC 7269549 at 0.32 arcsec, 0.74 arcsec binary), gaia-dr3-2127687048516119168 (KIC 9459487 at 0.79 arcsec, the TPF host itself) |

Two audit facts a referee will want up front:

1. **All eight entered the survivor set through the v1.5 relaxation.** Every one has `stage2_pass_v14 = False` with `stage2_fails_v14 = tls_sde` (TLS SDE 4.0-8.8; the four batch-A objects had 11.6-22.2). The v1.5 removal of `tls_sde_min` (post-hoc, 2026-10-17 ~15:55Z) raised the survivor count from 11 to 23; the twelve admitted that way include all eight here, and all eight are artefacts. The registered v1.4 cut was doing work even if its statistic is uncalibrated.
2. **The inverted-flux null was already level with or above the real signal at stage 1 for five of eight** (inverted / real BLS SDE: kic-005609180 11.09 / 10.93; kic-008620445 12.26 / 10.74; kic-005950408 12.12 / 9.83; gaia-dr3-2127687048516119168 13.31 / 10.01; kic-005866176 10.55 / 10.01; the other three 8.84 / 12.91, 8.40 / 11.37, 8.84 / 9.49). Those inverted runs were rejected only on duty-cycle / SNR cuts because the inverted winners were 7 h boxes at P < 1 d. Batch A had real / inverted 89 / 14, 46 / 13, 101 / 20, 36 / 10. A rule "real SDE must exceed the inverted SDE by a margin" (post-hoc) separates the two batches completely.

## Summary table (referee-sized)

Pipeline numbers from `runs/poolA_main/vet/<fname>/vet.json` (stage 2, rules v1.7; box depth, `snr_emp_pass2`, `event_weighted_depth_ppm`, `event_chi2_const/event_dof`), `runs/poolA_main/results.jsonl` (stage-1 BLS SDE real / inverted), `runs/poolA_main/pixel_v16/<fname>/pixel.json`; red-team and re-derivation numbers from `runs/poolA_main/adversarial/<fname>/`. "Residual" is the depth at the fixed ephemeris after the artefact cadences are removed, as measured by the red team (or, for the two provisional objects, the synthesis spot check); it is a post-hoc statistic in every case.

| rank | fname | P (d) | box depth (ppm) / dur (h) | T14max (h) | N events (usable) | SNR_emp s2 | BLS SDE real / inverted (s1) | TLS SDE | event-wtd depth vs box; chi2/dof | re-derivation blind | pixel v1.6 | red team | residual after artefact removal | identity / novelty (2026-10-17) | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | kic-005866176 | 29.4330 | 481 / 1.95 | 5.7 | 31 | 7.5 | 10.01 / 10.55 | 8.5 | 250 vs 481; 67/30 | not top peak (3:2 alias wins); collapses to SDE 3.5-4.0 under 0.5-2 d detrending | inconclusive (16/17 edge/outside) | KILL: two post-gap PSF/SAP anti-correlated excursions (ch 17) on an off-stamp star | 272 +- 48 ppm PSF, 276 +- 55 SAP (SNR 5.7 < 7.1); BLS 3.59 vs inverted 5.94 | KIC 5866176, not a Kepler target; host KIC 5866192 flat; nothing in KOI/TCE/EB/PS/TOI/SIMBAD | rejected |
| 2 | kic-005609180 | 21.1063 | 656 / 3.0 | 5.0 | 61 | 12.3 | 10.93 / 11.09 | 8.8 | 318 vs 656; 87/60 | 21.1 d SDE 1.66 (rank > 40) with 1 d post-gap mask; inverted 8.31 | not run (none warranted) | KILL: three channel-17 post-gap ramps; flux columns disagree by up to 100 ppt | 201-220 +- 65 ppm PSF only; SAP 97-113 +- 63; absent in season 2 | KIC 5609180, never a Kepler target; host KIC 5609183 clean; NEA/SIMBAD/B&S 2025/Moquin 2026 null | rejected |
| 3 | kic-004474607 | 32.9427 | 2141 / 1.95 | 5.6 | 17 | 11.6 | 12.91 / 8.84 | 5.7 | 1025 vs 2141; 45/16 | best period 32.98 d = Q13 downlink spacing; blind 1-50 d SDE 6.7 < inverted 7.7 | inconclusive; Q13-only dchi2 183 for an off-stamp source | KILL: two Q13 segment-start ramps (ch 49) carry the signal; control box at every segment start reproduces it | 884 +- 191 ppm PSF / 603 +- 166 SAP at fixed ephemeris; masked local SDE 2.2 | KIC 4474607; TPF host KIC 4474637 = KOI-3829.01 FP EB (3.887 d) does not leak (13 +- 21 ppm) | rejected |
| 4 | gaia-dr3-2127687048516119168 | 32.8687 | 351 / 1.5 | 5.4 | 38 | 9.1 | 10.01 / 13.31 | 6.2 | 137 vs 351; 101/37 | top peak 16.40 d (P/2), SDE 8.1 vs inverted 6.0; at P/2 115 +- 17 ppm still 67% two events | inconclusive (1 on-target quarter, Q3; sum dchi2 11) | **none** (provisional) | 63 +- 40 ppm PSF, 38 +- 36 SAP without epochs 31/32 (80% of depth; Q13 starts 1182.74 / 1215.72; PSF/SAP anti-correlated at 1215.76-1215.80) | **KIC 9459487 (Kp 15.25, Q1-Q17, the TPF host) at 0.79 arcsec; RUWE 1.94; pool row kic = 0**; no KOI/TCE within 180 arcsec | rejected (provisional) |
| 5 | kic-008620445 | 28.0717 | 1476 / 0.75 | 5.2 | 17 (of 44 epochs) | 11.3 | 10.74 / 12.26 | 7.6 | 243 vs 1476; 27/16 | top peak 40.277 d (SDE 9.7 vs inverted 8.1), 28.07 d not top-1 nor an alias; per-year 5121/38/2750/502 ppm | inconclusive (13/17 edge/outside; dchi2 289 for an alternative over 2 quarters) | **none** (provisional) | 329 ppm PSF / 182 SAP without epochs 2 and 35 (76% of depth: 2.4% step at the Q2 attitude tweak BKJD 200.29-200.36; Q12 post-gap ramp at 1126.65; 926.37 d = 33 x 28.072 = 23 x 40.277); events.csv 243 +- 115 ppm | KIC 8620445, never a Kepler target; host KIC 8620443; no KOI/TCE | rejected (provisional) |
| 6 | gaia-dr3-2102601072256824704 | 32.5571 | 191 / 3.375 | 5.6 | 43 | 8.2 | 9.49 / 8.84 | 4.0 | 118 vs 191; 108/42 | trimmed (0.3 d edges) full search SDE 1.4, rank 20104; untrimmed inverted 10.74 > real 10.60 | inconclusive; scaled-image hypothesis preferred (dchi2 86) | KILL: three segment-edge events (Q2 safe mode, Q6 start, Q15 gap) + target identity | 62 +- 20 ppm PSF (SNR 3.1), 50 +- 21 SAP | **KIC 7269549 (Kp 14.76, Q1-Q17) at 0.32 arcsec, unresolved 0.74 arcsec twin pair (RUWE 3.9); pool row kic = 0** | rejected |
| 7 | kic-005950408 | 38.6067 | 4104 / 0.75 | 6.0 | 7 (of 31 epochs) | 10.6 | 9.83 / 12.12 | 7.1 | 177 vs 4104; 1.8/6 (wrong events) | not recovered (rank 641); rank 114,626 without channel 25 | inconclusive (1 quarter) | KILL: channel-25 bleed-column slosh from saturated host KIC 5950394 (Kp 8.45); six cadences carry everything | 53 ppm outlier-free; seasons 0/1/2 -242 / +287 / +710 ppm | KIC 5950408; no KOI/TCE/EB/TOI within 3 arcmin; host giant flat | rejected |
| 8 | kic-006515744 | 36.2075 | 12,763 / 0.75 | 5.7 | 1 | 13.8 | 11.37 / 8.40 | 5.4 | 357 vs 12,763; dof 0 | 36.2 d below spectrum median after 5-MAD clip (SDE 7.4); inverted 9.7 vs 12.7 real | 0 of 8 quarters usable | KILL: four single-cadence PSF excursions (Q2 attitude tweak BKJD 200.33; +89 sigma pixel event BKJD 1370.09) | nothing (10 of 15 epochs consistent with zero) | KIC 6515744, not a Kepler target; nothing in catalogues | rejected (already stage-2 FAIL under v1.7 n_events and v1.4 tls_sde) |

Duration physics, for the record: six of eight have `duration_ratio` (box duration / T14max at b = 0) of 0.11-0.35; the 0.75 h boxes at 28-39 d (kic-008620445, kic-005950408, kic-006515744) need b > 1 on their 0.85-0.91 Rsun hosts, i.e. no planetary geometry produces them. This was in `vet.json` for every object before the red teams ran.

## Why this ranking

The order is "least dead first": how much of a residual at the fixed ephemeris survives removal of the artefact cadences, and whether that residual appears in more than one flux product. None of it changes the status.

1. **kic-005866176** is the only object with a residual that is (a) present in PSF, SAP and PSF_FLAT flux, (b) consistent across the three clean seasons (299 / 334 / 325 ppm) and (c) not reached by 1500 random phases. The red team still kills it, correctly: the ephemeris was selected on artefact-dominated data, the residual fails the registered `snr_emp_pass2_min = 7.1` (SNR 5.7), it is not blind-recoverable, and the star is 0.2-0.4 px outside the stamp in every quarter so pixel localisation is impossible. Kepler cannot lift it.
2. **kic-005609180**'s residual is PSF-only (SAP 1.5-1.8 sigma where CROWDSAP demands >= 95% of the PSF depth), absent in one of four seasons, carried by three events, and matched by an equal-amplitude brightening 3 h earlier.
3. **kic-004474607**'s 884 +- 191 ppm residual looks large but is a post-hoc measurement at a period that equals the Q13 downlink spacing; the masked light curve has no peak there (local SDE 2.2) and its best real peak is below the inverted null.
4. **gaia-dr3-2127687048516119168** and 5. **kic-008620445** have no red-team record, so they are placed below the three documented kills with residuals, but above the fully dead objects because the formal step is missing, not because anything survives: their residuals are 63 +- 40 and ~250-330 +- 115 ppm, their periods are gap-comb teeth, and the Gaia-only object is a Kepler target.
6-8. **gaia-dr3-2102601072256824704**, **kic-005950408**, **kic-006515744** have nothing left and, for the first, a second independent kill (target identity).

## Cross-object findings that matter more than any single object

- **The same cadences recur across unrelated stars and channels.** BKJD 200.33 (Q2 attitude tweak) drives kic-006515744 and kic-008620445. The Q13 segment starts at BKJD 1182.737 and 1215.718 drive kic-004474607 (channel 49), gaia-dr3-2127687048516119168 (different module) and kic-008620445 epoch 37. Their 32.98 d spacing, and the monthly cadence generally, is why 349 of 438 real first-screen passes sit in the 24-40 d class and why five of these eight have P = 28-39 d. A pool-wide bad-cadence list derived from the KBONUS sources themselves would remove these at stage 1.
- **PSF-vs-SAP anti-correlation is a clean artefact discriminator.** In kic-005866176, gaia-dr3-2127687048516119168, kic-005609180 and kic-005950408 the PSF flux dips while SAP flux rises (or vice versa) at the carrying cadences. No loss of light does that.
- **events.csv excluded the carrying epochs** for kic-005950408 and kic-008620445 (n_in >= 2 selection), so `event_chi2_const` tested the wrong events; `event_weighted_depth_ppm` nevertheless disagreed with the box depth by a factor 2-23 for all eight. That disagreement is already in the pipeline output and should be a rule.
- **Two "background dwarfs" are Kepler targets.** Both Gaia-only rows (kic = 0, empty kic_sep) are within 0.8 arcsec of a KIC star observed Q1-Q17; one is a 0.74 arcsec near-twin binary. The pool builder needs a positional cross-match against keplerstellar, not only the KBONUS kic column.

## Caveats that travel with every object

- **Post-hoc rules.** Stage 2 (v1.3), v1.4, v1.5, the v1.6 pixel rule and v1.7 were all registered after results had been seen (`pipeline/prereg_poolA.json`, `revision_log`). Every stage-2 and pixel number above is post-hoc relative to the stage-1 registration. The kills themselves do not depend on any threshold: they rest on cadence-level anatomy and on the inverted-flux null, which was registered.
- **Red-team measurements are own-code products**, not pipeline products with a run id; they are reproducible from the scripts and manifests in each `adversarial/<fname>/` directory. The two synthesis spot checks are likewise own code (`synthesis_spotcheck/*.py`, `MANIFEST.sha256`), read-only on the data.
- **Human inspection** of `fold_events.png` / `diffimages.png` is recorded by the red teams for six objects and has not been done by a human for the two provisional ones.
- **Novelty is moot** for artefacts, but all catalogue statements are dated 2026-10-17 (`data/known/manifest.json`, 15:46Z) and the Kepler Bonus team's own candidate list is unpublished.

## Per-survivor next steps

Listed in the ranked JSON for each object. In brief: retire all eight; write formal red-team records for kic-008620445 and gaia-dr3-2127687048516119168 (narrow BLS plus inverted null with 0.3-0.5 d excised after every segment start and the AttitudeTweak cadence masked; depths with and without the carrying epochs in PSF and SAP; host-TPF QUALITY flags at BKJD 200.33 for kic-008620445); re-classify the two Gaia-only pool rows as Kepler targets; run no pixel or TRICERATOPS work on any of them.

## Pipeline recommendations (post-hoc; to be entered in `revision_log` as v1.8 and passed by the positive controls before any pool rerun)

1. Mask 0.3-1.0 d after every gap / segment start and the cadence after any Kepler AttitudeTweak / NoFinePoint flag before detrending; require the masked BLS to still peak at the same period and to exceed the inverted-flux null.
2. Report and gate on real-vs-inverted SDE directly (all eight here have inverted / real >= 0.68; batch A <= 0.28).
3. Include every epoch with in-box data in `events.csv`; require event-weighted and box depth to agree within 2 sigma; enforce `event_chi2_const` with a red-noise-aware threshold.
4. Per-quarter PSF/SAP depth sign agreement and SAP/PSF ratio consistent with CROWDSAP; per-quarter RED_CHI2 / raw-MAD screen; saturated-TPF-host flag (Kp < ~10) with bleed-column distance.
5. Pool hygiene: positional cross-match of every row against keplerstellar (flag < 2 arcsec), RUWE / ipd_frac_multi_peak flags, off-stamp / PSFFRAC < 0.7 flag as non-localisable.
6. Reinstate a calibrated narrow-window detection statistic in place of the dropped `tls_sde_min`, or at minimum report v1.4 and v1.5 verdicts side by side as the revision_log already requires.

## Before anything goes public (CLAUDE.md checklist, status for this batch)

Not applicable: no object in this batch passes stage 2 in substance, none is pixel-localised, none survives red team or re-derivation. Nothing from batch B should appear in any public statement except, if desired, as worked examples of the false-positive classes above, with the post-hoc labels intact.

Evidence paths: `/home/user/ZyushgPlanet/runs/poolA_main/vet/<fname>/{vet.json,events.csv,fold_events.png}`, `/home/user/ZyushgPlanet/runs/poolA_main/pixel_v16/<fname>/pixel.json`, `/home/user/ZyushgPlanet/runs/poolA_main/results.jsonl`, `/home/user/ZyushgPlanet/runs/poolA_main/vet_summary.csv`, `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/<fname>/{redteam.md,redteam_verdict.json,rederive/,redteam*/,synthesis_spotcheck/}`, `/home/user/ZyushgPlanet/pipeline/prereg_poolA.json` (revision_log), `/home/user/ZyushgPlanet/data/known/manifest.json`.