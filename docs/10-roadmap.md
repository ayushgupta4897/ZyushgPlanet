# 10 · Where the search stands, and the road from here

_Last updated 2026-10-17 (UTC). Numbers trace to `runs/poolA_main/` (Modal, ranks 101–6465), `runs/poolA_pilot/` (local, ranks 1–100), `runs/controls*/`, and `pipeline/prereg_poolA.json` (revision log v1.0–v1.5)._

## 1. What was done today

| Step | Result |
|---|---|
| Pool definition (frozen, hashed) | 6,465 Kepler Bonus background dwarfs: never Kepler targets, not in Kuber's 1,800, PSF fraction > 0.5, crowding > 0.7, ≥ 8 quarters, Gaia-classified dwarfs, G < 17.5 (`targets/poolA_kbonus_bg_dwarfs_ranked.csv`). "Never searched" means: no Kepler-pipeline search and no published transit search; the Kepler Bonus team itself is mining the same data (`docs/04` §3). |
| Pre-registration | First-screen thresholds frozen before any pool star was searched (v1.0–v1.2). Stage-2 criteria were registered *after* a 208-star partial look at stage-1 results (v1.3), one criterion was relaxed after a control failed it (v1.4), and a post-hoc revision replaced an ill-defined statistic (v1.5). Every step is in the revision log with its commit; stage-2 results are post-hoc relative to stage 1. |
| Positive controls | KIC 5689121 (Kuber's candidate) passes stage 1, stage 2 and pixel localisation (13/13 quarters on target). Kepler-138 c and d pass stage 1 (c as the first signal, d as the second); c passes every stage-2 criterion except the intended known-object veto. The first version of the screen *failed* the controls (three statistic bugs); `runs/controls_v1.0_failed/` keeps that run. One positive control is a check, not a validation. |
| Stage 1 screen | Modal, 6,365 stars: 438 real vs 95 inverted passes, 0 errors, 31 skipped (no usable cadences), mean 125 s/star, wall time 2 h 17 m (`modal_status.json`). Local pilot, 100 stars: 10 vs 2. Pool total: **448 real vs 97 inverted** over 6,465 (6.9 % vs 1.5 %). |
| What the null means | Inverting the flux turns dips into bumps, so the inverted rate bounds *sign-symmetric* noise and variability only. By best-peak period: 0.5–2 d 24 vs 26; 2–10 d 33 vs 36; 10–20 d 12 vs 4; **20–40 d 379 vs 31**. The excess is a dip-shaped systematic (post-gap thermal ramps at Kepler's monthly downlinks), not a count of astrophysical signals. At stage 1 the SDE ≥ 9 and SNR ≥ 7.1 cuts alone pass 1,732 real vs 145 inverted stars; the duty-cycle and single-event cuts do most of the suppression. |
| Stage 2 vetting | 448 survivors re-derived with transits masked, per-event local depths with gap proximity, V-shape proxy, variability, duration plausibility, narrow-window TLS, ephemeris matching → **11 survive under the registered v1.4 rules, 23 under v1.5** (v1.5 drops the narrow-window TLS SDE, adds empirical SNR ≥ 7.1 at stage 2 and a two-sided secondary test). Per-quartile depths are computed and shown but are not a criterion. |
| Pixel localisation | Registered verdict rules (v1.5): quarters with the target at the stamp edge, < 3 events or < 2 significant pixels are inconclusive; on-target needs ≥ 5 conclusive quarters with the deficit within 2σ of the target and closer to it than to any other Gaia source in ≥ 70 % of them. Being re-run on both controls and all 23 survivors with the hardened estimator (the first pass, kept under `pixel/`, used no edge test and is descriptive only). |
| Adversarial pass | Independent from-scratch re-derivation (no shared code) + red-team review per survivor: **running on all 23**. |

Estimator definitions: stage 1 odd/even and secondary use the BLS weighted-fit depths (`depth_odd/even/phased`) with 60 random off-transit boxes; stage 2 uses unweighted box means on the masked re-detrend with 80 boxes; "SNR" below is `snr_emp_pass2` (stage-2 empirical), "depth" is `depth_pass2_ppm`, "events" is `n_events`; TLS formal SDE/SNR are in each `vet.json`.

## 2. The survivors (23 under v1.5; the 11 v1.4 survivors are marked ★; two v1.4 survivors drop under v1.5 and are listed at the end)

| fname | Host (Gaia) | P (d) | depth (ppm) | dur (h) | SNR | events (% positive) | R_p if real | Pixel (first pass, descriptive) | Pixel (v1.5 rules) | Adversarial |
|---|---|---|---|---|---|---|---|---|---|---|
| ★ kic-003831969 | G dwarf, G 16.5, 1,969 pc | 15.3134 | 997 | 3.8 | 23.8 | 86 (92 %) | ~2.8 R⊕ | 2.9″, 9/12 q, SNR 3.1 | pending | pending |
| ★ kic-010213563 | K dwarf, G 17.1, 1,673 pc | 36.7671 | 2,242 | 5.6 | 20.6 | 36 (100 %) | ~4.1 R⊕ | 2.3″, 8/13, SNR 5.7 | pending | pending |
| ★ kic-006699887 | K dwarf, G 16.2, 859 pc | 15.1136 | 1,023 | 3.8 | 20.0 | 86 (100 %) | ~2.3 R⊕ | 1.1″, 12/13, SNR 7.1 | pending | pending |
| ★ kic-011138423 | early-M, G 17.4, 894 pc | 4.8945 | 1,316 | 2.0 | 15.2 | 253 (87 %) | ~2.1 R⊕ | 3.1″, 10/12, SNR 3.7 | pending | pending |
| kic-006515744 | — | 36.2075 | 12,763 | 1.1 | 13.8 | — (100 %) | ~10.5 R⊕ | — | pending | pending |
| kic-005609180 | — | 21.1063 | 656 | 3.4 | 12.3 | — (75 %) | ~2.6 R⊕ | — | pending | pending |
| kic-008955698 | — | 28.3811 | 1,653 | 5.6 | 12.0 | — (86 %) | ~2.3 R⊕ | — | pending | pending |
| kic-004474607 | — | 32.9427 | 2,141 | 2.3 | 11.6 | — (88 %) | ~4.4 R⊕ | — | pending | pending |
| kic-008620445 | — | 28.0717 | 1,476 | 0.9 | 11.3 | — (71 %) | ~3.6 R⊕ | — | pending | pending |
| kic-005950408 | — | 38.6067 | 4,104 | 0.9 | 10.6 | — (86 %) | ~6.4 R⊕ | — | pending | pending |
| ★ kic-006448637 | K dwarf, G 15.7, 1,006 pc | 35.5873 | 846 | 4.5 | 10.0 | 26 (100 %) | ~2.4 R⊕ | inconclusive (target outside stamp in 3 q; 2–3 events/q) | pending | pending |
| ★ kic-006367225 | G dwarf, G 17.2, 4,803 pc | 12.9715 | 3,034 | 2.2 | 9.9 | 100 (88 %) | ~5.3 R⊕ | another source closer in 12/14 q | pending | pending |
| ★ kic-003963519 | K dwarf, G 17.3, 4,365 pc | 6.8156 | 1,124 | 2.4 | 9.6 | 148 (75 %) | ~2.8 R⊕ | 3.6″, 8/13, SNR 3.1 | pending | pending |
| kic-008636232 | — | 36.3024 | 1,369 | 3.4 | 9.2 | — (77 %) | ~3.4 R⊕ | — | pending | pending |
| gaia-dr3-2127687048516119168 | — | 32.8687 | 351 | 1.8 | 9.1 | — (76 %) | ~1.7 R⊕ | — | pending | pending |
| ★ kic-011285493 | G dwarf, G 16.8, 2,858 pc | 10.7910 | 796 | 4.5 | 8.7 | 96 (72 %) | ~2.7 R⊕ | 7.9″, SNR 2.5 | pending | pending |
| kic-010782709 | — | 31.1781 | 468 | 2.6 | 8.4 | — (71 %) | ~2.3 R⊕ | — | pending | pending |
| kic-011236638 | — | 26.3823 | 323 | 3.5 | 8.3 | — (83 %) | ~1.7 R⊕ | — | pending | pending |
| kic-006184850 | — | 34.7932 | 3,622 | 1.4 | 8.2 | — (88 %) | ~4.7 R⊕ | — | pending | pending |
| gaia-dr3-2102601072256824704 | — | 32.5571 | 191 | 3.8 | 8.2 | — (72 %) | ~1.3 R⊕ | — | pending | pending |
| ★ kic-009655224 | G dwarf, G 17.1, 2,235 pc | 6.7035 | 638 | 2.4 | 8.2 | 125 (75 %) | ~2.3 R⊕ | 4.6″, 6 q only | pending | pending |
| kic-009019826 | — | 36.3047 | 729 | 4.3 | 7.9 | — (80 %) | ~1.7 R⊕ | — | pending | pending |
| kic-005866176 | — | 29.4330 | 481 | 2.2 | 7.5 | — (81 %) | ~2.3 R⊕ | — | pending | pending |
| ★ kic-009520080 (drops under v1.5: −4.0σ brightening at phase 0.5) | K dwarf, G 17.5, 1,657 pc | 2.3770 | 1,569 | 1.1 | 28.3 | 468 (86 %); quartile depths 1957/1601/1091/1566 | ~3.0 R⊕ | 1.7″, 13/17, SNR 5.2 | pending | pending |
| ★ kic-001577049 (drops under v1.5: stage-2 empirical SNR 1.9) | K dwarf, G 17.3, 2,630 pc | 0.5379 | 1,324 | 1.1 | 1.9 | 771 (75 %) | ~3.1 R⊕ | deficit 28″ away | pending | pending |

Host columns for the twelve v1.5-only rows are filled from the target list once their pixel maps are in. Fourteen of the 23 lie at 21–39 d, inside the systematic band identified in §1; the per-event gap test and the red team are what separate a real 30-day transit from a post-gap ramp. Nothing in this table is a planet. The honest label for anything that clears the adversarial pass, the Gaia/RUWE check and a documented TRICERATOPS run is **planet candidate** (`CLAUDE.md`, Language); until then it is a *signal*.

## 3. What a survivor still needs before anyone says "candidate" in public

1. Adversarial pass clean: independent re-derivation agrees on period, depth and shape; red team finds no killer.
2. Pixel verdict on-target under the registered rules, with the deficit-position error quoted.
3. Deep single-target run (`pipeline/run.py`): full-range TLS, pass-2 fit with Kepler-band limb darkening, per-season dilution (SAP vs PSF depth vs the target's flux share), Gaia neighbours within 20″ with flux ratios.
4. Host star: Gaia DR3 parameters, RUWE, non-single-star flags; radius from colour + parallax, not a guess; implied planet radius and transit-derived stellar density consistency.
5. False-positive probability: TRICERATOPS with the resolved neighbours and their flux fractions passed in, ≥ 4 seeds, scenario breakdown.
6. Novelty with dates (`data/known/manifest.json`): KOI, DR25 TCE, Kepler EB catalogue, confirmed planets, Bienias & Szabó 2025, the Kepler Bonus team's own candidates (Moquin et al.), SIMBAD.
7. Evidence pack: one README with the ten key numbers, the plots, the JSONs, SHA-256 manifest, Zenodo deposit.

## 4. The honest answer to "how close are we"

- **To a candidate**: close for the handful of on-target, high-SNR, short-to-medium-period signals (the 15-day K dwarfs and the 4.9-day M dwarf are the natural first cases) if the adversarial pass and TRICERATOPS agree. That would be what Kuber has, found with a registered, control-tested pipeline over 3.6× more stars and with the limits of the null control stated rather than glossed.
- **To a validated planet**: far, and for the same structural reason as Kuber: hosts at 900–2,000 pc with G 16–17.5 cannot get the high-resolution imaging or spectroscopy that pushes the bound-companion scenario below 1.5 %. The realistic end state for Pool A finds is a refereed candidate note (AJ/MNRAS-class, with a professional co-author) plus an evidence pack, not a Kepler-NNNN name.
- **To a validatable find**: that needs Pool B. The verified landscape ranks faint nearby M dwarfs in TESS (T 13.5–16, < 100 pc) as the regime where a solo pipeline can still find something *and* ground-based 1-m telescopes can follow it up; the fresh-sector window (S108–109 cutouts already on TESScut, no TOIs yet from S103–107) is the time-critical bonus.

## 5. Next milestones

| # | Milestone | Depends on |
|---|---|---|
| 1 | Pixel verdicts under registered rules + adversarial verdicts for all 23; candidate-grade list | running |
| 2 | Deep runs + TRICERATOPS for candidate-grade signals; evidence packs under `candidates/<fname>/` | 1 |
| 3 | Draft candidate note; identify a professional collaborator (Kepler Bonus authors are the natural first contact; `docs/06`) | 2 |
| 4 | Pool B target list: TESS M dwarfs T 13.5–16 within 100 pc, north first; light-curve supply check (TESS-SPOC / QLP / TGLC / own TESScut photometry) | `docs/04b` §9 |
| 5 | Pool B pre-registration, controls (known TOIs around faint M dwarfs), pilot | 4 |
| 6 | Modal batch for Pool B; same triage/vetting/pixel stack with TESS-specific vetting (LEO-Vetter, centroid offsets, scattered light) | 5 |
| 7 | Rotate the Modal token; `modal app stop zyushg-planet-screen` when no further batches are planned | user |
