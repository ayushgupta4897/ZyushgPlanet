# 10 · Where the search stands, and the road from here

_Last updated 2026-10-17 (UTC). Numbers trace to `runs/poolA_main/`, `runs/poolA_pilot/`, `runs/controls/`._

## 1. What was done today

| Step | Result |
|---|---|
| Pool definition (frozen, hashed) | 6,465 Kepler Bonus background dwarfs: never Kepler targets, not in Kuber's 1,800, PSF fraction > 0.5, crowding > 0.7, ≥ 8 quarters, Gaia-classified dwarfs, G < 17.5 (`targets/poolA_kbonus_bg_dwarfs_ranked.csv`) |
| Pre-registration | `pipeline/prereg_poolA.json` v1.0 → v1.2 with a dated revision log; every change made before pool stars were judged |
| Positive controls | KIC 5689121 (Kuber's candidate) and Kepler-138 c + d pass stage 1, stage 2 and pixel localisation; the neighbour and all inverted runs fail. The first version of the screen *failed* the controls (three statistic bugs); `runs/controls_v1.0_failed/` keeps that run |
| Stage 1 screen | 6,465 stars on Modal (6,365) + locally (100); **438 real** vs **95 inverted** first-screen passes; 0 errors; mean 125 s/star |
| Stage 2 vetting | 438 survivors re-derived with transits masked, per-event local depths, gap proximity, per-quartile consistency, V-shape, variability, duration plausibility, TLS, ephemeris matching → **11 survive** |
| Pixel localisation | All 11 run on the originating target-pixel files, per quarter |
| Adversarial pass | Independent from-scratch re-derivation (no shared code) + red-team review per survivor: **in progress** |

The stage-1 null rate (95/6,365 = 1.5 % inverted passes against 6.9 % real passes) is the honest statement of what the registered thresholds mean: roughly one in five first-screen passes is noise or a systematic, and the rest still have to be explained one by one. The dominant failure classes at stage 2 were (i) events piling up next to Kepler's monthly downlink gaps and quarter edges, with 24–40 day "periods", and (ii) periodic variables and eclipsing binaries that pass in both polarities.

## 2. The eleven stage-2 survivors

| fname | Host (Gaia) | P (d) | depth (ppm) | dur (h) | SNR | events | R_p if real | Pixel localisation | Adversarial verdict |
|---|---|---|---|---|---|---|---|---|---|
| kic-006699887 | K dwarf, G 16.2, 859 pc | 15.1136 | 1,023 | 3.8 | 20.0 | 86 (100 % positive) | ~2.3 R⊕ | **on target**: 1.1″, 12/13 quarters, peak SNR 7.1 | pending |
| kic-009520080 | K dwarf, G 17.5, 1,657 pc | 2.3770 | 1,569 | 1.1 | 28.3 | 468 | ~3.0 R⊕ | **on target**: 1.7″, 13/17, SNR 5.2 (secondary −4σ: brightening at phase 0.5, needs explaining) | pending |
| kic-010213563 | K dwarf, G 17.1, 1,673 pc | 36.7671 | 2,242 | 5.6 | 20.6 | 36 (100 %) | ~4.1 R⊕ | on target: 2.3″, 8/13, SNR 5.7 | pending |
| kic-003831969 | G dwarf, G 16.5, 1,969 pc | 15.3134 | 997 | 3.8 | 23.8 | 86 | ~2.8 R⊕ | plausible: 2.9″, 9/12, SNR 3.1 | pending |
| kic-011138423 | early-M, G 17.4, 894 pc | 4.8945 | 1,316 | 2.0 | 15.2 | 253 | ~2.1 R⊕ | plausible: 3.1″, 10/12, SNR 3.7 | pending |
| kic-003963519 | K dwarf, G 17.3, 4,365 pc | 6.8156 | 1,124 | 2.4 | 9.6 | 148 | ~2.8 R⊕ | plausible: 3.6″, 8/13, SNR 3.1 | pending |
| kic-006367225 | G dwarf, G 17.2, 4,803 pc | 12.9715 | 3,034 | 2.2 | 9.9 | 100 | ~5.3 R⊕ | weak: another source closer in 12/14 quarters | pending |
| kic-009655224 | G dwarf, G 17.1, 2,235 pc | 6.7035 | 638 | 2.4 | 8.2 | 125 | ~2.3 R⊕ | weak: 6 quarters usable, 4.6″ | pending |
| kic-011285493 | G dwarf, G 16.8, 2,858 pc | 10.7910 | 796 | 4.5 | 8.7 | 96 | ~2.7 R⊕ | weak: 7.9″, SNR 2.5 | pending |
| kic-006448637 | K dwarf, G 15.7, 1,006 pc | 35.5873 | 846 | 4.5 | 10.0 | 26 | ~2.4 R⊕ | **off target**: no coherent deficit | pending |
| kic-001577049 | K dwarf, G 17.3, 2,630 pc | 0.5379 | 1,324 | 1.1 | 1.9 | 771 | ~3.1 R⊕ | **off target**: deficit 28″ away, 1,510 ppm in aperture | pending |

Reading guide: "on target" means the per-quarter missing-light maps put the deficit closer to this star than to any other Gaia source and a point source at the target fits better than a scaled copy of the stamp (the test that put Kuber's signal on his star in 13/13 quarters). Nothing here is a planet. The honest label for the best three, if the adversarial pass agrees, is **planet candidate**.

## 3. What a survivor still needs before anyone says "candidate" in public

1. Adversarial pass clean: independent re-derivation agrees on period, depth and shape; red team finds no killer.
2. Deep single-target run (`pipeline/run.py`): full-range TLS, pass-2 fit with Kepler-band limb darkening, per-season dilution (SAP vs PSF depth vs the target's flux share), Gaia neighbours within 20″ with flux ratios.
3. Host star: Gaia DR3 parameters, RUWE, non-single-star flags; radius from colour + parallax, not a guess; implied planet radius and transit-derived stellar density consistency.
4. False-positive probability: TRICERATOPS with the resolved neighbours and their flux fractions passed in, ≥ 4 seeds, scenario breakdown.
5. Novelty with dates: KOI, DR25 TCE, Kepler EB catalogue, confirmed planets, Bienias & Szabó 2025, the Kepler Bonus team's own candidates (Moquin et al.), SIMBAD.
6. Evidence pack: one README with the ten key numbers, the plots, the JSONs, SHA-256 manifest, Zenodo deposit.

## 4. The honest answer to "how close are we"

- **To a candidate**: close. If two or three of the on-target survivors clear the adversarial pass and TRICERATOPS, we will have what Kuber has, found with a registered, control-tested pipeline over 3.6× more stars, with the false-alarm rate measured rather than guessed.
- **To a validated planet**: far, and for the same structural reason as Kuber: hosts at 900–1,700 pc with G 16–17.5 cannot get the high-resolution imaging or spectroscopy that pushes the bound-companion scenario below 1.5 %. The realistic end state for Pool A finds is a refereed candidate note (AJ/MNRAS-class, with a professional co-author) plus an evidence pack, not a Kepler-NNNN name.
- **To a validatable find**: that needs Pool B. The verified landscape ranks faint nearby M dwarfs in TESS (T 13.5–16, < 100 pc) as the regime where a solo pipeline can still find something *and* ground-based 1-m telescopes can follow it up; the fresh-sector window (S108–109 cutouts already on TESScut, no TOIs yet from S103–107) is the time-critical bonus.

## 5. Next milestones

| # | Milestone | Depends on |
|---|---|---|
| 1 | Adversarial verdicts for all 11; candidate-grade list | running |
| 2 | Deep runs + TRICERATOPS for candidate-grade signals; evidence packs under `candidates/<fname>/` | 1 |
| 3 | Draft candidate note; identify a professional collaborator (Kepler Bonus authors are the natural first contact; `docs/06`) | 2 |
| 4 | Pool B target list: TESS M dwarfs T 13.5–16 within 100 pc, north first; light-curve supply check (TESS-SPOC / QLP / TGLC / own TESScut photometry) | `docs/04b` §9 |
| 5 | Pool B pre-registration, controls (known TOIs around faint M dwarfs), pilot | 4 |
| 6 | Modal batch for Pool B; same triage/vetting/pixel stack with TESS-specific vetting (LEO-Vetter, centroid offsets, scattered light) | 5 |
| 7 | Rotate the Modal token; `modal app stop zyushg-planet-screen` when no further batches are planned | user |
