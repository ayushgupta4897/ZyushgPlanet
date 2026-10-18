# 10 · Where the search stands, and the road from here

_Last updated 2026-10-18 (UTC). Numbers trace to `runs/poolA_main/` (Modal, ranks 101–6465), `runs/poolA_pilot/` (local, ranks 1–100), `runs/controls*/`, `pipeline/prereg_poolA.json` (revision log v1.0–v1.8), `runs/poolA_main/adversarial/` (re-derivations, red-team verdicts, `synthesis_batch{A,B,C}.md`), `report/data.json` and the evidence packs in `candidates/<fname>/`._

## 1. What was done

**Funnel (Pool A, pre-registration v1.8):** 6,465 targets → **448 real vs 97 inverted (2 − f)** first-screen passes → **22** stage-2 passes → **5** localise on target in pixel tests → **4** survive independent re-derivation and red team → **3 candidates**.

| Step | Result |
|---|---|
| Pool definition (frozen, hashed) | 6,465 Kepler Bonus background dwarfs: never Kepler targets, not in Kuber's 1,800, PSF fraction > 0.5, crowding > 0.7, ≥ 8 quarters, Gaia-classified dwarfs, G < 17.5 (`targets/poolA_kbonus_bg_dwarfs_ranked.csv`). "Never searched" means no Kepler-pipeline search and no published transit search; the Kepler Bonus team itself is mining the same data (`docs/04` §3). |
| Pre-registration | First-screen thresholds frozen before any pool star was searched (v1.0–v1.2). Stage-2 criteria registered *after* a 208-star partial look at stage-1 results (v1.3); one criterion relaxed after a control failed it (v1.4); post-hoc revisions replaced an ill-defined statistic (v1.5), replaced a pixel rule that failed its control (v1.6), added a minimum-event rule (v1.7) and a season-stacked pixel test (v1.8). Every step is in the `revision_log` with its motivation; stage-2 and pixel results are post-hoc relative to stage 1. |
| Positive controls | KIC 5689121 (Kuber's candidate) passes stage 1, stage 2 and both pixel tests. Kepler-138 c and d pass stage 1; c passes every stage-2 criterion except the intended known-object veto. The first version of the screen *failed* the controls (three statistic bugs); `runs/controls_v1.0_failed/` keeps that run. Season-stacked pixel test (v1.8): KIC 5689121 and kic-006699887 on target; kic-006367225 and kic-001577049 offset, all as required. |
| Stage 1 screen | Modal, 6,365 stars: 438 real vs 95 inverted passes, 0 errors, 31 skipped (no usable cadences), wall time 2 h 17 m (`runs/poolA_main/modal_status.json`). Local pilot, 100 stars: 10 vs 2. Pool total: **448 real vs 97 inverted** (6.9 % vs 1.5 %). |
| What the null means | Inverting the flux turns dips into bumps, so the inverted rate bounds *sign-symmetric* noise only. By best-peak period, real vs inverted: 0.5–2 d 24 vs 26; 2–10 d 33 vs 36; 10–20 d 12 vs 4; **20–40 d 379 vs 31**. The excess is a dip-shaped systematic (post-gap thermal ramps at Kepler's monthly downlinks), not a count of astrophysical signals. |
| Stage 2 vetting | Masked re-detrend, per-event depths with gap proximity, V-shape proxy, variability, duration plausibility, TLS, ephemeris matching → **11 pass under v1.4, 23 under v1.5, 22 under v1.7** (v1.7 removes kic-006515744, which passed with one usable event). All three counts are reported. |
| Pixel localisation | Per-quarter χ² model comparison (v1.6, `runs/poolA_main/pixel_v16/`): the missing-light map is fitted by a point source at the target, at every other Gaia source within 2 px, and by a scaled copy of the stamp; on-target needs Σ Δχ² ≥ 25 against the null and against the best alternative. Season-stacked version (v1.8, `runs/poolA_main/pixel_season/`) for signals with few events per quarter. **5 of the 22 localise on target**: kic-003831969, kic-003963519, kic-006699887, kic-010213563, kic-011138423. The earlier centroid rules (v1.5) failed their positive control and are not used as a gate. |
| Adversarial pass | Independent from-scratch re-derivation (no shared code) plus red-team review for every on-target or inconclusive survivor: **16 reviewed, 7 survive, 1 weak, 8 killed**. On the funnel path, kic-010213563 drops (the dip is real but its light comes from an unmodelled G 18.35 Gaia source 0.97″ away), leaving 4. |
| False-positive probability | See §1.1. |
| Transit fits and stellar density | `scripts/transit_fit.py` (batman, 15× supersampled 29.4-min cadence, free Kipping limb darkening, emcee 40 × 6000). Transit-implied stellar densities match Gaia DR3 FLAME within 0.5 σ: KIC 6699887 1.74 (+0.85/−1.2) vs 2.24 g/cc (−0.3 σ), KIC 3831969 1.53 (+0.69/−1.1) vs 1.16 (+0.3 σ), KIC 11138423 4.2 ± 2.8 vs 3.0 (+0.4 σ). Chains are ~7 autocorrelation times after burn-in: fine for medians, to be lengthened for publication. |
| Human review | All stage-2, pixel, deep-run and transit-fit plots for the three candidates inspected by a human on 2026-10-17 (`runs/poolA_main/human_review.json`). |
| Evidence and archiving | Evidence packs with SHA-256 manifests in `candidates/<fname>/`; data, code and evidence on Zenodo, DOI 10.5281/zenodo.23442225; RNAAS note submitted 2026-10-18 (AAS82486, under consideration). |

### 1.1 False-positive probabilities

`pipeline/fpp_triceratops.py`: TRICERATOPS 1.1.0 in Kepler mode, 8 seeds × 10⁶ draws, the originating TPF for the WCS, a 3×3-pixel aperture on the target standing in for the Kepler Bonus PSF fit; a nearby star is dropped only when its own Kepler Bonus light curve excludes the eclipse it would need at > 5 σ. Three settings, all reported (`runs/poolA_main/fpp*/`):

| Candidate | Baseline (TIC parameters, no companion constraint) | + Gaia DR3 contrast curve | **+ MOLUSC companions (Gaia RUWE + Gaia contrast curve), Gaia FLAME host** |
|---|---|---|---|
| KIC 6699887 | 2.4 % ± 2.1 % | 1.7 % ± 1.8 % | **0.63 % ± 1.0 %** (median 0.26 %; seeds 0.15–3.4 %) |
| KIC 3831969 | 2.7 % ± 0.3 % | 2.1 % ± 0.2 % | **0.71 % ± 0.5 %** (median 0.53 %) |
| KIC 11138423 | 3.6 % ± 0.2 % | not run | **0.70 % ± 0.2 %** (median 0.66 %) |

NFPP is 0 in every run. MOLUSC (Wood et al. 2021, unmodified; compatibility shims documented in `scripts/molusc_run.py`) used a multiplicity fraction of 0.46 (the solar-type value, conservative for K/M dwarfs) and Gaia's G ≤ 20 completeness; it lowers the companion prior from 25–31 % to 15–19 %. The remaining FPP is a planet on an unresolved bound companion or a diluted planet on the target; eclipsing-binary scenarios are < 10⁻⁴. The quoted FPPs include no multiplicity boost.

Other signals: the 47.69-d second signal on KIC 6699887 has its own FPP 2.4 % ± 0.5 % (baseline setting, `runs/poolA_main/fpp_posthoc/`), NFPP 0. KIC 3963519: FPP 39 % ± 3 %, driven by an unmodelled G ≈ 18 neighbour, so it stays a signal. kic-009655224: FPP 5.8 % ± 0.2 % with an inconclusive pixel test.

## 2. Disposition of every stage-2 survivor

The 23 v1.5 survivors plus the two v1.4 survivors that drop under v1.5. Full numbers per row are in `runs/poolA_main/vet/<fname>/vet.json`, `runs/poolA_main/pixel_v16/summary.md` and `runs/poolA_main/adversarial/<fname>/`.

| fname | P (d) | Pixel (v1.6) | Outcome |
|---|---|---|---|
| kic-006699887 | 15.1136 | on-target | **Candidate** (plus a second on-target 47.69-d signal) |
| kic-003831969 | 15.3134 | on-target | **Candidate** |
| kic-011138423 | 4.8945 | on-target | **Candidate** |
| kic-003963519 | 6.8156 | on-target | Survives red team; FPP ≈ 39 % → signal |
| kic-010213563 | 36.7671 | on-target | Red team weak: deficit sits on an unmodelled G 18.35 source 0.97″ away; host ambiguous |
| kic-009655224 | 6.7035 | inconclusive | Survives red team; FPP 5.8 %, pixel inconclusive → signal |
| kic-006448637 | 35.5873 | inconclusive | Survives red team; localised only by red-team code (saturated stamp host) → signal |
| kic-011285493 | 10.7910 | inconclusive | Killed: eclipse of KOI-7431.01 / KIC 11285625, 5.6′ away on the same CCD channel |
| kic-006367225 | 12.9715 | offset | Killed: background eclipsing binary on a G 20.5 neighbour 3.6″ away |
| kic-005950408 | 38.6067 | inconclusive | Killed: bleed columns of a saturated Kp 8.45 giant |
| kic-005609180 | 21.1063 | inconclusive | Killed: post-gap PSF-model ramps |
| kic-004474607 | 32.9427 | inconclusive | Killed: Q13 post-downlink thermal ramps |
| kic-005866176 | 29.4330 | inconclusive | Killed: post-gap artefact; target outside the stamp |
| gaia-dr3-2102601072256824704 | 32.5571 | inconclusive | Killed: segment-edge ramps; also a target-list leak (0.32″ from KIC 7269549) |
| kic-008620445 | 28.0717 | inconclusive | Rejected (provisional, spot-check): Q2 attitude-tweak step and a post-gap ramp |
| gaia-dr3-2127687048516119168 | 32.8687 | inconclusive | Rejected (provisional, spot-check): Q13 segment starts; also a target-list leak (0.79″ from KIC 9459487) |
| kic-006515744 | 36.2075 | inconclusive | Killed: single-cadence PSF artefacts; fails v1.7 `n_events_min` (not among the 22) |
| kic-008955698 | 28.3811 | offset | Not red-teamed |
| kic-008636232 | 36.3024 | offset | Not red-teamed |
| kic-010782709 | 31.1781 | offset | Not red-teamed |
| kic-011236638 | 26.3823 | offset | Not red-teamed |
| kic-006184850 | 34.7932 | offset | Not red-teamed |
| kic-009019826 | 36.3047 | offset | Not red-teamed |
| kic-009520080 (v1.4 only) | 2.3770 | on-target | Survives red team but fails the v1.5 two-sided secondary test (unexplained brightening at phase 0.5) → signal |
| kic-001577049 (v1.4 only) | 0.5379 | offset | Fails the v1.5 empirical-SNR floor |

Fourteen of the 23 lie at 21–39 d, inside the systematic band of §1; every re-derivation that disagreed with the pipeline pointed at post-gap ramps or single-cadence outliers there.

## 3. The three candidates and what each still needs

| Candidate | Host (Gaia DR3) | P (d) | R_p (R⊕) | FPP (MOLUSC) | Open issues |
|---|---|---|---|---|---|
| KIC 6699887 | K dwarf, G 16.2, 859 pc | 15.1136 | 2.76 (+0.17/−0.14) | 0.63 % | Large seed-to-seed FPP scatter; second 47.69-d signal on target (season-stacked pixel test, post-hoc) |
| KIC 3831969 | G dwarf, G 16.5, 1,969 pc | 15.3134 | 3.47 (+0.38/−0.32) | 0.71 % | ~2 σ season dependence of depth |
| KIC 11138423 | early M, G 17.4, 894 pc | 4.8945 | 2.83 (+0.26/−0.24) | 0.70 % | Faintest host; SAP/PSF depth ratio explained by the red team, not by a model |

What all three still need before anyone may call them validated:

1. **High-resolution imaging** (speckle or adaptive optics) of each G 16–17 host, to replace the statistical bound-companion constraint inside ~0.5″. This needs a professional collaborator with access to a large telescope.
2. **Radial velocities are impractical** at this brightness; no mass measurement is planned.
3. **A refereed methods paper** describing the pipeline, nulls, controls and pixel tests. ExoFOP lists community candidates only after refereed publication, and accepts an RNAAS as that publication only if the methods were refereed first (`docs/06`).

## 4. Lessons the red team forced

All are post-hoc observations; none is applied retroactively to Pool A results. They belong in the next pre-registration.

1. **Ephemeris matching must cover the whole CCD channel, per season**, not a 3′ cone: a 12.5 %-deep eclipsing binary 5.6′ away was imprinted channel-wide on Kepler Bonus PSF photometry (kic-011285493).
2. **Cross-match "never a Kepler target" by position, not KIC id**: two Gaia sources in the pool were Kepler targets or their TPF hosts (0.32″ and 0.79″).
3. **List unmodelled Gaia neighbours within 5″ for every survivor**: Gaia sources absent from the Kepler Bonus catalogue are absorbed into the target's PSF flux and produced two false hosts (kic-010213563, kic-006367225).
4. **Add a segment-edge veto at stage 2** (events within 0.3–0.5 d of a gap): the 21–40 d band is dominated by post-gap ramps.
5. **Require a margin between real and inverted SDE** (real ≥ k × inverted): several rejected signals had inverted SDE level with or above the real one and passed only because the inverted winners failed the duty-cycle cut.
6. **Build a pool-wide bad-cadence list**: the same spacecraft cadences (Q2 attitude tweak, Q13 segment starts) recur across unrelated stars and channels; PSF flux dipping while SAP flux rises at the same cadence is a clean artefact flag.
7. **Run every pixel rule on positive and offset controls first**: the v1.5 centroid rules failed the positive control because background sources sit at stamp edges by construction.

## 5. Next milestones

| # | Milestone | Status / depends on |
|---|---|---|
| 1 | RNAAS decision on the three-candidate note | Submitted 2026-10-18 (AAS82486), under consideration |
| 2 | Speckle/AO imaging of the three hosts with a professional collaborator (Kepler Bonus authors are the natural first contact; `docs/06`) | Collaborator |
| 3 | Refereed methods paper (pipeline, nulls, controls, pixel tests, lessons of §4) | 2 helps; prerequisite for an ExoFOP community-candidate listing |
| 4 | ExoFOP community-candidate listing | 3 |
| 5 | Extend to further pools: Pool B (faint TESS M dwarfs, T 13.5–16, < 100 pc; `docs/04b` §9) and remaining Kepler Bonus sources, with a new pre-registration that registers §4 before any search | Pre-registration, controls, pilot |
| 6 | Shut down the Modal deployment and rotate its token when no further batches are planned | Owner |
