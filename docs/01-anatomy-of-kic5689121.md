# 01 · Anatomy of the KIC 5689121 claim

> Source of truth for this chapter: Kuber Mehta's own evidence repository
> (https://github.com/Kuberwastaken/planethunt, Zenodo doi:10.5281/zenodo.23325261),
> his site (https://kuber.studio/planethunt/), and direct queries we ran against the
> NASA Exoplanet Archive and MAST on 2026-10-17. Where the tweet and the repo disagree,
> the repo is quoted.

## 1. The claim vs. the record

| | The tweet (2026-10-11) | The repository / RNAAS draft (2026-10-13) |
|---|---|---|
| Wording | "Found a new planet with Codex this weekend" | "**No confirmed planet is claimed.** ... a well-vetted **planet candidate**" |
| Host | "a red dwarf in Lyra ~300 pc away hiding in old Kepler data" | KIC 5689121 = TIC 120686465 = Gaia DR3 2103792461826460288; Teff ≈ 3633 K, R ≈ 0.51 R☉, Kp 16.33, d ≈ 298 pc |
| Signal | "a dip that repeats every 13.77 days" | P = 13.769448 ± 0.000065 d, T0 = BJD_TDB 2455630.2454, depth 1.4–2.4 ppt (method-dependent), duration 2.1–2.4 h, Rp ≈ 2.5–3 R⊕ |
| Novelty | "went undiscovered" | "Novelty is provisional." The star appears in Barber et al. 2022 (AJ 164, 88) Table 5 with a *Kepler = Y* flag (its Kepler light curve was good enough for their Notch transit search). Whether they looked at this signal is unresolved. |
| Status | planet | TRICERATOPS FPP = 0.069 (validation requires < 0.015). RNAAS manuscript AAS82222 **declined** on the day of submission. |

**Bottom line:** a real, repeating, on-target-looking transit signal in a dataset almost nobody searches, found in a weekend with an AI agent, but *not* a validated planet, and probably not validatable with any telescope an amateur can reach.

## 2. Why the star was "hiding": it was never a Kepler target

We checked on 2026-10-17:

- NASA Exoplanet Archive `cumulative` (KOI) table: **no row** for `kepid = 5689121`.
- MAST: **no** `kplr005689121` light-curve or target-pixel products exist. The star was never on the Kepler target list.
- MAST HLSP **KBONUS-BKG** (Kepler Bonus, Martínez-Palomera, Hedges & Dotson 2023, AJ 166, 265, arXiv:2310.17733): products exist for `Gaia DR3 2103792461826460288`. These are PSF-photometry light curves extracted from the *target pixel file of the neighbour* KIC 5689116 (7.1″ west, ≈2.1× brighter in Gaia G), which *was* a Kepler target.
- TESS: 10 FFI sectors (14, 26, 40, 41, 53, 54, 74, 75, 80, 81) exist for TIC 120686465, but the star is far too faint (T ≈ 15) for a 2 ppt dip to be recovered in a 21″-pixel blend. Kuber's fixed-ephemeris folds give −0.5 to −2.1 ppt (apertures) and +0.45 ± 1.8 ppt (TGLC): no recovery, no rejection.

So the "undiscovered" part is explained by *coverage*, not by a mistake of the Kepler team: the official pipeline only searched the ~200,000 targeted stars. Kepler Bonus added ~600,000 background sources in 2023, and nobody has systematically searched them for planets.

## 3. The Kepler-169 e coincidence (and why it is not contamination)

Kepler-169 (KIC 5689351, KOI-505) lies 4.39′ away and has a planet, Kepler-169 e, with **P = 13.76710058 d** (NASA Exoplanet Archive). The candidate's period is 13.769448 d. That is suspicious: Kepler's DR25 Robovetter has an explicit *ephemeris-match* test (Coughlin et al. 2014) because electronic cross-talk, column anomalies and PRF wings put copies of one star's eclipses into other stars' light curves.

Kuber's check (we agree with the arithmetic): the two ephemerides are **~2.9 d out of phase** and the period difference (0.0023 d) only drifts that offset from 2.80 to 3.04 d across the four-year baseline, so *no candidate transit ever overlaps a Kepler-169 e transit* (nearest pair 67 h apart). Folding the candidate's light curve on all five Kepler-169 ephemerides gives null depths. A period coincidence with no epoch coincidence is not contamination.

**Lesson:** ephemeris matching against every known KOI/TCE/EB within a few arcminutes *and* on the same CCD module/output (including harmonics) must be an automatic step in your pipeline.

## 4. The search, step by step (from `analysis/editorial-audit-20261013/README.md`)

### 4.1 Sample construction

| Cumulative cut | Sources |
|---|---:|
| Kepler Bonus v1.0 **M-star catalogue** (catalogue-selected potential M dwarfs) | 29,800 |
| catalogue minimum `PSFFRAC` > 0.5 (the source's share of the PSF flux) | 9,204 |
| catalogue minimum `CROWDSAP` > 0.7 (aperture crowding metric) | 6,477 |
| ≥ 8 available PSF quarters | 2,495 |
| background only: `kic != tpforg` (source is not the star the pixels were downloaded for) | 1,836 |
| exclude Kepler primary targets | **1,800** |

Sorted by ascending Gaia G. The Codex "pilot" searched ranks 1–64 on 10 Oct; the extension searched ranks 65–1800 on 11 Oct (1,672 actually searched; 63 skipped for < 6 usable quarters, 1 download timeout).

### 4.2 Search and first screen

- Reserve the last `max(2, ceil(0.25·N))` usable quarters as a **hold-out** before searching.
- Clean discovery and reserved subsets independently.
- **BLS** over **0.5–25 d**, up to two peaks per star (code imported from `paraschopra/tess-planet-candidates` and `eumemic/exoplanets`, both MIT).
- First-screen pass requires: BLS SNR ≥ 9; ≥ 3 good events in each subset; held-out SNR ≥ 4; odd/even difference < 3σ; secondary significance < 3; max single-event contribution < 0.5. Stellar density assumed 5 ρ☉ for candidate generation.

### 4.3 Null controls (the most important table in the whole repo)

Every star was also run **inverted** (`f → 2 − f`, so dips become bumps) through the identical pipeline:

| Stage | Stars | Real passes (signals / stars) | Inverted passes |
|---|---:|---:|---:|
| Initial 64 | 64 | 2 / 2 | 0 / 0 |
| Ranks 65–576 | 497 | 10 / 9 | 10 / 9 |
| Ranks 577–1800 | 1,175 | 19 / 18 | 22 / 19 |
| **Total** | **1,736** | **31 / 29** | **32 / 28** |

The first screen passes *noise as often as data*. Even at the candidate's own strength (P > 1.5 d, BLS SNR ≥ 26.7, held-out SNR ≥ 15.5), **two inverted stars survive** (KIC 8701308 at 4.25 d, KIC 5283912 at 8.78 d). Kuber's own conclusion, corrected on 13 Oct: "The case for KIC 5689121 rests on the vetting (shape, density, localization, dilution pattern), not on detection strength."

**Lesson:** a BLS SNR of 27 in PSF photometry of a Kp 16 star is *not* a 27σ detection. Inverted and time-scrambled controls are mandatory, and they must be run on the whole sample, with thresholds fixed *before* you look.

### 4.4 Vetting stack that was applied to the survivor

1. **Which star dims?** The neighbour KIC 5689116's own KBONUS PSF light curve at the candidate ephemeris: 24 ± 41 ppm (it would need ≈ 1100 ppm to host the signal). Target: 2348 ± 128 ppm.
2. **Dilution pattern.** Per quarter, the aperture (SAP) depth tracks the target's share of aperture flux (χ² = 15.4/12, p = 0.22); constant depth is rejected (p = 0.005); a contaminant-hosted signal is rejected (p ≈ 1e−36).
3. **Raw-pixel centroids** from the 13 calibrated TPFs of KIC 5689116: in-transit shift detected; target favoured over the neighbour by Δχ² = 934 and over every other Gaia source within 30″ by Δχ² ≥ 4300. (The KBONUS `CENTROID_COLUMN/ROW` columns show *no* response; they are model-derived and useless for this test.)
4. **Odd/even** depths agree; **no secondary eclipse**; transit-derived stellar density 4.4 (+2.3/−2.8) ρ☉ is consistent with an M dwarf.
5. **Catalogue novelty:** no KOI, DR24/DR25 TCE, TOI or CTOI for the source or its neighbour; no ExoFOP follow-up uploads.
6. **TRICERATOPS 1.1** (driven offline without its MAST initialiser, Kepler exposure settings, Gaia-built background population): FPP = 0.069 (runs 0.058–0.089). Scenario shares: planet on target 0.75; planet on target with bound companion 0.17; **planet on an unresolved bound companion 0.068**; everything else ≲ 0.01. A *hypothetical* Keck AO contrast curve only brings FPP to ≈ 0.027, because companions inside ≈ 0.1″ (≈ 30 au at 298 pc) remain allowed. NFPP was reported as 0 only because resolved neighbours were not passed in.
7. **MCMC transit fit** (PyTransit/batman + emcee, exposure-time integrated for 29.4-min cadence).

### 4.5 Where it could not go further

- **No high-resolution imaging** exists (archival CFHT image saturates both stars). Keck/NIRC2 AO is the obvious ask; the 2027A NASA-Keck call closed 2026-09-10 and requires a US-institution PI.
- **TESS cannot test it** (above). Spitzer, SPHEREx, ZTF, PTF, NEOWISE all miss the transit windows or are too shallow.
- **No RV mass is realistic** for a Kp 16.3 M dwarf.
- Realistic end state, in Kuber's words: *statistically validated, if imaging and Gaia-RUWE companion limits squeeze the bound-companion scenario below ~1%*. Otherwise it stays a candidate forever.

## 5. The publication attempt

- RNAAS manuscript AAS82222, "A 13.77-day transit candidate on the background M dwarf KIC 5689121 from Kepler Bonus photometry", sole author (independent researcher), AI assistance (Codex, ChatGPT, Claude) disclosed. Submitted 2026-10-13 IST, **declined the same day**, before a scientific editor was assigned.
- His post-decision audit addresses four topics: (1) how the sample was really constructed; (2) the *other* 30 first-screen signals and the inverted controls; (3) the TESS non-detection; (4) the Barber et al. 2022 catalogue match (prior search of the same light curve).
- Zenodo DOI and a public GitHub evidence pack with SHA-256 manifests, `claim-sources.json` tracing every number in the note to a file, `inventory.json`, pinned requirements, and an LLM entry point (`planethunt.md`).

## 6. Timeline

| Date (2026) | Event |
|---|---|
| 9 Oct | TESS search of 175 nearby K/M dwarfs begins (Codex). Two leads; one is probably an eclipsing binary (FPP ≈ 0.57). |
| 10 Oct | Kepler Bonus background-M-dwarf pilot (64 stars) → KIC 5689121 candidate. Evidence archives released. |
| 11 Oct | Tweet ("Found a new planet"). Extension to 1,736 stars with inverted controls. Validation pack: centroids, dilution, hold-out, TRICERATOPS. |
| 12 Oct | Archival follow-up: Spitzer, 10 TESS sectors, TGLC, ZTF, PTF, NEOWISE, CFHT, SPHEREx. Nothing recovers or rejects the transit. |
| 13 Oct | RNAAS submission; declined. Post-decision audit. |
| 14 Oct | (Context) Scientific American profiles Pavel Rabtsevich's AI-assisted TESS candidate TIC 4206066, which the TESS team expects to become an official TOI. |

## 7. What to copy, what to avoid

**Copy**
- Searching a dataset the official pipelines never searched (Kepler Bonus background stars).
- Freezing the target list *before* looking at photometry; ranking by brightness.
- Reserving hold-out quarters before the search; matched inverted controls on every star; a positive control (Kepler-138 c/d recovered).
- Pixel-level "which star dims" tests using the original TPFs, dilution-vs-quarter consistency, ephemeris matching against Kepler-169, TRICERATOPS, and a public, hash-manifested evidence pack with every number traceable to a file.
- Candidate language throughout the repository and the manuscript.

**Avoid**
- Announcing "found a new planet" publicly for a candidate with FPP ≈ 6.9 %. The repository's own documents, written over the following days, use candidate language throughout.
- Choosing a regime (Kp 16, 300 pc) where *no* amateur-accessible follow-up can ever validate the result.
- Post-hoc screens (the 12.5–14.5 d exclusion and the SNR ≥ 26.7 cut were chosen after seeing the candidate).
- Submitting a journal note for an unvalidated candidate three days after the search, before independent review by someone outside the project.
- Reporting a detection SNR that the null controls show is reachable by noise.
