# 11. TESS feasibility: faint, nearby northern M dwarfs

Status: exploration, 2026-10-19. No search has been run and nothing here is a signal. Every number comes from
`scripts/tess_feasibility.py` and the files it wrote in `runs/tess_feasibility/`.

## Why look

Docs/02 §5 item 2 ranked TESS M dwarfs at T ≈ 13.5–16 within 100 pc as the second-best open ground. The motive now is
follow-up. The three Kepler Bonus candidates sit 860–1,970 pc away. A TESS candidate around an M dwarf within 100 pc has
a transit several ppt deep, which a 1-m telescope can see from the ground.

That doesn't shortcut the paper. TOIs come from the TESS team's own pipelines, and a CTOI needs a refereed paper.

## Pool (`pool_summary.json`)

- **Selection:** Gaia DR3 stars with parallax > 10 mas (at > 10 σ), 1.8 < BP−RP < 5, M_G > 8 and northern ecliptic
  latitude. That gives 46,254 stars at G 13.8–17.5.
- **Brightness cut:** with T from the Stassun et al. (2019) colour relation, 33,689 have T 13.5–16, of which 5,348 are
  within 60 pc.
- **Sub-groups:** 5,378 sit at |b| < 10°, 4,976 at T > 15.5, and 29,812 have RUWE < 1.4.
- **Known objects:** 27 TOI rows fall within 10″ of a pool star. About 13 of those are genuine pool members (their
  catalogue Tmag agrees within 0.6 mag); the rest are brighter neighbours.

## Light curves (`coverage_summary.json`, 120 random pool stars)

| Product | Stars with any data | Median sectors |
|---|---|---|
| TARS (Boyle, Bouma & Mann 2026; CPM; S1–96, T < 16) | 95 % | 3 |
| TESS-SPOC FFI | 78 % | 2 |
| QLP | 73 % | 2 |
| SPOC 2-min | 18 % | 0 |
| TGLC (Han & Brandt 2023; PSF; S1–55) | 50 % | 0.5 |

TGLC files exist for 220 of the 314 sector-pointings that tess-point predicts in S1–55.

TARS stores time and flux only, with 30-min to 200-s cadence depending on sector. It is about 0.1 MB per star and
sector, so the whole pool is a few tens of GB.

## Who has already searched these stars

- **NEMESIS I** (Feliz et al. 2021, AJ 161, 247): 100-pc M dwarfs, sectors 1–5 (south).
- **NEMESIS II** (poster, Zenodo 12973447): about 191,000 FFI light curves of 100-pc M dwarfs in S1–26. That includes the
  northern S14–26, but there is no published candidate list.
- **QLP faint-star search:** it is producing TOIs in this exact range. Examples are TOI-7281 (T 14.65), TOI-7415 and
  TOI-7437, all from `qlp-s8x-faintsearch`. SPOC FFI multi-sector runs (`spoc-s14-s60`) add more.
- **Eschen & Kunimoto (2024, MNRAS):** nine M-dwarf candidates from S1–63 FFIs.
- **zacharynapier/tess-mdwarf-search** (GitHub, read 2026-10-19): 15,512 M dwarfs at T 13–15.5 and d < 100 pc, all sky.
  It requires |b| > 10°, ≥ 4 TESS-SPOC sectors and no 2-min data. Its results are "Pending" and it says "a Research Note
  is under review".
  - **Overlap with our sample:** only 19 of 120 stars meet a reproduction of that selection, which over-counts the overlap
    because the Teff, log g and contamination cuts are not reproduced.
  - **Complement:** 65 of 120 (54 %) fall outside it, have no 2-min data and have ≥ 2 TARS sectors. Those stars are
    fainter (mean T 14.93) and have fewer sectors (median 3).

## Pipeline controls (`controls.json`)

The Kepler pipeline's cleaning failed on TESS. With one global noise estimate, 200-s sectors swamped the BLS. The fix,
now in `clean()`:

- biweight detrending with a 0.5 d window
- trim 0.3 d after and 0.1 d before each gap > 0.5 d
- 10-min bins
- per-segment robust noise, so each sector is weighted by its own scatter
- upward-only 4σ clip

With these settings, on the 13 genuine pool TOIs:

- **In-window signals (P 0.5–15 d):** recovered as the top BLS peak on every star TARS covers, 9 of 9. TGLC recovers 5 of
  5 where it has data.
- **Out-of-window periods:** found at aliases. The 32.07 d and 34.93 d signals came back at P/3, and the 0.164 d signal at
  4P.
- **Second planets:** not recovered (TOI-2096.01, TOI-6022.02). The search would need iterative masking.
- **Depth:** TARS keeps 30–96 % of the catalogue depth on the signals not flagged as off-target, median 70 %. The CPM eats
  into transits; the low end is TOI-1696.01 after 10-min binning.
- **Off-target false positives:** two FPs nearly vanish in TARS, which is useful because they are likely on neighbours:
  - TOI-1632.01: 7,618 vs 117,724 ppm
  - TOI-1474.01: 539 vs 15,868 ppm
- **No data:** TOI-6838 (12 pc) has neither TARS nor TGLC data.

## Sensitivity and yield (`inject.json`, `yield.json`)

Boxes of 1–4 R⊕ were injected into 40 random pool stars' TARS light curves. The stars have a median of 4 sectors and a
median 10-min scatter of 7,394 ppm. Periods were drawn log-uniform in 0.5–10 d, with b < 0.8.

**Noise floor.** On the same stars without injection, the top-peak BLS SDE has percentiles 50/90/95 = 6.0 / 9.7 / 10.1.
So 15 % of stars exceed SDE 9 and 7.5 % exceed SDE 10 on noise alone. A first screen on SDE alone would pass thousands of
noise signals, so the inverted null and per-transit tests must do the real work, as in pool A.

**Recovery (period match and SDE ≥ 10):**

| R⊕ | 1 | 1.5 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| all | 0.03 | 0.08 | 0.15 | 0.48 | 0.75 |
| T 13.5–14.5 | 0.05 | 0.16 | 0.21 | 0.47 | 0.68 |
| T 14.5–15.5 | 0 | 0 | 0.18 | 0.82 | 1.0 |
| T 15.5–16 | 0 | 0 | 0 | 0.10 | 0.60 |

**Expected yield.** With rough M-dwarf occurrence rates (`OCC`, uncertain by ~2×) and geometric transit probabilities, the
estimate is 7.4 detections per 1,000 stars at SDE ≥ 10:

- about 248 for the whole pool
- about 39 within 60 pc
- about 134 in the complement of the closest public search

These are upper bounds and are not novelty-filtered:

- the injection ignores the CPM's depth loss (median 70 % retained on the controls)
- the 40-star recovery fractions carry ±0.1–0.15
- known TOIs and other teams' finds are not removed

A realistic expectation is a few new, small candidates per 1,000 stars, mostly 2–4 R⊕ around the smallest and nearest stars.

## Verdict

**Feasible, with a narrower claim than docs/02 made.** The pool is not untouched. Professional and amateur searches cover
part of it, and one public project overlaps closely.

**The open part:**
- northern stars outside the closest public search's selection (< 4 TESS-SPOC sectors, |b| < 10°, or T > 15.5)
- all TARS sectors through S96 stacked together
- an inverted-flux null and pixel tests carried through every stage

**Before any pool star is searched:**
1. A pre-registration (`pipeline/prereg_tess_v1.json`) fixing the period range, durations, SDE/SNR cuts and the
   inverted-null rule.
2. Iterative masking for multi-planet systems.
3. TESS pixel localisation from TESScut difference images (21″ pixels make neighbours the main false positive).
4. A novelty cross-match against live ExoFOP TOI/CTOI, the NEA, NEMESIS, Eschen & Kunimoto and the public searches above.
5. The 13 pool TOIs as positive controls, which must pass first.

**Suggested first run:** the 5,348 stars within 60 pc. They are the best for ground follow-up, and the run is small:
about 2–3 GB of TARS data and about 100 CPU-hours for real plus inverted flux at about 34 s per search.
