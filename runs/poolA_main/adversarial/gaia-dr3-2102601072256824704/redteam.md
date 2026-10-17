# Red-team review: gaia-dr3-2102601072256824704 (run poolA_main, 2026-10-17)

**Verdict: KILL.** The 32.556-d first-screen signal is an assembly of three detrending
artefacts at Kepler segment edges (safe-mode entry, post-roll quarter start, post-anomaly
gap exit). With those edges handled, nothing remains above the inverted-flux null. On top
of that, the target is not a background dwarf: it is the Kepler target KIC 7269549 itself,
an unresolved 0.74-arcsec pair of G dwarfs, already searched by the Kepler pipeline.

All numbers below come from `runs/poolA_main/adversarial/gaia-dr3-2102601072256824704/redteam/redteam_check.json`
and `redteam_check2.json` (scripts `redteam_check.py`, `redteam_check2.py` in the same
directory; wotan biweight detrend, window 1.0 d, break_tolerance 0.5 d; astropy
BoxLeastSquares; nothing imported from `pipeline/`), from the pipeline's
`runs/poolA_main/vet/gaia-dr3-2102601072256824704/{vet.json,events.csv}`, from the
independent re-derivation `.../rederive/results.json`, and from the FITS headers of
`data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2102601072256824704_kepler_v1.0_lc.fits`.
Catalogue dates: `data/known/manifest.json` fetched 2026-10-17T15:46 UTC; Gaia DR3 cone
search via astroquery on 2026-10-17.

## 1. What the "signal" is made of (independent third detrend)

Pipeline ephemeris: P = 32.55706 d, t0 = 148.8235 BKJD, duration 3.375 h, stage-2 depth
191 +- 32 ppm (vet.json, rules v1.7, stage2_pass = true; failed v1.4 on tls_sde, TLS SDE 4.0).

My biweight detrend, PSF flux, no edge trimming (`PSF_edge0`): depth 147 +- 22 ppm
(empirical SNR 6.8), odd 189 / even 109 ppm, 43 events. Per-event local depths: 40 events
scatter within +-300 ppm; three stand at **992, 1044 and 969 ppm** (epochs n = 1, 12, 39),
each 0.04-0.14 d from a data-segment edge. Rule applied before looking at the numbers:
|local depth| > 5x median |local depth| flags exactly these three and nothing else. Excluding
them, the depth is **66 ppm**.

Same detrend with wotan `edge_cutoff = 0.3 d` (`PSF_edge0.3`): depth **62 +- 20 ppm (SNR
3.1)**, odd 81 / even 45 ppm, no outlier events.

SAP flux: 82 +- 22 ppm untrimmed, 50 +- 21 ppm trimmed; odd/even 148/20 -> 53/48.

Per-channel depth (season), PSF, untrimmed -> trimmed:
ch22 (Q1,5,9,13,17) 65 -> 62; **ch38 (Q2,6,10,14) 301 -> 91**; ch62 (Q3,7,11,15) 158 -> 78;
ch46 (Q4,8,12,16) 26 -> 10 ppm. The untrimmed signal lives on the channels that host the
three edge events (Q2, Q6 on ch38; Q15 on ch62); after trimming all four seasons sit at the
~50 ppm noise floor with ch46 consistent with zero.

The three events (BKJD, UTC, context from the time array; raw PSF slope from a linear fit
+-0.6 d):
- n = 1, tc 181.42, 2009-07-01, Q2: **0.095 d before the 2.21-d gap starting 181.515
  (2009-07-02)**, i.e. the last ~2 h before the first Kepler safe mode; raw flux ramping
  **-5680 ppm/d**, local range 4548 ppm.
- n = 12, tc 539.53, 2010-06-25, Q6: **0.059 d after the first cadence of Q6** (1.84-d
  Q5->Q6 roll/downlink gap); raw flux ramping **+5552 ppm/d**, range 3795 ppm.
- n = 39, tc 1418.53, 2012-11-20, Q15: **0.069 d after a 5.66-d gap** ending 1418.461
  (2012-11-19; far longer than a routine ~1-d downlink); raw flux ramping **+3730 ppm/d**,
  range 3832 ppm.
Over a 3.4-h window those ramps move the flux by 500-800 ppm; a 1-d one-sided smoother
cannot follow them, so the residual is a fake 1000-2000 ppm "dip". The pipeline's own
events.csv records the same three events (1114, 1016, 988 ppm at gap_dist 0.13, 0.038,
0.088 d) and event_chi2_const = 108.5 for 42 dof, which should already have failed the
star; the stage-2 gap-proximity rule evidently does not treat a segment start/end as a gap
when the gap itself lies outside the window.

Odd/even: the three artefacts are n = 1, 39 (odd) and n = 12 (even); they alone explain
the 2.6-sigma odd/even asymmetry (189/109 -> 81/45 after trimming).

## 2. Periodogram with the edges handled

Narrow BLS 28-37 d, likelihood objective, SDE in units of the periodogram scatter:
- untrimmed, real: top 32.5554 d, SDE 9.2; at the pipeline period SDE 8.6, depth 125 ppm.
  Inverted: SDE -0.2 at that period.
- **edge-trimmed (0.3 d), real: SDE 2.35 at the pipeline period, depth 72 +- 14 ppm (formal)**;
  the top real peak elsewhere is 28.15 d at SDE 7.2 against an inverted-null top of 34.55 d
  at SDE 7.0, i.e. noise.
This agrees with the re-derivation's full 1-50 d search (edge-trimmed SDE 1.4 at 32.556 d,
rank 20104; inverted-null max 10.74 > real max 10.60 even on the untrimmed data). By the
repo's own rule the real detection was never distinguishable from the inverted null.

## 3. Stellar variability / rotation aliases
Lomb-Scargle 0.2-60 d on the detrended flux: top peak 1.4724 d, amplitude ~10 ppm (same for
both polarities); power near 32.56, 16.28 and 10.85 d is < 1e-4 of the scatter. No rotation
or pulsation alias can produce the signal. (The 1.47-d LS peak is close to the inverted
BLS top at 1.4956 d, both at the ~10 ppm level; irrelevant to the 32.6-d question.)

## 4. The target is not a background dwarf; it is a 0.74" binary Kepler target
- `keplerstellar`: KIC 7269549 at RA 287.59302, Dec 42.877159, **0.32 arcsec** from the
  Gaia source; Kp 14.761, Teff 5802 K, R 0.99 Rsun, observed Q1-Q17. It is the TPF host
  (TPFORG 7269549). Not a KOI, not a DR25 TCE, not in the Kepler EB catalogue (Kirk+2016,
  J/AJ/151/68); no KOI/TCE/confirmed planet within 3 arcmin.
- Gaia DR3 cone search (12"): **2102601072256824704 (G 15.276, BP-RP 0.923, plx 0.697 +-
  0.074 mas, RUWE 3.87, ipd_frac_multi_peak 70) and 2102601072259527296 at 0.74" (G 15.762,
  BP-RP 0.921, plx 0.718 +- 0.053, pmra/pmdec 0.14/3.67 vs 0.28/3.56 mas/yr)**. Same
  distance, same colour, common proper motion: a physical pair of near-twin G dwarfs. Their
  combined G = 14.75 matches Kp 14.76. Primary flux fraction 0.61.
- KBONUS serves neither `kic-007269549` nor the companion (both URLs 404 on 2026-10-17),
  so the single PSF light curve we searched is the unresolved pair fitted as one source.
  Any signal would be diluted by 1/0.61 (primary) or 1/0.39 (companion) and could not be
  assigned to either star from this photometry; and the Kepler pipeline already searched
  this star for 17 quarters with no TCE. The pool's premise ("never searched") fails here;
  the ranked CSV carries kic = 0 and an empty kic_sep for this row, which is how it got in.
- The pipeline's rstar_used = 0.87 Rsun is derived from the primary's G alone; fine for
  duration arithmetic but the system is a binary.

## 5. Ephemeris-match contamination
Channels: 22 (season 3), 38 (season 0), 62 (season 1), 46 (season 2). KOIs/DR25 TCEs with
P within 0.2% of 32.556 d or of its 2x, 3x, 1/2, 1/3 harmonics: 32 KOIs and ~70 TCEs
field-wide; the closest on sky are KIC 7430876 (32.612 d, 0.42 deg, epoch offset -4.3 d),
KIC 7515212 (16.258 d, 0.27 deg, -9.2 d) and KIC 6850504 (10.854 d, 0.54 deg, -10.3 d).
None matches in epoch (a match needs |offset| < ~0.15 d); KIC 7887124 (32.496 d, offset
0.14 d) is 2.25 deg away and its 0.18% period difference drifts 2.5 d across the baseline.
No ephemeris-match contamination. Moot anyway given section 1.

## 6. Duration / density / implied radius (if it were real)
For R* 0.87-0.99 Rsun, M* ~1 Msun, P = 32.557 d: a/R* 43-49, central T14 5.1-5.7 h;
observed 3.4-3.6 h implies b ~ 0.75, allowed. Observed 147 ppm would mean Rp 1.5 Re on the
primary or 1.7 Re on the companion; the edge-trimmed 62 ppm would mean ~1.0-1.1 Re. Nothing
unphysical, so this test does not kill it; it is simply not needed.

## 7. Checks that could not be run
- No pixel.json exists and `pipeline.pixel_vet` was not run. It should not be: there is no
  coherent signal to localise, and at 0.74" the pair is unresolvable in 4" pixels anyway.
- Bienias & Szabo 2025 and the Kepler Bonus team's own candidate lists were not queried
  live; with the signal gone and the target being KIC 7269549 (no KOI/TCE), the novelty
  question is moot. Record as unchecked.

## 8. Recommendations (post-hoc; label as such in any pre-registration update)
1. Add segment-edge distance (distance to the first/last cadence of each contiguous segment,
   not only to gaps inside the window) to the per-event table, and require survival with
   events < 0.25 d from an edge removed, or trim 0.25-0.3 d at every edge before detrending.
2. Enforce the event_chi2_const rule: 108.5/42 dof is a 6-sigma rejection of a constant
   depth that v1.5-v1.7 let through.
3. Pool hygiene: flag any target whose TPF host KIC lies < 2" away (this row: 0.32") and any
   with a Gaia neighbour < 2" at comparable brightness, or RUWE > 2 / ipd_frac_multi_peak >
   30. Such sources are the Kepler target itself, already searched, and KBONUS PSF flux is a
   blend. Re-check the other batch-A survivors for all three failure modes.
