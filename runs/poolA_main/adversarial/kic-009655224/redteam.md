# Red-team review: kic-009655224 (Gaia DR3 2128120393543836928)

Run id: `poolA_main/adversarial/kic-009655224/redteam`, 2026-10-17. Candidate language only: this is a
*signal* in a KBONUS-BKG background-dwarf PSF light curve. Nothing here is a planet, and this review does
not promote it past "signal"; it records the attempts made to kill it and their outcomes.

Inputs read: `runs/poolA_main/vet/kic-009655224/{vet.json,fold_events.png,events.csv}`,
`runs/poolA_main/pixel/kic-009655224/{pixel.json,diffimages.png}`, the independent re-derivation
(`rederive/rederive_result.json` and its supplements), `targets/poolA_kbonus_bg_dwarfs_ranked.csv` row,
`data/known/` (fetched 2026-10-17T15:46 UTC, see `data/known/manifest.json`), the ten cached Kepler
target-pixel files of the TPF host KIC 9655236 (Q3,4,5,7,8,9,11,12,13,14), and the two neighbour KBONUS
light curves the re-derivation downloaded (read from its scratch `untrusted_dl/`, never executed).
All new numbers below trace to the scripts and logs in `redteam/` (SHA-256 in `redteam/MANIFEST.sha256`).

**Verdict: survives.** None of the tests below produced a kill or a >3 sigma anomaly. The residual
weaknesses are listed in section 9 and must be carried into any further write-up.

## 1. Signal as inherited (for reference)
P = 6.70346 d (quote as 6.7036 +/- 0.0002), t0 = 265.1832 BKJD, duration 2.4 h, depth 638 +/- 78 ppm
(pipeline, masked) / 614 +/- 72 ppm (re-derivation), empirical SNR 8.2 / 8.5, BLS SDE 15.6 vs inverted
7.1, 125 events, odd/even 0.33 sigma, secondary -29 ppm (-0.37 sigma). Host from the pool table: G 17.10,
BP-RP 0.98, Teff 5170 K, Rstar 0.83 Rsun, 2.2 kpc, RUWE 1.04. Kepler channels by season (TPF headers):
season 1 (Q3,7,11,15) ch 23 mod 8.3; season 2 (Q4,8,12) ch 39 mod 12.3; season 3 (Q5,9,13) ch 63 mod
18.3; season 0 (Q6,10,14) ch 47 mod 14.3; skygroup 39. The target's KBONUS light curve has Q3,4,5,6,8,9,
10,12,14,15; the host TPF exists for Q3,4,5,7,8,9,11,12,13,14, so Q6/Q10/Q15 of the target come from a
different TPF and Q7/Q11/Q13 pixel tests use the host TPF only.

## 2. Is it already known? (`redteam/catalog_check.log`, Gaia/SIMBAD queries 2026-10-17)
- KIC 9655224, 9655236, 9655226: no row in KOI cumulative, DR25 TCE, or the Kirk+2016 Kepler EB catalogue.
- Within 3 arcmin: one KOI (K01714.01, KIC 9655129, P 2.744 d, 4% deep EB false positive, 145" away)
  and 17 DR25 TCEs on 11 stars; none at P, P/2, 2P, P/3 or 3P.
- Field-wide period coincidences (|dP/P| < 2e-3) with epoch check: at P the nearest are KIC 8757619
  (P 6.70842, 39 ppm, 89' away, epoch offset 8.6 h) and KIC 8314525 (P 6.70746, 819 ppm, 163', -6.6 h);
  at P/2 the EB KIC 9719636 (P 3.35157, 5.5% deep, Kp 17.3, 52' away, epoch offset 4.5 h; its neighbour
  KIC 9719634 carries the same period at 340 ppm). All epoch offsets exceed the 2.4 h duration, the
  period differences drift by 0.7-0.9 d over the 180 cycles, and for the 3.35-d EB a crosstalk ghost would
  put equal events at phase 0.5 of our fold, where the data show -29 ppm. No ephemeris match on any
  channel. (The MAST KIC sky-group query failed; the epoch test does not need it.)
- Gaia DR3 within 15": target (ruwe 1.04, ipd_frac_multi_peak 0, non_single_star 0, astrometric excess
  noise 0.09 mas), the G=15.9 star 2128120393543835648 = KIC 9655226 at 5.27" (plx 0.92 mas, foreground
  K dwarf), the host G=16.4 at 9.36" (plx 0.23), and a G=20.3 source at 11.8". No source within 5".
- SIMBAD 30": only BOKS 26915 (the host) and UCAC4 682-069738 (the G=15.9 star); nothing at the target.
- Bienias & Szabo 2025 report hot-Jupiter candidates with P 2.5-4.8 d and depths of per cent; Moquin et
  al. 2026 BOI-1.01 is P 3.37 d on a G=17.8 M dwarf. Neither matches this star or period. Their full
  target lists were not available locally; the exclusion is by period/depth regime, not by KIC lookup.

## 3. Pixel-level origin (own analysis: `redteam/pixel_check.{py,log,json}`, `redteam/pixel_fit.{py,log,json}`)
Method, independent of `pipeline/pixel_vet.py`: for each of the 10 host TPFs, each transit event is
differenced locally (per-pixel linear trend from the two flanking 0.25-d windows, evaluated at the
in-transit cadences), events are averaged, errors are the event-to-event scatter. Results in e/s of light
lost, so they can be compared with the 1.07 +/- 0.14 e/s the target's PSF light curve loses:
- Unclipped whole-stamp deficit, weighted over quarters: **1.32 +/- 0.34 e/s** (the pipeline's
  `deficit_total_e_s` of 3.5-4.6 e/s is a positive-clipped sum, pixel_vet.py line 221, and is not a
  measure of light lost).
- Deficit in the single pixel nearest each source (10 quarters): **target 0.49 +/- 0.08 e/s** (6.5 sigma);
  G=15.9 star 0.19 +/- 0.09; host 0.08 +/- 0.07.
- Joint three-source Gaussian-PRF amplitude fit (sigma 0.6 px; 322 pixels): **target 1.05 +/- 0.15 e/s**,
  G=15.9 star 0.39 +/- 0.18, host -0.06 +/- 0.14. Summed chi2: null 444, target-only 375, neighbour-only
  421, host-only 434, all three 353. Target-only beats neighbour-only by dchi2 = 46; with sigma 0.8 and 1.0
  px the ordering is unchanged (target 1.16 +/- 0.21 and 1.20 +/- 0.29 e/s).
- The pipeline's own `pixel.json`, summed over the 10 quarters, agrees: chi2 for a Gaussian at the target
  299.5 vs 328.6 at the G=15.9 star, 340.9 at the host, 324.8 scaled-image, 348.3 null. Its "deficit
  centroid closer to the target in only 2 of 6 quarters" statistic is a significance-weighted centroid
  of 2-3 pixels with the target sitting 0.3-1.3 px from the stamp edge in most quarters (Q8, Q11, Q12,
  Q14 margins 0.75, 0.29, 0.70, 0.36 px), which biases the centroid inward, toward the neighbours; it
  should not be read as evidence against the target.
- Total light lost agrees across all products: PSF 1.07 +/- 0.14, SAP 1.04 +/- 0.12 (re-derivation),
  whole stamp 1.32 +/- 0.34, PRF-fit target amplitude 1.05 +/- 0.15 e/s.
Reading: the light is lost at the target's position to the ~1-pixel (4") precision Kepler allows, and the
amount lost is what the target's own light curve says. A source blended within ~1 px of the target is not
excluded by any Kepler pixel test.

## 4. The two bright neighbours' own light curves (`redteam/neighbour_allq.log`)
Folded at the target ephemeris using every quarter they have (the re-derivation used only overlapping
quarters): KIC 9655226 (G 15.9, 5.3", Q3,7,10,11,13): PSF 70 +/- 37 ppm = 0.43 +/- 0.23 e/s, SAP
116 +/- 33 ppm = 0.76 +/- 0.22 e/s. KIC 9655236 (host, G 16.4, 9.4", Q3-17): PSF -40 +/- 32 ppm =
-0.15 +/- 0.12 e/s. For either to be the origin its own deficit would have to exceed 1.07 e/s (plus
whatever leaks into the target), i.e. >= 2 e/s; the G=15.9 star is >= 4 sigma below that and the host is
excluded outright. The 0.4-0.8 e/s "echo" in the G=15.9 star is the size expected from its 4-pixel
aperture and PSF-fit sharing with a target 1.3 px away, and matches the 0.39 +/- 0.18 e/s neighbour term
in my PRF fit; it is consistent with leakage but is a 2-sigma item that a true Kepler-PRF fit should close.

## 5. Eclipsing-binary diagnostics (`redteam/lc_checks.log`, own 1.0-d masked running-median detrend)
- Global depth 648 +/- 89 ppm; SAP 300 +/- 34 ppm (ratio 0.46; expected 0.48 from the raw FLUX/SAP level).
- Secondary: deepest box at any phase 0.1-0.9 is 150 ppm at phase 0.16 (1.7 sigma); phase 0.5 -29 ppm.
- Odd/even 0.33 sigma (pipeline), 0.39 (re-derivation): a 2P ephemeris with alternating eclipses is out.
- Shape: trapezoid grid on the folded curve prefers a flat bottom; a fully V-shaped (ingress = half
  duration) model is worse by dchi2 = 9.0, a 0.3-fraction ingress by 5.6. Expected ingress fraction for
  k = 0.025, b ~ 0.65 is ~0.04, consistent with the box preference.
- Duration/density: for Ms 0.80-0.90, Rs 0.83, a/Rs = 16.7-17.4 and central T14 = 3.0-3.15 h; the
  observed 2.33-2.44 h needs b ~ 0.65 (duration_ratio 0.76 in vet.json). Consistent with a dwarf host;
  implied radius 2.25 Re for depth 614 ppm if the host is the target. An undiluted stellar-mass companion
  on this host would give a much deeper, V-shaped, secondary-bearing signal; a grazing EB is disfavoured
  by the shape at ~3 sigma only.

## 6. Variability and aliases (`redteam/lc_checks.log`)
Lomb-Scargle 0.2-100 d on per-quarter-normalised undetrended flux: peak at 37.0 d (~900 ppm half
amplitude, quarter-scale systematics); power at 6.703 / 3.352 / 13.407 / 2.234 d ranks 9018 / 191008 /
1651 / 17413 of 200,000 frequencies. No rotation signal at P or its harmonics. Pipeline LS: 0.56 d, 70 ppm.
Blind-BLS alternatives (re-derivation): 1.434 d (dlogL 14.9) and P/2 (14.4) against 27.3 at P.

## 7. Kepler systematics (`redteam/lc_checks.log`, `redteam/flag_check.log`)
- Gap proximity: events with gap distance > 0.5 / 1 / 2 d give 547 +/- 74 / 525 +/- 76 / 540 +/- 81 ppm;
  Spearman(depth, gap distance) = -0.09 (p 0.34).
- Time trend of per-event depths: -204 +/- 173 ppm per 1000 d (-1.2 sigma).
- Per-season depths (same channel each season): season 0 (ch 47) 716 +/- 124, season 1 (ch 23) 368 +/-
  112, season 2 (ch 39) 790 +/- 176, season 3 (ch 63) 653 +/- 181 ppm; chi2 6.3 for 3 dof (p ~ 0.10).
  Per quarter: Q3 453, Q4 690, Q5 966, Q6 845, Q8 955, Q9 330, Q10 745, Q12 775, Q14 563, Q15 281 ppm
  (errors 156-347). The signal is present in every season; season 1 is 2 sigma shallow.
- Quality flags on the same cadences (host TPF QUALITY, 11 flag types with >= 1 occurrence): in-transit
  cadences are not closer than random to desaturations (z = -1.1), coarse point (0.4), Earth points, safe
  modes or attitude tweaks; the largest is reaction-wheel zero crossings at z = +2.6, unremarkable among
  11 tests and moot because SAP_QUALITY != 0 cadences are excluded before the search.
- Event dominance: 5 of 124 events are > 3 sigma from the mean depth (expected ~0.3; consistent with the
  20% error underestimate the re-derivation noted); removing the 5 highest-SNR events leaves 429 +/- 75
  ppm (5.7 sigma); median(in) - median(out) = 644 ppm equals the mean-based 648, so the signal is not
  carried by tails.

## 8. Dilution consistency
CROWDSAP 0.72-0.85 and FLFRCSAP 0.85-0.96 in the per-quarter headers are mutually inconsistent with the
PSF/SAP flux ratio (re-derivation), so no CROWDSAP-based dilution correction is applied. The absolute
light-lost agreement in section 3 is the usable statement: a resolved contaminator would lose more light
in the stamp/SAP than in the PSF fit, and none of the three products disagrees.

## 9. What survives and what does not (red-team summary)
Kill attempts that failed: known-object and ephemeris contamination; pixel localisation (target preferred
over both bright neighbours at dchi2 >= 29 in two independent codes, total light lost matches); the
neighbours' own light curves; secondary eclipse at any phase; odd/even; V-shape; duration vs host
density; stellar variability; gap/flag/trend correlations; single-event dominance.
Residual weaknesses, none individually above ~2.6 sigma:
1. Season-1 (ch 23) depth 368 +/- 112 ppm vs ~700 elsewhere (p ~ 0.1); Q15 281 +/- 162 is the "late
   weakening" the re-derivation flagged.
2. A 0.4 +/- 0.2 e/s echo of the dip in the G=15.9 neighbour's PSF flux and in my PRF fit, consistent
   with leakage but not demonstrated to be leakage (needs the real Kepler PRF, not a Gaussian).
3. Per-event errors underestimated by ~20%; the random-phase SNR (8.2-8.5) already absorbs this.
4. Kepler pixels cannot exclude a blended eclipsing source within ~4" of the target; Gaia lists nothing
   within 5" but is incomplete for G > 20 that close to a G = 17 star. This is the standard residual that
   only a TRICERATOPS-style FPP/NFPP calculation can quantify.
5. Host is G 17.1 at 2.2 kpc: ground-based confirmation is impractical; the claim can only ever be
   statistical.

What would settle it: (a) difference-image fit with the Kepler PRF model (prf files for ch 23/39/47/63)
instead of a Gaussian, reporting target vs neighbour amplitudes per season; (b) TRICERATOPS with the Gaia
DR3 neighbours and the KBONUS PSF-flux depth, reporting FPP and NFPP; (c) a per-season depth comparison
after a PRF-based re-extraction to test whether the season-1 shortfall is photometric; (d) a KIC lookup in
the Bienias & Szabo 2025 and Moquin et al. 2026 target lists once they are obtained.

Status after this review: a signal that passed stage 2, pixel-level checks (pipeline and this independent
test), an independent re-derivation and this red team. It is not a candidate until item (b) is run and
documented, and it is never "discovered".
