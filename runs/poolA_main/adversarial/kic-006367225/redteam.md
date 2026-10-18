# Red-team review: kic-006367225 (KIC 6367225, Gaia DR3 2077645628644907776, Kp 17.2)

Date 2026-10-17. Run id `poolA_main`. Candidate language only; nothing here is a planet.
Verdict: **KILL** — the periodic 12.972-d flux deficit is real but is not on this star. It is localised on the
unmodelled Gaia DR3 2077645628637721728 (G = 20.49) 3.6 arcsec north of the target, where its implied depth
is ~4 %, V-shaped, consistent with a background eclipsing binary. Every number below traces to a file under
`runs/poolA_main/adversarial/kic-006367225/redteam/` (own code, written for this review), to
`runs/poolA_main/vet/kic-006367225/vet.json`, `runs/poolA_main/pixel/kic-006367225/pixel.json`, or to the
independent re-derivation JSONs in `../rederive/`.

## 1. What was on the table
- Pipeline (vet.json, rules v1.7): P = 12.971528 d, t0 = BKJD 142.7530, box 2.25 h, pass-2 depth 3034 +/- 305 ppm
  (empirical SNR 9.85, box noise), 100 events, 88 % positive, odd/even 0.90 sigma, secondary -0.28 sigma, TLS SDE 13.4,
  Rp 5.28 Re (box) / 4.73 Re (TLS) for Rstar 0.879 Rsun. stage2_pass = true, known_matches = [].
- Pipeline pixel.json summary: deficit closer to the target than to any other Gaia source in only **2 of 14** quarters;
  median nearest-other-source-to-deficit 1.35 arcsec; "target Gaussian beats scaled image" in 11 quarters (that test
  compares the target against the bright TPF star, not against the faint neighbour).
- Independent re-derivation (../rederive/): recovered blind, but depth 1972 +/- 158 ppm (not 3034), SAP/PSF = 0.34,
  per-season PSF depths 1129-3925 ppm, raw PSF flux assigned to the target 0.25x-1.9x expected depending on module,
  G = 20.5 neighbour unmodelled by KBONUS and nearest to the pipeline's deficit in 4 quarters.

## 2. Own pixel-level localisation (the kill) — `redteam/pixel_own.py`, `pixel_own.json`, `pixel_own.png`
Method (no code shared with `pipeline/pixel_vet.py`): 17 KIC 6367220 long-cadence TPFs read as data
(`/root/.cache/lightkurve/.../kplr006367220-*_lpd-targ.fits.gz`, QUALITY == 0). For each event at the pipeline
ephemeris: in-transit |dt| < 1.125 h, local out-of-transit 2.5-8 h either side (kills slow trends); per-event
difference image (OOT mean minus in-transit mean) with per-pixel variance from the OOT scatter; inverse-variance stack
per quarter; events within 1 d of a detected one-cadence step in the summed stamp flux are dropped.
WCS systematics are calibrated per quarter by fitting a Gaussian to the mean image of the G = 13.8 star KIC 6367220
and shifting all Gaia positions by the measured offset (|offset| <= 0.12 px in every quarter). The stacked deficit is
fitted with a Gaussian of the stellar width (0.51-0.58 px) with free position; position error from a 150-draw
bootstrap over events; fixed-position chi2 for the target, the G = 20.5 star, KIC 6367220 and the G = 20.9 star.

Per quarter (peak SNR; distance of fitted deficit to target / to G20.5 in px; fixed-position chi2 target / G20.5):
```
Q2  mod 8  SNR  9.2  0.87 / 0.02   82 / 15
Q4  mod18  SNR  8.2  2.81 / 2.41   80 / 32   (fit poorly constrained, boot std 6 px)
Q6  mod 8  SNR 10.7  1.35 / 0.51  124 / 27
Q8  mod18  SNR  7.3  0.80 / 0.12   78 / 20
Q10 mod 8  SNR 15.2  0.91 / 0.08  222 / 15
Q12 mod18  SNR 12.0  1.87 / 1.01  153 / 39
Q14 mod 8  SNR 14.5  0.75 / 0.20  207 / 41
Q16 mod18  SNR  8.6  0.99 / 0.12   92 / 28
Q1,3,5,7,9,11,13,15 (mod 14/12): peak SNR 2.1-3.4, individually uninformative; all prefer G20.5 weakly.
```
Combined over the 13 quarters with peak SNR > 3, in sky coordinates:
- deficit offset from target = (-1.33, +2.85) +/- 0.27 arcsec; the G = 20.5 star is at (-1.00, +3.43) arcsec.
- distance to target **3.15 arcsec = 11.6 sigma**; distance to G20.5 star 0.66 arcsec = 2.4 sigma.
- chi2 of the per-quarter positions: at target 147.2 / 26 dof; at G20.5 18.0 / 26 dof.
- summed fixed-position chi2: target 1199.9, **G20.5 351.2**, KIC 6367220 1365.8, G20.9 1369.0, null 1474.5.
  Delta chi2 (target minus G20.5) = **849**.
- The offset direction is north (toward the G20.5 star), not toward or away from the bright star (SW) and not
  toward the stamp interior (edge truncation would bias toward lower columns), so neither a PSF-wing nor an
  edge-truncation systematic produces it. The same direction is seen in three seasons (modules 8, 18 and 12/14),
  i.e. it is fixed on the sky, not on the detector.

Controls — `redteam/pixel_tests.py`, `pixel_tests.json` (same code path, same TPFs):
- Null, ephemeris shifted by +0.5 P: median peak SNR 1.93, 1/17 quarters above 3, summed chi2 target 33.3 vs G20.5
  33.5 vs null 35.0 (no deficit anywhere). Null at +0.25 P: the same (chi2 39.2 / 38.6 / 39.4).
- Injection of a 2.8 e/s Gaussian deficit (the measured absolute deficit, section 3) at the **target** at a different
  epoch (+0.37 P): recovered 0.47 +/- 0.17 arcsec from the target, 3.87 arcsec from G20.5; chi2 target 434 vs G20.5 809
  (Delta 375 favouring the target).
- Same injection at the **G20.5 star**: recovered 3.73 arcsec from the target, 1.25 arcsec from G20.5; chi2 G20.5 417
  vs target 866 (Delta 449 favouring G20.5).
So the method separates the two stars (0.9 px apart) in both directions at the real signal's amplitude; the real
signal (Delta chi2 849 for G20.5, 0.66 arcsec from it) is unambiguously on the neighbour.

## 3. Flux bookkeeping confirms the pixels — `redteam/lc_checks.py`, `lc_checks.json`
Per-quarter box depths at the pipeline ephemeris from the KBONUS per-quarter extensions (masked 1-d running median,
raw un-stitched medians), grouped by module (inverse-variance means):
```
module  quarters        PSF depth ppm   PSF abs e/s   SAP depth ppm   SAP abs e/s   raw PSF e/s (expected 1401)
 8      2,6,10,14        973 +/- 125    2.30 +/- 0.29   1026 +/- 107   2.81 +/- 0.29   2246-2604 (1.6-1.9x)
12      3,7,11,15        700 +/- 276    1.11 +/- 0.43    395 +/- 122   0.81 +/- 0.25   1479-1679 (1.1-1.2x)
14      1,5,9,13,17     1691 +/- 396    1.63 +/- 0.38    236 +/-  95   0.90 +/- 0.36    863-1017 (0.6-0.7x)
18      4,8,12,16       6055 +/- 806    2.50 +/- 0.33    405 +/-  54   2.90 +/- 0.39    351- 445 (0.25-0.32x)
chi2(const, 3 dof): PSF ppm 43.1; PSF e/s 8.6; SAP ppm 35.0; SAP e/s 41.9
```
- The fractional PSF depth is not a property of the star: it tracks 1/(flux the KBONUS PSF fit assigns to the
  target), which varies 0.25x-1.9x between modules. The absolute PSF deficit is ~2 e/s in every season.
- The pipeline's 3034 ppm headline is reproduced (own masked box on the stitched flux: 3058 +/- 146 ppm, n_in 415),
  but it is dominated by the module-18 quarters where 16 of 24 events exceed 5000 ppm (events_by_module in
  lc_checks.json: weighted event depth 7371 ppm on module 18 vs 1252-2272 ppm elsewhere). The re-derivation's
  1972 ppm is the same signal with a different normalisation. Neither is an on-target depth.
- Maximum absolute deficit captured, 2.8-2.9 e/s (modules 8 and 18, 5-6 px apertures): on the Kp 17.2 target that
  would be ~2100 ppm; on the G = 20.49 neighbour (Kepler-band flux ~70 e/s, 3.3 mag fainter) it is **~4 %**.
  Flux alone allows both; the pixels (section 2) choose the neighbour.
- SAP absolute deficit varies 0.8-2.9 e/s with aperture (2 px on module 12 vs 6 px on module 18), i.e. the
  eclipsing source sits at the aperture edge, where the G20.5 star is (it lies further from the aperture centre
  than the target in every season).

## 4. Shape and duration
- Re-derivation trapezoid: T14 = 4.40 h, T23/T14 = 0.27 (V-shaped). For the nominal host (0.879 Rsun, logg 4.655 ->
  M = 1.27 Msun, a/Rs = 28.6) a k = 0.044 planet gives T14 <= 3.6 h and T23/T14 = 0.92 at b = 0. T23/T14 = 0.27 with
  k = 0.044 needs b >= 0.95 (grazing); at lower b it needs k >= 0.085-0.57, i.e. a stellar companion. Grazing
  geometries shorten T14 (1.8 h at b = 0.9), the opposite of the observed 4.4 h. A V-shaped, ~4 %, 12.97-d eclipse
  on a faint star is a textbook background EB; no secondary is required (PSF limit 158 ppm, in the neighbour's units
  ~0.3 %, so a faint or eccentric companion is allowed).

## 5. Known objects, ephemeris contamination, novelty (catalogue download date 2026-10-17, data/known/manifest.json)
- Within 3 arcmin: no KOI, no Kepler EB (Kirk+2016), no confirmed planet. DR25 TCEs on KIC 6367260 (0.97', Kp 12.3,
  eight TCEs, P = 0.83, 11.32, 46-179 d) and KIC 6367365 (2.2', P = 1.1, 41, 57, 86 d): no period within 0.2 % of
  12.9715 d or its 1/3, 1/2, 2x, 3x harmonics. KIC 6367220 and 6367225 are in none of KOI/TCE/EB.
- Field-wide period matches (tolerance 0.2 %): nearest is KOI-1798.01 = KIC 6867766 (Kepler EB, morph 0.03,
  P = 12.964718 d, depth 3733 ppm, Kp 14.4) at 39.2 arcmin. Excluded as the source: our fold at its exact ephemeris
  gives 213 +/- 158 ppm (lc_checks.json `koi1798_fold`), and the period difference 0.0068 d is 4.5x our period
  uncertainty (0.77 d of drift over 113 cycles). Same module/channel is therefore irrelevant. Kepler-1283 b
  (P = 12.946 d, 131') and Kepler-772 b (12.992 d, 223') are excluded by the same drift argument.
- Our channels per season: 48 (mod 14), 24 (mod 8), 40 (mod 12), 64 (mod 18) (KBONUS per-quarter headers).
- Bienias & Szabo 2025 Appendix A (Zenodo 14618891, PDF sha256 c9edf3a3..., read as data): seven hosts
  (KIC 4459924, 4543171, 8175131, 8293539, 9040849, 9674230, 11033410), none of our identifiers; their P 2.5-4.8 d.
  Moquin et al. 2026 BOI-1.01: P = 3.37 d, not ours. SIMBAD 15 arcsec: only KIC 6367220 (plain star); 60 arcsec:
  KIC 6367275, 6367260, nothing variable. No prior report of this signal found as of 2026-10-17; also none expected,
  since the eclipsing source is a G = 20.5 star with no KBONUS light curve.

## 6. Systematics checked and not implicated
- Gaps: 6/100 events within 1 d of a gap; mean depth near-gap 2611 vs far 3156 ppm; corr(depth, gap distance) = -0.004.
- Steps: the +6 %/+12 % one-cadence jumps in Q2 and Q6 of the KBONUS flux are absent from the summed TPF flux
  (step detector found none in Q2/Q6, four in Q15), so they are KBONUS photometry artefacts; events near steps were
  dropped from the pixel stacks anyway.
- Rotation/variability: LS top peak 0.7647 d, 189 ppm; power at P, P/2, 2P, P/3 is <= 0.17 of the peak. No alias.
- Odd/even (0.8-0.9 sigma) and secondary (-0.3 to 1.5 sigma) are clean, as expected for a faint EB with P = 12.97 d
  whose secondary is below the diluted detection limit.

## 7. Verdict and what would overturn it
KILL. Most likely explanation: a background eclipsing binary on Gaia DR3 2077645628637721728 (G = 20.49, 3.6 arcsec N
of the target), unmodelled by KBONUS-BKG, whose ~4 % V-shaped eclipse leaks ~2-3 e/s into the target's PSF fit and
aperture and is reported as a 1000-7000 ppm "transit" depending on how much flux the PSF model assigns to the target
that season. Nothing in this review supports a signal on KIC 6367225 itself.
To overturn: a difference-image or PRF forward-model fit (both stars free) placing the deficit within ~1 arcsec of
the target at > 3 sigma in at least two seasons, with the fixed-position chi2 preferring the target. The injection
test above shows that outcome would have been detected had it been true.

## 8. Files (this review)
- `redteam/lc_checks.json`  sha256 4928b339f5b0a5c3...
- `redteam/lc_checks.py`  sha256 1bb08b0ff7145728...
- `redteam/pixel_own.json`  sha256 4dc421a8cdd577ac...
- `redteam/pixel_own.png`  sha256 dd0883be47791eae...
- `redteam/pixel_own.py`  sha256 50377ef1e2086114...
- `redteam/pixel_tests.json`  sha256 fae6c828b4b6b008...
- `redteam/pixel_tests.py`  sha256 6dc99bebfbc46c7c...
- `../rederive/rederive_results.json`, `../rederive/rederive_extra.json` (independent re-derivation, read after own numbers existed)
