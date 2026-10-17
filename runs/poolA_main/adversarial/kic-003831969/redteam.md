# Red-team review: kic-003831969 (Gaia DR3 2100154070473374976, TPF host KIC 3831977)

Date: 2026-10-17. Run id: poolA_main. Reviewer role: adversarial (try to kill the signal).
Candidate language only: this is a stage-2 *signal* with pixel localisation; it is not a planet.

Inputs read: `runs/poolA_main/vet/kic-003831969/{vet.json,events.csv,fold_events.png}`,
`runs/poolA_main/pixel/kic-003831969/{pixel.json,diffimages.png}`, the independent re-derivation
(`adversarial/kic-003831969/rederive/`), `targets/poolA_kbonus_bg_dwarfs_ranked.csv` row 2404,
`data/known/*` (downloaded 2026-10-17 per file mtimes), KBONUS source catalogue, the 17 cached
TPFs of KIC 3831977 (`~/.cache/lightkurve/.../kplr003831977_*`), and the KBONUS light curve of
KIC 3831977 (downloaded today to scratch).

My own code and outputs are in `adversarial/kic-003831969/redteam_work/`
(`lc_tests.py -> lc_tests.json`, `pix_free_sigfix.py -> pix_free_sigfix.{json,png}`,
`season2_map.py`, `shape_test.py`). Ephemeris used throughout: P = 15.31335 d, t0 = BKJD 144.28486,
box duration 3.75 h (pipeline values from vet.json).

## Verdict: SURVIVES (with stated caveats)

I could not kill it. The signal is on (or within ~1 arcsec of) the target, is not on the TPF host,
matches nothing known, and has a flat-bottomed shape with no secondary and no odd/even difference.
What remains open is the generic unresolved-blend space (a faint EB or hierarchical triple inside
~1-2 arcsec of a G=16.5 star at 2 kpc), which only a statistical validation (TRICERATOPS) plus
high-resolution imaging can address. Two mild (2-3 sigma) tensions are recorded below; neither
kills a 3-4 Earth-radius signal, but both belong in any write-up.

## 1. Is it on the TPF host or a catalogued neighbour? (blended EB)

**Host KIC 3831977 (G=15.0, 7.9 arcsec E, Kp 15.03) is flat at the ephemeris.** KBONUS PSF flux:
-44 +/- 22 ppm (random-phase box error); SAP flux: +20 +/- 21 ppm; per-season values all within
+/-2.5 sigma of zero and of mixed sign (lc_tests.json, `host3831977`). A dip on the host leaking
into the target's PSF fit would need the host to dim by >= 1000 ppm x F_t/F_h / leak-fraction;
even for a generous 20% leak that is ~1200 ppm on a Kp 15 star with 17 quarters, which the Kepler
pipeline would have flagged (no DR25 TCE on KIC 3831977) and which is directly excluded by the
host fold. Host origin: killed.

**Free-position difference-image fit (my own, `pix_free_sigfix.py`).** For every quarter I built
per-event differenced images (in-transit core minus a per-pixel linear baseline from the local
out-of-transit window), inverse-variance co-added them by season (same pixel grid), and fitted a
2-D Gaussian (sigma fixed at 0.6 px from the free fits) with free centre. Results (pixel = 3.98"):

| season | quarters | target px | fitted deficit centre | offset | chi2 free / at target / null (npix) | chi2 at best catalogued neighbour |
|---|---|---|---|---|---|---|
| 0 | 2,6,10,14 | (0.71,1.02) | (0.79,1.05) +/- (0.08,0.09) | 0.33" | 20.2 / 21.2 / 129 (25) | 123 (S11, G19.2) |
| 1 | 3,7,11,15 | (1.09,0.64) | (1.06,0.57) +/- (0.12,0.08) | 0.32" | 27.0 / 28.1 / 110 (24) | 106 (N11, G19.9) |
| 3 | 5,9,13,17 | (1.08,0.99) | (1.12,1.08) +/- (0.11,0.10) | 0.40" | 39.3 / 40.2 / 120 (22) | 116 (SSW13, G20.5) |
| 2 | 4,8,12,16 | (0.51,0.46) | degenerate (target in corner pixel; fit runs off-stamp) | - | 103 / 166 / 197 (20) | 187 (S11) |

In three of four seasons the deficit is centred on the target to 0.3-0.4" (formal +/-0.4"), and
every Gaia source within 25" (host G15.0 at 8"; G19.2 at 11" S; G19.9 at 11" N; G20.5 at 13" SSW;
G19.2 at 16" NNW; G17.3 at 22" NNW; G20.2 at 26" E) is rejected by delta-chi2 > 65 in each season.
Because the pixel grid rotates by 90 deg between seasons, an off-stamp source would have to sit
within ~0.5 px (2") of the target in *every* orientation, i.e. be unresolved. The pipeline's own
per-quarter numbers (pixel.json: Gaussian-at-target beats scaled-image in 12/12 quarters; median
target-to-deficit 2.9") are weaker only because they are per-quarter (peak SNR ~3) and the
Gaussian-vs-scaled-image comparison is the wrong question when the host dominates the image.

Season 2 caveat: the target sits in the corner pixel (0.51,0.46) of a 4x5 stamp, so half its PSF
is off-stamp. The deficit is in the bottom row (pixels [0,0],[0,1],[0,2] at SNR 5.8, 7.6, 3.1) and
is elongated along the stamp edge towards +x; the remaining chi2 is a +/- dipole on the host's
bright pixels ([1,2] +5.1, [2,2] -5.6 sigma; a pointing-jitter signature, not a deficit). Season 2
alone cannot distinguish "target" from "just below the stamp", but the other three seasons fix the
sky position, so this is a geometry limitation, not evidence against the target.

**Deficit amplitude vs expectation.** Integrated Gaussian deficit in the stamp: 2.8, 2.6, 2.9 e/s
(seasons 0,1,3; errors ~0.6 e/s) vs expected 3217 e/s (target PSF flux) x per-season depth x
psffrac (0.83-0.84) = 2.9, 2.2, 3.2 e/s. The missing light is exactly what a ~1000 ppm dip on a
star of the target's brightness produces; there is no room for a much brighter hidden source. (It
cannot exclude a *fainter* unresolved EB with a deeper eclipse; see section 7.)

**SAP vs PSF (re-derivation + mine).** SAP depth 680 +/- 45 ppm vs PSF 988 +/- 57 ppm (ratio
0.69 +/- 0.06). Per season: 0.69, 0.63, 0.72, 0.81 (+/-0.10-0.15) against KBONUS CROWDSAP for the
target of 0.94, 0.88, 0.93, 0.87. If the dip is on the target and the PSF depth is the true depth,
the SAP depth should equal PSF depth x CROWDSAP; the measured ratio is ~0.2 lower (2-3 sigma
overall). Possible benign causes: uncatalogued light/background in the 4-8 px aperture, or the PSF
depth being ~20% over-estimated for a star whose PSF is truncated by the stamp edge. Recorded as a
mild tension, not a kill. **Resolved in the addendum below**: with the measured SAP and PSF flux
levels in place of the header CROWDSAP the expected ratio is 0.71 and the tension disappears.

## 2. Stellar variability / rotation aliases

vet.json: Lomb-Scargle peak at 0.306 d with 41 ppm half-amplitude (4% of depth), no period match.
The inverted-flux null tops out at SDE 9.5 (re-derivation) vs 52.6 for the real data. The BLS peaks
2 and 3 are the P/2 and 2P aliases of the same signal, not independent periodicities. No kill.

## 3. Kepler systematics (gaps, quarter edges, safe modes, rolling band, thermal ramps)

- Channels per season (TPF headers): 12 (mod 4.4), 16 (mod 6.4), 76 (mod 22.4), 72 (mod 20.4).
  The signal is present in all four channels/seasons (section 6 table), which rules out a
  single-channel electronic artefact (rolling band, column anomaly) as the source.
- Events near gaps: 10% within 1 d, 16% within 2 d (vet.json), consistent with random placement
  for a 15.3-d ephemeris in a 1460-d baseline with ~40 gaps. Single-event dominance 0.076; 91% of
  86 events individually positive; removing any one event changes nothing.
- 15.31 d is not a Kepler cadence (31.0-d downlink, 93-d quarters, 372.5-d orbit, 3-d momentum
  dumps) nor a harmonic of one. Events drift through the quarter (93/15.31 = 6.07).
- Per-event depth scatter is 1.5x the white-noise expectation (chi2 195/85 pipeline, 220/85
  re-derivation, 35/16 per quarter in my test). This is ordinary red noise at the 4-h timescale
  for a G=16.5 PSF light curve and is also why the pipeline's empirical SNR (23.8) is optimistic
  relative to the random-phase estimate (19.6). Not a kill, but quote the conservative SNR.

## 4. Ephemeris-match contamination from known objects

Local catalogues (data/known, downloaded 2026-10-17): KOI cumulative (9564 rows), DR25 TCE (34032),
Kepler EB catalogue Kirk+2016 (2876), confirmed planets in the field, TESS EB catalogue.
- Nothing of any kind within 3 arcmin of the target (KOI, TCE, EB, confirmed planet, TESS EB).
- Period match at any harmonic 1/3, 1/2, 2/3, 1, 3/2, 2, 3 within 0.2%: 40 KOIs, 70 TCEs, 8 EBs
  across the field, the nearest at 81.5 arcmin (K07684.01, P=10.224 d = 2/3 P, phase 0.865: not a
  match). The closest period twin, K00475.02 (P=15.31296 d, 570 arcmin away), has our t0 at phase
  0.204 of its ephemeris (3.1 d off): not a match. The only phase-coincident objects (K01430.02
  P=22.93 d phase 0.018; TCE 5726667 P=22.94 d phase 0.998; TCE 9596300 P=46.01 d phase 0.995)
  are 3/2 or 3x harmonics 9-11 deg away, where only every third of our events would coincide, and
  they are not on modules 4, 6, 20 or 22. No known source anywhere in the field matches both
  period and epoch, so same-module video crosstalk from a catalogued EB/KOI is excluded.
- Also required for crosstalk: a deficit of ~3 e/s at 1e-3 to 1e-4 coupling implies a source
  eclipse of 3e3-3e4 e/s, i.e. a bright (Kp < 11) deep EB with this ephemeris, of which none exists.

## 5. Already known?

- KIC 3831969 / KIC 3831977: not in KOI, DR25 TCE, Kirk+2016 EB, confirmed-planet lists.
- Bienias & Szabo 2025 (A&A 694, A209; Appendix B from Zenodo 10.5281/zenodo.14618752, fetched
  today): their seven hot-Jupiter candidates are KIC 4459924, 9674230, 4543171, 9040849, 8175131,
  8293539, 11033410 (P = 2.5-4.8 d). Neither our KIC, the host, nor the Gaia ID appears.
- Kepler Bonus team: no published KBONUS-BKG planet-candidate list exists (the 2021 EXBA-mask paper
  reported one candidate and 69 EBs from custom-aperture regions; our star is in a normal target
  TPF, not an EXBA mask). A 2026 workshop poster (Moquin et al., "BOI-1.01") is a single unrelated
  object.
- SIMBAD cone (30"): only UCAC4 645-065182 (a positional catalogue entry, no variability type).
- TOI/CTOI (ExoFOP lists downloaded 2026-10-17: 8151 TOIs, 5138 CTOIs; `redteam_work/toi_check`
  below): target is TIC 120425498 (Tmag 16.01, 0.03" from the Gaia position), host TIC 120425499
  (Tmag 14.49); neither, nor the four faint neighbour TICs (1715992116-9), is a TOI or CTOI, and
  there is no TOI/CTOI within 3 arcmin. (TESS could not reach 1000 ppm at T=16 anyway.)

## 6. Depth consistency across seasons / quarters (my per-season test, PSF flux)

Season (Q1,5,9,13,17): 1199 +/- 135 ppm; (Q2,6,10,14): 1089 +/- 111; (Q3,7,11,15): 836 +/- 100;
(Q4,8,12,16): 772 +/- 102. chi2 about the weighted mean 939 ppm = 9.3 for 3 dof (p = 0.03) with
white errors, 4.1/3 after inflating errors by the measured per-event excess (x1.5). SAP flux shows
the same ordering (830, 688, 606, 625). A 2-sigma hint that the two "shallow" seasons (channels 76
and 72, where psffrac is lowest, 0.77-0.83) under-recover the dip, i.e. that the recovered depth
depends slightly on how much of the target's PSF is on the stamp. Consistent with a genuine dip on
the target whose PSF photometry is imperfect at the stamp edge; it does not look like a blend
signature (a blend from a fixed off-target position would also be season-dependent, but the
pixel fits already place the deficit on the target in all three usable seasons).

## 7. Shape, secondary, odd/even, duration, density, implied radius

- Shape (`shape_test.py`, trapezoid with 29.4-min cadence supersampling, masked re-detrend):
  depth 1030 ppm, T14 = 4.06 h, T23/T14 = 0.80; delta-chi2 = 2 for T23/T14 in [0.60, 0.99];
  a V-shape (T23/T14 <= 0.4) is disfavoured by delta-chi2 >= 8 and a pure V by 22. The
  re-derivation's T23/T14 = 0.65 did not include cadence smearing and should not be read as a
  grazing geometry.
- Secondary at phase 0.5: -4 +/- 48 ppm (re-derivation), +1 ppm (pipeline), -27 ppm (mine):
  3-sigma limit ~140 ppm = 14% of the primary. Odd/even 0.3-0.5 sigma. Deepest non-primary box
  anywhere in phase 120 ppm. A twin-EB at 2P is excluded; an EB with a secondary < 14% of the
  primary (cool companion) is not.
- Duration vs host: for the catalogue host (0.83 Rsun, assume 0.85 Msun) a/R* = 29.6 and a central
  transit lasts 4.08 h; measured 4.06 h (box-smeared trapezoid) to 4.29 +/- 0.25 h (re-derivation),
  i.e. b ~ 0 on the assumed host. At 1968 pc and this line of sight A_G ~ 0.4-0.8 mag is plausible,
  which would make the host bluer and larger (~1.0-1.1 Rsun, T14(b=0) = 4.7-5.0 h) and the
  duration comfortably central-to-moderate b. Either way: consistent with a dwarf host, no
  giant/subgiant required.
- Implied radius: rp/rs = 0.030-0.033 -> 2.7-3.0 Re on 0.83 Rsun, 3.6-4.0 Re on 1.0-1.1 Rsun.
  A sub-Neptune at 15.3 d around a K/G dwarf: ordinary for Kepler. Nothing exotic needed.
- Gaia: RUWE 1.057, parallax S/N 12.6, no sign of a resolved companion.

## 8. What would still kill it (and what to do next)

1. An unresolved EB/hierarchical triple within ~1-2": run TRICERATOPS with the Gaia DR3 field and
   the measured depth/duration; the secondary limit (< 14% of primary) and odd/even (< 1 sigma)
   are the inputs. High-resolution imaging (Keck/Gemini AO or speckle) of a G=16.5 target is
   feasible if anyone cares enough.
2. A planet-sized companion on a different star inside the Kepler PSF but not in Gaia (G > 20.5
   within 1"): this is the same TRICERATOPS space.
3. The season-dependence of the recovered depth (section 6) should be re-measured with a proper
   PRF model fit rather than a Gaussian before any depth or radius is quoted to better than ~25%.
   (The SAP/PSF vs CROWDSAP discrepancy is resolved in the addendum: header CROWDSAP understates
   the host's contamination of the SAP aperture.)
4. Quote the conservative empirical SNR (19.6) and the per-event red-noise-inflated depth error
   (~70-80 ppm), not the pipeline's 23.8 / 42 ppm.

## Checks run (summary)

- Host KIC 3831977 KBONUS PSF and SAP flux folded at the ephemeris (flat).
- Free-centre Gaussian deficit fit on season-stacked, per-event-differenced TPF images, with chi2
  at target and at every Gaia neighbour; integrated deficit vs expected.
- Target PSF and SAP depth per season and per quarter; SAP/PSF vs CROWDSAP.
- Trapezoid shape fit with cadence supersampling (V-shape test).
- Ephemeris/harmonic match vs KOI, DR25 TCE, Kirk+2016 EB, confirmed planets, TESS EBs; 3-arcmin
  cone in all; channel/module per season from TPF headers.
- Novelty: KOI/TCE/EB/confirmed, Bienias & Szabo 2025 Appendix B, KBONUS team literature, SIMBAD.
- Duration/density/implied radius for the catalogue host and an extinction-corrected host.
- Re-read of pipeline vet.json, events.csv, fold_events.png, pixel.json, diffimages.png and the
  independent re-derivation.

## Addendum (2026-10-17, after the session was interrupted and resumed)

**Traceability of the pixel numbers.** In `pix_free_sigfix.json` the field `fit_sigma_px` (0.7) is the
unused initial guess; the Gaussian sigma is fixed at 0.6 px (`SIGFIX` in `pix_free_sigfix.py`). The
"integrated Gaussian deficit" quoted in section 1 (2.8, 2.6, 2.9 e/s for seasons 0, 1, 3) is
`fit_amp` x 2 pi x 0.6^2 (fit_amp = 1.22, 1.14, 1.28). The JSON field `deficit_sum_e_s` (3.1, 1.6, 1.7)
is a different quantity, the raw sum of the whole difference image including noise pixels, and is
consistent with the Gaussian integral within its ~0.65 e/s error.

**CROWDSAP in the target table is the per-quarter minimum.** The ranked target CSV gives
crowdsap 0.763, flfrcsap 0.734, psffrac 0.754; those are the minima over quarters (Q12, Q4, Q16).
The per-quarter LIGHTCURVE_Qn headers give CROWDSAP 0.76-0.97 (season means 0.94, 0.88, 0.93, 0.87),
FLFRCSAP 0.73-0.97, PSFFRAC 0.75-0.85. The independent re-derivation's "~0.76 expectation" for the
SAP/PSF depth ratio therefore used the most contaminated quarter, not the typical one.

**SAP/PSF depth ratio from measured flux levels (`sap_psf_fluxlevel.py` -> `.json`).** Rather than
trust the header, use the light curves themselves: median SAP flux 4051 e/s vs PSF flux 3217 e/s
(mean FLFRCSAP 0.898). If the dip is on the target, SAP depth / PSF depth = FLFRCSAP x F_psf / F_sap
= 0.714 (0.794 if the aperture captured all the target light). Measured 0.689 +/- 0.060: 0.4 sigma
(1.8 sigma) from expectation. Per season: expected 0.83, 0.70, 0.71, 0.61 vs measured
0.69 +/- 0.13, 0.63 +/- 0.10, 0.73 +/- 0.13, 0.81 +/- 0.15 (all within 1.3 sigma). The "tension" in
section 1 was an artefact of header CROWDSAP (0.87-0.94) understating how much of the 4x brighter
host 8" away lands in the 4-8 pixel SAP aperture; the real SAP contamination is ~26-29%. This also
adds an argument independent of the pixel fits: an eclipse on a *resolved* neighbour would have to
contribute more to the SAP aperture than to the target's PSF-fit flux, giving a ratio > F_psf/F_sap
= 0.79 and typically >> 1; the measured 0.69 +/- 0.06 is below that bound. Only an unresolved source
(inside the PSF core) evades it, the same residual space as before.

**Net effect on the verdict**: unchanged (survives). One of the two caveats (SAP/PSF) is withdrawn; the
other (season dependence of the PSF depth, chi2 9.3/3 white, 4.1/3 red-noise inflated) stands and is
a photometry-systematics caveat for the quoted depth, not evidence against the signal.

Files added in this addendum: `redteam_work/sap_psf_fluxlevel.py`, `redteam_work/sap_psf_fluxlevel.json`,
`redteam_work/toi_check.txt`.
