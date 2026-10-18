# Red-team review: kic-004474607 (Gaia DR3 2052438633081525632, TPF host KIC 4474637)

Date: 2026-10-17. Run id: poolA_main. Reviewer role: adversarial (try to kill the signal).
Candidate language only; nothing here is a planet. This is a stage-2 *signal* whose pixel
localisation (pixel_v15, pixel_v16) was already "inconclusive".

Inputs read: `runs/poolA_main/vet/kic-004474607/{vet.json,events.csv,fold_events.png}`,
`runs/poolA_main/pixel_v15/kic-004474607/pixel.json`, `runs/poolA_main/pixel_v16/kic-004474607/{pixel.json,diffimages.png}`
(no `runs/poolA_main/pixel/kic-004474607/` exists), the stage-1 record in `runs/poolA_main/results.jsonl`,
the independent re-derivation (`adversarial/kic-004474607/rederive/{summary.json,eb_crosscheck.json,events_gallery.png}`),
`targets/poolA_kbonus_bg_dwarfs_ranked.csv` rank 6100, `pipeline/prereg_poolA.json` (thresholds only),
`data/known/*` (manifest fetched 2026-10-17T15:46Z), the KBONUS light curve
`data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-004474607_kepler_v1.0_lc.fits` (headers and
LIGHTCURVE_STITCHED), and cached ExoFOP TOI/CTOI lists and Bienias & Szabo 2025 text from earlier
red teams (scratchpad `verify-landscape/downloads/{toi,ctoi}_20261017.csv`, `verify-kbonus-dataset/bienias.txt`).

My own code and outputs: `adversarial/kic-004474607/redteam_work/` (`lc_tests.py -> lc_tests.json`,
`lc_tests.log`, `q13_events_raw.png`, `MANIFEST.sha256`). No code under `pipeline/` was read or imported.
Ephemeris used: P = 32.94271710658947 d, t0 = BKJD 754.5746921888925, box 1.95 h (vet.json).

## Verdict: KILL

The 32.94 d signal is a Kepler post-downlink thermal-settling artefact in Q13, not an astrophysical
dip. Every statistic that made it pass the first screen and stage 2 comes from two "events" that are
the first 4 cadences after a data gap. With those cadences masked, the period is not detected and the
inverted-flux null wins.

## 1. The two deep events are post-gap ramps (Kepler systematics)

`q13_events_raw.png` shows raw per-quarter-normalised PSF and SAP flux around epochs 13 and 14:

| epoch | tc (BKJD) | nearest segment start (BKJD) | tc - start | local depth (ppm, raw-flux linear baseline) |
|---|---|---|---|---|
| 13 | 1182.830 | 1182.737 (Q13 start, after 0.715 d gap) | 0.093 d | 8629 +/- 1303 |
| 14 | 1215.773 | 1215.718 (after 0.736 d monthly-downlink gap) | 0.055 d | 8130 +/- 1777 |

(`lc_tests.json: events_pipeline_ephemeris`). In both cases the "in-transit" cadences are the first
four after the gap; the raw flux there is 5,000-10,000 ppm below the following day, then rises to
+5,000 ppm (PSF) / +15,000-20,000 ppm (SAP) within 0.3-0.5 d and decays over ~2 d. That is the
classic thermal ramp after a downlink/pointing gap. The pipeline's own `events.csv` lists these two
events at 8061 and 8101 ppm with `gap_dist` 0.093 and 0.055 d; the other 15 events are 2.9-23.8 d
from any segment start.

**Control (same recipe, no ephemeris):** a 1.95 h box placed 0.14 d after *every* segment start of
the light curve (37 starts; `lc_tests.json: control_box_0p14d_after_every_segment_start`) gives
"depths" of 5947 +/- 1454, 4341 +/- 2137 and 3572 +/- 858 ppm at the three Q13 monthly segment starts
(BKJD 1182.74, 1215.72, 1245.86), i.e. *every* Q13 segment start produces a 3.5-6 ppt "transit" at this
offset, including the one (1245.86) that is not on the ephemeris. Outside Q13 the same control gives
+3942 (Q11, BKJD 1096.65) and +2387 ppm (Q15) at two other starts and noise elsewhere. Q13 is on
channel 49 (module 15, output 1; same as Q9 and Q17), the season with the lowest PSFFRAC (0.61-0.62)
and only 2-3 SAP pixels.

**The period is the downlink spacing.** The two affected segment starts are separated by
1215.718 - 1182.737 = 32.981 d. The re-derivation's best multi-duration period is 32.98173 d, i.e. the
"period" of the signal is the interval between two consecutive Q13 monthly downlinks. The pipeline's
32.9427 d and TLS 32.9392 d are teeth of the same alias comb that still put both post-gap cadence
groups inside the box (the two ephemerides coincide only at epochs 11-13, as the re-derivation noted).

## 2. Without the post-gap cadences there is no detection

My own re-detrend (per-segment 1 d running median with the 32.94 d transits masked, the pipeline's
recipe) reproduces the pipeline's headline depth exactly: box depth at the pipeline ephemeris
**2162 +/- 190 ppm** (vs `depth_pass2` 2141 ppm; `lc_tests.json: box_depth_psf_all`). So the headline
number is reproducible; it is simply dominated by 8 post-gap cadences out of 64. Masking 0.5 d after
every segment start (856 of 27,819 cadences):

| statistic | all cadences | post-gap (0.5 d) masked |
|---|---|---|
| PSF box depth at pipeline ephemeris | 2162 +/- 190 ppm | 884 +/- 191 ppm |
| SAP box depth at pipeline ephemeris | 1581 +/- 164 ppm | 603 +/- 166 ppm |
| per-event weighted stack (17 / 15 events) | 1090 +/- 182 ppm | 862 +/- 185 ppm |
| BLS 25-40 d (transit-unmasked blind re-detrend, durations 1.5-4.5 h, 15.4k periods): top peak | 32.939 d, SDE 10.3 | 27.10 d, SDE 6.0 |
| BLS local SDE at 32.94 +/- 0.05 d | 10.3 | **2.2** |
| inverted (2-f) BLS 25-40 d top peak | 35.06 d, SDE 7.5 | 35.06 d, SDE 7.2 |

(`lc_tests.json: bls_25_40d`, errors from 2000 random-phase boxes). With the post-gap cadences
removed, 32.94 d is a 2-sigma tooth of the periodogram, the best real peak (SDE 6.0) is below the
inverted null (7.2), and nothing in the 25-40 d window passes the pre-registered SDE >= 9 /
SNR >= 7.1. The residual 884 +/- 191 ppm box at the fixed ephemeris is a post-hoc measurement at a
period that was selected by the artefact; it is not an independent detection (the masked BLS does not
pick it out), it is 4.6 sigma on an estimator whose grid-maximum noise is of that order, and the
pipeline's own quartile depths (490-1270 ppm) and the re-derivation's 609 +/- 181 ppm without epochs
13/14 agree it is sub-threshold. The 2012 per-year depth (2097 +/- 296 ppm vs 494 and 830 in
2011/2013, re-derivation) and the per-event chi2 (45.3/16 pipeline, 70/16 re-derivation) are the same
two events.

The inverted-flux top peak (35.06 d, 1.5 h, "depth" 4.8 ppt) is the *positive* side of the same ramps
(the +5-15 ppt overshoot 0.3-0.5 d after the Q13/Q15 gaps), which is further evidence that the long-
period periodogram of this light curve is set by gap systematics in both polarities.

## 3. Pixel localisation is contaminated by the same ramp

pixel_v15/v16 (8 TPF quarters) are "inconclusive": 0 conclusive quarters, target at the stamp edge or
outside in 6 of 8 (target_px ~ (0.4,-0.4) in season 3, (1.0,0.15) in season 2). v16's chi2 comparison
prefers Gaia DR3 2052438568667567488 (G 17.5, 28.7 arcsec away, outside the stamp) over the target by
dchi2 = 183; that figure is Q13 alone (chi2 best-alternative 1735 vs target 1915; the other three
quarters sum to ~ -3). In Q13 the "deficit" map has peak SNR 21.8, depth-in-aperture 18,574 ppm and
total missing light 101 e/s, centred 4.2 arcsec from the stamp's flux centroid (the bright star KIC
4474645 / KOI-1657, G 14.2), i.e. the whole stamp dims in the post-gap cadences, so the brightest
pixels dominate the difference image. This is the ramp, not an eclipsing source. The pixel evidence
therefore neither localises the signal nor points at a real alternative host.

## 4. Blended / background EB

- TPF host KIC 4474637 = KOI-3829.01 (FALSE POSITIVE, `koi_fpflag_ss=1`; Kirk+2016 EB, P = 3.8867254 d,
  depth 37 %, Kp 16.81, 18.3-20 arcsec). KIC 4474645 = KOI-1657.01 (FALSE POSITIVE, P = 3.886732 d,
  Kp 14.2, 24.6 arcsec; flagged `co`/`ec`, i.e. contaminated by the same EB). Both are 3.8867 d
  objects; 32.94/3.8867 = 8.476 is not a low-order rational (re-derivation: 17/2 drifts 1.4 d over
  the baseline). The re-derivation's fold of the target PSF flux on the EB ephemeris gives 13 +/- 21
  ppm (0.004 % of the EB depth): the EB does not leak into the target's PSF flux, and excluding EB
  eclipse cadences leaves the 32.9 d box depth unchanged (1009 +/- 223 ppm). Epoch 13 is at EB phase
  0.60 (out of eclipse), epoch 14 at 0.08 (in primary eclipse); their depths are the same, so the EB
  is not what makes those two events deep either.
- SAP/PSF depth ratio: with post-gap cadences masked, 603/884 = 0.68 +/- 0.24 (CROWDSAP 0.70-0.94
  predicts ~0.7-0.9 for a target-hosted signal, > 1 for a PSF-removed contaminant): non-discriminating
  at this SNR, and moot once the signal is gone.
- Odd/even (1.0 sigma pipeline, 0.4 re-derivation) and secondary (-0.07 / -1.5 sigma) are null, as
  expected for a two-event artefact at a 33 d period.

## 5. Stellar variability / rotation aliases

Lomb-Scargle (vet.json): best period 1.18 d, half-amplitude 68 ppm, LS power 0.0007, amplitude/depth
0.03, no match to 32.94 d or harmonics. The stage-1 second and third peaks (0.598, 0.622 d, 93 ppm)
also appear as the *inverted* top peaks (0.622, 0.598 d), i.e. low-amplitude red noise at both
polarities. Not the explanation of the 32.9 d peak, which is explained in section 1.

## 6. Ephemeris-match contamination and "already known?"

- 3 arcmin cone (KOI, DR25 TCE, Kirk+2016 EB, confirmed planets, TESS EB; catalogues dated
  2026-10-17): KOI-3829.01 / TCE 4474637-01 (3.8867 d), KOI-1657.01 / TCE 4474645-01 (3.8867 d) and
  six DR25 TCEs on KIC 4474484 (0.81-349 d, 89 arcsec; none within 0.2 % of 32.94 d or its 1/3-3x
  harmonics). No ephemeris match.
- Whole-field period match (0.2 % tolerance, harmonics 1/3-3): 30 KOIs and ~80 DR25 TCEs land near
  the harmonics (33 d is a crowded period for DR25 TCEs), but the closest in phase, TCE 4663658-01
  (32.9614 d, 15 arcmin, 2.05 h), is 2.17 d off in epoch; none is within one duration in phase.
  Channels used: 77 (Q8/12/16), 49 (Q9/13/17), 33 (Q11/15); no same-channel object with this
  period+epoch is needed because the signal is explained by gap timing, not by a contaminant.
- Not previously reported: KIC 4474607 / Gaia DR3 2052438633081525632 is absent from the KOI, DR25
  TCE, Kirk+2016 EB, confirmed-planet and TESS-EB lists; not among the Bienias & Szabo 2025 candidates
  (KIC 4459924, 4543171, 8175131, 8293539, 9040849, 9674230 in the cached text); no TOI/CTOI within
  3 arcmin (ExoFOP lists of 2026-10-17; 8151 TOIs, 5138 CTOIs); SIMBAD within 40 arcsec returns only
  KOI-3829 (EB*, 18.3 arcsec) and KOI-1657 (24.6 arcsec), nothing at the target. Novelty is moot for an
  artefact, but is recorded for completeness.

## 7. Duration / density / implied radius

Catalogue host: Teff 5435 K, Rstar 0.87 Rsun, BP-RP 0.92, d = 4.18 kpc, RUWE 0.93. For
Mstar ~ 0.95 Msun, a = 0.198 AU, T14(b=0) = 5.16 h (vet.json: 5.59 h). The observed 1.95 h box needs
b = 0.93 (grazing), which is allowed but would make the pipeline's `core_to_box_depth_ratio` (0.93) and
flat-looking fold mildly surprising. Implied radius 3.2 Re for 1100 ppm, 4.4 Re for 2141 ppm, 2.3 Re
for the 600-ppm post-mask residual. Not a kill by itself; listed only because the write-up asked.

## 8. Why the pipeline passed it (for the audit trail; no thresholds changed)

Stage 1: SDE 12.9 (real) vs 8.8 (inverted), SNR 9.9, max_single_event_frac 0.39 (< 0.5). Stage 2 v1.7:
gap_frac_1d 0.18 (< 0.5), event reduced chi2 2.83 (< 3.0), snr_emp_pass2 11.6, TLS SNR 8.7. All passes
are driven by the same 8 cadences. The existing rules measure gap proximity per event (`gap_dist`) but
do not exclude post-gap events from the depth/SNR stack, and the chi2 ceiling of 3.0 was set from
controls (1.95, 1.50) with a gap systematic at 25, so a two-event artefact among 17 slips under it.
Any new cut (e.g. "re-stack with events < 0.5 d after a segment start removed and require SNR >= 7.1
and a masked BLS peak at the same period") must be logged in `prereg_poolA.json` `revision_log` and
every affected result labelled post-hoc. I did not change anything under `pipeline/`.

## What would settle it (if anyone wanted to resurrect it)

- Nothing cheap: the masked BLS (section 2) already shows no peak at 32.94 d. A full 1-40 d search of
  the post-gap-masked light curve with the pre-registered thresholds and its inverted null would
  formally close it; I expect no pass.
- The ~600-900 ppm residual at the fixed ephemeris would need > 10 additional clean transits
  (i.e. data that do not exist for Kepler) to become a 7-sigma, independently-selected detection.

Online sources consulted: SIMBAD (astroquery, 2026-10-17). Everything else is local and dated above.
