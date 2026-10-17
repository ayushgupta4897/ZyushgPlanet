# Red-team review: kic-004474607 — verdict: KILL

Full notes: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-004474607/redteam.md`; verdict JSON: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-004474607/redteam_verdict.json`; my code and outputs: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-004474607/redteam_work/` (`lc_tests.py`, `lc_tests.json`, `lc_tests.log`, `q13_events_raw.png`, `MANIFEST.sha256`). Candidate language only; nothing here is a planet. No code under `pipeline/` was read or imported.

## The kill

**The two deep events are post-gap thermal ramps.** Epoch 13 (BKJD 1182.830) and epoch 14 (BKJD 1215.773) are the first four cadences after segment starts at BKJD 1182.737 and 1215.718 (0.093 and 0.055 d). In `q13_events_raw.png` the raw PSF and SAP flux there sit 5-10 ppt below the following day, then overshoot by +5 ppt (PSF) / +15-20 ppt (SAP) within 0.3-0.5 d and decay over ~2 d: the classic Kepler post-downlink settling ramp. Local depths 8629 +/- 1303 and 8130 +/- 1777 ppm vs 0-2 ppt for the other 15 events (`lc_tests.json: events_pipeline_ephemeris`), matching the pipeline's own `events.csv` (8061/8101 ppm, gap_dist 0.093/0.055 d).

**The period is the downlink spacing.** 1215.718 - 1182.737 = 32.981 d is the re-derivation's best period (32.98173 d); the pipeline's 32.9427 d is a tooth of the same alias comb that still captures both post-gap cadence groups.

**Control without an ephemeris.** A 1.95 h box 0.14 d after *every* segment start (37 starts) gives 5947 +/- 1454, 4341 +/- 2137 and 3572 +/- 858 ppm at the three Q13 monthly starts (BKJD 1182.74, 1215.72, 1245.86), the third of which is not on the ephemeris; Q13 is channel 49, PSFFRAC 0.61-0.62, 2-3 SAP pixels.

**Without the post-gap cadences there is no detection** (0.5 d masked after each segment start, 856 of 27,819 cadences):

| statistic | all cadences | post-gap masked |
|---|---|---|
| PSF box depth at pipeline ephemeris | 2162 +/- 190 ppm (reproduces depth_pass2 2141) | 884 +/- 191 ppm |
| SAP box depth | 1581 +/- 164 ppm | 603 +/- 166 ppm |
| per-event weighted stack | 1090 +/- 182 (17 ev) | 862 +/- 185 (15 ev) |
| BLS 25-40 d top real peak | 32.939 d, SDE 10.3 | 27.10 d, SDE 6.0 |
| BLS local SDE at 32.94 d | 10.3 | **2.2** |
| inverted (2-f) top peak | 35.06 d, SDE 7.5 | 35.06 d, SDE 7.2 |

The real peak is below the inverted null and below every pre-registered threshold (SDE 9, SNR 7.1). The residual 600-900 ppm at the fixed ephemeris is a post-hoc measurement at an artefact-selected period, not an independent detection. The inverted top peak (35 d, +4.8 ppt, 1.5 h) is the positive overshoot of the same ramps.

**Pixel evidence is the same ramp.** pixel_v15/v16 are inconclusive (0 conclusive quarters; target at the stamp edge in 6/8). v16's dchi2 = 183 preference for Gaia DR3 2052438568667567488 (G 17.5, 28.7 arcsec, off-stamp) is Q13 alone (1735 vs 1915; other quarters sum to ~ -3), where the "deficit" is 18,574 ppm in aperture with peak SNR 21.8 centred 4.2 arcsec from the stamp's flux centroid: the whole stamp dims in the post-gap cadences.

## What did not kill it (recorded)
- EB host KIC 4474637 (KOI-3829.01 FP, 3.8867 d, 37 %, 18-20 arcsec) and KIC 4474645 (KOI-1657.01 FP, same period, 24.6 arcsec): no leak into target PSF flux (13 +/- 21 ppm), 32.94 d is not a low-order alias, and epochs 13/14 are at EB phases 0.60 and 0.08 with equal depth.
- No KOI/TCE/EB/confirmed/TESS-EB/TOI/CTOI matches period and epoch within 3 arcmin or in the field (closest phase twin TCE 4663658-01 at 32.9614 d is 2.17 d off); target absent from all catalogues, Bienias & Szabo 2025 and SIMBAD (only KOI-3829 at 18.3 arcsec and KOI-1657 at 24.6 arcsec within 40 arcsec).
- Odd/even and secondary null; LS variability 68 ppm at 1.18 d.
- Duration 1.95 h needs b = 0.93 on the 0.87 Rsun host (T14(b=0) = 5.16 h); Rp 2.3-4.4 Re depending on the depth believed. Not a kill by itself.

## Why the pipeline passed it (audit; no thresholds changed)
max_single_event_frac 0.39 (< 0.5), gap_frac_1d 0.18 (< 0.5), event reduced chi2 2.83 (< 3.0), SNR_pass2 11.6: all driven by 8 cadences out of 64. The rules record per-event gap proximity but do not exclude post-gap events from the stack. Any new cut (re-stack with events < 0.5 d after a segment start removed, require the masked BLS to still peak at the same period) must be logged in `prereg_poolA.json` `revision_log` and labelled post-hoc.

## What would settle it
A full 1-40 d pre-registered search plus inverted null on the post-gap-masked light curve, to record the non-detection formally; nothing else is worth spending on.