# kic-005866176: red-team review (2026-10-17)

**Verdict: kill.** Full notes: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005866176/redteam.md`; verdict JSON `redteam_verdict.json`; scripts, JSON, plots, logs and `MANIFEST.sha256` under `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005866176/redteam/`. Nothing from `pipeline/` was imported. Language: first-screen signal only.

## The kill: two post-gap PSF-fit excursions make the detection

Per-event local depths at the pipeline ephemeris (`redteam/lc_checks.json`): BKJD 231.44 (Q2, ch 17, 0.045 d after a gap) **3251 +- 607 ppm**; BKJD 937.83 (Q10, ch 17, 0.027 d after a gap) **3568 +- 229 ppm**; the other 29 events: weighted mean 275 ppm, chi2 37.0/28. With the two events the chi2 is 258/30 (the pipeline's own vet.json already shows event_weighted_depth 250 ppm vs box depth 481 and chi2 67/30).

Raw-flux anatomy (`redteam/event_rawflux.log`, 3 h bins): at 231.48 the target PSF `FLUX` is -7830 ppm recovering to +2700 over 2 d while the target `SAP_FLUX` is **+18,590 ppm** decaying to +1,700; TPF column 0 (the only column the target's wing reaches, and the whole KBONUS SAP aperture) drifts -3,261 -> -21,954 ppm while the host's columns drift -450 -> -1,430; `POS_CORR1` drifts -0.023 -> -0.042 px. At 937.85: PSF -1,523 -> +2,500 ppm step, SAP +11,152 -> 0, column 0 +22,938 -> +10,000, POS_CORR1 0.040 -> 0.030. Quarter box depths: Q2 PSF +1459 +- 135 vs SAP **-1067 +- 147** ppm; Q10 PSF +1040 +- 134 vs SAP **-992 +- 147**. A loss of light cannot raise the SAP flux by 1-2 % while lowering the PSF flux: these are post-downlink focus/pointing settling events seen through a PSF model for a star centred **outside the stamp** (Gaia position at x = -0.22 to -0.39 px in all 11 quarters; host at x = 3.0; PSFFRAC 0.52-0.69).

Post-gap ramps are generic here: first-0.15 d means after single gaps reach -3,909, -2,438, -1,962, -678, -609 ppm; of 1,500 random ephemerides at the same period 4.7 % hit >= 2 post-gap windows and that subset has box-depth scatter 122 ppm with tails to ~430 ppm (52 ppm otherwise). Survey-wide, 349/438 real first-screen passes sit in the 24-40 d monthly-gap class.

## What is left after excising 0.3 d after every gap (2.2 % of data)

| statistic | all data | excised |
|---|---|---|
| 1.95 h box, PSF (400 random-phase boxes) | 514 +- 50 ppm | **272 +- 48 ppm (SNR 5.7 < 7.1 floor)** |
| same, SAP | -8 +- 56 | 276 +- 55 |
| own BLS 20-40 d: SDE at 29.43 d / top | 6.65 (top, 29.4318 d) | **3.59** / top 30.89 d SDE 5.50 |
| inverted-flux null, top SDE | 5.69 | 5.94 |
| odd/even | 363/640 (3.1 sigma) | 362/185 (1.9 sigma) |
| per season PSF S0/S1/S2/S3 (ch 17/57/65/25) | 907+-76 / 298+-86 / 335+-180 / 326+-99 | 168+-74 / 299+-94 / 334+-181 / 325+-92 |

A ~270-330 ppm feature at the fixed ephemeris survives in PSF, SAP and PSF_FLAT (360 +- 52) flux and 0/1500 random phases reach it, so a small dip is not excluded; but the ephemeris was chosen on artefact-dominated data, it is not blind-recoverable, and it fails the registered rules.

## Pixel level and host

Season-stacked event-differenced deficit maps on the 11 host TPFs (`redteam/pixel_checks.json`): SNR at the pixel nearest the target 2.7 (S0, excised) / 2.5 (S1) / 1.2 (S2) / 0.75 (S3); max |SNR| anywhere 2.7; in S1 a PRF template at the host fits as well as at the target; grid optima scatter. Localisation is impossible for an off-stamp star. The host KIC 5866192's own KBONUS light curve at the ephemeris: 19 +- 35 ppm (PSF), 14 +- 25 ppm (SAP), < 2.1 e/s at 2 sigma: the host is not the source either.

## Known-object and other checks

No KOI/TCE/Kepler-EB/PS for KIC 5866176 or 5866192; within 3 arcmin only unrelated objects (KOI-506.01 at 1.58 d, EB 5866138 at 342 d, long-period TCEs); 176 harmonic period matches field-wide with only three far-off, incoherent phase coincidences; five cumulative KOIs at 29.38-29.49 d misaligned by 7-21 d; 0 TOIs in a 3 arcmin box; SIMBAD within 20 arcsec shows only the host; not in the Bienias & Szabo 2025 list as transcribed by the sibling review; six other pool stars peak near 29.43 d with unrelated phases. Lomb-Scargle: no rotation alias (power near P 0.004). Quality flags 8.7 % in transit vs 8.7 % overall. Duration 1.95-3.0 h needs b = 0.93-0.82 for R* 0.95 Rsun. Implied radius 2.3 Re (481 ppm) or 1.7 Re (272 ppm), with uncalibrated season-dependent dilution.

## Recommendation

Retire the signal; log as a post-gap PSF-fit false positive (edge-of-stamp background star, channel 17). Register (post-hoc, v1.8) a post-gap 0.3 d excision before the stage-2 depth, a PSF/SAP sign-agreement test per quarter, and an off-stamp flag; these will matter for the 24-40 d class.