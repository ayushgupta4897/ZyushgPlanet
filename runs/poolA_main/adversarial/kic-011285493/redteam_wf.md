# kic-011285493 — red-team review (2026-10-17)

Candidate language only. This is a transit-like **signal** in the KBONUS-BKG PSF light curve of KIC 11285493 (Gaia DR3 2131988887763678720; G = 16.84, BP-RP = 0.94, GSP-Phot Teff 5286 K, parallax 0.350 +- 0.041 mas = 2.86 kpc, RUWE 1.00; originating TPF KIC 11285496, Kp 13.66, 9.2 arcsec away). Nothing here is a planet. Run ids: first screen `runs/poolA_main/triage_all.csv` (prereg v1.7), stage-2 `runs/poolA_main/vet/kic-011285493/vet.json`, pixel `runs/poolA_main/pixel/kic-011285493/pixel.json`, independent re-derivation `runs/poolA_main/adversarial/kic-011285493/rederive/results.json`. KBONUS FITS SHA-256 `2044c692119d311a…e89` (full hashes in `redteam/MANIFEST.sha256`). Known-object catalogues: `data/known/manifest.json`, fetched 2026-10-17 15:46 UTC. My scripts, logs and outputs: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-011285493/redteam/`; notes: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-011285493/redteam.md`.

**Verdict: KILL.** The 10.79-d, ~800-ppm dip is an *ephemeris-match false positive*: it is the primary eclipse of the bright detached eclipsing binary **KIC 11285625 = KOI-7431.01 (Kp 10.14)**, 5.63 arcmin away on the **same CCD channel in every season**, whose Kepler-EB-catalogue period is 10.790444 d and whose primary minimum coincides with our t0 to 0.6 h. The eclipse is imprinted on KBONUS light curves channel-wide by the KBONUS extraction (it appears in the KBONUS PSF flux of KIC 11285496 and of a second pool star, kic-011285117, but is absent from the official Kepler SAP/PDCSAP of KIC 11285496 and from the target's own raw pixels in the season where the KBONUS dip is strongest). The pipeline missed it because the ephemeris match is confined to a 3-arcmin cone.

## 0. What the signal is (evidence pack)
- Pipeline: P = 10.791014 d, t0 = 131.52143 BKJD, 4.5-h box depth 796 +- 92 ppm (box_noise 91 ppm, empirical SNR 8.75), 96 events, odd/even 1.35 sigma, secondary 160 ppm (1.76 sigma), TLS SDE 11.1 / SNR 8.0, Rp ~ 2.5-2.65 Re for R* = 0.86 Rsun, duration_ratio 1.18, event chi2/dof 2.05, `known_matches: []`.
- Re-derivation (blind): P = 10.790158 d, t0 = 843.73888 BKJD, 801 +- 116 ppm at 4.5 h (SNR 6.9), BLS SDE 19.3, trapezoid 891 ppm / T14 5.9 h / T23/T14 0.41, 90 events, per-event chi2/dof 3.0, Y4 254 +- 199 ppm, pixel localisation inconclusive.
- First-screen inverted null: real SDE 13.2 / SNR_emp 7.86 passes; inverted top peak SDE 9.5 at 0.577 d fails on snr_emp and duty cycle.

## 1. The ephemeris match (the kill)
| object | source | period (d) | epoch (BKJD) | depth | notes |
|---|---|---|---|---|---|
| KIC 11285625 | Kirk+2016 EB catalogue | **10.790444 +- 0.000032** | BJD0 = **142.33451** | morph 0.23 | Kp 10.143, Teff 6882 K |
| K07431.01 | KOI cumulative | 5.395256 (= P/2) | 131.542553 | 104,242 ppm, 7.84 h | FALSE POSITIVE, fpflag_ss = 1 |
| KIC 11285625.01 | DR25 TCE | 5.39528 | 131.540 | 98,050 ppm | MES 399 |
| this signal (pipeline) | vet.json | 10.791014 | 131.52143 | 796 ppm | dP/P = +5.3e-5; t0 - KOI t0 = **-0.5 h** |
| this signal (re-derivation) | results.json | 10.790158 | 843.73888 | 801 ppm | dP/P = -2.6e-5; (t0 - BJD0)/P = 65.0024 -> **+0.6 h** |

- Separation 5.63 arcmin: outside the pipeline's 3-arcmin cone. Chance of such a period+epoch match among the ~22 TCEs within 12 arcmin ~1e-5. By the Coughlin et al. (2014) standard this is an ephemeris-match FP regardless of mechanism.
- Physical association excluded: EB Gaia DR3 parallax 1.865 +- 0.017 mas (536 pc) vs 2.86 kpc for the target.

### 1a. Same CCD channel every season (`redteam/geometry.log`, `eb_tpf.log`)
| season | quarters | our TPF | EB TPF | target col/row | EB col/row |
|---|---|---|---|---|---|
| 0 | Q2, 6, 10, 14 | ch 79, mod 23.3 | ch 79, mod 23.3 | 643 / 943 | 726 / 922 |
| 1 | Q3, 7, 11, 15 | ch 51, mod 15.3 | ch 51, mod 15.3 | 643 / 939 | 726 / 918 |
| 2 | Q4, 8, 12, 16 | ch 7, mod 3.3 | ch 7, mod 3.3 | 640 / 938 | 722 / 917 |
| 3 | Q1, 5, 9, 13, 17 | ch 35, mod 11.3 | ch 35, mod 11.3 | 639 / 943 | 722 / 922 |

### 1b. The EB's own eclipse profile matches the dip, including its suppressed secondary
Official Kepler SAP of KIC 11285625 folded at the EB ephemeris: primary 0.1251, secondary 0.0869 at phase 0.4995, ratio 0.695, width 6.8 h, rounded V. Fitting that profile to the target's masked-detrended PSF flux (`redteam/eb_template.log`): primary -> **952 +- 136 ppm (random-phase error), SNR 7.0**; secondary -> 247 ppm (2.6 sigma), secondary/primary scale 0.37 (same suppressed ratio, ~0.15, in the other contaminated sources below; the pipeline's 160 ppm and the re-derivation's 109 ppm 'secondary' are this residual EB secondary). **Amplitude depends on the CCD**: 825 ppm (s0, ch 79), 1634 ppm (s1, ch 51), 391 ppm (s3, ch 35). T14-box depths by season: 525 +- 162, 1151 +- 146, 368 +- 134 ppm (s1 vs s3: 3.9 sigma); per quarter Q7 2265 +- 305, Q11 1563 +- 278 vs Q9 -194 +- 242 ppm. An astrophysical transit cannot know which CCD the star is on.

### 1c. The same clock is in other KBONUS light curves on the channel (`redteam/neighbours_fold.log`)
| KBONUS source | Kp | host TPF (ch s0/s1) | PSF flux (e/s) | primary dip | e/s | secondary | per season (ppm) |
|---|---|---|---|---|---|---|---|
| kic-011285493 (this) | 16.85 | KIC 11285496 (79/51) | 1,584 | **952 +- 136 ppm** (7.0 sigma) | 1.5 | 247 | 825 / 1634 / – / 391 |
| kic-011285496 (TPF host) | 13.66 | itself (79/51) | 41,422 | **65 +- 16 ppm** (4.2 sigma) | 2.7 | 9 | 59 / 99 / – / 38 |
| kic-011285117 (14.6' from EB) | 15.96 | KIC 11285120 (79/51) | 3,394 | **582 +- 94 ppm** (6.2 sigma) | 2.0 | 91 | 1402 / 126 / – / 452 |
| kic-011285661 (2.1' from EB) | – | KIC 11285663 | 18,182 | 108 +- 60 ppm | 2.0 | 203 | 235 / -14 / – / 183 |
| kic-011285855 (6.0' from EB) | 16.35 | KIC 11285862 (79/51) | 4,083 | 49 +- 51 ppm | 0.2 | 13 | -26 / 119 / – / 54 |

Additive contamination of 1.5-2.7 e/s in unrelated sources of very different brightness; ~1.4e-5 of the ~1.4e5 e/s the EB loses in eclipse.

### 1d. Where the leak happens: KBONUS processing, not the pixels
- Official Kepler photometry of KIC 11285496 at the EB ephemeris: SAP **-14 +- 13 ppm**, PDCSAP **-7 +- 10 ppm** (`redteam/official_sap.log`); the KBONUS PSF flux of the same star shows +65 +- 16 ppm = +2.7 e/s. The dip is introduced between the calibrated pixels and the KBONUS light curve (inference: the channel-wide time-dependent PRF/perturbation model of psfmachine, fitted on bright channel stars including a Kp 10.1 EB).
- Target's own raw TPF pixels (`redteam/pixel_events.log`, quarters with >= 4 events): r < 1.5 px deficit **0.48 +- 0.20 e/s** overall; in season 1 (KBONUS says 1634 ppm = 1.8 e/s expected) **-0.06 +- 0.45 e/s**, 4 sigma short. The light is not missing from the target. The pipeline's pixel test (median peak SNR 2.5, 3/8 quarters) was rightly inconclusive; the re-derivation's 'same 1.3 e/s in SAP and PSF' is the same additive leak seen through two apertures.

## 2. Other tests (consistent with the artefact)
- Odd/even 586 vs 778 ppm (0.6-1.35 sigma elsewhere) — a copied eclipse repeats identically.
- Secondary: see 1b; nothing deeper than 298 ppm at any other phase.
- Lomb-Scargle: 0.7155 d, 88 ppm; nothing at P, P/2, 2P, P/3.
- Gaps: 2 of 94 events within 0.5 d of a gap (epoch -61 after the Q2 safe mode drives them); far-from-gap depth 625 ppm; without the 3 strongest events 576 ppm.
- Known? KIC 11285493 is in no KOI/TCE/EB/confirmed table nor SIMBAD; the EB is the known object. Bienias & Szabo 2025 (hot Jupiters 2.5-4.8 d) and Moquin BOI-1.01 (3.37 d) are unrelated.
- Host parameters: GSP-Phot R = 0.868 Rsun at distance_gspphot 1948 pc and A_G 0.004, both inconsistent with the parallax (2858 pc, 2.2 sigma) and expected extinction; at the parallax distance M_G = 4.56 implies R ~ 1.3-1.6 Rsun, so `rstar_used = 0.86` would have under-reported Rp (~4-5 Re) and t14_max (5.5-7 h). Moot here.

## 3. Recommendations (nothing changed under `pipeline/`)
1. Extend the ephemeris match from a 3-arcmin cone to the whole channel/module, harmonics 1/3..3, |dP/P| < 1e-3, epoch within one duration; use the Kepler EB-catalogue period for koi_fpflag_ss = 1 objects.
2. For every survivor, fold the KBONUS light curves of the TPF host and of other same-channel KBONUS sources at the survivor's ephemeris.
3. Compare the TPF host's KBONUS flux with its official SAP/PDCSAP at the ephemeris.
4. kic-011285117 (stage-1 survivor at 39.89 d) carries this EB's clock at 6.2 sigma; re-examine it with that in mind.
5. Treat season-dependent depth as a hard flag.

## 4. Bottom line
Real in the KBONUS PSF light curve, correctly recovered blind, but not light missing from KIC 11285493: it is the 12 %-deep primary eclipse of KOI-7431.01 / KIC 11285625 (P = 10.790444 d), a Kp 10.1 eclipsing binary 5.6 arcmin away on the same CCD channel, leaking at ~2 e/s into KBONUS light curves across the channel. Period match 2.6e-5, epoch match 0.6 h, same channel, season-dependent amplitude, residual EB secondary. **Kill; do not carry forward.**