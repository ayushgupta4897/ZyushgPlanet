# kic-011285493 — red-team review (2026-10-17)

Candidate language only. This is a transit-like **signal** in the KBONUS-BKG PSF light curve of
KIC 11285493 (Gaia DR3 2131988887763678720; G = 16.84, BP-RP = 0.94, GSP-Phot Teff 5286 K,
parallax 0.350 +- 0.041 mas = 2.86 kpc, RUWE 1.00; originating TPF KIC 11285496, Kp 13.66, 9.2 arcsec away).
Nothing here is a planet. Run ids: first screen `runs/poolA_main/triage_all.csv` (prereg v1.7),
stage-2 `runs/poolA_main/vet/kic-011285493/vet.json`, pixel `runs/poolA_main/pixel/kic-011285493/pixel.json`,
independent re-derivation `runs/poolA_main/adversarial/kic-011285493/rederive/results.json`.
KBONUS FITS SHA-256 `2044c692119d311a…e89` (full hashes in `redteam/MANIFEST.sha256`).
Known-object catalogues: `data/known/manifest.json`, fetched 2026-10-17 15:46 UTC.
My scripts, logs and outputs: `runs/poolA_main/adversarial/kic-011285493/redteam/` (scratch copy in
`/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic011285493/`).

**Verdict: KILL.** The 10.79-d, ~800-ppm dip is an *ephemeris-match false positive*: it is the primary
eclipse of the bright detached eclipsing binary **KIC 11285625 = KOI-7431.01 (Kp 10.14)**, 5.63 arcmin away
on the **same CCD channel in every season**, whose Kepler-EB-catalogue period is 10.790444 d and whose
primary minimum coincides with our t0 to 0.6 h. The eclipse is imprinted on KBONUS light curves channel-wide
by the KBONUS extraction (it appears in the KBONUS PSF flux of KIC 11285496 and of a second pool star,
kic-011285117, but is absent from the official Kepler SAP/PDCSAP of KIC 11285496 and from the target's own
raw pixels in the season where the KBONUS dip is strongest). The pipeline missed it because the
ephemeris match is confined to a 3-arcmin cone.

## 0. What the signal is (numbers from the evidence pack)
- Pipeline: P = 10.791014 d, t0 = 131.52143 BKJD, 4.5-h box depth 796 +- 92 ppm (box_noise 91 ppm, empirical
  SNR 8.75), 96 events, odd/even 1.35 sigma, secondary 160 ppm (1.76 sigma), TLS SDE 11.1 / SNR 8.0,
  Rp ~ 2.5-2.65 Re for R* = 0.86 Rsun, duration_ratio 1.18, event chi2/dof 2.05, `known_matches: []`.
- Re-derivation (blind): P = 10.790158 d, t0 = 843.73888 BKJD, 801 +- 116 ppm at 4.5 h (SNR 6.9), BLS SDE 19.3,
  trapezoid 891 ppm / T14 5.9 h / T23/T14 0.41, 90 events, per-event chi2/dof 3.0, Y4 depth 254 +- 199 ppm,
  pixel localisation inconclusive.
- First-screen inverted null (`triage_all.csv`): real SDE 13.2 / SNR_emp 7.86 passes; inverted top peak
  SDE 9.5 at 0.577 d fails on snr_emp and duty cycle (reported together, as required).

## 1. The ephemeris match (the kill)
`redteam/ephem_match.py`, `redteam/ephem_match.log`; catalogue rows quoted from `data/known/`.

| object | source | period (d) | epoch (BKJD) | depth | notes |
|---|---|---|---|---|---|
| KIC 11285625 | Kirk+2016 EB catalogue (`J_AJ_151_68_catalog.csv`) | **10.790444 +- 0.000032** | BJD0 2454975.33451 = **142.33451** | morph 0.23 (detached) | Kp 10.143, Teff 6882 K |
| K07431.01 | KOI cumulative | 5.395256 (= P/2) | 131.542553 | 104,242 ppm, T14 7.84 h | FALSE POSITIVE, `koi_fpflag_ss = 1` |
| KIC 11285625.01 | DR25 TCE | 5.39528 | 131.540 | 98,050 ppm | MES 399 |
| this signal (pipeline) | vet.json | 10.791014 | 131.52143 | 796 ppm | dP/P = +5.3e-5; t0 - KOI t0 = -0.021 d = **-0.5 h** |
| this signal (re-derivation) | results.json | 10.790158 | 843.73888 | 801 ppm | dP/P = -2.6e-5; (t0 - BJD0)/P = 65.0024 -> **+0.6 h** after the EB primary |

- Separation 5.63 arcmin (337 arcsec): outside the pipeline's 3-arcmin cone, hence `known_matches: []`.
- Chance of coincidence: among the ~22 DR25 TCEs (9 stars) within 12 arcmin, the probability that one lands
  within |dP/P| < 3e-5 of a harmonic of our period *and* within +-0.6 h in epoch on a 10.8-d period is
  of order 1e-6 per object, ~1e-5 in total. This is not a coincidence; by the Coughlin et al. (2014) standard
  (period and epoch match with a brighter, deeper object on the same channel) it is an ephemeris-match FP
  regardless of mechanism.
- Physical association excluded: the EB's Gaia DR3 parallax is 1.865 +- 0.017 mas (536 pc; source 2131979370116194432,
  G 10.08) vs 2.86 kpc for the target.

### 1a. Same CCD channel every season (`redteam/geometry.log`, `redteam/eb_tpf.log`, `redteam/tpf_headers.log`)
| season | quarters | our TPF (KIC 11285496) | EB TPF (KIC 11285625) | target CCD col/row | EB CCD col/row |
|---|---|---|---|---|---|
| 0 | Q2, 6, 10, 14 | ch 79, mod 23.3 | ch 79, mod 23.3 | 643 / 943 | 726 / 922 |
| 1 | Q3, 7, 11, 15 | ch 51, mod 15.3 | ch 51, mod 15.3 | 643 / 939 | 726 / 918 |
| 2 | Q4, 8, 12, 16 | ch 7, mod 3.3 | ch 7, mod 3.3 | 640 / 938 | 722 / 917 |
| 3 | Q1, 5, 9, 13, 17 | ch 35, mod 11.3 | ch 35, mod 11.3 | 639 / 943 | 722 / 922 |

The EB sits 83 columns and 20 rows from the target on the same output in all four seasons (so this is not the
classic same-column anomaly; see 1d for the actual path). KBONUS has no season-2 quarters for this target.

### 1b. The EB's own eclipse profile matches the dip, including its suppressed-secondary signature
`redteam/eb_lc.py`, `redteam/eb_template.py` (official Kepler SAP of KIC 11285625, 15 quarters, folded at the
EB-catalogue ephemeris; profile in `redteam/eb_profile.csv`): primary depth 0.1251, secondary 0.0869 at
phase 0.4995 (circular), ratio 0.695, primary width (>10 % of depth) 6.8 h, rounded-V shape.

Fitting that profile as a template to the target's masked-detrended PSF flux (`redteam/eb_template.log`):
- primary-template amplitude -> **952 +- 136 ppm (random-phase error), SNR 7.0**; the EB's primary shape fits
  the dip as well as any trapezoid (my own trapezoid grid, `redteam/lc_battery.log`, cannot distinguish
  T23/T14 = 0 … 1, delta chi2 < 2.2; T14 = 5.25-7 h).
- secondary-template amplitude -> 247 ppm (2.6 sigma); ratio of secondary to primary *scale* = 0.37 where a
  perfect photometric copy would give 1.0 (3.1 sigma low). The same suppressed ratio (~0.15) is seen in the two
  other contaminated KBONUS sources below, so it is a property of the contamination path, not evidence for a
  transit. (The pipeline's 160 ppm and the re-derivation's 109 ppm "secondary" are this residual EB secondary.)
- **amplitude depends on the CCD**: 825 ppm (season 0, ch 79), 1634 ppm (season 1, ch 51), 391 ppm
  (season 3, ch 35). Box depths by season (T14 = 5.9 h box, white errors): 525 +- 162, 1151 +- 146, 368 +- 134 ppm;
  season 1 vs season 3 differ by 3.9 sigma; per quarter Q7 = 2265 +- 305, Q11 = 1563 +- 278 vs Q9 = -194 +- 242 ppm.
  An astrophysical transit cannot know which CCD the star is on. This is the "over-dispersed per-event depths"
  and the "quartile 2 vs 3" tension already noted by the pipeline and the re-derivation, now explained.

### 1c. The same clock is in other KBONUS light curves on the channel (`redteam/neighbours_fold.log`)
EB-primary template fitted to other KBONUS-BKG PSF light curves at the EB ephemeris:

| KBONUS source | Kp | host TPF (channel s0/s1) | PSF flux (e/s) | primary dip | = e/s | secondary | per season (ppm) |
|---|---|---|---|---|---|---|---|
| kic-011285493 (this signal) | 16.85 | KIC 11285496 (79/51) | 1,584 | **952 +- 136 ppm** (7.0 sigma) | 1.5 | 247 ppm | 825 / 1634 / – / 391 |
| kic-011285496 (the Kp 13.6 neighbour) | 13.66 | itself (79/51) | 41,422 | **65 +- 16 ppm** (4.2 sigma) | 2.7 | 9 ppm | 59 / 99 / – / 38 |
| kic-011285117 (pool star 14.6' from the EB) | 15.96 | KIC 11285120 (79/51) | 3,394 | **582 +- 94 ppm** (6.2 sigma) | 2.0 | 91 ppm | 1402 / 126 / – / 452 |
| kic-011285661 (2.1' from the EB) | – | KIC 11285663 | 18,182 | 108 +- 60 ppm (1.8 sigma) | 2.0 | 203 ppm | 235 / -14 / – / 183 |
| kic-011285855 (6.0' from the EB) | 16.35 | KIC 11285862 (79/51) | 4,083 | 49 +- 51 ppm | 0.2 | 13 ppm | -26 / 119 / – / 54 |

The contamination is *additive*, 1.5-2.7 e/s in three unrelated sources of very different brightness
(65 ppm on a 41,000 e/s star, ~950 ppm on a 1,600 e/s star), with the same suppressed secondary and a
season-dependent amplitude. A 10.4 % eclipse of a 1.1e6 e/s star removes ~1.4e5 e/s; 2 e/s is 1.4e-5 of that.

### 1d. Where the leak happens: KBONUS processing, not the pixels (`redteam/official_sap.log`, `redteam/pixel_events.log`)
- Official Kepler photometry of KIC 11285496 (same star as row 2 above) at the EB ephemeris: SAP
  **-14 +- 13 ppm** (-0.6 e/s), PDCSAP **-7 +- 10 ppm**; per season 4 / -25 / -169 / -14 ppm (SAP). The KBONUS PSF
  flux of the same star shows +65 +- 16 ppm = +2.7 e/s. The dip is therefore introduced between the calibrated
  pixels and the KBONUS light curve, i.e. by the KBONUS/psfmachine extraction (its channel-wide, time-dependent
  PRF/perturbation model is fitted on bright stars of the channel; a Kp 10.1 eclipsing binary losing 12 % of
  its light every 10.79 d is a natural way to imprint that clock on every source on the channel). I did not
  reverse-engineer psfmachine; this is an inference from the SAP-vs-KBONUS contrast and the channel-wide
  occurrence.
- Target's own raw TPF pixels (event-differenced, per-event errors, quarters with >= 4 events): deficit in
  r < 1.5 px around the Gaia position **0.48 +- 0.20 e/s** overall (chi2/dof 2.6 across quarters), far pixels
  0.26 +- 0.12 e/s, whole stamp 2.8 +- 1.2 e/s. In season 1, where the KBONUS flux says 1634 ppm (= 1.8 e/s
  expected inside r < 1.5 px for a target-hosted dip), the target pixels give **-0.06 +- 0.45 e/s** (Q3, 7, 11, 15),
  4 sigma short. The light is not missing from the target. (The pipeline's pixel test, median peak SNR 2.5,
  3 of 8 quarters closer to the target than any other source, was rightly inconclusive; the SAP/PSF "same 1.3 e/s"
  result of the re-derivation is the same additive leak seen through two apertures, as that author suspected.)

## 2. Other tests (all consistent with the artefact; none rescues the signal)
- Odd/even: 586 vs 778 ppm (my T14 box), 0.6-1.35 sigma elsewhere. Expected: the contamination copies the
  same primary eclipse every 10.79 d.
- Secondary: see 1b. Nothing deeper than 298 ppm at any other phase (look-elsewhere over ~70 phases).
- Lomb-Scargle (both polarities, out-of-eclipse, detrended): 0.7155 d, 88 ppm half-amplitude; nothing at P, P/2, 2P, P/3.
- Gaps and quarter edges: 2 of 94 events within 0.5 d of a gap (depth 2226 ppm, driven by epoch -61 at BKJD 185.54,
  1.8 d after the Q2 safe mode, as the re-derivation found); far-from-gap depth 625 ppm. Dropping the three
  strongest events leaves 576 ppm. The signal is not a gap artefact, it is an EB artefact.
- Is it already known? KIC 11285493 is in no KOI, DR25 TCE, Kepler EB, confirmed-planet (Kepler field) table
  and not in SIMBAD (only 2MASS J19003107+4901037 = KIC 11285496 within 40 arcsec). The EB itself is the known
  object. Bienias & Szabo 2025's seven candidates are hot Jupiters at 2.5-4.8 d (docs/02, docs/04), not this.
- Host parameters (Gaia DR3 + astrophysical_parameters, `redteam/gaia_simbad.log`): GSP-Phot gives
  R = 0.868 Rsun, M_FLAME 0.87 Msun, L 0.53 Lsun but at distance_gspphot = 1948 pc and A_G = 0.004, both at odds
  with the parallax (2858 pc, 2.2 sigma) and with ~0.3-0.5 mag of extinction expected at 2.9 kpc in this field.
  At the parallax distance M_G = 4.56 implies L ~ 1.1-1.8 Lsun and R ~ 1.3-1.6 Rsun, so the pipeline's
  `rstar_used = 0.86` would have under-reported Rp (2.65 -> ~4-5 Re) and t14_max (3.8 -> 5.5-7 h). Moot here,
  but the target list's Rstar column should not be trusted for sub-giant-ish stars at kpc distances.

## 3. What would have caught it, and what to do (recommendations only; I changed nothing under `pipeline/`)
1. Extend `vet_batch`'s ephemeris match from a 3-arcmin cone to **the whole channel/module** (Coughlin et al.
   2014 matched over the entire field), at harmonics 1/3 … 3, with |dP/P| < 1e-3 and epoch within one duration;
   use the Kepler EB catalogue period (not the KOI half-period) for EBs with `koi_fpflag_ss = 1`.
2. For every survivor, fold the KBONUS light curves of the *originating TPF star* and of two or three other
   KBONUS sources on the same channel at the survivor's ephemeris: a channel-wide artefact shows up at once.
3. Compare the KBONUS flux of the TPF host with its official SAP/PDCSAP at the ephemeris, as done in 1d.
4. **kic-011285117** (stage-1 survivor at 39.89 d) demonstrably carries this EB's clock at 6.2 sigma in its
   KBONUS flux; its own 39.9-d signal is unrelated in period but the light curve is contaminated and must be
   re-examined with that in mind.
5. Treat season-dependent depth (here 3.9 sigma between seasons 1 and 3) as a hard flag: it is the fingerprint
   of anything that depends on the CCD rather than on the star.

## 4. Bottom line
The signal is real in the KBONUS PSF light curve and was correctly recovered blind, but it is not light missing
from KIC 11285493. It is the 12 %-deep primary eclipse of KOI-7431.01 / KIC 11285625 (P = 10.790444 d), a Kp 10.1
eclipsing binary 5.6 arcmin away on the same CCD channel, leaking at the ~2 e/s level into KBONUS light curves
across the channel — present in the KBONUS flux of KIC 11285496 and kic-011285117, absent from the official
Kepler photometry of KIC 11285496 and from the target's own pixels. Period match 2.6e-5, epoch match 0.6 h,
same channel, season-dependent amplitude, residual EB secondary. **Kill; do not carry forward.**
