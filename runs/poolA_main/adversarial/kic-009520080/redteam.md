# Red-team review: kic-009520080 (KIC 9520080, Gaia DR3 2127366472162960768)

Date 2026-10-17. Role: try to kill the 2.377 d signal. Inputs read: `runs/poolA_main/vet/kic-009520080/{vet.json,events.csv,fold_events.png}`, `runs/poolA_main/pixel/kic-009520080/{pixel.json,diffimages.png}`, `runs/poolA_main/pixel_v15/` and `pixel_v16/` summaries, the independent re-derivation under `adversarial/kic-009520080/rederive/`, `pipeline/prereg_poolA.json`, `data/known/` (manifest fetched 2026-10-17T15:46Z). Own tests (no `pipeline/` imports) are in `adversarial/kic-009520080/redteam_out/` (scripts + JSON/log outputs); scratch copies under the session scratchpad `redteam-kic009520080/`.

Language: this is a stage-2 **signal**. Nothing below calls it a planet.

## Verdict: SURVIVES (with two open items that must travel with it)

I could not find a non-planet explanation or an instrumental origin. The signal is on the target's pixel, not on either bright neighbour, in all 17 quarters and all four CCD channels; no catalogued object shares its ephemeris; the inverted null is clean; odd/even, timing and depth are consistent across seasons. The two open items are (1) a marginal (2.2-2.8 sigma, estimator-dependent) brightening at phase 0.5 that no planet produces and that makes the signal a formal stage-2 FAIL under the post-hoc v1.5 rule (c), and (2) the generic unresolved-blend risk (a G > 21 background EB inside the ~4 arcsec pixel, or a bound companion EB) that Kepler pixels cannot address and that only TRICERATOPS-style statistics can bound.

## 0. Governance note (must be stated wherever this signal is reported)
`vet.json`: `stage2_pass_v14 = true`, `stage2_pass = false` under `poolA_kbonus_background_dwarfs_v1.5`, `stage2_fails = ["secondary_2sided"]`. `prereg_poolA.json` revision_log v1.5(c) records that the two-sided secondary rule was added POST-HOC *because of this star's* -4 sigma phase-0.5 brightening. So by the currently registered rules this signal does not pass stage 2; it passes the pre-v1.5 rules. Both facts have to be reported together.

## 1. Evidence-pack facts (pipeline numbers, `vet.json` unless stated)
P = 2.3770144 d, t0 = BKJD 132.15454, box duration 1.125 h, depth pass-2 (masked re-detrend) 1569 +- 62 ppm, empirical SNR 28.3, 468 events with local depths, TLS SDE 22.2 / SNR 24.7 (narrow window), odd/even 1.13 sigma, secondary -222 ppm (-4.0 sigma, i.e. brightening), LS 0.286 d at 92 ppm (not harmonically related; 2.377/0.286 = 8.3), Rp(box) 3.03 Re for Rstar 0.70, known_matches []. Pixel v1.6 chi2 verdict `on-target` (sum dchi2 null-target 917, best-alternative-target 335; `pixel_v16/kic-009520080/pixel.json`); centroid v1.5 verdict `inconclusive` because the target lies at/over the stamp edge in 12/17 quarters.

## 2. Tests and results

### 2.1 Blended / background eclipsing binary
**TPF host and bright neighbours (independent per-pixel test, `redteam_out/pixel_depths.py` -> `pixel_depths.json/.log`).** On the 17 cached TPFs of KIC 9520088 (Kp 14.78, Teff 6150 K, R 1.11 Rsun from the TPF headers; this IS Gaia DR3 2127366575242177536, G 14.77, 10.4 arcsec from the target), I built per-pixel light curves, per-quarter normalised, transit-masked 1-d running median, box depth at the ephemeris:
- Target's pixel: deficit present in **every quarter**, SNR 3.2-10.7 per quarter, combined **1.02 +- 0.04 e/s (24.5 sigma)**, i.e. 1000-1840 ppm of that pixel's flux (median ~1400 ppm). Expected in-stamp missing light for 1570 ppm on the target: 1537 e/s (KBONUS `psf_mean_flux`) x 1570e-6 x psffrac 0.61 = 1.5 e/s, spread over the edge pixel and its neighbours (row sums 0.8-1.9 e/s). Quantitatively consistent with the whole dip being on the target.
- Host's brightest pixel: **-0.33 +- 0.51 e/s** combined over the five brightest (host-dominated) pixels, median **-11 ppm** of their flux. A blend from the host would put ~1500 ppm into every host-dominated pixel; observed < ~30 ppm (2 sigma). Excluded by a factor > 50.
- KBONUS PSF light curves of the neighbours (downloaded to scratch, same detrend): **KIC 9520088 (G 14.8): depth at the ephemeris -6 +- 12 ppm** (scatter 327 ppm, 1163 in-transit cadences); **KIC 9520089 (G 13.9, 19 arcsec): -190 +- 210 ppm** (noisy star, 8300 ppm scatter). Neither can supply a 1500 ppm dip to a nearly clean (crowdsap 0.98) target aperture.
- The only other Gaia source near the target is DR3 2127366472160644352, **G 20.88 at 3.46 arcsec** (0.87 px), parallax 0.08 +- 1.2 mas (background). In seasons 1-3 it falls in a different pixel from the target; the deficit follows the target's pixel (mean 1.01 e/s) with only PSF-spill in the G 20.9 pixel (mean 0.24 e/s over 13 quarters where the pixels differ). A source at the G 20.9 position would give the reverse ratio. The pipeline chi2 fits agree (target preferred over this source in 15/17 quarters, summed dchi2 ~ 220). It would need a ~3.8 % eclipse (if 100 % of its light leaked) to mimic the signal; not supported by the pixel data.
- Gaia DR3 cone to 15 arcsec (TAP, 2026-10-17) returns exactly the four sources above the pipeline already used: 17.49, 20.88 (3.5"), 14.77 (10.4"), 19.54 (14.3"). Target flags: RUWE 1.15, ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.012, non_single_star 0, astrometric_excess_noise_sig 4.0 (unremarkable at G 17.5). No resolved companion.
- Secondary eclipse: none; the phase-0.5 box is a *brightening* (see 2.3). Odd/even 1.13 sigma (pipeline), 1.14 (re-derivation), odd 1516 / even 1667 ppm in my box (`lc_tests2.log`). Odd and even transit centres agree: best offsets +1.8 and +0.4 min (difference 1.4 min, well inside the ~5 min resolution), so no e cos w signature of a period-doubled eccentric EB.
- Depth vs aperture: SAP/PSF = 1490/1590 = 0.94 (mine), 0.94 +- 0.05 (re-derivation); crowdsap predicts ~0.98. No sign of a brighter contaminant in the SAP aperture.
- Per-season depth (seasons 0-3 = channels 62/46/22/38): 1500 +- 144, 1518 +- 97, 1552 +- 132, 1789 +- 126 ppm (`lc_tests.log`): same depth on four different detectors (chi2 ~ 3.9/3).

**What is NOT excluded:** a background EB fainter than Gaia's limit (G > ~21) inside the target's pixel (needs an eclipse of >= 4 % at G 21, 10 % at G 22, 24 % at G 23, 60 % at G 24 if fully blended), or a bound companion EB (hierarchical triple). Kepler's 4 arcsec pixels cannot test either; this is the TRICERATOPS / follow-up item.

### 2.2 Stellar variability / rotation aliases
LS on the out-of-transit residuals (both polarities implicitly, since LS is sign-blind): best 0.2863 d at ~125 ppm (pipeline 92 ppm); power at P, P/2 and 2P is 4e-5, 1.5e-5, 5.5e-6 (vs 1e-3 at the peak). No periodic variability at the transit period or its harmonics. The 1.1 h, 1500 ppm box is not a sinusoid alias.

### 2.3 The phase-0.5 brightening (open item 1)
Own detrend (`lc_tests2.py`): full phase scan with a transit-duration box (177 boxes, phase 0.06-0.94):
- mean: **+192 ppm**, scan std 56 ppm, **z = 2.84**, rank 1/177 (most positive box);
- 3-sigma-clipped mean: **+126 ppm, z = 2.24**, rank 2/177;
- median: **+51 ppm, z = 1.32**, rank 20/177;
- per-epoch medians averaged (robust to single-cadence spikes): **+237 +- 69 ppm vs random-phase control +46 +- 69 -> z = 2.7**;
- 13 cadences above +3 sigma in the box (expected 7.9 from the light curve's own heavy positive tail, 0.7 % > 3 sigma vs 0.14 % Gaussian), at **13 different epochs** (Q1-Q14), so it is not one flare; excluding those epochs leaves +110 ppm.
- Per season: +187 +- 143, +218 +- 90, +200 +- 120, +161 +- 136 ppm: consistent across four channels, i.e. not a channel-specific systematic, but individually 1.1-2.4 sigma.
- Re-derivation: -189 +- 58 ppm (their sign = depth), narrow (~1 h), no cos/sin phase curve above 40 ppm.
Reading: a 2.2-2.8 sigma excess at a pre-specified phase, with no physical mechanism in a planet scenario and none obvious in a blend scenario either (an EB's secondary is a dip). With a pre-specified phase the one-sided chance probability is ~0.3-1 %; across the ~11-30 stars that reached this stage, one such fluctuation is not surprising. It does not generate false transits, so it is not a kill reason, but it is unexplained and it makes the signal a stage-2 FAIL under the v1.5 rule. Carry it as an open item; do not drop it, do not over-read it.

### 2.4 Kepler systematics
- Channel/season from the TPF headers (`redteam_out/tpf_headers.json`): season 0 = mod 18.2 / ch 62, season 1 = 14.2 / ch 46, season 2 = 8.2 / ch 22, season 3 = 12.2 / ch 38. Same depth on all four (2.1), so not rolling band, not a single-channel electronic artefact.
- Gaps: `gap_frac_1d = 0.09`, `gap_frac_2d = 0.18`; weighted depth by gap distance 0-0.5 d: 1709, 0.5-1: 1510, 1-2: 1977, 2-5: 1510, > 5 d: 1467 ppm (`events.csv`). No gap-driven excess.
- 2.377 d is not a momentum-dump (~3 d), cadence (116.3 cadences per period, not integer) or quarter-length alias.
- Per-event outliers: `events.csv` reduced chi2 912/467 = 1.95 is inflated by one event, n = 342 (BKJD 945.09, Q10): local depth -5686 +- 491 ppm (z = -14.7). The raw flux shows several isolated +6000 to +8000 ppm single-cadence spikes within +-1 d of that epoch in both PSF and SAP flux (`lc_tests.log`); the formal per-event error (491 ppm for 2 cadences of 1900 ppm scatter) is too small. Dropping that one event gives ~ 695/466 = 1.5, in line with the positive controls (1.95, 1.50). This single event also explains the pipeline's low third time-quartile (weighted 1091 ppm; the quartile *median* is 1708 ppm) and the re-derivation's chi2 discrepancy. Depth is constant in time: per year 1713, 1500, 1323, 1454 (chi2 4.7/3, re-derivation), per season as above.

### 2.5 Ephemeris-match contamination from any known object
Field-wide period match at 1:1, 1:2, 2:1, 1:3, 3:1 against KOI (9564 rows), DR25 TCE (34032), Kepler EB catalogue (Kirk+16, 2876), confirmed planets (2839); `data/known/manifest.json` fetched 2026-10-17T15:46Z:
- Within 3 arcmin: **no** KOI, TCE, EB or confirmed planet. Nearest KOI/TCE of any period is K07174.01 / KIC 9458266 at 502 arcsec (P 1.5849 d, not a harmonic within 0.1 %).
- Closest fractional period deviations anywhere in the field: KOI 1:1 2.2e-4 (KIC 9898447, 4 deg away, 32 ppm); TCE 1:1 1.0e-4 (KIC 9471846); 2:1 5.9e-6 (KIC 12252424 = Kepler-113 c, 748 ppm, 4.8 deg away, 17382 arcsec); EB 1:1 4.0e-4 (KIC 11099351, 6.7 deg). A contaminating ephemeris must match to ~1e-6 fractional to stay coherent over 614 cycles (our period is known to 3e-7 d); nothing does except Kepler-113 c at 2:1, which is 748 ppm on a Kp 13.5 star 4.8 degrees away: its epoch phase relative to ours is 0.48 (so its transits would fall at our phase ~0.5, not 0), it is far too shallow and far too distant to leak, and the odd/even test at 2P (1.13 sigma) shows no alternation. Not a match.
- TPF host KIC 9520088 and KIC 9520089: no KOI, no TCE, not in the EB catalogue; both were Kepler targets, KIC 9520080 never was (`kepler_target_kepids.csv`).

### 2.6 Is it already known?
KOI/TCE/EB/confirmed: no (above). Bienias & Szabo 2025 background candidates (Appendix B list, 7 objects: KIC 4459924, 4543171, 8175131, 8293539, 9040849, 9674230, 11033410): not this star or TPF. TOI/CTOI lists (ExoFOP CSVs downloaded 2026-10-17 by a sibling red-team run): nothing within 5 arcmin. SIMBAD (2026-10-17): within 30 arcsec only 2MASS J19160408+4610016 and UCAC4 681-067998 (the bright neighbours), no variable/EB classification; 22 objects within 3 arcmin, none flagged as EB/planet host. KBONUS team's own candidate list: none exists in `data/kbonus-catalogs/` to check against; not verified beyond the catalogue flags (`psf_fraction_flag` 011111111111111111 for the target).

### 2.7 Duration / density / implied radius
Target: G 17.49, BP-RP 1.25, Teff 4546 K, Rstar 0.70, d 1657 pc, RUWE 1.15 (`targets/poolA_kbonus_bg_dwarfs_ranked.csv` rank 3064). For M = 0.70 Msun, R = 0.70 Rsun: a = 0.0310 AU, a/Rs = 9.5, T14(b = 0) = 1.99 h; the observed 1.125-1.15 h needs **b ~ 0.85**, Rp(box) 3.0 Re, Teq ~ 1040 K. Depth vs box width (0.5 h: 1903, 1.0 h: 1625, 1.125 h: 1590, 2.0 h: 993 ppm) is the smearing of a ~1.1-1.3 h event by 29.4 min cadence, so the true depth is ~1900-2100 ppm and Rp ~ 3.2-3.4 Re at b ~ 0.85; ingress/egress shape is unconstrained (2.3 cadences per transit). A duration consistent with a high-impact-parameter transit of the K dwarf is also consistent with a blend on a denser star, so this is not discriminating; it simply does not contradict the host.

## 3. What would settle it
1. Ground-based transit photometry at the predicted times (depth ~2000 ppm, 1.1 h, G 17.5: a 1-2 m telescope job) to confirm the event on the target at sub-arcsec resolution and separate it from the G 20.9 star at 3.5 arcsec.
2. TRICERATOPS (or equivalent) false-positive probability with the Gaia neighbour list and a contrast curve; the G > 21 background-EB and bound-companion-EB scenarios are the ones that remain.
3. High-resolution imaging (AO or speckle) of the target to bound companions inside 1 arcsec.
4. A longer, independent look at the phase-0.5 brightening (e.g. a second detrending family; a per-quarter PSF-vs-SAP comparison at phase 0.5): if it is real at ~200 ppm it is not a planet signature and the signal should be downgraded.
5. Keep the v1.5 `secondary_2sided` fail attached to any report and state that the rule was introduced after this star was seen.

## 4. Files written
- `runs/poolA_main/adversarial/kic-009520080/redteam.md` (this file)
- `runs/poolA_main/adversarial/kic-009520080/redteam_out/`: `lc_tests2.py`, `lc_tests2.json`, `lc_tests2.log` (phase scan, odd/even timing, SAP/PSF), `pixel_depths.py`, `pixel_depths.json`, `pixel_depths.log` (per-pixel depths on the 17 TPFs), `neighbours.py`, `neighbour_catalog_rows.json` (KBONUS catalogue rows), `tpf_headers.json` (channel/module/output per quarter). Downloaded neighbour light curves and the first-pass `lc_tests.py/.log` remain in the scratchpad (`redteam-kic009520080/`).
