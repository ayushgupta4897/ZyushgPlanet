# Red-team review: kic-005950408 (KIC 5950408, Gaia DR3 2102123605035883392, TPF host KIC 5950394)

Date: 2026-10-17. Run id: poolA_main. Reviewer role: adversarial (try to kill the signal).
Candidate language only: this was a stage-2 *first-screen signal*; it is not a planet and, after this review, not a signal.

Inputs read: `runs/poolA_main/vet/kic-005950408/{vet.json,events.csv,fold_events.png}` (no `pixel/`
directory exists for this star), `runs/poolA_main/{triage_all,survivors,survivors_stage2,vet_summary}.csv`
rows for kic-005950408, `pipeline/prereg_poolA.json` revision_log, the independent re-derivation
(`adversarial/kic-005950408/rederive/`), `targets/poolA_kbonus_bg_dwarfs_ranked.csv` row 2739,
`data/known/*` (manifest fetched 2026-10-17T15:46Z), the KBONUS light curve
(`data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005950408_kepler_v1.0_lc.fits`, all 18 per-quarter
extensions and headers), and the Kepler TPFs of the host KIC 5950394 for Q1, Q5 and Q13 (downloaded today
to scratch via lightkurve; `kplr005950394-{2009166043257,2010174085026,2012179063303}_lpd-targ.fits.gz`).

My own code and outputs are in `adversarial/kic-005950408/redteam_work/`:
`rt_lc_tests.py -> rt_lc_tests.{json,png}` (per-cadence dissection, dropout census, random-ephemeris null,
Lomb-Scargle, duration geometry) and `rt_pixels.py -> rt_pixels.{json,png}` (host-TPF pixel check at the
three carrying epochs). No code from `pipeline/` was imported. Ephemeris used throughout: P = 38.60671 d,
t0 = BKJD 151.46748, box 0.75 h (vet.json); the 0.9-h box of the task text was run as well.

## Verdict: KILL

The 38.6-d / ~4100-ppm "signal" is not a transit of anything. It is the sum of six cadences, all on CCD
channel 25 (module 9, output 1; season 3), that sit inside multi-hour 2-6 % excursions of a light curve
whose season-3 photometry is dominated by the saturated wings and bleed columns of the Kp 8.45 red giant
KIC 5950394 (BD+41 3275) 30.6 arcsec away, the star whose target-pixel file the "target" lives in. The
target's 4-pixel KBONUS aperture lies 1.2-2.2 px from the giant's bleed columns in every season-3 quarter,
collects 8,500-10,900 e/s against the target's own ~5,100-5,400 e/s, and the PSF fit in those quarters has
reduced chi2 5e5-2e6. In the pixels the whole stamp loses no light at any of the three "events"; the
aperture content merely sloshes by +/-10 % as the giant's wing moves. Remove channel 25 and nothing is
left (random-ephemeris maximum over 30,000 trials 621 ppm; the pipeline period ranks 114,626th of 271,529
in the re-derivation). Independently of the pixel evidence, the depth/duration combination is unphysical
for the host (0.75 h at P = 38.6 d on a 0.91 Rsun dwarf needs b = 1.05, where a 4100-ppm depth is
impossible), the pipeline's own per-event statistics already said the events are consistent with zero
(event_weighted_depth 177 ppm vs box depth 4104 ppm), and the inverted-flux null on this star reached a
higher SDE than the real signal at stage 1.

## 1. What carries the depth (own code, `rt_lc_tests.json`)

At the vet.json ephemeris with the 0.75-h box there are 43 quality-0 cadences in the box (pipeline
n_in_transit = 43). Their local-detrended PSF flux (1-d running median, own code) gives:

| statistic (box 0.75 h) | PSF flux | SAP flux |
|---|---|---|
| mean in-box depth | 4980 ppm | 2976 ppm |
| median in-box depth | 468 ppm | 323 ppm |
| mean depth after dropping the 6 cadences below -1 % | **53 ppm** | - |
| in-box standard deviation | 13,065 ppm | - |
| formal error of the mean (std / sqrt 43) | 1992 ppm | - |
| out-of-box MAD (per-point) | 938 ppm | - |

With the 0.9-h box (52 cadences) the mean is 4095 ppm, i.e. I reproduce the pipeline's 4104 ppm exactly,
and the median is 445 ppm, the outlier-free mean 17 ppm. The pipeline's own `depth_pass2_formal_err_ppm`
is 1578 ppm: the headline depth is a 2.6-sigma number whose "empirical SNR" of 10.6 comes from dividing
by a 388-ppm box noise that does not describe the in-box scatter (13,000 ppm).

Per epoch (0.9-h box, 31 epochs with data): epoch 0 = -36,072 ppm (2 cadences), epoch 10 = -52,116 (1),
epoch 18 = -12,925 (1), epoch 29 = -37,491 (2); the other 27 epochs lie between -2,128 and +2,384 ppm.
Per season: season 3 (channel 25) 13,559 ppm from 15 cadences; seasons 0/1/2 = -242 / +287 / +710 ppm
(10 / 17 / 10 cadences). Everything is in season 3.

The six carrying cadences (all channel 25; "SAP" = KBONUS 4-px aperture flux, same detrend):

| cadence | BKJD | quarter | PSF ppm | SAP ppm | RED_CHI2 | distance to quarter end |
|---|---|---|---|---|---|---|
| 2081 | 151.457 | Q1 | -40,851 | -29,813 | 5.5e5 | 13.5 d |
| 2082 | 151.477 | Q1 | -31,294 | -24,537 | 5.5e5 | 13.5 d |
| 20975 | 537.530 | Q5 | -52,116 | -26,372 | 8.2e5 | **0.04 d (last cadences of Q5)** |
| 36091 | 846.401 | Q9 | -12,925 | -12,570 | 1.9e6 | 37.9 d |
| 56873 | 1271.055 | Q13 | -40,005 | -13,620 | 5.4e5 | **2.0 d, on the 43 % flux collapse at BKJD ~1273** |
| 56874 | 1271.075 | Q13 | -34,978 | -11,107 | 5.4e5 | 2.0 d |

None is an isolated dip of transit width. Neighbouring cadences (same file):
- epoch 0: -88 min -15,777 ppm, -59 min -29,234, 0 -40,851, +29 -31,294, +88 -31,065, +118 -30,615,
  +177 -27,913 ppm. A >= 5-hour, 3 % trough inside the Q1 segment whose raw PSF MAD is 68,748 ppm and
  which swings by +/-15 % over days (top panel of `rt_lc_tests.png`).
- epoch 10: -294 min -25,083, -177 -32,566, -88 -41,389, 0 -52,116, +59 -58,178 ppm, then the quarter
  ends (next cadence 46.6 h later, Q6, channel 17, flat). A monotonic 6 % ramp into the Q5 roll.
- epoch 29: -177 min +13,885, -59 -25,113, 0 -40,005, +29 -34,978, +59 -16,869, +118 +249,
  +147 +24,683 ppm: +/-4 % cadence-to-cadence oscillation on the shoulder of the Q13 collapse.

Dropout census (quality-0 cadences, local-detrended PSF flux): channel 25 has 13,169 cadences of which
425 are below -1 %, 48 below -3 % and 374 above +1 %. Channels 17, 57 and 65 together have 36,640
cadences with **zero** below -1 % and 6 above +1 %. The light curve has two regimes, and every carrying
cadence is in the bad one.

Depth x duration is conserved between stages: stage 1 (triage_all.csv) found 1972 ppm with a 1.5-h box,
stage 2 found 4104 ppm with a 0.75-h box (2958 vs 3078 ppm h). A transit's box depth does not double when
the box is halved; a fixed handful of outlier cadences does exactly that.

## 2. Pixel level: the events are not at the target (own code, `rt_pixels.json`, `rt_pixels.png`)

The host TPF (KIC 5950394, Kp 8.445, Teff 4893 K, R 16-19 Rsun in keplerstellar; 30.6 arcsec from the
target) has 22-24 saturated pixels and two bleed columns carrying 5.66-5.85e6 e/s. From the KBONUS
per-quarter headers (ROW, COLUMN) and the TPF WCS keywords the target sits at stamp x = 8.15, 9.23 and 8.19
px in Q1, Q5 and Q13, i.e. 1.15, 2.23 and 2.19 px from the nearest bleed column. The KBONUS 4-px aperture
sums to 8,531 / 10,901 / 8,832 e/s in Q1 / Q5 / Q13, against the target's own PSF flux of ~5,100-5,400 e/s
in seasons 0-1 (PSF_FLUX_NOVA column; the SAP level in seasons 0-2 is 3,300-3,800 e/s), so at least 40-50 %
of the light in the target aperture in season 3 is the giant's wing. This is why the PSF fit fails there
(median RED_CHI2 1.4e6-2.3e6 vs 440-5,800 on the other three channels) and why the raw per-quarter PSF MAD
is 21,000-89,000 ppm in Q1/Q5/Q9/Q13 vs 850-3,400 ppm elsewhere.

At each carrying epoch (in-box cadences minus the local out-of-box median, per pixel):
- Q1 / epoch 0: whole-stamp flux +129 and +309 ppm on the two in-box cadences (flat); bleed columns
  +289 / +459 ppm; the 4-px aperture -9,805 / -4,742 ppm while drifting from +85,000 to -20,000 ppm over
  the surrounding 28 h. Only one pixel in the whole stamp is below -3 sigma, at (row 11, col 8), 9 rows
  from the target. No light is lost; light moves.
- Q5 / epoch 10: stamp +61 / +159 ppm; aperture -426 / -606 ppm, in the middle of a drift from +80,000 to
  -100,000 ppm over 28 h ending at the quarter roll. The significant difference-image pixels form a
  dipole: (row 9, cols 6-8) at -4.0 to -4.7 sigma at the foot of the bleed column, (row 7, cols 7-8) at
  +4.2 sigma, and (row 23, cols 7-8) at +4.4 / +5.1 sigma, i.e. the pixels next to the target
  *brightened*. That is a sub-pixel motion of the saturated wing, not a dip on a star.
- Q13 / epoch 29: stamp -1,280 / -1,029 ppm and bleed columns -1,126 / -846 ppm (the giant's own
  0.1 % wobble on the Q13 shoulder); aperture -23,396 / -19,730 ppm, i.e. the aperture amplifies the
  giant's wobble twenty-fold. The only pixels below -3 sigma are (10,6), (11,7), (13,6), (15,6), all on
  the bleed column 7-12 rows from the target. The difference-image sum over the stamp is -11,680 e/s,
  84 % of it in the bleed columns.

The missing light expected from a 4 % dip on the target is ~210 e/s; the aperture difference is indeed
-211 / -245 / -248 e/s, but it is produced by the giant's wing redistributing, as the flat stamp totals,
the dipole structure and the bleed-column pixels show. There is nothing here for `pipeline/pixel_vet.py`
to localise, and the KBONUS PSF photometry of this star is not usable in season 3.

## 3. Statistics of the box depth (random-ephemeris null)

30,000 random ephemerides (P log-uniform 20-40 d, random phase, 0.75-h box, same detrended PSF flux,
same mean-in-box depth statistic): robust std 414 ppm, which is the pipeline's "box noise" (388 ppm), but
the distribution is not Gaussian: p99 = 1,458 ppm (3.5 robust sigma), p99.9 = 2,417 ppm, maximum 4,166 ppm,
and 2 of 30,000 random ephemerides reach the pipeline's 4,104 ppm. A BLS search over 0.5-40 d tests far
more than 10^4 independent (period, phase) placements, so one 4,100-ppm box is expected from this light
curve by chance. Without the channel-25 quarters the same null has robust std 131 ppm, p99.9 = 416 ppm and
a maximum of 621 ppm over 30,000 trials: the heavy tail is channel 25 and nothing else.

The pipeline's own stage-1 inverted-flux null on this star (triage_all.csv): top inverted peak P = 0.527 d,
SDE 12.12, higher than the real-flux SDE 9.83; it was rejected only on snr_emp, odd_even and duty_cycle.
Real vs inverted pass count for this star is therefore 1 : 0 on the full rule set but 0 : 1 on SDE alone,
which is what a periodogram dominated by symmetric +/- excursions gives (374 cadences above +1 % and 425
below -1 % on channel 25). The second real signal at 35.15 d failed on odd/even at 8.6 sigma.

## 4. Pipeline-internal evidence that was already fatal

- `event_weighted_depth_ppm` = 177 and `depth_by_quartile_ppm` = 106 / 268 / 253 / -114 ppm against a
  headline `depth_pass2_ppm` of 4104: the pipeline's per-event depths are 20x smaller than its box depth.
  `events.csv` contains 7 events (epochs 5, 8, 19, 24, 26, 30, 34; all with 2 in-box cadences, depths
  -657 to +380 ppm, every one consistent with zero); the three epochs that carry the depth (0, 10, 29) are
  not among them, so `event_chi2_const` = 1.77 / 6 dof tested only the events that do not contain the
  signal. A box depth that the per-event depths cannot reproduce should be an automatic stage-2 failure.
- `depth_pass2_formal_err_ppm` = 1578: the depth is 2.6 sigma formally; `snr_emp_pass2` = 10.6 uses a
  388-ppm out-of-box box noise that ignores the 13,000-ppm in-box scatter.
- `duration_ratio` = 0.12 (0.75 h of a 6.07-h maximum): section 6 shows this is unphysical for the host.
- `stage2_pass_v14` = False (TLS SDE 7.08 < 9). The v1.5 revision (prereg revision_log) dropped the TLS SDE
  rule as uncalibrated, after which the star passed under v1.7. I do not argue to reinstate that rule; the
  rules that should have caught this are the event-vs-box depth consistency and a per-channel scatter or
  RED_CHI2 screen (see section 8).
- Stage 1 `real_max_single_event_frac` = 0.34 with the 1.5-h box, i.e. one event carried a third of the
  power already at the first screen.
- Odd/even: pipeline 0.05 sigma (4082 vs 4122 ppm). With only three carrying epochs (0 and 10 even, 29
  odd; 18 even) that agreement is a coincidence of two outliers in each parity, not a test passed; the
  re-derivation finds 4.9 sigma at a t0 shifted by 0.002 d, which shows how fragile the statistic is here.

## 5. Blended / background eclipsing binary, ephemeris contamination, variability

- Ephemeris match (local KOI cumulative, DR25 TCE, Kirk+2016 EB, confirmed Kepler-field planets;
  harmonics 1/4, 1/3, 1/2, 1, 2, 3, 4 at 0.2 % tolerance): 29 KOIs, 82 TCEs and 6 EBs match a period
  harmonic somewhere in the field; none matches the epoch (closest |delta phase| = 0.086 for K00784.01 at
  P/2 and 0.095 for the KIC 5873209 TCE at P, i.e. 3.3-3.7 d off; the box is 0.75 h). Nothing within
  3 arcmin in any of the four catalogues or in the TESS EB catalogue. Live NASA Exoplanet Archive TAP
  (2026-10-17): no KOI or DR25 TCE for KIC 5950394, 5950408, 5950377 or 5950466; no TOI within a
  0.14 x 0.10 deg box. The signal exists only on channel 25, so a contaminating KOI/EB would have to be on
  module 9 output 1 *and* match period and epoch; none does.
- The only plausible contaminator, the Kp 8.45 giant itself, is flat at epochs 0 and 10 (stamp totals
  within +/-300 ppm) and wobbles by 0.1 % at epoch 29; it is not an EB (not in Kirk+2016) and has no TCE
  in 17 quarters of saturated-aperture photometry. Its light is not dipping; it is moving across the
  target's aperture.
- SAP/PSF depth ratio: 0.46-0.60 (mine: 2976 / 4980 = 0.60 with the 0.75-h box; re-derivation 0.46 with
  0.9 h) against >= 1 for any astrophysical dip in the aperture. Consistent only with the PSF model
  misbehaving, as the RED_CHI2 values say directly.
- V-shape / secondary / odd-even are meaningless for a signal made of six cadences in three quarters; for
  completeness: secondary 305 ppm / 0.8 sigma (pipeline), -207 ppm / -0.6 sigma (re-derivation).
- Lomb-Scargle (own, 0.1-100 d, both polarities implicit): per-quarter-normalised PSF flux peaks at
  52.7 d with 2,548 ppm half-amplitude, driven entirely by the channel-25 swings (the other channels have
  no excursion beyond 1 %); locally detrended flux peaks at 3.70 d with 248 ppm. Neither is a harmonic of
  38.6 d (ratios 0.73 and 10.4), and neither is astrophysical. The pipeline's 0.826 d / 110 ppm is likewise
  unrelated. No rotation alias.
- Kepler systematics: of the three carrying epochs, epoch 10 is the final two cadences of Q5 (0.04 d from
  the quarter roll) and epoch 29 is on the Q13 end-of-quarter flux collapse (43 % in PSF flux at
  BKJD ~1273); epoch 0 is in the thermally unsettled first weeks of Q1. All three are classic quarter-edge /
  thermal-ramp loci, amplified here by the saturated neighbour.

## 6. Duration, density and implied radius

Host (target table): G = 16.30, BP-RP = 0.86, Teff 5610 K, d = 2420 pc, Rstar 0.91 Rsun, RUWE 0.94.
With M = 1.0 Msun, a/Rs = 52.7 and the b = 0 duration for a k = 0.064 body is 5.96 h (pipeline t14_max
6.07 h). A 0.75-h (or 0.9-h) total duration requires b = 1.052-1.056, within 0.01 of the grazing limit
1 + k = 1.064, where the covered area is a sliver and the depth would be a small fraction of k^2, not
4,100 ppm. There is no stellar radius and mass that makes a 0.75-h, 4,100-ppm box at P = 38.6 d a planetary
transit (a central transit this short would need R*/M^(1/3) eight times smaller than the Sun's). The
nominal Rp of 6.4 Re (box) or 3.9 Re (TLS) is therefore not a physical quantity for this signal.

## 7. Is it already known?

Not in KOI, DR25 TCE, Kirk+2016 EB or the confirmed-planet list (local copies, 2026-10-17; live TAP
confirms for the four KICs within 1 arcmin). KIC 5950408 was never a Kepler target (absent from
keplerstellar; observed only as background in the KIC 5950394 TPF). SIMBAD (astroquery, today): nothing
within 30 arcsec of the target; within 1 arcmin only BD+41 3275 (= KIC 5950394, 30.6"), Gaia DR3
2102123403176985216 (34.9"), 2MASS J19142664+4114497 (41.3"), UCAC4 657-070385 (56.3") and
2MASS J19143364+4115477 (57.1"). Bienias & Szabo 2025 (A&A 694, A209): their seven hot-Jupiter candidates
(KIC 4459924, 9674230, 4543171, 9040849, 8175131, 8293539, 11033410, P = 2.5-4.8 d, as recorded today from
the Zenodo Appendix B by the kic-003831969 red team) do not include this star or its host. Moquin et al.
2026 (BOI-1.01, P = 3.37 d) is unrelated. Novelty is irrelevant: there is no signal to be novel.

## 8. Recommendations (for the maintainers; any rule change is post-hoc and must be logged)

1. Drop kic-005950408 from every survivor list; mark the v1.7 stage-2 pass as a false positive of the
   "saturated-neighbour / failed-PSF-fit" class.
2. Add, with a dated revision_log entry, a per-quarter data-quality screen before BLS: e.g. exclude quarters
   whose median RED_CHI2 exceeds ~1e5 or whose raw PSF MAD exceeds ~10x the star's median, and require a
   signal to survive without them. On this star that removes Q0/1/5/9/13/17 and the signal.
3. Require stage-2 consistency between `event_weighted_depth_ppm` and `depth_pass2_ppm` (e.g. ratio within
   a factor 2 at 3 sigma) and include *every* epoch with in-box data in events.csv, not only those with
   n_in >= 2 and n_out >= 6; the carrying epochs here were exactly the ones excluded.
4. Report `depth_pass2_formal_err_ppm` next to `snr_emp_pass2`; a 2.6-sigma depth with a "10.6-sigma" SNR
   should not pass.
5. Flag any pool star whose TPF host is brighter than Kp ~10 (saturation) for a mandatory
   bleed-column-distance check in the season(s) where PSFFRAC/SAP levels jump.

## Checks run

- Read vet.json, events.csv, fold_events.png, the four poolA_main summary CSV rows, the prereg revision log
  and the independent re-derivation (results JSON, log, fold figure).
- Own per-cadence dissection at the pipeline ephemeris (0.75-h and 0.9-h boxes): PSF and SAP local
  depths, mean vs median vs outlier-free depth, in-box scatter and formal error, per-epoch and per-season
  depths, RED_CHI2, centroid offsets, FLATTEN_MASK, gap and quarter-edge distances, +/-4-cadence
  neighbourhoods of each carrying cadence.
- Dropout census per quarter and per channel (|excursion| > 1 % and > 3 %).
- Random-ephemeris null for the mean-in-box depth statistic (30,000 trials, P 20-40 d, 0.75-h box), with
  and without the channel-25 quarters.
- Lomb-Scargle 0.1-100 d on per-quarter-normalised and on locally detrended PSF flux; harmonic ratios to P.
- Duration/density geometry for the catalogue host; implied impact parameter and radius.
- Host-TPF pixel check at epochs 0, 10, 29 (Q1, Q5, Q13 downloaded via lightkurve): saturation and bleed
  columns, target position from KBONUS headers, per-pixel difference images and SNR maps, whole-stamp /
  bleed-column / KBONUS-aperture / non-bleed flux series cadence by cadence, Kepler QUALITY flags.
- Known-object checks: 3-arcmin cone and period/harmonic + epoch match against KOI cumulative, DR25 TCE,
  Kirk+2016 EB, confirmed Kepler-field planets, TESS EB; live NEA TAP for the four KICs within 1 arcmin and
  for TOIs; SIMBAD 1-arcmin cone; Bienias & Szabo 2025 list and Moquin et al. 2026.
- Real vs inverted pass accounting from the pipeline's own stage-1 record.

## What would settle it (nothing is recommended)

The only thing that could resurrect a 38.6-d signal on this star is a re-reduction of seasons 0, 1 and 2
alone showing a coherent dip at the ephemeris. Those seasons already give -242 / +287 / +710 ppm in the
box (10 / 17 / 10 cadences; re-derivation dataset B: 477 +/- 122 ppm with the period ranking 114,626th), so
the answer is already in hand: no.
