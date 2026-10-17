# Red-team review: kic-006699887 (KIC 6699887, Gaia DR3 2077554781499350144)

Date 2026-10-17. Run id poolA_main. Role: try to kill the 15.11-d signal. Candidate language only; nothing here is a planet.

Evidence read: `runs/poolA_main/vet/kic-006699887/{vet.json,fold_events.png,events.csv}`, `runs/poolA_main/pixel/kic-006699887/{pixel.json,diffimages.png}`, `runs/poolA_main/results.jsonl` (stage-1 real + inverted), the independent re-derivation under `adversarial/kic-006699887/rederive/`, target row 452 of `targets/poolA_kbonus_bg_dwarfs_ranked.csv`, the 17 cached TPFs of the originating target KIC 6699878 (`/root/.cache/lightkurve/mastDownload/Kepler/kplr006699878_lc_Q011111111111111111/`), and the KBONUS source catalogue.

Own code (no `pipeline/` imports), outputs and figure: `runs/poolA_main/adversarial/kic-006699887/redteam/` — `pixel_edge_test.py` (+`_ephem.py`), `lc_systematics.py`, `pixel_edge_test.json`, `pixel_edge_test_p47.json`, `lc_systematics.json`, `second_signal.json`, `phase_scan_clean.json`, `redteam_summary.png`. Scratch copies in `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic006699887/` (the host's KBONUS light curve was downloaded there, `dl/`).

## Verdict: survives

I could not kill it. Every attack below either fails outright or ends in the generic residual (an unresolved blend inside ~1.5 arcsec that only a TRICERATOPS-style calculation can price). Two findings go beyond the evidence pack: (i) the pixel deficit is pinned to the target even when the source is allowed to sit *beyond* the stamp edge the target lives on, and (ii) the pipeline's stage-1 "second signal" at 23.845 d (rejected for a 10-sigma odd/even mismatch) is a real, independent 47.69-d transit-like signal on the same star that fell outside the 0.5-40 d search window and was picked up at its P/2 alias.

## 1. Geometry that makes this star fragile

- The target (G 16.22, Kp 16.26, KBONUS PSF flux 4629 e/s) is a background source in the postage stamp of KIC 6699878 (G 15.09, Kp 15.12, PSF flux 12632 e/s, Teff 4516 K, R 0.60-0.64 Rsun: a K dwarf like the target), 9.70 arcsec away (Gaia DR3 TAP, 2026-10-17). KIC 6699887 was never a Kepler target (`data/known/kepler_target_kepids.csv`: absent).
- TPF headers: skygroup 63; season 0 = Q2/6/10/14 channel 23 (mod 8 out 3), season 1 = Q3/7/11/15 ch 39 (mod 12.3), season 2 = Q4/8/12/16 ch 63 (mod 18.3), season 3 = Q1/5/9/13/17 ch 47 (mod 14.3). Stamps are 4x5, 4x6, 5x5, 4x6 px.
- In every season the target sits in the **first column** of the stamp: x = 0.84 (s0), 0.18 (s1), 0.68 (s2), -0.05 (s3) px. KBONUS psffrac = 0.73, npixsap = 3 (the SAP aperture is columns 0-1 only). So a significance-weighted deficit centroid (what `pixel_vet.py` reports, median 1.1 arcsec) cannot by construction move beyond the edge; the pipeline's test never asked whether the source could be *outside* the stamp on the far side of the target. That was my main attack.

## 2. Pixel-level: can the deficit be beyond the stamp edge, or on the host? (own code)

`pixel_edge_test.py`: per quarter, event-differenced per-pixel deficit maps (in-transit core vs a local linear out-of-transit baseline, inverse-variance stacked over events), Gaussian PSF width fitted on the out-of-transit image (sigma 0.55-0.73 px), then a single-Gaussian deficit fit with the position free over -3 to nx+2 px (i.e. allowed up to 3 px beyond the edge), and chi2 profiles along the edge-normal. Season stacks (identical geometry):

| season | target x | free-fit position (px) | offset from target | 3-sigma x range | deficit at target (e/s) | deficit at host (e/s) | chi2 null / at host / at target / free |
|---|---|---|---|---|---|---|---|
| 0 (Q2,6,10,14) | 0.84 | (0.60, 2.40) | 1.0" | [0.5, 0.7] | 4.1 | 0.1 | 250 / 250 / 36 / 17 |
| 1 (Q3,7,11,15) | 0.18 | (-0.20, 2.90) | 1.7" | [-0.6, 0.0] | 6.0 | 0.9 | 406 / 402 / 133 / 99 |
| 2 (Q4,8,12,16) | 0.68 | (0.80, 2.00) | 0.7" | [0.8, 0.9] | 6.0 | 0.5 | 434 / 433 / 141 / 118 |
| 3 (Q1,5,9,13) | -0.05 | (-0.50, 2.00) | unconstrained | [-3, -0.1] | 5.2 | 0.3 | 338 / 337 / 55 / 45 |

- Seasons 0 and 2, where the target centre is inside the stamp, pin the source to within 0.2-0.3 px (about 1 arcsec) of the Gaia position of the target and exclude any position beyond the edge (x < -0.5) at far more than 3 sigma. Season 3 (target centre just outside the edge) is, as expected, unable to distinguish target from beyond-edge; season 1 is intermediate. The beyond-the-edge hypothesis is therefore dead on seasons 0 and 2 alone.
- The host KIC 6699878 is excluded in every season: the deficit evaluated at the host position is 0.1-0.9 e/s (a 400-ppm dip on the host would be 5 e/s) and chi2 at the host equals the null.
- Deficit amplitude in e/s (4.1-6.0 per season; Gaussian-extrapolated total) is consistent with the light-curve depth on a 4629 e/s star (1.08e-3 x 4629 = 5.0 e/s). It is *not* consistent with a dip on a much brighter or much fainter unresolved star unless the dip is correspondingly shallower/deeper, which is the usual degeneracy.
- Pipeline note: `pixel.json` reports Q2, Q13, Q14 as "no significant pixels" (Q13 chi2_null = 0.1). My maps for those quarters show the deficit at the target at the same level as other quarters (Q13: chi2 null 116 -> 47 at target; Q2: 73 -> 32; Q14: 78 -> 28). Those pipeline quarters look like an error-inflation glitch in `pixel_vet.py`, not a data problem; worth checking because the same glitch could hide a real off-target deficit on another star.

## 3. Host star's own light curve

KBONUS light curve of KIC 6699878 (kic-006699878, downloaded 2026-10-17), same masked 1-d running-median detrend, same ephemeris and 3.75-h box: PSF depth **-1.8 +- 15.8 ppm**, SAP depth 17 +- 17 ppm (scatter 330 ppm). If the dip were on the host and leaking into the target's PSF flux, the host itself would show at least 1079 x 4629/12632 = 395 ppm. Excluded at >20 sigma. (Second signal, section 8: host -21 +- 24 ppm.)

## 4. Light-curve attacks (own detrend: per-quarter median, 1-d running median with transits masked; 60,092 cadences)

- Depth, 3.75-h box: PSF 1079 +- 33 ppm; SAP 959 +- 34 ppm; SAP/PSF = 0.89 vs CROWDSAP 0.885 (expected for an on-target signal; a contaminant would push SAP above PSF).
- 86 events; 77/86 individually > 2 sigma, 0 negative; weighted mean 1074 ppm; leave-one-out range 1060-1091 ppm (no single event drives it). Event chi2 153/85 with formal errors; the pipeline's empirical box noise (51 vs 33 ppm formal, x1.5) brings this to ~68/85, so the excess scatter is red noise at the 4-h timescale, not depth variability.
- Gaps / downlinks / ramps: 8 events within 1 d of a >0.5-d gap, depth 1086 +- 103 vs 1073 +- 32 far from gaps; 5 events within 2 d after a gap (thermal ramp window) 1226 +- 132 vs 1066 +- 31 (1.2 sigma); observed fraction of events in the 2-d post-gap window 5.8% vs 9.6% expected for random phasing; only 2 events within 0.5 d of a quarter edge. 15.11 d is not commensurate with the ~31-d downlink or ~93-d roll cadence, and the signal is present in all 17 quarters on four different channels: not rolling band, not a channel-specific artefact.
- Per quarter (17) chi2 38/16 with formal errors (~17/16 with the x1.5 empirical scaling). Per season PSF depths 998 +- 69 / 1045 +- 64 / 1278 +- 66 / 1036 +- 64 ppm (chi2 11/3 formal, ~5/3 scaled); SAP 809 / 1037 / 1195 / 837. Season 2 is 20% deeper in both PSF and SAP and is not driven by one event (season-2 event median 1260 ppm). I read this as season-dependent flux assignment of an edge-of-stamp star by the PSF fit (crowding with the host changes with the roll), which makes the **absolute depth, hence Rp, uncertain at the +-15-20% level**; it does not discriminate on-target from blend.
- Secondary (clean scan, reference = all out-of-transit cadences with both signals masked): phase 0.5, 3.75-h box: **-26 +- 29 ppm**; 12-h box: -10 +- 17 ppm. Box depth at all 702 phases outside the transit: z-scores have std 1.02, max +3.07 (90 +- 29 ppm at -158 h), min -2.76: pure noise. No secondary at any phase, so no P/2 EB alias (depth at phase 0.5 would need to be ~1000 ppm) and no eccentric-EB secondary.
- Shape: pipeline core/box ratio 1.09 and the rederivation's trapezoid give a flat bottom; at 30-min cadence T23/T14 for a 2.4 Re object is unmeasurable, so V-shape is not a usable discriminator here, only a non-kill.
- Variability: Lomb-Scargle on quarter-normalised out-of-transit flux, 0.1-60 d: top peaks 43, 42, 48, 29, 50 d (long-term systematics, 180-225 ppm amplitude); power at P and 2P = 0.002 / 0.001 (nothing). The pipeline's LS period 0.687 d (39 ppm) and the stage-1 inverted-null peak at 0.657 d (38 ppm) are a sub-50-ppm short-period ripple, irrelevant to a 1000-ppm 3.75-h box.
- Inverted flux (2-f), masked detrend, BLS 14.5-15.8 d: best 15.48 d at 154 ppm, nothing at P (the rederivation's un-masked inverted peak at P is the running-median shoulder artefact it describes). Stage-1 inverted null for this star peaked at 0.657 d, SDE 14.1, 38 ppm (fails snr_emp and duty cycle); real s1 SDE 89, SNR_emp 22.3, 86 transits.

## 5. Ephemeris-match / contamination from known objects

- Local catalogues (`data/known/`, 9564 KOIs, 34,032 DR25 TCEs, 2876 Kepler EBs, 2839 confirmed Kepler-field planets): nothing within 3 arcmin except the 0.2012-d contact EB KIC 6699679 at 2.73 arcmin (unrelated period). Period matches at ratios 1:1, 2:1, 1:2, 3:1, 1:3, 3:2, 2:3, 4:1, 1:4, 5:1, 1:5, 6:1, 1:6 within 0.2%: the only ones with epoch offsets < 5 h are TCE 9714572 (6P, 254 arcmin away, impossible for the same module), TCE 11974090 (2P, 504 arcmin), TCE 5396936 (P/5, 273 arcmin) and KOI-176.01 / Kepler-483 b (2P, 113 arcmin, -5 h, and a 2P contaminant would fail odd/even). Nothing at 15.1136 d itself closer than 0.04% (TCE 8574270, epoch off by 10-119 h, 189 arcmin). No KOI/TCE/EB shares this ephemeris; direct-PRF, column-anomaly and video-crosstalk contamination all require an ephemeris match and there is none.
- Already known? Not a KOI, not a DR25 TCE (the host KIC 6699878 has no TCE either, consistent with the signal being outside its optimal aperture: pipeline depth_in_aperture ~22 ppm), not in the Kepler EB catalogue, not a confirmed planet host, not one of the 7 Bienias & Szabo 2025 hosts (KIC 4459924, 4543171, 8175131, 8293539, 9040849, 9674230, 11033410), not Moquin et al. 2026 BOI-1.01 (P 3.37 d, G 17.8 M dwarf). SIMBAD within 30 arcsec returns only UCAC4 661-077711 (= the host, no classification); VizieR Gaia DR3 variability tables and the EB catalogue return nothing. Catalogue download dates are those recorded in `data/known/`. The KBonus team's TraSHDUMP pipeline is active and unpublished, so "unreported" cannot be asserted, only "not in any public list on 2026-10-17".

## 6. Host plausibility, duration, radius

- Gaia DR3: parallax 1.164 +- 0.032 mas (859 pc), BP-RP 1.35, M_G 6.55: a main-sequence K dwarf (not a giant, so no "dip too shallow for a giant" angle). RUWE 1.087, astrometric_excess_noise 0.16 mas (2.3 sigma), ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.015, non_single_star 0: no sign of a resolved or semi-resolved companion down to ~0.3-0.5 arcsec at modest contrast.
- Duration: intrinsic T14 ~3.5 h (3.75-4.0 h observed includes ~0.5 h of 30-min smearing). For b = 0 and k = 0.034: rho* = 3.3 (3.5 h), 2.7 (3.75 h), 2.2 (4.0 h) g/cc vs 3.6 g/cc for the catalogue R 0.66 / M 0.73. Consistent with a near-central transit and/or R* 10% larger; the host KIC 6699878 (3.8-4.6 g/cc) would fit equally well, so duration does not discriminate, only the pixels and the host light curve do (sections 2-3).
- Radius: 2.4 Re for 1079-1141 ppm on 0.66 Rsun, with the +-15-20% depth systematic above, i.e. ~2.2-2.7 Re. A 1000-ppm dip is far too shallow for a stellar companion on the target itself; an on-target EB would need a secondary with R ~ 0.03 Rsun, which does not exist, so the only EB routes are blends.

## 7. What is left (why "survives" is not "validated")

- An unresolved eclipsing binary within ~1-1.5 arcsec of the target (Gaia-unresolved if dG >~ 3) would need an eclipse depth >= 1.1e-3 x 10^(0.4 dG): 1.7% at dG = 3, 11% at dG = 5, 43% at dG = 6.5. EBs do that. The flat bottom and the absence of any secondary at 30 ppm reduce but do not remove this; it needs a TRICERATOPS-style FPP/NFPP with the Kepler PRF and the Gaia-based contrast curve, which this stage does not provide.
- Depth/radius carry a +-15-20% systematic from season-dependent PSF flux assignment of an edge-of-stamp star.

## 8. Bonus finding: a second transit-like signal on the same star (P = 47.69 d)

The stage-1 second search (`real.s2` in results.jsonl) found 23.845 d, 4.5 h, 442 ppm, with odd/even mismatch 10.4 sigma (odd 101 / even 771 ppm) and was dropped. With the 15.11-d transits removed entirely and a BLS over 40-60 d: **P2 = 47.6906 d, t0 = 135.7185 BKJD (BJD 2454968.7185), duration 4.5 h, depth 877 +- 49 ppm**, 26 events, 23/26 individually > 2 sigma, event chi2 26/25, largest single-event share 7%, SAP depth 781 +- 54 ppm, secondary at phase 0.5 126 +- 53 ppm (2.4 sigma, noise-compatible), **host KIC 6699878: -21 +- 24 ppm**. Pixel localisation with the same free-position Gaussian test (`pixel_edge_test_p47.json`): seasons 0 and 2 pin the deficit to within 0.2 px of the target (free fits (0.60, 2.30) and (0.80, 2.00); 3-sigma x ranges [0.4, 0.9] and [0.5, 1.0] px), deficit at the host -1.1 to 0.7 e/s. Period ratio 47.6906/15.1136 = 3.1555 (just wide of 3:1, as many Kepler pairs are); duration 4.0-4.5 h implies b ~ 0.4-0.6 for the catalogue density; Rp ~ 2.1 Re. This is the pipeline's s2 at its P/2 alias: the 0.5-40 d window clipped the true period and the "odd/even failure" was the alias. For the red team this matters because (a) two independent, flat, secondary-free, pixel-localised signals on one star make the blended-EB explanation for the 15.11-d signal much less likely (both would have to come from the same unresolved contaminant), and (b) it is a concrete pipeline note: real signals with P > 40 d will surface as odd/even failures at P/2 and could be discarded for the wrong reason. It is post-hoc (found while reviewing), unvetted beyond the checks listed, and is not being claimed here.

## 9. Pipeline notes (not part of the verdict)

1. `pixel_vet.py` quarters Q2, Q13, Q14 report no significant pixels while my own maps show the usual deficit at the target; check the per-pixel error estimate.
2. The deficit-centroid statistic cannot move beyond the stamp edge; for edge sources add a free-position PRF fit (as here) or at least report that the test is one-sided.
3. Second-signal searches should run to at least 2x the primary window, or the odd/even cut will discard real long-period signals at their P/2 alias.
4. The stage-1 inverted null on an un-masked detrend will show a peak at the transit period for strong detections (running-median shoulders); the rederivation already noted this.

## 10. Numbers to quote

Pipeline (vet.json, poolA_main): P 15.1136115 d, t0 BJD 2454967.72786, 3.75-h box depth 1023 +- 33 ppm (formal), empirical SNR 20.0, TLS SDE 13.7, 86 events, odd/even 0.03 sigma, secondary -37 ppm (-0.73 sigma), Rp 2.30 Re. Stage-1: BLS SDE 89.2, SNR_emp 22.3, inverted null SDE 14.1 at 0.657 d (fails). Red team (this file): host-star PSF depth -1.8 +- 15.8 ppm; free-position pixel deficit within 1.0" (s0) and 0.7" (s2) of the target with beyond-edge positions excluded at >3 sigma; secondary -26 +- 29 ppm; second signal 47.6906 d, 877 +- 49 ppm, host -21 +- 24 ppm.
