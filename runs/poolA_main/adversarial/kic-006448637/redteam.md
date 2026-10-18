# kic-006448637 — red-team review (2026-10-17)

Candidate language only. This is a transit-like **signal** in the KBONUS-BKG PSF light curve of
KIC 6448637 (Gaia DR3 2077617724235317504; G = 15.73, BP-RP = 1.10, Teff 5024 K, R* ~ 0.77 Rsun,
1006 pc, RUWE 1.02; originating TPF KIC 6448574). Nothing here is a planet and nothing here is
"discovered". Run ids: stage-2 `runs/poolA_main/vet/kic-006448637/vet.json` (rules v1.5),
pixel `runs/poolA_main/pixel/kic-006448637/pixel.json`, independent re-derivation
`runs/poolA_main/adversarial/kic-006448637/rederive/`. KBONUS FITS SHA-256 prefix `fbfdb14841d7a20a`.
Known-object catalogues: `data/known/manifest.json`, fetched 2026-10-17 15:46 UTC.

**Verdict: survives** (passes everything I could throw at it; not yet a candidate — see "what is still
missing"). My job was to kill it and I could not. The one formal failure in the evidence pack, the
pipeline pixel test (0 of 13 quarters), is a failure of that test on this TPF, not of the signal:
a localised re-analysis of the same TPF pixels puts the missing light on the target in 12 of 12
quarters where the target is on the stamp.

## 0. What the signal is (numbers from the evidence pack)
- P = 35.58726 d (pipeline BLS) / 35.58777 +- 0.00045 d (re-derivation timing fit); t0 = 186.6862 BKJD.
- Depth 846 +- 52 ppm formal, box_noise 85 ppm -> empirical SNR 10.0 (pipeline, 4.5-h box);
  897 +- 45 ppm, SNR 20 (re-derivation, 4.13-h box). The pipeline's error is the conservative one.
- 26 events with in-transit cadences (24 with >= 50% coverage); all positive; event chi2 35/25.
- Odd/even 790 vs 911 ppm (0.7 sigma pipeline; 2.0 sigma re-derivation; 0.4 sigma without epoch 28).
- Secondary at phase 0.5: 24 ppm (0.3 sigma) / 44 +- 44 ppm; nothing deeper than ~100 ppm at any phase 0.2-0.8.
- Inverted-flux null: SDE 9.9 vs 50.6 for the real signal (re-derivation blind search).
- TLS (narrow window): SDE 11.6, SNR 10.7, Rp ~ 2.3 Re. LS variability 32 ppm at 0.74 d.

## 1. Blended / background eclipsing binary
### 1a. Pipeline pixel test is structurally insensitive here
`pixel.json`: deficit centroid 70-115 arcsec from the target, median peak SNR 3.6, in-aperture depth
87 ppm, "0 of 13 quarters closer to target than any other source". Look at `diffimages.png`: the TPF
belongs to KIC 6448574, a **Kp 8.75 saturated star** with a 40-80-row bleed column; the stamp is
52x13 to 79x12 pixels and the target sits in the bottom 1-4 rows, 1.7 px beside the bleed column
(`tpf_local.log`: target (x,y) = (5.2, 1.9) in Q4, (5.3, 1.1) in Q5, (3.3, 18.3) in Q2; host at
(6.9, 25.0) etc.). The pipeline's significance-weighted centroid over the whole stamp is set by
bleed-column pixels whose event-differenced noise is hundreds of e/s; the target's whole eclipse
deficit is ~5 e/s (850 ppm of a ~7,000 e/s star). "Depth in aperture" is the depth in the
*saturated host's* aperture (0.7 ppm expected). The test cannot see the target at all on this TPF.
This should be fed back to `pipeline/pixel_vet.py` as a saturated-host / edge-of-stamp mode; I did not
modify the pipeline.

### 1b. Localised re-analysis of the same TPF (my code, `redteam/tpf_local.py`, `redteam/tpf_radius.py`)
Per quarter, box photometry on raw TPF pixels (r <= 1.6 px around the Gaia position, bleed column
excluded) with a 1-d running median masked at the ephemeris, plus a per-event out-minus-in deficit map.
Target-box depth (ppm) in the 12 quarters where the target is on the stamp (Q0, Q1, Q3, Q6, Q7, Q9: target
off the stamp or < 2 events):

| Q | season | depth +- err (ppm) | deficit (e/s) | nbr box (G=16.0 star) |
|---|---|---|---|---|
| 2 | 0 | 952 +- 197 | 3.0 +- 6.0 | -0.5 +- 1.0 |
| 4 | 2 | 633 +- 187 | 5.1 +- 2.3 | -1.6 +- 3.9 |
| 5 | 3 | 669 +- 167 | 4.1 +- 1.6 | -0.1 +- 1.2 |
| 8 | 2 | 561 +- 172 | 4.7 +- 2.5 | +1.0 +- 4.2 |
| 10 | 0 | 954 +- 157 | 6.2 +- 4.9 | off stamp |
| 11 | 1 | 511 +- 142 (610 incl. bleed px) | 3.3 +- 1.4 | off stamp |
| 12 | 2 | 754 +- 162 | 5.1 +- 2.8 | off stamp |
| 13 | 3 | 1134 +- 342 | 5.2 +- 0.9 | off stamp |
| 14 | 0 | 783 +- 146 | 6.5 +- 2.9 | off stamp |
| 15 | 1 | 1228 +- 183 | 7.4 +- 1.6 | off stamp |
| 16 | 2 | 698 +- 160 | 4.3 +- 2.1 | off stamp |
| 17 | 3 | 379 +- 382 | 4.0 +- 7.5 | off stamp |

Every quarter shows a positive deficit at the target of 3-7 e/s, i.e. the same ~850 ppm seen in the
KBONUS PSF and SAP fluxes, from raw pixels with none of the KBONUS machinery. Radial profile
(`tpf_radius.log`): the deficit is already complete at r <= 1.0-1.6 px (3-8 px) and does not grow
to r <= 3.5 px; the 2.5-4.5 px annulus is consistent with zero in every quarter -> a point source
at the target, not an extended or offset contribution. Deficit centroid relative to the target,
averaged over 12 quarters: **dx = 0.00 +- 0.06 px, dy = +0.21 +- 0.07 px (1 px = 3.98 arcsec)**,
i.e. within ~0.9 arcsec of the Gaia position. The G = 16.03 neighbour is at (dx, dy) = (-2.2, -1.4)
px; the saturated host at dy = +23 px. The 0.2-px offset toward +y is at the level of Kepler WCS /
PSF-truncation errors (the target's PSF is cut by the stamp edge in most quarters; psffrac = 0.50).

### 1c. Neighbours' own KBONUS light curves at this ephemeris (same box, own detrend)
| source | G | sep | PSF flux (e/s) | depth at ephemeris | required if it were the source |
|---|---|---|---|---|---|
| KIC 6448637 (target) | 15.73 | 0 | 7218 | 812 +- 74 (PSF), 811 +- 75 (SAP) | 0.085 % |
| KIC 6448650 | 16.03 | 10.4" (2.6 px) | 5237 | **-75 +- 121** (PSF); 438 +- 325 (SAP) | 0.11 % (1120 ppm) -> excluded ~10 sigma |
| KIC 6448629 | 18.42 | 22" (5.6 px) | 743 | 971 +- 694 | 1.0 % -> excluded ~13 sigma |
| KIC 6448666 | 10.9 | 30" (7.6 px) | 452,500 | 232 +- 199 | 0.0014 % (not testable; see below) |

KIC 6448666 (Kp 11.0, A/F star, DR25 TCE at 0.554 d, 16-49 ppm) is 7.6 px away; 35.587/0.554 = 64.2 is not
a harmonic, and the pixel deficit is compact at the target, 7.6 px from it. Remaining Gaia sources on the
stamp (G 19.4-21.1 at 1.9-5.8 px) would need 2.4-12 % eclipses sitting 1.9-5.8 px from where the deficit
is; excluded by the centroid (> 25 sigma in the nearest case). SIMBAD (astroquery, 2026-10-17): no object
within 40" except TYC 3143-1410-1 (= the Kp 11 star at 30"); no entry for KIC 6448637.

### 1d. Shape, secondary, odd/even, dilution
- Trapezoid chi2 grid on the detrended PSF flux with 29.4-min cadence integration (this review):
  box (T23/T14 = 1) is the best fit; delta chi2 = +0.8 at 0.7, +5.9 at 0.5, +15.8 at 0.3, **+24 for a pure V**.
  V-shape excluded at ~4.9 sigma; T23/T14 >= ~0.55 at 2 sigma. Not a grazing EB. (The re-derivation's
  T23/T14 = 0.70 is not a measurement: delta chi2 < 1 over 0.7-1.0.)
- Secondary: 3-sigma limit ~130 ppm at phase 0.5; nothing at any phase. A blended EB would need a
  primary/secondary eclipse-depth ratio > 7 or be at 2P with equal eclipses (odd/even 0.4-2.0 sigma —
  the classic 2P degeneracy is **not** excluded by Kepler data alone).
- Depth by season (PSF flux, this review): 857 +- 76, 857 +- 81, 728 +- 65, 853 +- 101 ppm for
  seasons 0-3 — consistent (season 2 low by 1.3 sigma). A contaminant at a different position relative
  to the pixel grid each season would be expected to modulate the depth; it does not.
- SAP/PSF depth ratio 0.96 +- 0.07 (re-derivation) with crowdsap 0.94: no extra dilution in the aperture sum.
- Implied Rp (undiluted): Rp/R* = 0.029-0.030 -> **2.4-2.5 Re** for R* = 0.77 Rsun.

What this does not exclude: an unresolved companion/background EB within ~1-2 arcsec (Gaia-incomplete for
delta G > ~3 there; such a source would need a >= 2.5 % eclipse), or an EB at 2P with equal eclipses.
Those are TRICERATOPS / high-resolution-imaging questions, not Kepler-pixel questions.

## 2. Stellar variability / rotation aliases
LS on the out-of-transit detrended flux: strongest peak 0.734-0.743 d with power 8e-4 (32 ppm
half-amplitude, 4 % of the depth). LS is polarity-agnostic; nothing at 35.6 d or its harmonics. No.

## 3. Kepler systematics
- Event timing vs gaps (`vet/events.csv`, re-derivation `events.csv`): of 24 counted events, 2 are within
  0.5 d of a data gap (epoch 23: 0.13 d; epoch 28: 0.40 d after the Q12/Q13 gap). Expected ~0.7 by chance.
  Depth with all events 812 ppm; without epoch 28, 784; without 23 and 28, 794 ppm (this review, PSF flux,
  same box). The signal does not live in gap edges.
- Epoch 28 (BKJD 1183.14, Q13 start = thermal-ramp territory, noisiest quarter): local depth 2268 +- 349 ppm,
  a 4.2-sigma pull; it drives the per-year depth drift (chi2 13/4 -> 4.5/4 without it) and the odd/even
  tension. It should be marked in the evidence pack as a known thermal-ramp-adjacent outlier; the pipeline
  keeps it (gap_frac_1d counts it only as "near a gap"). Dropping it is a post-hoc choice; both numbers
  should be quoted.
- Fold figure (`fold_events.png`): a positive ~3,000-5,000-ppm bump at -8 h in the binned fold comes from
  one high point (epoch 39, Q17; per-event depth 4380 +- 7900 ppm). It is a brightening, not a dip, and
  does not affect the in-transit statistic.
- Edge-of-stamp and bleed adjacency: half of the target's PSF falls off the stamp (psffrac 0.50) and the
  core is 1.7 px from the saturated host's bleed column. Per-quarter depth shows no trend with the
  target's row position (Q5 at row 1.1: 740 ppm; Q2 at row 18.3: 897 ppm; Q15 at row 3.7: 1186; Q11 at
  row 3.7: 592), and excluding bleed-column pixels does not change the deficit (Q11: 511 vs 610;
  Q15: 1228 vs 1254 ppm). Rolling band / monthly thermal cycles have no 35.6-d periodicity, and the
  inverted null shows no comparable box signal.

## 4. Ephemeris-match contamination from a known object on the same channel
TPF KIC 6448574 channel by season (headers): season 0 mod 8 out 3 ch 23; season 1 mod 12 out 3 ch 39;
season 2 mod 18 out 3 ch 63; season 3 mod 14 out 3 ch 47.

**Kepler-1715 b = KOI-308.01 (KIC 6291837, Kp 12.35, P = 35.5969684 d, depth 684 ppm, T14 6.69 h,
CONFIRMED)** lies 18.6 arcmin from the target and is on the **same channel** (its TPF headers: Q2 ch 23,
Q4 ch 63 — identical to the host TPF in those seasons; by the field-rotation symmetry the same holds in
seasons 1 and 3). The periods agree to 2.7e-4. This is the one thing a referee would seize on, so:
- Period: 35.5970 vs 35.58777 +- 0.00045 d (24-event timing fit) -> **20 sigma apart**; 0.0097 d/epoch.
- Timing: our mid-transits are **0.84 d (epoch 0) to 1.20 d (epoch 36) before** the nearest Kepler-1715 b
  transit; minimum |offset| over all 40 epochs is 0.835 d = 4.5 x our 4.5-h duration. The two event
  trains never overlap. Column-anomaly / video-crosstalk contamination preserves timing; this does not
  match. Coincidence in period only.
- Other matches within |dP/P| < 5e-4 and harmonics (KOI/TCE/EB tables, catalogue date 2026-10-17):
  Kepler-1716 b (KIC 7700622, 35.5835 d, 3 deg away, 9.5 sigma, phase offset 2.3 d);
  TCE KIC 8447096 (35.5935 d, 12.7 sigma); P/2 KOIs on KIC 7812167/7812175/7812179 (17.794 d, three KIC ids
  at one position 4.5 deg away, phase offset 0.04 P = 1.4 d); 2P TCEs 7971540 / 9083354 (71.17 d, 2.4-4.7 deg away, phase
  offsets 0.14-0.35 P). None is on-channel-and-on-time. Nearest KOI/EB to the target is KIC 6448768
  (16.4868 d, 127"), no harmonic relation (ratio 2.158).

## 5. Is it already known?
Not in: KOI cumulative (9,564 rows), DR25 TCE (34,032), Kirk+2016 Kepler EB catalogue (2,876), confirmed
planets in the Kepler field (2,839) — none for KIC 6448637, 6448574, 6448650, 6448629; nothing within
3 arcmin except KIC 6448768 (16.49 d). Bienias & Szabo 2025 (A&A 694, A209) appendix B: their 7
hot-Jupiter candidates (P 2.5-4.8 d) use source apertures KIC 11033410, 4543171, 8175131, 8293539,
9040849, 9674230, 4459924 — not this TPF. Moquin et al. 2026 BOI-1.01 is a 3.37-d signal on a G = 17.8
M dwarf — not this. SIMBAD: no entry. "Not in the catalogue" is still not "unreported": the KBonus team's
own search (TraSHDUMP) is unpublished and could contain it.

## 6. Duration / density / implied radius
P = 35.587 d, K dwarf M* ~ 0.75-0.85 Msun, R* = 0.768 Rsun: a/R* = 54-56, T14(b=0) = 5.0-5.2 h,
rho* = 1.7-1.9 rhosun. Observed T14 4.1-4.8 h -> rho*(b=0) upper bound 2.1-3.3 rhosun -> needs
b ~ 0.3-0.6 on this host: consistent (a giant or an F star would be inconsistent; a fainter background
dwarf would not be). Rp = 2.4-2.5 Re undiluted. Nothing here discriminates against a background dwarf
EB, but nothing is inconsistent with the Gaia host.

## 7. Verdict and what is still missing
**survives.** Every test I could run on Kepler data points to a ~850-ppm, 4-5-h, flat-bottomed, periodic
(35.588 d) loss of light from within ~1 arcsec of KIC 6448637, with no secondary, no odd/even difference
beyond the epoch-28 outlier, no season dependence, no gap/edge dependence, no stellar-variability alias, no
neighbour carrying it, no on-channel known object sharing the ephemeris, and no catalogue entry.

Most likely explanation: an on-target transit-like signal of a ~2.5-Re body at 35.59 d on a K dwarf, **or**
an unresolved (< ~1-2 arcsec) faint eclipsing binary / EB at 2P that Kepler pixels cannot separate.
Before "candidate" language:
1. `pipeline/pixel_vet.py` must be rerun on this star in a mode that handles a saturated-host /
   edge-of-stamp TPF (local box around the Gaia position, bleed column excluded) so the localisation is a
   pipeline product, not a red-team scratch result. My scripts and logs are in `redteam/`.
2. TRICERATOPS with the Gaia neighbour list (needs FPP < 0.015, NFPP < 0.001 only in a refereed paper).
3. Quote both depth sets (with / without epoch 28) and the pipeline's conservative SNR = 10.0.
4. Live cross-match of ExoFOP TOI/CTOI, VSX and Gaia NSS, and a check against the KBonus team's list
   when it appears.
Weak points to keep saying out loud: SNR 10 (per-event SNR ~3-4); 2P-EB degeneracy; Gaia incompleteness
inside 1-2 arcsec; half the PSF is off the stamp and 1.7 px from a bleed column.

## Files
- `redteam.md` (this note)
- `redteam/tpf_local.py`, `redteam/tpf_local.log`, `redteam/tpf_local.json`: per-quarter target-box /
  neighbour-box / host-aperture depths and 7x7 deficit-SNR cutouts on the KIC 6448574 TPF.
- `redteam/tpf_radius.py`, `redteam/tpf_radius.log`: radial deficit profile and deficit centroid per quarter.
- Neighbour-light-curve, catalogue-match, per-season, gap-exclusion, LS, trapezoid-grid and density
  numbers above were produced inline in this session (astropy / pandas, no `pipeline/` code);
  their raw outputs are reproduced in the appendix below.

## Appendix: raw outputs of the inline checks (verbatim)
```
neighbour KBONUS LCs at P=35.58726, T0=186.68617, 4.5-h box, per-quarter normalised, 1-d masked running median:
kic-006448629 FLUX    G=18.4 medflux 743.1   depth  971.0 +- 694.2 ppm n_in 40  nev 6   scatter 4362 ppm
kic-006448629 SAP     medflux 170.9          depth -282.6 +- 792.2 ppm n_in 93  nev 14
kic-006448650 FLUX    G=16.0 medflux 5237.1  depth  -75.4 +- 120.6 ppm n_in 40  nev 6   scatter 693 ppm
kic-006448650 SAP     medflux 1127.0         depth  437.5 +- 325.2 ppm n_in 190 nev 26
kic-006448666 FLUX    G=10.9 medflux 452532  depth  232.3 +- 198.8 ppm n_in 215 nev 27
kic-006448666 SAP     medflux 141766         depth   88.3 +- 197.8 ppm n_in 298 nev 37
kic-006448637 FLUX    medflux 7217.7         depth  812.2 +-  73.7 ppm n_in 180 nev 26  scatter 521 ppm
kic-006448637 SAP     medflux 6216.7         depth  810.6 +-  75.0 ppm n_in 182 nev 26

per-quarter box depth (ppm) +- scatter/sqrt(n_in), n_in, season  [PSF FLUX]:
Q2 s0 897+-110 (20)  Q4 s2 712+-146 (13)  Q5 s3 740+-134 (12)  Q8 s2 627+-131 (15)  Q10 s0 908+-180 (7)
Q11 s1 592+-106 (21) Q12 s2 762+-127 (16) Q13 s3 1011+-190 (21) Q14 s0 770+-131 (13) Q15 s1 1186+-123 (17)
Q16 s2 791+-121 (18) Q17 s3 570+-325 (7)
per-season: s0 857+-76 (40)  s1 857+-81 (38)  s2 728+-65 (62)  s3 853+-101 (40)
all 812 ppm (n_in 180); no ep28 784 (174); no ep23,28 794 (168)
LS (detrended, out-of-transit) top periods: 0.7342 d (power 0.00079), 0.7424-0.7426 d (0.00072-0.00076)

trapezoid grid (N=1132 pts within +-0.6 d, sigma 520 ppm, 9x supersampled 29.4-min cadence):
T23/T14=1.0 dchi2=0.00 T14=4.10 h depth=867   0.9: +0.78 (4.40 h)   0.8: +0.18 (4.70 h)   0.7: +0.78 (4.90 h)
0.6: +2.37   0.5: +5.92   0.4: +10.99   0.3: +15.80   0.2: +19.44   0.1: +22.78   0.0: +24.15 (T14 5.7 h, 1226 ppm)

Kepler-1715 b (KOI-308.01) vs ours, offset = ours - nearest KOI transit (d):
epoch 0: -0.835  4: -0.875  8: -0.915  12: -0.956  16: -0.996  20: -1.036  24: -1.076  28: -1.116  32: -1.156  36: -1.197
min |offset| over epochs 0-39 = 0.835 d.  KOI row: P 35.5969684, t0 187.52765 BKJD, dur 6.691 h, depth 684.5 ppm, Kp 12.351,
fpflags 0/0/0/0, sep 1118 arcsec.  TPF headers KIC 6291837: Q2 mod 8 out 3 ch 23; Q4 mod 18 out 3 ch 63.

density: M*=0.75/0.80/0.85 -> a/R* 53.9/55.0/56.2, T14(b=0) 5.19/5.08/4.98 h, rho* 1.66/1.77/1.88 rhosun;
Rp = 2.44 Re (846 ppm) / 2.51 Re (897 ppm).  T14 = 4.13/4.50/4.84 h -> rho*(b=0) = 3.31/2.56/2.05 rhosun.
T23/T14 for k=0.03: b=0 0.94, 0.5 0.92, 0.7 0.89, 0.8 0.85, 0.9 0.72, 0.95 0.49.
```
