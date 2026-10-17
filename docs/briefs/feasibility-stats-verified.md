# Feasibility statistics for a solo archival M-dwarf transit search (Kepler / Kepler Bonus / K2 / TESS FFI)

Date: 2026-10-17. Verified in `/home/user/ZyushgPlanet/.venv` (Python 3.13, astropy 8.0.2, lightkurve 2.6.0, astroquery 0.4.11, wotan, scipy 1.18, ticgen 1.0.4). Original scratch work: `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/feasibility-stats/` (scripts `*.py`, outputs `*.out`, catalogues in `dl/`, light curves in `lkcache/`). Verification log and new measurements: `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/verify-feasibility-stats/VERIFICATION_LOG.md` and `kbonus_flux_vs_G.txt`. Kuber's repo (`scratchpad/untrusted/planethunt`) was read, never executed. Every number below is marked **[verified]** (re-derived or re-queried by the checker), **[from agent output]** (taken from the research agent's saved output and not re-run) or **[unverified]**.

## 0. Bottom line

| Question | Answer |
|---|---|
| Transiting 1-4 Re planets with P<50 d per 1000 M dwarfs (DC15, 0.5 Msun/0.5 Rsun) | 76.5 (25.8 with Rp 2-4 Re; 23.7 are 1-1.5 Re at P<10 d). DC15's 0.5-1 Re row adds ~25 more sub-Earths per 1000 stars that only Kepler-depth data at Kp<15 can reach **[verified against DC15 Table 5]** |
| Kuber's pool (1,836 eligible KBONUS background M dwarfs; median G 18.47; 36 brighter than G 16, 104 brighter than G 17) | model 111 detectable 1-4 Re planets at P<50 d (47 with Rp>2 Re) if DC15 rates apply; with the magnitude-resolved Kepler realism factor (0.45 at Kp<16, **0.15 at Kp>16**) the expectation is **~15-20**, not ~40; his null-dominated first screen yielded 1 candidate + ~5 unresolved leads from 1,736 stars **[funnel and bins verified; model from agent output]** |
| Kepler primary cool dwarfs (4,821 DR25 stars Teff<4000, logg>4) | model 356 detectable vs 129 catalogued KOIs (P<50 d, 1-4 Re) -> ~230 residual near threshold (SNR 7-10), already Robovetted; hardest ground **[counts verified]** |
| TESS FFI, GCNS M dwarfs within 100 pc (119k-132k with T<~15) | model 1,000-2,500 detectable (2-8 sectors) vs 258 TOIs; BUT these stars are **not unsearched**: TESS-SPOC FFI targets include all d<=100 pc stars with T<=16 and crowding>=0.5, and SPOC has run TPS on them since May 2021. The residual is where SPOC's single-sector search (Sectors 1-55), the 10k-per-CCD cap, crowding<0.5, or d>100 pc leave gaps; eumemic is already there with 645 counted candidates **[selection rules verified]** |
| K2 cool dwarfs | K2TARGETS has EPIC Teff/logg: 27,666 unique long-cadence cool dwarfs with Kp<16.5 (not the brief's 13,891 proxy); model yield roughly doubles to ~1,400 detectable before discount; known 134-150; EPIC parameters for faint cool stars are photometric and unreliable **[counts verified]** |
| Multiple testing | 0.5-25 d BLS over 4 yr of Kepler = ~6e7 effective trials; white noise gives max BLS-SNR 5.8-6.1 (SDE 4.6-5.5); SNR>=9 or SDE>=9.1 would be <1 FA per 10k stars in white noise, but real data are red: Kuber's inverted null passed his SNR>=9 first screen on 1.8% of stars (32/1,736) and 2/1,736 inverted stars beat his candidate's scores (per-star 1.2e-3, 90% CI 2.0e-4 to 3.6e-3; expected ~2 equivalent events in his search, P(>=1) = 0.86). Detection statistics alone cannot establish a candidate; localisation/dilution evidence must **[verified]** |
| Validation end state | Kp 16 / 300 pc (Kuber-type): RV K ~3.3 m/s on V~16.6 is out of reach; speckle only at the Gemini faint limit (R>16 needs 30-45 min; 0.1 arcsec = 30 au); 1-m photometry gives ~3.6 sigma per 2,400-ppm transit; **ExoFOP now accepts Kepler/K2 candidates but only after a refereed paper (or an RNAAS citing a published methodology)** -> realistic end state = refereed note + ExoFOP 'Other Discovery Data Source' listing, otherwise Zenodo/GitHub evidence pack. Tmag 12 / 30 pc (Rabtsevich-type): 1.5-ppt transit, J~10, K ~2.3 m/s, 0.1 arcsec = 3 au -> SG1 1-m photometry, speckle and 1 m/s RVs are standard; publication is still the gate for ExoFOP **[policy verified 2026-10-17]** |
| Compute | BLS (astropy, 1 core): Kepler 4-yr 46k points, 0.5-25 d Ofir grid OS=2 (274k periods x 6 durations) = 124 s; OS=1 = 54 s; 1 quarter = 3.8 s; TESS 2 sectors 2.8 s; 10 sparse sectors 54-178 s **[from agent output]**. Kuber: 18.7 s/star median incl. download (64 stars) **[verified]**. Never call `autoperiod` on multi-year data: default spacing gives 1.0e8 periods for a 1-h minimum duration **[verified]** |

## 1. Occurrence rates and transiting planets per 1000 M dwarfs

### 1.1 Published rates **[all abstracts verified 2026-10-17]**

| Study | Sample | Headline numbers | Source |
|---|---|---|---|
| Dressing & Charbonneau 2015 (DC15) | 2,543 Kepler M dwarfs (Huber+2014 radii), Q1-Q16; 156 planet candidates | 2.5 +/- 0.2 planets per M dwarf with 1-4 Re and P<200 d; 0.56 (+0.06/-0.05) Earth-size (1-1.5 Re) and 0.46 (+0.07/-0.05) super-Earths (1.5-2 Re) per star for P<50 d | [arXiv:1501.01623](https://arxiv.org/abs/1501.01623), abstract, Sect. 7.1, Table 5 |
| Hsu, Ford & Terrien 2020 | 1,746 Kepler DR25 M dwarfs, 89 planet candidates | 4.2 (+0.6/-0.6) [Dirichlet prior, their baseline] or 8.4 (+1.2/-1.1) [uniform prior] planets per star, 0.5-4 Re, 0.5-256 d | [arXiv:2002.02573](https://arxiv.org/abs/2002.02573) |
| Ment & Charbonneau 2023 | 363 mid-to-late M dwarfs within 15 pc (0.1-0.3 Msun, median 0.17), TESS S1-42 | 0.61 (+0.24/-0.19) planets >0.5 Re per star for P 0.4-7 d; 1-sigma upper limit 0.07 per star for Rp>1.5 Re at P 0.5-7 d; terrestrials:sub-Neptunes = 14:1 | [arXiv:2302.04242](https://arxiv.org/abs/2302.04242) |
| Hardegree-Ullman et al. 2019 | Kepler M3-M5 V (spectroscopic) | 1.19 (+0.70/-0.49) planets per star, 0.5-2.5 Re, P 0.5-10 d; M5 V: 3.07 (+5.49/-2.49) | [arXiv:1905.05900](https://arxiv.org/abs/1905.05900) |
| Bergsten et al. 2022 (FGK comparison) | Kepler FGK 0.56-1.63 Msun | eta_Earth = 15 (+6/-4)%; super-Earth/sub-Neptune transition period scales as M*^(1.7+/-0.2); sub-Neptunes ~2x super-Earths in the HZ | [arXiv:2209.04047](https://arxiv.org/abs/2209.04047) |
| Fressin et al. 2013 (false positives) | Kepler candidates | global FP rate 9.4%; Earth-size (0.8-1.25 Re) 12.3%; small Neptunes (2-4 Re) 6.7%; giants 17.7% | [arXiv:1301.0842](https://arxiv.org/abs/1301.0842) |

Implication: for 1-4 Re planets inside 50 d, M dwarfs host ~1.6 planets per star (DC15) versus ~0.5-0.8 for FGK dwarfs, and their small radii make each planet 2-4x easier to detect. Ment & Charbonneau's <0.07 per star for Rp>1.5 Re at P<7 d warns that the 2-4 Re row below is over-generous for the reddest (mid/late-M) stars that dominate colour-selected faint samples.

### 1.2 Model table **[DC15 Table 5 values verified line by line; model arithmetic from `yield_model.py`]**

DC15 Table 5 cumulative rates (Huber-radii block) converted to a piecewise-constant density in ln P; M* = 0.5 Msun, R* = 0.5 Rsun; transit probability 1.08 R*/a (DC15's Kipping 2013 eccentricity factor, Sect. 7).

| Rp | P | planets per star | transiting per 1000 stars | mean p_transit |
|---|---|---|---|---|
| 1.0-1.5 Re | <10 d | 0.212 | 23.7 | 0.111 |
| 1.0-1.5 Re | 10-50 d | 0.350 | 7.5 | 0.021 |
| 1.0-1.5 Re | 50-200 d | 0.327 | 2.3 | 0.007 |
| 1.5-2.0 Re | <10 d | 0.108 | 12.0 | 0.111 |
| 1.5-2.0 Re | 10-50 d | 0.351 | 7.5 | 0.021 |
| 1.5-2.0 Re | 50-200 d | 0.239 | 1.9 | 0.008 |
| 2.0-4.0 Re | <10 d | 0.151 | 16.8 | 0.111 |
| 2.0-4.0 Re | 10-50 d | 0.423 | 9.0 | 0.021 |
| 2.0-4.0 Re | 50-200 d | 0.315 | 2.8 | 0.009 |

Totals: 2.48 planets per star (1-4 Re, 0.5-200 d; matches DC15's 2.5 +/- 0.2), 83 transiting per 1000 stars; P<50 d and 1-4 Re: 1.59 per star, 76.5 transiting per 1000; P<50 d and Rp>2 Re: 0.57 per star, 25.8 per 1000. Not included: DC15's 0.5-1.0 Re row (0.514 per star at P<50 d, 0.597 at P<200 d; injection recovery only 36-54%).

Geometry for M* = 0.5 Msun, R* = 0.5 Rsun **[formula verified]**:

| P (d) | a (AU) | p_transit | mean T14 (h) |
|---|---|---|---|
| 1 | 0.0155 | 0.162 | 0.90 |
| 3 | 0.0323 | 0.078 | 1.30 |
| 5 | 0.0454 | 0.055 | 1.54 |
| 10 | 0.0721 | 0.035 | 1.94 |
| 13.77 | 0.0892 | 0.028 | 2.15 |
| 30 | 0.150 | 0.017 | 2.79 |
| 50 | 0.211 | 0.012 | 3.31 |
| 100 | 0.335 | 0.0075 | 4.17 |
| 200 | 0.531 | 0.0047 | 5.25 |

For Kuber's actual star (M* 0.505, R* 0.507): a = 0.0895 AU, a/R* = 38.0 (his fit 39.7), p_transit = 0.026, mean T14 2.15 h (his 2.1-2.4 h). The candidate is geometrically ordinary: one expects ~1 transiting 2-4 Re planet with P 10-50 d per ~110 M dwarfs.

Caveats: the 1-1.5 Re / 5.5-18.2 d DC15 cell is an upper limit (<26.5%), absorbed via the cumulative table; Hsu 2020 would roughly double every number; the mid/late-M deficit of >1.5 Re planets (Ment & Charbonneau) cuts the 2-4 Re row for the reddest stars.

## 2. Detectability: noise versus magnitude and minimum detectable radius

### 2.1 Kepler official pipeline **[verified]**

Christiansen et al. 2012 define rms CDPP so that a transit of depth equal to the CDPP for its duration has SNR 1 ([arXiv:1208.0595](https://arxiv.org/abs/1208.0595), Sect. 2). Their Table 2 (Q3, 6-h rms CDPP, 0.5-mag bins, median for dwarfs with log g>4): Kp 9: 32.4 ppm; 10: 38.5; 11: 37.0; 12: 41.7; 13: 54.5; 14: 80.8; 15: 139.1 ppm (10th percentiles 12.0-100.3 ppm).

Cool dwarfs are noisier than the all-dwarf medians. NASA Exoplanet Archive TAP table `Q1_Q17_DR25_KS` (the DR25 stellar table; note `q1_q17_dr25_stellar` does not exist), columns `rrmscdpp03p0/06p0/12p0`, `teff<4000 and logg>4`: 4,821 stars (3,876 with Kp<16.5; 4 lack CDPP):

| Kp bin | n | median 3-h | median 6-h | median 12-h | mean 6-h |
|---|---|---|---|---|---|
| 10-12 | 68 | 62 | 78 | 91 | 99 |
| 12-13 | 62 | 89 | 102 | 120 | 221 |
| 13-14 | 182 | 115 | 104 | 92 | 237 |
| 14-15 | 585 | 166 | 132 | 113 | 172 |
| 15-16 | 2,461 | 300 | 234 | 191 | 264 |
| 16-17 | 1,151 | 549 | 417 | 324 | 469 |
| 17-20 | 294 | 1,196 | 912 | 690 | 2,182 |

(Means in the bright bins are 2x the medians because of active stars; a median-based model is optimistic there. The 3h/6h ratio implies CDPP ~ T^-0.16 to T^-0.3, far from white-noise T^-0.5, so scaling 6-h CDPP to a 2-h M-dwarf transit costs ~20-40%, not 73%.)

### 2.2 KBONUS PSF photometry **[paper quotes verified; target noise re-measured independently]**

Paper (Martinez-Palomera, Hedges & Dotson 2023, [arXiv:2310.17733](https://arxiv.org/html/2310.17733v1); MAST HLSP https://archive.stsci.edu/hlsp/kbonus-bkg): 606,900 light curves from Gaia DR3 inputs, 406,548 new background sources; "comparable quality as those computed by the Kepler pipeline, with CDPP values <100 ppm for sources G<16" (abstract; the body says "a large number of sources brighter than G=16 have CDPP values under 100 ppm"); "PSF photometry performs up to 40% better, in CDPP value, compared to aperture photometry for sources fainter than G=13.25" (Summary); M-star catalogue selected from Gaia/WISE/2MASS colours plus Gaia Teff cuts (Sect. III.2; "almost 27,500 new light curves in the M-type dwarf region" plus ~3,700 Kepler-targeted M dwarfs). Header keywords: PSFFRAC = fraction of the expected PSF captured in the TPF; CROWDSAP = target flux / total flux in the aperture.

Measured on the stitched KBONUS-BKG light curve of Gaia DR3 2103792461826460288 = KIC 5689121 (`hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits`, 8.5 MB, 50,182 cadences, quarters 1-4, 6-8, 10-12, 14-16; header KEPMAG 16.332, GMAG 16.282, TPFORG 5689116; extensions LIGHTCURVE_STITCHED + per-quarter LIGHTCURVE_Qn/APERTURE_Qn; columns TIME, CADENCENO, QUARTER, FLUX, FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR, PSF_FLAT_FLUX, PSF_FLAT_FLUX_ERR, SAP_QUALITY, FLATTEN_MASK). Per-quarter normalisation, wotan biweight (1-d window), MAD of running boxcar means:

| flux column | median flux_err | point-to-point | 2.3-h CDPP | 3-h | 6-h | 12-h |
|---|---|---|---|---|---|---|
| `FLUX` (PSF) | 812 ppm | 940 ppm | 571 | 534 | 405 | 285 |
| `SAP_FLUX` | 624 | 760 | 454 | 420 | 315 | 215 |
| `PSF_FLAT_FLUX` | 823 | 939 | 543 | 507 | 374 | 251 |

So at G 16.3 the background-source PSF light curve has a 6-h CDPP of ~400 ppm, 4x the paper's "<100 ppm for G<16" headline (which describes the catalogue as a whole, dominated by Kepler targets) and 1.15x the DR25 cool-dwarf median at the same magnitude. A fold at Kuber's ephemeris (P 13.769448, T0 BJD 2455630.245431) gives 76 transits with data; box depth 1,850-1,970 ppm depending on window definitions; white-noise SNR 26-34 (his "26.7 BLS SNR"; his red-noise-corrected aperture local SNR is 13.7).

Noise versus magnitude from 13 eligible background M dwarfs (12 sampled across G 13.7-18.5 plus the target), PSF flux, 6-h CDPP vs the DR25 cool-dwarf median at Kp~G (`kbonus_noise_vs_G.csv`) **[CSV values re-read; two anomalous stars re-measured]**:

| KIC | G | PSFFRAC | CROWDSAP | quarters | 2.3-h CDPP | 6-h CDPP | ratio to Kepler median | G - Kp implied by raw flux |
|---|---|---|---|---|---|---|---|---|
| 9178867 | 13.74 | 0.93 | 0.82 | 18 | 139 | 90 | 0.89 | -0.01 |
| 11241109 | 13.99 | 0.72 | 1.00 | 9 | 93 | 60 | 0.58 | -0.03 |
| 6545415 | 14.46 | 0.78 | 0.71 | 18 | 343 | 254 | 2.03 | +0.18 |
| 8474897 | 14.58 | 0.88 | 0.77 | 18 | 733 | 644 | 4.97 (SAP 5,648 ppm: variable) | -0.08 |
| 11390946 | 15.49 | 0.79 | 0.74 | 13 | 292 | 185 | 0.84 | -0.19 |
| 8025611 | 15.51 | 0.50 | 0.85 | 11 | 273 | 174 | 0.79 | -0.15 |
| 5689121 (target) | 16.28 | 0.67 | 0.75 | 13 | 571 | 405 | 1.15 | -0.03 |
| 10198500 | 16.49 | 0.77 | 0.83 | 18 | 200 | 145 | 0.36 | **+2.26** |
| 4570168 | 16.50 | 0.66 | 0.94 | 14 | 532 | 336 | 0.83 | -0.06 |
| 3658797 | 17.50 | 0.70 | 0.81 | 14 | 1,007 | 737 | 0.99 | +0.92 |
| 8250618 | 17.50 | 0.87 | 0.89 | 18 | 935 | 598 | 0.81 | +0.17 |
| Gaia DR3 2126972915714414464 | 18.50 | 0.73 | 0.74 | 8 | 1,483 | 976 | 0.85 | +0.94 |
| 8565344 | 18.50 | 0.70 | 0.71 | 17 | 415 | 270 | 0.23 | **+2.30** |

Median ratio 0.84. Lessons: (i) typical KBONUS background M dwarfs are as quiet as Kepler targets of the same magnitude, so a 4-yr baseline makes 1.5-2 Re planets at P<20 d detectable even at G 17-18; (ii) the two "impossibly quiet" stars are measured contamination cases: the Kepler zero point Kp + 2.5 log10(median FLUX in e/s) = 25.31 +/- 0.20 (10 files with valid KEPMAG) implies Kp 14.24 for KIC 10198500 (26,899 e/s at G 16.49) and Kp 16.20 for KIC 8565344 (4,393 e/s at G 18.50): the PSF flux is ~8x too bright and belongs to the originating targets (KIC 10198505 / 8565347); both files also lack KEPMAG (1e+20). Two more stars are 2.3x too bright (+0.9 mag). The target KIC 5689121 passes (-0.03). **Mandatory pre-filter: flag G - Kp_implied > 0.5 and 6-h CDPP < 0.5x or > 3x the DR25 cool-dwarf median at G; a transit in a flagged light curve is a dilution puzzle, not a candidate.**

### 2.3 TESS FFI noise and products **[model verified; products corrected]**

Pre-launch model (Stassun et al. 2018, [arXiv:1706.00495](https://arxiv.org/abs/1706.00495), as coded in `ticgen` 1.0.4 `get_oneHourNoiseLnsigma`, `uv pip install ticgen`): 1-hour noise 94 ppm at T 8, 132 at T 9, 203 at T 10, 336 at T 11, 595 at T 12, 1,122 at T 13, 2,285 at T 14, 5,190 at T 15, 13,957 at T 16.

Measured **[from agent output]** on TESS-SPOC FFI light curves (200-s PDCSAP, biweight 0.75 d, MAD of 1-h means): Kepler-138 (T 11.52, S82) 506 ppm (model 595); Kepler-169 (T 13.56, S81, CROWDSAP 0.48) 2,205 ppm (model ~1,800); Kepler-42 (T 13.74, S82, CROWDSAP 0.61) 2,408 ppm (model ~2,000). The model is accurate to ~20% for isolated stars and ~25% optimistic for crowded Kepler-field stars. **[verified]** TARS light curve of TIC 120686465 itself (S81, TESSMAG 15.21, 200 s, 11,033 points): 1-h noise 5,029 ppm (model ~5,900).

Products for faint stars **[verified 2026-10-17]**: `lk.search_lightcurve('TIC 120686465')` returns 10 HLSP light curves by author **TARS** (Sectors 14, 26 at 1800 s; 40, 41, 53, 54 at 600 s; 74, 75, 80, 81 at 200 s). TARS = TESS All-Sky Rotation Survey (Boyle, Bouma & Mann; DOI 10.17909/2sff-fn29; https://archive.stsci.edu/hlsp/tars): CPM-corrected FFI light curves for ~7.5 million TIC stars with T<16 and parallax>2 mas (d<500 pc), Sectors 1-96, 0.14 MB per star-sector (lightkurve cannot parse them; read with astropy.io.fits: columns TIME, FLUX). No TESS-SPOC or QLP product exists for this star. For M dwarfs T - G is about -1.0 to -1.1 (this star: -1.07), so GCNS G should be converted with T ~ G - 1.0.

Mission state **[verified]**: TESS is in Sector 109 (Orbit 229); Sector 110 begins 2026-10-31 ([tess.mit.edu/observations](https://tess.mit.edu/observations/)); 200-s FFIs since Sector 56 (EM2). TESS-SPOC FFI light curves ([Caldwell et al. 2020](https://arxiv.org/abs/2011.05495); https://archive.stsci.edu/hlsp/tess-spoc): targets are (1) all 2-min targets, (2) "include" targets with H<=10 OR distance<=100 pc, crowding>=0.5 and T<=16, (3) field stars in T order with T<=13.5, logg>=3.5, crowding>=0.8, capped at 10,000 per CCD (152,359-160,000 per sector); "Beginning May 2021, the TESS-SPOC pipeline processing was extended to run planet search on TESS-SPOC FFI targets" (single-sector DV from Sector 36, multi-sector from Sector 56; no multi-sector search of Prime/EM1 FFI data). QLP produces light curves for T<13.5 (Kunimoto et al. 2021/2022 HLSP notes, [arXiv:2110.05542](https://arxiv.org/abs/2110.05542), [arXiv:2211.04386](https://arxiv.org/abs/2211.04386)) and "starting in Sector 41 ... select fainter M dwarfs"; the TESS Faint-star Search vetted 10.5<T<13.5 (Kunimoto et al. 2022, ApJS 259, 33, [arXiv:2112.02176](https://arxiv.org/abs/2112.02176)). **Consequence: a GCNS M dwarf within 100 pc with T<=16 is, cap permitting, a SPOC-searched star; "never searched" applies to d>100 pc faint stars, crowding<0.5 stars, cap-dropped stars in dense fields, and periods beyond single-sector reach in Sectors 1-55.** (Chopra's GUIDE.md records that his "never searched by SPOC" claim for T<13.5 stars was false.)

### 2.4 K2 **[EVEREST quote verified; star counts corrected]**

EVEREST 2.0: "For most K2 campaigns, we recover the original Kepler precision to at least Kp = 14, and to at least Kp = 15 for campaigns 1, 5, and 6" ([arXiv:1702.05488](https://arxiv.org/abs/1702.05488)). Model assumption: Kepler CDPP x1.0 to Kp 14.5, degrading linearly to x2.0 at Kp 17; 80-d campaigns, 95% duty cycle. Star counts: the TAP table `K2TARGETS` has EPIC parameters `k2_teff`, `k2_logg`, `k2_rad`, `k2_mass` (populated for 236,785 of 392,156 long-cadence rows, `k2_lcflag=1`); `k2_teff<4000 and k2_logg>4` gives 47,749 campaign-rows and **27,666 unique EPICs with Kp<16.5** (rows by Kp: 13: 2,323; 14: 6,881; 15: 11,002; 16: 12,229). The brief's colour/PM proxy (J-K 0.78-1.05, Kp-J>2.3, PM>20 mas/yr: 13,891 at Kp<16.5; colour-only 24,623) undercounts by ~2x. EPIC Teff/logg for faint cool stars are photometric (Huber et al. 2016) and include misclassified giants/reddened dwarfs, so 27,666 is an upper bound on genuine M dwarfs.

### 2.5 Noise table and minimum detectable radius (SNR >= 7.1, >= 3 transits, 0.5 Rsun / 0.5 Msun host) **[from agent output `min_radius_table.csv`; inputs verified]**

| mag | Kepler 6-h CDPP | KBONUS 6-h (0.84x) | K2 6-h | TESS 1-h (Tmag=mag) |
|---|---|---|---|---|
| 12 | 90 | 76 | 90 | 595 |
| 13 | 103 | 87 | 103 | 1,122 |
| 14 | 118 | 100 | 118 | 2,285 |
| 15 | 183 | 155 | 220 | 5,190 |
| 16 | 325 | 275 | 520 | 13,957 |
| 17 | 664 | 561 | 1,328 | 48,659 |
| 18 | 912 | 770 | 1,823 | 250,067 |

Minimum detectable planet radius (Re) at P = 5 / 13.77 / 30 d:

| mag | Kepler 4 yr | KBONUS 13 q | K2 80 d | TESS 2 sectors | TESS 4 sectors | TESS 13 sectors |
|---|---|---|---|---|---|---|
| 12 | 0.38 / 0.48 / 0.57 | 0.37 / 0.47 / 0.55 | 0.78 / 0.98 / - | 1.80 / 2.13 / - | 1.51 / 1.79 / 2.04 | 1.13 / 1.33 / 1.52 |
| 13 | 0.41 / 0.51 / 0.61 | 0.40 / 0.50 / 0.59 | 0.83 / 1.05 / - | 2.47 / 2.93 / - | 2.08 / 2.46 / 2.80 | 1.55 / 1.83 / 2.09 |
| 14 | 0.43 / 0.55 / 0.65 | 0.42 / 0.53 / 0.63 | 0.89 / 1.12 / - | 3.53 / 4.17 / - | 2.96 / 3.51 / 4.00 | 2.21 / 2.61 / 2.98 |
| 15 | 0.54 / 0.68 / 0.81 | 0.53 / 0.66 / 0.79 | 1.22 / 1.52 / - | 5.31 / 6.29 / - | 4.47 / 5.29 / 6.02 | 3.33 / 3.94 / 4.49 |
| 16 | 0.72 / 0.90 / 1.08 | 0.70 / 0.88 / 1.05 | 1.87 / 2.35 / - | 8.71 / 10.3 / - | 7.33 / 8.67 / 9.88 | 5.46 / 6.46 / 7.36 |
| 17 | 1.03 / 1.29 / 1.54 | 1.01 / 1.26 / 1.50 | 2.99 / 3.75 / - | 16.3 / 19.3 / - | 13.7 / 16.2 / 18.4 | 10.2 / 12.1 / 13.7 |
| 18 | 1.21 / 1.51 / 1.80 | 1.18 / 1.48 / 1.76 | 3.50 / 4.39 / - | 36.9 / 43.7 / - | 31.0 / 36.7 / 41.8 | 23.1 / 27.4 / 31.1 |

("-" = fewer than 3 transits.) Kepler-length data are the only archive in which Earths and sub-Earths around Kp 15-17 M dwarfs are detectable; TESS FFIs reach 2-4 Re planets only for T<~13 (2-4 sectors) and sub-Neptunes at T 14-15 only in the continuous viewing zones. Kuber's Rp ~2.7 Re at P 13.8 d around a G 16.3 star is 3x above the KBONUS threshold (0.88 Re): the lack of brighter, easier detections in his pool is a sample-size effect, not a sensitivity effect. The SNR here is a CDPP-based matched-filter SNR (Kepler MES-like); astropy BoxLeastSquares `objective='snr'` with red-noise-inflated `dy` can differ by 2x (Kuber: 26.7 vs 13.7), so always state which statistic a threshold refers to.

## 3. Yield: how many unreported detectable planets remain per dataset

Method (`yield_model.py`): N_stars x [occurrence density x p_transit] integrated over 1-4 Re and P<50 d, weighted by a detection probability: STEP (P_det = 1 if SNR >= 7.1 and >= 3 transits; upper envelope) or REALISTIC (logistic, 50% at SNR 9, 17% at 7.1, 92% at 12, mimicking DR25 gamma-CDF completeness). "Known" = catalogued planets around Teff<4000 hosts in the same P/Rp window from the Exoplanet Archive (`cumulative`, `toi`, `k2pandc`).

Known-planet census around cool hosts (Teff<4000), 2026-10-17 **[re-queried via TAP; verified]**:
- KOI `cumulative`, CANDIDATE+CONFIRMED: 214 KOIs on 143 hosts; P<50 d and 1-4 Re: 129; P<50 d and Rp>2: 49; by Kp: <14: 8, 14-15: 24, 15-16: 104, 16-17: 19, >17: 2. Median depth 870 ppm, median model SNR 28.
- `toi` excluding FP/FA: 581 TOIs (PC 307, CP 219, KP 28, APC 27; excluded FP 100, FA 13); 1-4 Re, P<50 d, within 100 pc: 258; by Tmag for Rp<4 (any period/distance): <10: 28, 10-11: 57, 11-12: 114, 12-13: 78, 13-14: 28, 14-15: 7, >15: 0. Median TOI depth 4,700 ppm.
- `k2pandc`: 223 unique cool-host planets (candidate+confirmed); P<50 d, 1-4 Re: 134 (default parameter row; 51 lack a radius there) to 150 (any row).

Star counts **[verified unless noted]**:
- Kepler primaries: 4,821 DR25 cool dwarfs (Teff<4000, logg>4); 3,876 with Kp<16.5.
- KBONUS: joining the M-star catalogue (29,800 rows, `hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits`, 2.4 MB) to the 606,900-row source catalogue (159 MB) reproduces Kuber's funnel: 29,800 -> PSFFRAC>0.5: 9,204 -> CROWDSAP>0.7: 6,477 -> >=8 PSF quarters: 2,495 -> kic != tpforg: 1,836 (-> 1,800 after primary-target exclusion). KIC 5689121 is rank 52 by ascending G. Magnitude distribution of the 1,836: G<14: 2; 14-15: 10; 15-16: 24; 16-16.5: 26; 16.5-17: 42; 17-18: 341; 18-19: 1,391 (median 18.47; 104 with G<17). The whole M-star catalogue has 2,040 Kepler targets (kic == tpforg) and 27,760 background sources; 1,246 background sources are brighter than G 17 **[from agent output]**.
- TESS: Gaia GCNS (Smart et al. 2021; Gaia archive table `external.gaiaedr3_gcns_main_1`) M-dwarf candidates (BP-RP>1.84, M_G>8): within 100 pc 213,665; G<13: 17,547; G<14: 46,679; G<15: 87,768; G<16: 131,614 **[G<16 from agent output]**; G<12.7: 11,610; G<13.7: 36,700; G<14.7: 74,331; G<15.7: 118,970. Within 30 pc: 5,025 total, 4,161 with G<15.7. With T ~ G - 1.0 the "T<15" pool is ~131,614, ~10% above the brief's 118,970.
- K2: 27,666 unique EPIC cool dwarfs with Kp<16.5 (EPIC Teff/logg; upper bound).

### 3.1 Results (REALISTIC mode; STEP totals in brackets) **[from agent output `yield_model_realistic.out`; inputs verified]**

Kepler primary cool dwarfs (official PDC photometry, 4 yr):

| Kp | n stars | transiting 1-4 Re P<50 | detectable 1-4 Re | detectable 2-4 Re | known KOIs | residual |
|---|---|---|---|---|---|---|
| <12 | 68 | 5.2 | 5.2 | 1.8 | 0 | 5.2 |
| 12-13 | 62 | 4.7 | 4.7 | 1.6 | 2 | 2.7 |
| 13-14 | 183 | 14.0 | 14.0 | 4.7 | 6 | 8.0 |
| 14-15 | 585 | 44.7 | 44.7 | 15.1 | 22 | 22.7 |
| 15-16 | 2,464 | 188.5 | 187.0 | 63.6 | 84 | 103 |
| 16-17 | 1,151 | 88.0 | 83.1 | 29.7 | 13 | 70.1 |
| >17 | 294 | 22.5 | 17.3 | 7.4 | 2 | 15.3 |
| total | 4,821 | 368 | 356 [362] | 124 | 129 | ~227 |

The known/model ratio is the empirical discount between "occurrence x idealised completeness" and "planets that survived detection, vetting and cataloguing" in the best-searched dataset. It is **magnitude dependent: 0.45 for Kp<16 (114/255.6) and 0.15 for Kp>16 (15/100.4); 0.36 overall**. Use the bin-appropriate value, not 0.36 everywhere. (The ratio also absorbs KOI-host Teff crossmatch issues and the use of CDPP medians that hide active stars; DC15's own 156 candidates from 2,543 stars = 61 per 1000 stars suggests completeness ~0.6-0.7 for their brighter sample.)

KBONUS eligible background M dwarfs (1,836; PSF photometry; 13 quarters; noise 0.84x Kepler):

| G | n stars | transiting | detectable 1-4 Re | detectable 2-4 Re |
|---|---|---|---|---|
| 13-14 | 2 | 0.2 | 0.2 | 0.1 |
| 14-15 | 10 | 0.8 | 0.8 | 0.3 |
| 15-16 | 24 | 1.8 | 1.8 | 0.6 |
| 16-17 | 68 | 5.2 | 4.9 | 1.8 |
| >17 | 1,732 | 132.5 | 103.8 | 43.8 |
| total | 1,836 | 140 | 111.5 [121] | 46.6 |

With 1.5x pessimistic noise the total drops to 80-94. Applying the magnitude-resolved realism factors (0.45 for G<16, 0.15 for G>16): **~15-20 planets, ~7 with Rp>2 Re**, of which only ~2-3 lie among the 104 stars brighter than G 17 where vetting is tractable. Kuber searched 1,736 of these stars and reported 1 candidate plus unresolved leads (KIC 6629201 at 2.638 d, G 17.9; KIC 12164160 at 19.42 d; KIC 8743279 at 21.96 d; KIC 8803342 at 5.01 d; the last three at G 18.4-18.8 with 0.3-1% depths) **[verified in kbonus-extension README]**. The remaining gap is explained by (a) a first screen that passes noise as often as data (Sect. 4), (b) the >=3-events-in-both-subsets holdout rule (removes P > ~20-25 d and few-quarter stars), (c) uncertain M-dwarf status of colour-selected G 18-19 sources, (d) PSF-flux contamination (Sect. 2.2), and (e) an expected astrophysical FP fraction >=10% (Fressin 2013), much higher for background sources. The relaxed funnel (PSFFRAC>0.5, CROWDSAP>0.7, background, >=4 quarters) adds only ~500 stars (2,379 total; 133 detectable [147]), so KBONUS background M dwarfs are a small pool: a few to ~10 vettable new candidates in total, essentially all fainter than G 17.

TESS FFIs, GCNS M dwarfs within 100 pc (dilution 0.85, 3-transit minimum; known = TOIs within 100 pc, 1-4 Re, P<50 d):

| Tmag | n stars | transiting | detectable (2 sec) | detectable (4 sec) | detectable (8 sec) | known TOIs | residual (4 sec) |
|---|---|---|---|---|---|---|---|
| <12 | 11,610 | 888 | 379 | 542 | 694 | 170 | 372 |
| 12-13 | 25,090 | 1,919 | 425 | 663 | 961 | 58 | 605 |
| 13-14 | 37,631 | 2,878 | 198 | 396 | 675 | 24 | 372 |
| 14-15 | 44,639 | 3,414 | 26 | 82 | 204 | 6 | 76 |
| total | 118,970 | 9,100 | 1,029 [1,297] | 1,684 [2,098] | 2,534 [3,058] | 258 | ~1,425 |

With a 0.36-0.45 realism factor the residual is ~370-1,000 unreported detectable planets, ~70% at T 12-14. **Caveat (corrected): these stars are inside the TESS-SPOC d<=100 pc inclusion rule and have been TPS-searched since May 2021 (single-sector; multi-sector only from Sector 56), so "unreported" must be checked against TESS-SPOC DV/TCE products and QLP faint-M outputs, not only the TOI/CTOI lists.** eumemic is already mining this pool (9,979 K/M dwarfs; 22,448 M dwarfs at 11.5<T<=12.5 and 84,863 at 12.5<T<=13.5; 645 counted candidates as of 2026-10-16). Differentiators: vetting depth, the 13.5<T<15 band, P 20-50 d with >=6 sectors (eumemic's 20-100 d search covered only 3,661 K/M dwarfs), crowding<0.5 and d>100 pc stars, and periods that single-sector SPOC searches of Sectors 1-55 could not reach.

K2 (cool dwarfs; the agent's run used the 13,891-star proxy):

| Kp | n stars | transiting | detectable | detectable 2-4 Re | known (approx. split of 132) | residual |
|---|---|---|---|---|---|---|
| <14 | 2,316 | 177 | 154 | 54 | 40 | 114 |
| 14-15 | 3,843 | 294 | 247 | 89 | 40 | 207 |
| 15-16 | 4,805 | 368 | 216 | 104 | 45 | 171 |
| 16-16.5 | 2,927 | 224 | 74 | 48 | 7 | 67 |
| total | 13,891 | 1,063 | 691 [754] | 295 | 132 | ~560 |

With the EPIC-based count (27,666 unique cool dwarfs Kp<16.5, 70% of them at Kp 15-16.5) the model detectable count roughly doubles to ~1,300-1,400 and the known count is 134-150; with a 0.15-0.45 realism factor the remaining pool is ~200-500, mostly P<10 d around Kp 14-16 stars. K2 M dwarfs have been searched by several groups (Dressing+2017/2019, Hardegree-Ullman+2020 "Scaling K2", Zink+2021, Kruse+2019), so the true residual is smaller; K2 has no community candidate list beyond the Exoplanet Archive, which makes "unreported" harder to establish than for TESS or Kepler.

### 3.2 Assumptions list
1. Occurrence: DC15 Table 5 (Huber radii); density uniform in ln P inside each DC15 period interval and uniform in ln R inside each 0.5-Re radius sub-bin.
2. Host: 0.5 Msun, 0.5 Rsun for every star (redder GCNS/KBONUS samples are smaller: lower 2-4 Re occurrence, deeper transits).
3. Transit probability 1.08 R*/a; SNR = depth / CDPP(T14) x sqrt(N_transits) (DC15 Eq. 2); CDPP(T) = CDPP(6 h) x (6/T)^0.16 for Kepler/KBONUS/K2; white-noise scaling for TESS.
4. Baselines/duty: Kepler 1,460 d x 0.92; KBONUS 13 quarters (1,170 d) x 0.90; K2 80 d x 0.95; TESS 27.4 d per sector x 0.90; >=3 transits; P<50 d only.
5. Dilution 1.0 for Kepler/KBONUS/K2, 0.85 for TESS FFI apertures.
6. "Known" counts use Exoplanet Archive host Teff<4000 and exclude FP/FA (TESS) or FALSE POSITIVE (KOI); they are lower bounds on "already found" because SPOC/QLP TCEs and community candidates not promoted to KOI/TOI are not counted.
7. TESS sector count per star is a scenario parameter (2/4/8); the real distribution after 109 sectors is bimodal (most of the sky 2-5 sectors, CVZs 20-40). `tess-point` (`pip install tess-point`, `tess_stars2px.py`) gives the measured histogram for a target list **[tool not installed here; unverified]**.
8. No astrophysical false-positive term; add >=10% (Fressin 2013) for primaries and compute TRICERATOPS FPP per candidate for background sources.

## 4. The multiple-testing problem **[arithmetic verified; MC from agent output]**

Trial counts (`multiple_testing.py`; Ofir-style grid with frequency step dur_min x f / (OS x S), dur_min 1 h; 6 durations; ~P/2h independent phases per period):

| dataset | periods (OS=1 / OS=3) | effective independent trials | expected max SNR, white | SNR for <1 FA per 10k stars | P(SNR>=9) per star, white |
|---|---|---|---|---|---|
| Kepler 4 yr, 0.5-25 d | 137k / 411k | 6e7 | 5.99 | 6.96 | 7e-12 |
| Kepler 4 yr, 1-30 d | 119k / 358k | 7e7 | 6.02 | 6.99 | 8e-12 |
| TESS 2 sectors | 5k / 15k | 2e6 | 5.41 | 6.49 | 3e-13 |
| TESS 10 sectors over 5 yr | 177k / 532k | 8e7 | 6.03 | 7.00 | 9e-12 |
| K2 80 d | 7.5k / 22.5k | 3e6 | 5.48 | 6.55 | 4e-13 |

Measured (`bls_timing2.py`, astropy BoxLeastSquares, objective 'snr', Gaussian noise): max BLS SNR 5.38 (4.96-5.77, 8 realisations) for TESS 2 sectors; 5.43 (4.95-6.17) for one Kepler quarter; 5.82 (5.70-5.99, 3 realisations) for 4-yr Kepler at OS=1, 6.13 at OS=2; SDE ((max-mean)/std) 3.9-5.5. TLS FAQ (>10,000 white-noise runs): SDE 5.7 (FAP 10%), 6.1 (5%), 7.0 (1%), 8.3 (0.1%), 9.1 (0.01%), with "Real data often has partially correlated (red) noise. Then, the FAP estimates are too optimistic" ([TLS FAQ](https://transitleastsquares.readthedocs.io/en/latest/FAQ.html); Hippke & Heller 2019 [arXiv:1901.02015](https://arxiv.org/abs/1901.02015)).

In white noise, SDE >= 9.1 or BLS SNR >= 7 gives <1 false alarm per 10,000 stars and a single-star SNR>=9 event has probability ~1e-11. Real light curves are not white: Kuber's inverted (f -> 2-f) run passed his full first screen (BLS SNR>=9, >=3 events in both subsets, holdout SNR>=4, odd/even <3 sigma, no secondary, max single-event fraction <0.5) on 32 signals / 28 stars out of 1,736 versus 31 / 29 for the real data (stage table: 64 stars 2/2 vs 0/0; ranks 65-576 10/9 vs 10/9; ranks 577-1800 19/18 vs 22/19; binomial p = 1.00). The first screen's empirical per-star false-alarm rate is 1.8% (90% CI 1.35-2.47%), ~180 false passes per 10,000 stars, ten orders of magnitude above the white-noise expectation. Only an empirical null (inverted or scrambled light curves through the identical pipeline, including detrending, holdout split and vetting) calibrates this.

Kuber's own numbers (`analysis/editorial-audit-20261013/README.md`): at his candidate's strength (P>1.5 d, BLS SNR>=26.7, holdout SNR>=15.5) two inverted stars (KIC 8701308, KIC 5283912) still pass: per-star rate 1.2e-3 (90% CI 2.0e-4 to 3.6e-3); expected inverted-equivalent events in a 1,736-star search ~2 (90% CI 0.36-6.3); P(>=1) = 0.86. The detection statistics carry essentially no global significance, and the audit withdrew the "2/1,736 as an empirical false-alarm rate" framing. The case for KIC 5689121 rests on the astrophysical tests (per-quarter depth tracking the target's aperture share, chi2 15.4/12 vs p~1e-36 for a contaminant host; raw-TPF centroid delta chi2 = 934 favouring the target; the neighbour's own KBONUS light curve flat at 24 +/- 41 ppm; odd/even and secondary null; and, added here, the PSF flux level matching the star's own G) plus TRICERATOPS FPP 0.069, which is a candidate, not a validation (<0.015 needed). A 1,736-star inverted null can only bound per-star FARs down to 3/1,736 = 1.7e-3 (95%); <1 FA per 10,000 stars needs >=30,000 inverted stars with zero passes, which is why eumemic ran 798,769.

Chopra (https://github.com/paraschopra/tess-planet-candidates/blob/main/docs/GUIDE.md) **[verified]**: 743 real vs 752 flipped signals pass the search cuts (SDE>=7, SNR>=6, >=3 transits; p = 0.84); 31 vs 39 with no hard red flags; 10 vs 17 flag-free; strongest flag-free flipped SDE 8.4 vs real candidates TIC 370612437 (SDE 11.3, SNR 8.4) and TIC 252721481 (SDE 10.3, SNR 8.1); zero flag-free flipped SDE>=9 in 526 stars caps the per-star FAR at 5.7e-3 (95%): "Even zero fakes in 526 flipped stars only caps the rate at about 3 per survey (95%)"; the earlier 0.01 and 0.2 figures were withdrawn as "unsupported"; 38 scrambles per candidate (best fakes SDE 8.5 and 9.3) give p < 1/38 = 0.026 per candidate, uncorrected for 526 stars; TRICERATOPS "never completed"; ~1e6 trial periods per star; 1.5 h on 12 cores.

eumemic (https://github.com/eumemic/exoplanets README, 16 Oct 2026) **[verified]**: 798,769 inverted stars (two ~400k samples) scaled to 4,039,203 searched stars predict 96 false alarms among clean SNR>=9 signals (19 found before scaling); summing the README's columns, 2,462 including neighbour-within-15-arcsec and centroid-offset classes vs 9,938 real SNR>=9 signals (~25% FA fraction), and 2,756 vs 3,482 at SNR 7.5-9 (~79%). 645 counted candidates (587 full-sky + 58 earlier). Lesson: at SNR 7.5-9 most survivors are noise even after LEO-Vetter-style flags; SNR>=9 plus pixel-level tests is the floor; the inverted null must be >=10x the search sample to say anything at the <1e-3 level.

## 5. Validation feasibility by magnitude (`followup_feasibility.py`) **[masses/K re-derived; literature attributions corrected; ExoFOP policy verified]**

| | A: Kp 16, M1.5 V, 300 pc (KIC 5689121-like) | B: Tmag 12, mid-M, 30 pc (TIC 4206066-like) |
|---|---|---|
| approximate V / Ic / J | 16.6 / 15.2 / 13.7 (+/-0.3); T = 15.21 | 13.6 / 11.6 / 10.0 |
| depth, T14 | 2,378 ppm, 2.4 h (2.7 Re / 0.507 Rsun) | 1,509 ppm, 1.1 h (1.4 Re / 0.33 Rsun) |
| planet mass (Chen & Kipping 2017), RV K | 7.8 Me, 3.3 m/s | 2.5 Me, 2.3 m/s |
| 0.1 arcsec and 1 arcsec projected | 30 au, 300 au | 3 au, 30 au |
| photon noise per transit, 1-m / 0.4-m in Ic (40% throughput, 1,500 A, 3 arcsec aperture, dark sky) **[agent model]** | 435 / 1,088 ppm | 99 / 248 ppm |
| single-transit SNR with systematics floors (500 ppm 1-m, 800 ppm 0.4-m) | 3.6 (1-m), 1.8 (0.4-m) | 3.0 (1-m), 1.8 (0.4-m) |
| ground photometry verdict | marginal: 2-m class or 3-4 stacked 1-m transits; events every 13.8 d, 2.4-h window, Dec +41 -> ~5 usable/yr from one site | routine TFOP SG1 work; 0.3-0.4 m observers quote 1.3-2.6 ppt per 10 min ([arXiv:2107.03480](https://arxiv.org/pdf/2107.03480)) |
| speckle / AO | Gemini 'Alopeke/Zorro: "limiting magnitude ... r~18 in good observing conditions" ([Scott et al. 2021](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2021.716560/full)); typical targets "R magnitudes in the range of 6-15", exposures "30-45 min for fainter sources (R > 16)", contrast 4-5 mag near the 20-30 mas diffraction limit and 7-10 mag at 1 arcsec ([Howell & Furlan 2022](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.871163/full)); feasible but competitive, and Kuber's own TRICERATOPS test shows a Keck AO curve only reaches FPP ~0.027 because companions inside 30 au survive; CFHT archival image saturates both cores | standard; 0.02-1.2 arcsec = 0.6-36 au excluded, which drives TRICERATOPS FPP below 0.015 for nearby hosts |
| RVs | out of reach: K ~3.3 m/s on V 16.6 / J 13.7. MAROON-X's design goal was "1 m/s for mid M dwarfs as faint as V=16.5 in 30 minutes or less" assuming a 6.5-m telescope ([arXiv:1606.07140](https://arxiv.org/abs/1606.07140)); it is now on 8.1-m Gemini North with a ~30 cm/s instrumental floor ([arXiv:2210.06563](https://arxiv.org/abs/2210.06563)); a 3-sigma mass would still need ~30-50 h of 8-10 m time. CARMENES-class surveys target bright M dwarfs with 3-4 m/s jitter ([arXiv:1711.06576](https://arxiv.org/abs/1711.06576)) | hard but routine for V~13.6 M dwarfs with MAROON-X/ESPRESSO/CARMENES (K 2.3 m/s needs ~20-40 epochs at 1-1.5 m/s) |
| community catalogue route | **ExoFOP (policy of 2026-08-19): "reopened support for community planet candidates from TESS, Kepler, K2, and beyond"; candidates "must first be published in a peer-reviewed journal with online access" (or an RNAAS citing a previously published detection/vetting methodology), with approval via the "Published Candidate Upload Request"; naming TID.nn (TIC ID: 120686465.01), discovery data source "Kepler", >=2 of period/epoch/depth, supporting plots, TOI/candidate cross-check** (https://exofop.ipac.caltech.edu/tess/candidate_help.php). Kepler-field stars also have TESS data (Sectors 14,15,26,40,41,53-55,74,75,80-82) but TESS folds give -483 to -2,075 ppm, so no cTOI (TESS-discovered) claim | same ExoFOP policy applies: the cTOI listing follows publication; TFOP SG1 photometry, SG2 recon spectroscopy, SG3 high-resolution imaging, SG4 precise RVs ([tess.mit.edu/followup](https://tess.mit.edu/followup/)) |
| realistic end state | candidate with Zenodo DOI + GitHub evidence pack; a refereed note (RNAAS/AJ) only with new pixel-level or photometric evidence (AAS82222 was declined), after which an ExoFOP "Other Discovery Data Source" listing becomes possible; statistical validation not reachable (FPP 0.069 vs 0.015; NFPP undefined without resolved neighbours) | refereed paper or RNAAS citing a published methodology -> ExoFOP cTOI -> TOI promotion -> validated/confirmed planet with a TFOP co-author team; Rabtsevich's TIC 4206066 (FPP 0.03-0.04, flagged independently by LEO-Vetter and the TESS team) is on this track |

Rule of thumb: a candidate is "closable" by amateurs only if the host is T<~13 and within ~50 pc (SG1 photometry on 0.4-1 m telescopes, speckle in minutes, FPP<0.015 achievable with Gaia + contrast curves), and even then the ExoFOP listing now requires a refereed publication first. For Kp>15 Kepler-field stars the ceiling is a carefully audited candidate.

## 6. Time, compute and data volume

Measured BLS runtimes (1 CPU, astropy 8.0.2 BoxLeastSquares, 'snr' objective, 6 durations 1-6 h, Ofir log-frequency grid; `bls_timing2.out`) **[from agent output; CPUs shared, +/-30%]**:

| light curve | points | periods | wall time |
|---|---|---|---|
| Kepler 1 quarter (90 d), 0.5-25 d OS=3 | 2,867 | 25,349 | 3.8 s |
| TESS 2 sectors 10-min (54 d), 0.5-25 d OS=3 | 6,997 | 15,208 | 2.8 s |
| TESS 10 sectors over 5.2 yr, OS=1 / OS=3 | 34,974 | 177k / 532k | 54 s / 178 s |
| Kepler 4 yr (46k points), 1-30 d OS=1 / OS=2 | 45,980 | 119k / 238k | 54 s / 99 s |
| Kepler 4 yr, 0.5-25 d OS=2 | 45,980 | 274k | 124 s |

Warning **[verified in astropy source]**: `BoxLeastSquares.autoperiod` uses df = frequency_factor x min(duration) / baseline^2, which for a 1-h minimum duration and a 1,460-d baseline over 0.5-25 d is 1.0e8 trial periods (5e7 for 2 h) and never finishes; build your own Ofir/log-frequency grid (N ~ OS x S/dur_min x ln(f_max/f_min)).

Observed rates: Kuber's pilot median 18.7 s per star (mean 20.0 s; 64 stars; range 9.3-35.4 s; includes KBONUS download, cleaning, 0.5-25 d BLS with up to 2 peaks, holdout vetting; `evidence/codex/database_expansion/kbonus/pilot.log`) **[verified]** -> his 1,736-star search is ~9-10 CPU-h. Chopra: 526 stars, ~1e6 trial periods each, 1.5 h on 12 cores = ~120 core-s per star **[verified]**. eumemic: stack-slide ~9 s per star (coherent 17 s, semi-coherent 4 s), c8g.4xlarge Graviton4 2.6x faster than c6i at equal vCPU **[verified]**.

CPU-hour budget (single-core seconds per star -> hours):

| scenario | s/star | 1k stars | 10k stars | 100k stars |
|---|---|---|---|---|
| TESS 2-4 sectors, BLS OS=3 | 3-6 | 1-2 h | 8-17 h | 80-170 h |
| TESS 10 sectors, OS=1 | 54 | 15 h | 150 h | 1,500 h |
| Kepler/KBONUS 4 yr, OS=1 (137k periods) | 54 | 15 h | 150 h | 1,500 h |
| Kepler/KBONUS 4 yr, OS=2 (274k periods) | 124 | 34 h | 344 h | 3,440 h |
| TLS instead of BLS (typically 5-20x slower) | 300-2,000 | 80-560 h | 800-5,600 h | impractical |
| Kuber pipeline incl. download + holdout | 19 | 5 h | 53 h | 530 h |
| eumemic stack-slide (optimised, multi-sector) | 9 | 2.5 h | 25 h | 250 h |

On the 4-CPU sandbox: 1,836 KBONUS stars at OS=1 = ~7 h wall (14 h with the inverted run); 10k TESS 4-sector stars = ~4 h wall. On a 64-core cloud node (~$2-3/h spot) the 10k Kepler-length job is ~5 h; 100k Kepler-length stars ~2 days and ~$150. Always run the inverted null (doubles cost) and budget injection-recovery (x5-10 for 300 injections, as eumemic did).

Data volume: KBONUS stitched light curve 1.1-11.5 MB per star (mean 7.9 MB over 16 cached files; the MAST M-dwarf bundle is 134 GB for 29,800 stars, i.e. 4.5 MB average **[agent]**); Kuber's 1,836 eligible stars ~15 GB; Kepler PDC long-cadence 0.47 MB per quarter **[agent]** (4,821 cool dwarfs x 17 quarters = 38 GB; 2.3 GB for a 300-star pilot); TESS-SPOC FFI light curve 1.1 MB per sector **[agent]** (100k stars x 4 sectors = 440 GB; eumemic reports ~130 GB transferred / ~12 GB stored for 73,226 K/M-dwarf light curves); **TARS FFI light curves 0.14 MB per star-sector [verified]** (100k stars x 4 sectors = 56 GB, sectors 1-96, T<16); TESScut pixel cutouts 10-50 MB per sector, needed only for vetting. Download, not CPU, is the bottleneck: lightkurve `search_lightcurve` + `download` took ~30 s per KBONUS star (MAST latency), so a 2,000-star pool is ~17 h of serial downloads unless parallelised or fetched from the bulk bundle.

## 7. Practitioner conclusions for the playbook

1. The arithmetic favours many-star, modest-noise pools: 76 transiting 1-4 Re planets per 1000 M dwarfs (26 with Rp>2 Re) means a search needs >=500 stars to expect a handful of real transiting planets and >=2,000 to expect several detectable 2-4 Re ones. Kuber's pool was small (1,836) and faint (median G 18.5); with the faint-star realism factor (0.15) his expected harvest was ~15-20 planets and his effective harvest after a null-dominated screen was 1 candidate plus a few unlocalised leads.
2. KBONUS background M dwarfs are nearly exhausted as a "bright, vettable" pool: only 104 eligible stars are brighter than G 17 and Kuber searched 1,736 of 1,836. New ground means looser catalogue cuts (PSFFRAC/CROWDSAP/quarter count: +500 stars), different detrending (PSF vs SAP depths differ by 40%), the 2,040 Kepler-target M dwarfs' PSF light curves, or the 406,548 non-M background sources at a much lower occurrence rate.
3. The largest residual of detectable, unreported planets is in TESS FFIs at T 12-15, but the GCNS 100-pc stars are TESS-SPOC targets searched by TPS since 2021; differentiate on period (20-50 d with >=6 sectors; Sectors 1-55 had only single-sector SPOC searches), on the 13.5<T<16 band (TARS light curves exist for T<16), on crowding<0.5 or d>100 pc stars, and on vetting quality. K2 cool dwarfs (27,666 at Kp<16.5) are a ~2x larger pool than the brief assumed, mostly at Kp 15-16.5. Kepler primaries near threshold (~230) are already Robovetted.
4. Significance has to be measured, not computed: run the identical pipeline on inverted and >=38 scrambled copies of every star; report the maximum null score and the per-star pass rate at the candidate's own thresholds (Kuber: 2/1,736; Chopra: 0/526 flag-free at SDE>=9); require the real candidate to exceed the entire null distribution; size the null at >=10x the search sample to quote a <1e-3 per-star FAR; and state which SNR statistic you mean (CDPP/MES-like vs BLS 'snr' with inflated errors differ by 2x).
5. Pre-filter every light curve against its magnitude before searching: (a) G - Kp_implied from the raw PSF flux (zero point 25.31 +/- 0.2) must be within +/-0.5 mag, (b) 6-h CDPP within 0.5-3x the DR25 cool-dwarf median at its G (e.g. 145 ppm at G 16.5, 270 ppm at G 18.5 are too quiet), (c) a missing KEPMAG is a warning. Low-noise outliers are neighbour-dominated PSF fluxes; high outliers are variables that dominate false alarms.
6. Budget 1-2 CPU-min per Kepler-length star (BLS) and a few seconds per TESS star; the 4-CPU sandbox handles ~2,000 Kepler-length stars overnight including the inverted run; download latency dominates unless you use bulk bundles or TARS.
7. Decide the end state before searching: under ExoFOP's August-2026 policy every community candidate, TESS or Kepler, needs a refereed publication (or an RNAAS citing a published methodology) before it can be listed; T<13 within 50 pc hosts can then proceed to TFOP SG1/SG3 follow-up and FPP<0.015 validation, while Kp>15 Kepler-field hosts stop at an audited candidate.

## Files
Agent scratch (`/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/feasibility-stats/`): `yield_model.py` (+ `yield_model.out`, `yield_model_realistic.out`, `min_radius_table.csv`), `multiple_testing.py`, `bls_timing2.py` (+ `.out`), `kbonus_noise2.py` (+ `kbonus_target_noise.csv`), `kbonus_noise_from_cache.py` (+ `kbonus_noise_vs_G.csv`), `kbonus_funnel.py` (+ `kbonus_eligible_background_mdwarfs.csv`), `tap_cdpp2.py` (+ `dr25_cdpp_by_kp.csv`, `dl/dr25_cool_dwarfs.csv`), `known_bins.py`, `gaia_gcns3.py`, `tess_noise_measure.py`, `followup_feasibility.py`, catalogues/papers in `dl/`. Checker scratch (`/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/verify-feasibility-stats/`): `VERIFICATION_LOG.md`, `kbonus_flux_vs_G.txt`, `exofop_candidate_help.html`, `exofop_home.html`, `maroonx.txt`, `dl_tars/` (TARS S81 light curve of TIC 120686465).