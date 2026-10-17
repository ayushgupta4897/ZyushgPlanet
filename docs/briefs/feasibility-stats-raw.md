# Feasibility statistics for a solo archival M-dwarf transit search (Kepler / Kepler Bonus / K2 / TESS FFI)

Date: 2026-10-17. All sandbox commands were run in `/home/user/ZyushgPlanet/.venv` (Python 3.13, astropy 8.0.2, lightkurve 2.6.0, wotan, astroquery 0.4.11). Scratch files: `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/feasibility-stats/` (scripts `*.py`, outputs `*.out`, data `dl/`, tables `*.csv`). Kuber's repo was read (never executed) at `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/untrusted/planethunt`.

## 0. Bottom line (numbers justified below)

| Question | Answer |
|---|---|
| Transiting 1-4 Re planets with P<50 d per 1000 M dwarfs (DC15 rates, 0.5 Msun/0.5 Rsun) | 76 (of which 26 are 2-4 Re; 24 are 1-1.5 Re at P<10 d) |
| Kuber's actual pool (1,836 eligible KBONUS background M dwarfs; median G = 18.5; only 36 stars brighter than G 16, 104 brighter than G 17) | model: 110-120 detectable 1-4 Re planets (P<50 d) if the DC15 rates hold for these faint colour-selected stars, 44 of them 2-4 Re; empirically-discounted expectation 40 +/- a factor 2; his null-dominated first screen harvested 1 candidate + ~5 unresolved leads from 1,736 stars |
| Kepler primary cool dwarfs (4,821 DR25 stars Teff<4000, logg>4) | model 356 detectable vs 129 catalogued KOIs (P<50 d, 1-4 Re) -> the residual ~230 is near-threshold (SNR 7-10) Earths where DR25 completeness is 20-60% and where Kepler's TPS/Robovetter already looked; hardest ground |
| TESS FFI, GCNS M dwarfs within 100 pc (119k with T<~15) | model 1,000-2,500 detectable (2-8 sectors) vs 258 TOIs; empirically discounted ~400-900 unreported detectable planets, 70% at Tmag 12-14 where the official FFI searches (SPOC/QLP, T<13.5) are shallowest; eumemic is already harvesting this pool at 645 counted candidates |
| K2 (13.9k M-dwarf proxies Kp<16.5, 80-d campaigns) | model ~690 detectable vs 132 known cool-host planets; discounted ~200 remaining, most at P<10 d around Kp 14-16 stars |
| Multiple testing | a 0.5-25 d BLS over 4-yr Kepler data has ~6e7 effective trials; white noise alone gives max BLS-SNR ~5.8-6.1 (measured) and SDE ~4.6-5.5; SNR>=9 or SDE>=9.1 would be <1 false alarm per 10k stars in white noise, but real data are red: Kuber's inverted null passed his SNR>=9 first screen on 1.8% of stars (32/1,736), and 2/1,736 inverted stars exceeded his candidate's own scores (per-star FAR 1.2e-3, 90% CI 2e-4 to 3.6e-3). Detection statistics alone therefore do not establish his candidate (expected ~2 equivalent inverted events in his search); localisation/dilution evidence has to carry it |
| Validation end state | Kp 16 / 300 pc (Kuber-type): RV K ~3 m/s on a V~16.6 star is impossible; speckle is possible only at the Gemini faint limit (R>16 needs 30-45 min and 0.1 arcsec = 30 au leaves most companions unexcluded); a 1-m telescope gets ~3.6 sigma per 2,400 ppm transit at best -> realistic end state = well-documented candidate on Zenodo/GitHub (no CTOI route exists for Kepler-field stars), possibly an RNAAS note if pixel-level localisation is airtight. Tmag 12 / 30 pc (Rabtsevich-type): 1.5 ppt transit, J~10, 0.1 arcsec = 3 au -> SG1 1-m photometry, speckle and 1 m/s RVs are all standard -> CTOI -> TOI -> validated/confirmed planet is realistic |
| Compute | BLS (astropy, 1 core): Kepler 4-yr 46k points, 0.5-25 d Ofir grid OS=2 (274k periods x 6 durations) = 124 s; OS=1 = 54 s; 1 quarter = 3.8 s; TESS 2 sectors (7k pts) = 2.8 s; TESS 10 sparse sectors = 54-178 s. 10k Kepler-length stars ~ 150-350 CPU-h; 100k TESS 2-sector stars ~ 80 CPU-h. Kuber measured 18.7 s/star median (download + search + holdout) on 64 stars |

## 1. Occurrence rates and transiting planets per 1000 M dwarfs

### 1.1 Published rates used

| Study | Sample | Headline numbers | Source |
|---|---|---|---|
| Dressing & Charbonneau 2015 (DC15) | 2,543 Kepler M dwarfs (Huber+2014 radii), Q1-Q16 | 2.5 +/- 0.2 planets per M dwarf with 1-4 Re and P<200 d; 0.56 (+0.06/-0.05) Earth-size (1-1.5 Re) and 0.46 (+0.07/-0.05) super-Earths (1.5-2 Re) per star for P<50 d; HZ eta_Earth 0.16-0.24 | [arXiv:1501.01623](https://arxiv.org/abs/1501.01623), Sect. 7.1 and Table 5 (text extracted with pdftotext in the sandbox) |
| Hsu, Ford & Terrien 2020 | 1,746 Kepler DR25 M dwarfs, 89 planet candidates | 4.2 (+0.6/-0.6) [Dirichlet prior] or 8.4 (+1.2/-1.1) [uniform prior] planets per star, 0.5-4 Re, 0.5-256 d; prior-sensitive; fM/fFGK ratios in Table 2 run ~1-5 for 1-2.5 Re at P 2-4 d | [arXiv:2002.02573](https://arxiv.org/abs/2002.02573), Table 2 |
| Ment & Charbonneau 2023 | 363 mid-to-late M dwarfs within 15 pc (0.1-0.3 Msun), TESS S1-42 | 0.61 (+0.24/-0.19) terrestrial planets (>0.5 Re) per star for P 0.4-7 d; 1-sigma upper limit 0.07 planets per star with Rp>1.5 Re at P 0.5-7 d; terrestrials:sub-Neptunes = 14:1 | [arXiv:2302.04242](https://arxiv.org/abs/2302.04242) |
| Hardegree-Ullman et al. 2019 | Kepler M3-M5 V (spectroscopic) | 1.19 (+0.70/-0.49) planets per star, 0.5-2.5 Re, P 0.5-10 d; M5 V: 3.07 (+5.49/-2.49) | [arXiv:1905.05900](https://arxiv.org/abs/1905.05900) |
| Bergsten et al. 2022 (FGK comparison) | Kepler FGK 0.56-1.63 Msun | eta_Earth (HZ) = 15 (+6/-4)%; super-Earth/sub-Neptune transition period scales as M*^(1.7+/-0.2); sub-Neptunes ~2x super-Earths in the HZ | [arXiv:2209.04047](https://arxiv.org/abs/2209.04047) |

Implication: for 1-4 Re planets inside 50 d, M dwarfs host roughly 1.6 planets per star (DC15) versus ~0.5-0.8 for FGK dwarfs (Hsu 2020 Table 2 ratios of 1-5 per bin), and their small radii make each planet 2-4x easier to detect in depth. This is why every October-2026 search targeted M dwarfs.

### 1.2 Model table (my code, `yield_model.py`; DC15 Table 5 cumulative rates, Huber-radii block, converted to a piecewise-constant density in ln P; M* = 0.5 Msun, R* = 0.5 Rsun; transit probability 1.08 R*/a with the Kipping eccentricity factor DC15 used)

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

Totals: 2.48 planets per star (1-4 Re, 0.5-200 d), 83 transiting per 1000 stars; P<50 d and 1-4 Re: 1.59 planets per star, 76.5 transiting per 1000 stars; P<50 d and Rp>2 Re: 0.57 planets per star, 25.8 transiting per 1000 M dwarfs.

Geometry for M* = 0.5 Msun, R* = 0.5 Rsun:

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

The 2.15 h mean duration at 13.77 d matches Kuber's measured 2.1-2.4 h, and a/R* = 0.0892/(0.5 x 0.00465) = 38.4 matches his 39.7; his candidate is geometrically ordinary (a 2.8% transit-probability event, i.e. one expects ~1 transiting 2-4 Re planet with P 10-50 d per 110 M dwarfs).

Caveats: DC15's 1-1.5 Re / 5.5-18.2 d cell is an upper limit (<26.5%) which I absorbed via the cumulative Table 5; Hsu 2020 would roughly double every number; Ment & Charbonneau 2023 implies far fewer >1.5 Re planets around mid/late M dwarfs at P<7 d (<0.07 per star), so the 2-4 Re row is over-generous for the reddest stars in Gaia-selected samples.

## 2. Detectability: noise versus magnitude and minimum detectable radius

### 2.1 Kepler official pipeline (Christiansen et al. 2012 and DR25 stellar table)

Christiansen et al. 2012 define rms CDPP as the effective white noise seen by a transit pulse of a given duration: "A CDPP of 20 ppm for 3-hour transit duration indicates that a 3-hour transit of depth 20 ppm would be expected to have an SNR of 1" ([arXiv:1208.0595](https://arxiv.org/abs/1208.0595), Sect. 2; Table 2 extracted with pdftotext). Their Table 2 (Q3, 6-h rms CDPP, dwarfs with log g>4, median): Kp 9: 32.4 ppm; 10: 38.5; 11: 37.0; 12: 41.7; 13: 54.5; 14: 80.8; 15: 139.1 ppm. The lower envelope is "the minimum noise floor, with contributions from shot noise ... and the typical read noise".

Cool dwarfs are noisier than the all-dwarf medians (activity, fainter in Kp). From the NASA Exoplanet Archive TAP table `Q1_Q17_DR25_KS` (query `tap_cdpp2.py`, 2026-10-17; columns `rrmscdpp03p0/06p0/12p0`, Teff<4000, logg>4; 4,821 stars, 3,876 with Kp<16.5):

| Kp bin | n | median 3-h CDPP | median 6-h CDPP | median 12-h CDPP |
|---|---|---|---|---|
| 11-12 | 68 | 62 | 78 | 91 |
| 12-13 | 62 | 89 | 102 | 120 |
| 13-14 | 182 | 115 | 104 | 92 |
| 14-15 | 585 | 166 | 132 | 113 |
| 15-16 | 2,461 | 300 | 234 | 191 |
| 16-17 | 1,151 | 549 | 417 | 324 |
| 17-20 | 294 | 1,196 | 912 | 690 |

(Bright bins are inflated by active stars; the 3h/6h ratio implies CDPP ~ T^-0.16 to T^-0.3, far from white-noise T^-0.5, so scaling 6-h CDPP to a 2-h M-dwarf transit costs only ~20-40%, not 73%.)

### 2.2 KBONUS PSF photometry (paper claim vs sandbox measurement)

Paper claim: light curves "of comparable quality as those computed by the Kepler pipeline, with CDPP values <100 ppm for sources G<16", PSF photometry up to 40% better than aperture photometry for G>13.25, 606,900 light curves of which 406,548 are new background sources ([arXiv:2310.17733](https://arxiv.org/html/2310.17733v1), abstract, Sect. III.1.2, Fig. 8; MAST page https://archive.stsci.edu/hlsp/kbonus-bkg).

Sandbox measurement (`kbonus_noise2.py`, `kbonus_noise_from_cache.py`): downloaded the stitched KBONUS-BKG light curve of Gaia DR3 2103792461826460288 (= KIC 5689121; file `hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits`, 8.5 MB, 50,182 cadences, quarters 1-4, 6-8, 10-12, 14-16; header KEPMAG 16.33, GMAG 16.28, TPFORG 5689116), per-quarter normalised, wotan biweight (1-d window) detrend, MAD of running boxcar means:

| flux column | median flux_err | point-to-point | 2.3-h CDPP | 3-h | 6-h | 12-h |
|---|---|---|---|---|---|---|
| `flux` (PSF) | 812 ppm | 943 ppm | 571 | 534 | 405 | 285 |
| `sap_flux` | 624 | 762 | 454 | 420 | 315 | 215 |
| `psf_flat_flux` | 823 | 938 | 543 | 507 | 374 | 251 |

Per-quarter 2.3-h CDPP ranges 464-674 ppm (PSF). So at G 16.3 the background-source PSF light curve has a 6-h CDPP of ~400 ppm, 4x the paper's headline "<100 ppm for G<16" (which describes the whole catalogue, dominated by Kepler targets). A fold at Kuber's ephemeris (P 13.769448, T0 BJD 2455630.245431, 2.3-h box) gives 76 transits with data, box depth 1,856 ppm, white-noise SNR 33.7; the CDPP-based expectation (2,165 ppm / 571 ppm x sqrt 76) is 33.0. Both numbers are consistent with his "26.7 BLS SNR" once his red-noise-inflated errors are applied (his aperture local SNR 13.7 uses corrected uncertainties).

Noise versus magnitude from 13 eligible background M dwarfs (12 sampled across G 13.7-18.5 plus the target), PSF flux, 6-h CDPP, compared with the DR25 cool-dwarf median at the same magnitude (`kbonus_noise_vs_G.csv`):

| KIC | G | PSFFRAC | CROWDSAP | quarters | 2.3-h CDPP | 6-h CDPP | ratio to Kepler median at Kp~G |
|---|---|---|---|---|---|---|---|
| 9178867 | 13.74 | 0.93 | 0.82 | 18 | 139 | 90 | 0.89 |
| 11241109 | 13.99 | 0.72 | 1.00 | 9 | 93 | 60 | 0.58 |
| 6545415 | 14.46 | 0.78 | 0.71 | 18 | 343 | 254 | 2.03 |
| 8474897 | 14.58 | 0.88 | 0.77 | 18 | 733 | 644 | 4.97 |
| 11390946 | 15.49 | 0.79 | 0.74 | 13 | 292 | 185 | 0.84 |
| 8025611 | 15.51 | 0.50 | 0.85 | 11 | 273 | 174 | 0.79 |
| 5689121 (target) | 16.28 | 0.67 | 0.75 | 13 | 571 | 405 | 1.15 |
| 10198500 | 16.49 | 0.77 | 0.83 | 18 | 200 | 145 | 0.36 |
| 4570168 | 16.50 | 0.66 | 0.94 | 14 | 532 | 336 | 0.83 |
| 3658797 | 17.50 | 0.70 | 0.81 | 14 | 1,007 | 737 | 0.99 |
| 8250618 | 17.50 | 0.87 | 0.89 | 18 | 935 | 598 | 0.81 |
| (no KIC) Gaia DR3 2126972915714414464 | 18.50 | 0.73 | 0.74 | 8 | 1,483 | 976 | 0.85 |
| 8565344 | 18.50 | 0.70 | 0.71 | 17 | 415 | 270 | 0.23 |

Median ratio 0.84 (range 0.23-4.97). Two lessons: (i) typical KBONUS background M dwarfs are as quiet as Kepler targets of the same magnitude (the 1.5x pessimism one might assume is unwarranted), so a 4-yr baseline makes 1.5-2 Re planets at P<20 d detectable even at G 17-18; (ii) ratios of 0.23-0.36 are physically impossible for an isolated G 18.5 source (photon noise alone would be ~1,000 ppm per 6 h) and flag light curves whose PSF flux is actually dominated by the brighter neighbour. Any transit found in such a light curve is unlocalised by construction. A noise-versus-magnitude sanity test (flag stars with CDPP < 0.5x the Kepler median at their G) should be a mandatory pre-filter in the playbook.

### 2.3 TESS FFI noise

Pre-launch model (Stassun et al. 2018, [arXiv:1706.00495](https://arxiv.org/abs/1706.00495) Eq. on p.24, as implemented in `ticgen` 1.0.4, installed with `uv pip install ticgen`): 1-hour noise = 94 ppm at T 8, 132 at T 9, 203 at T 10, 336 at T 11, 595 at T 12, 1,122 at T 13, 2,285 at T 14, 5,190 at T 15, 13,957 at T 16.

Sandbox check on TESS-SPOC FFI light curves (`tess_noise_measure.py`, 200-s cadence, PDCSAP, biweight 0.75-d, MAD of 1-h means): Kepler-138 (T 11.52, S82): 1-h 506 ppm (model 595); Kepler-169 (T 13.56, S81, CROWDSAP 0.48): 2,205 ppm (model ~1,800); Kepler-42 (T 13.74, S82, CROWDSAP 0.61): 2,408 ppm (model ~2,000). The model is accurate to ~20% for isolated stars and ~25% optimistic for crowded Kepler-field stars. No TESS-SPOC product exists for TIC 120686465 itself (its light curves are TESS FFI cutouts/TGLC), consistent with Kuber's use of TESScut and TGLC.

Mission state: TESS is in Sector 109 (Sector 110 begins 2026-10-31), ~27-day sectors, 13-16 sectors per year ([tess.mit.edu/observations](https://tess.mit.edu/observations/)). SPOC FFI light curves cover up to 160,000 targets per sector ([Caldwell et al. 2020, arXiv:2011.05495](https://arxiv.org/abs/2011.05495)); QLP covers T<13.5 (Huang et al. 2020, [arXiv:2011.06459](https://arxiv.org/abs/2011.06459)). Stars with T>13.5 are the "never searched by the official pipelines" regime (this is where eumemic's 22,448 + 84,863 faint M dwarfs at 11.5<T<=13.5 sit; GUIDE.md records that Chopra's claim "never searched by SPOC" was false for his T<13.5 stars).

### 2.4 K2

EVEREST 2.0 light curves "recover the original Kepler precision to at least Kp = 14, and to at least Kp = 15 for campaigns 1, 5, and 6" ([arXiv:1702.05488](https://arxiv.org/abs/1702.05488)). Model assumption: Kepler CDPP x1.0 to Kp 14.5, degrading linearly to x2.0 at Kp 17; 80-d campaigns, 95% duty cycle. The K2TARGETS table has no Teff; an M-dwarf proxy (J-K 0.78-1.05, Kp-J>2.3, PM>20 mas/yr, long cadence) gives 2,316 stars Kp<14, 6,159 Kp<15, 10,964 Kp<16, 13,891 Kp<16.5 (colour-only without PM: 24,623 Kp<16.5) out of 392,156 K2 long-cadence targets.

### 2.5 Noise table and minimum detectable radius (SNR >= 7.1, >= 3 transits, 0.5 Rsun / 0.5 Msun host; `min_radius_table.csv`)

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

("-" = fewer than 3 transits.) Reading: Kepler-length data are the only archive in which sub-Earths and Earths around faint (Kp 15-17) M dwarfs are detectable; TESS FFIs reach 2-4 Re planets only for T<~13 (2-4 sectors) and sub-Neptunes at T 14-15 only in the continuous viewing zones. Kuber's Rp ~2.7 Re at P 13.8 d around a G 16.3 star is 3x above the KBONUS threshold (0.88 Re), so the lack of brighter, easier detections in his pool is a sample-size effect, not a sensitivity effect.

## 3. Yield: how many unreported detectable planets remain per dataset

Method (`yield_model.py`): for each dataset and magnitude bin, N_stars x [occurrence density x p_transit] integrated over 1-4 Re and P<50 d, weighted by a detection probability. Two modes: STEP (P_det = 1 if SNR >= 7.1 and >= 3 transits; an upper envelope) and REALISTIC (logistic P_det with 50% at SNR 9, 17% at 7.1, 92% at 12, mimicking the DR25 gamma-CDF completeness of Christiansen et al. 2020). "Known" = catalogued planets around Teff<4000 hosts in the same P/Rp window from the Exoplanet Archive TAP (`known_bins.py`, 2026-10-17: `cumulative`, `toi`, `k2pandc`).

Known-planet census around cool hosts (Teff<4000) pulled 2026-10-17:
- KOI cumulative, CANDIDATE+CONFIRMED: 214 KOIs on 143 hosts; by Rp x P: 1-1.5 Re (34/13/3 for <10/10-50/50-200 d), 1.5-2 (23/21/6), 2-4 (17/18/2), >4 (8/5/5); by Kp: <14: 8, 14-15: 24, 15-16: 104, 16-17: 19, >17: 2. P<50 d and Rp>2: 49; P<50 d and 1-2 Re: 94. Median depth 870 ppm, median model SNR 28.
- TOI (excluding FP/FA): 581 TOIs (307 PC, 219 CP, 28 KP, 27 APC); 1-4 Re and P<50 d within 100 pc: 295 (149 with Rp>2); by Tmag for Rp<4: <10: 23, 10-11: 46, 11-12: 113, 12-13: 73, 13-14: 27, 14-15: 7, >15: 0. Median TOI depth 4,700 ppm.
- K2: 223 unique cool-host planets (candidate+confirmed) on 183 hosts; P<50 d, 1-4 Re: 132; host Kp median 14.7.

Star counts:
- Kepler primaries: 4,821 DR25 dwarfs with Teff<4000 (Kp bins as in Sect. 2.1).
- KBONUS: the M-star catalogue (29,800 rows; `hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits`, 2.4 MB) joined to the 606,900-row source catalogue (159 MB) reproduces Kuber's funnel exactly (`kbonus_funnel.py`): 29,800 -> PSFFRAC>0.5: 9,204 -> CROWDSAP>0.7: 6,477 -> >=8 PSF quarters: 2,495 -> kic != tpforg: 1,836. KIC 5689121 is rank 52 by ascending G (G 16.28, PSFFRAC 0.67, CROWDSAP 0.75, 13 quarters). Magnitude distribution of the 1,836: G<14: 2; 14-15: 10; 15-16: 24; 16-16.5: 26; 16.5-17: 42; 17-18: 341; 18-19: 1,391 (median 18.47). The whole M-star catalogue has 2,040 Kepler targets (kic == tpforg) and 27,760 background sources; only 1,246 background sources are brighter than G 17.
- TESS: Gaia GCNS (Smart et al. 2021; table `external.gaiaedr3_gcns_main_1`) M-dwarf candidates (BP-RP>1.84, M_G>8; `gaia_gcns3.py`): within 100 pc 213,665 total; G<13: 17,547; G<14: 46,679; G<15: 87,768; G<16: 131,614. Approximate Tmag (T ~ G-0.7; Stassun+2019 gives T-G ~ -0.9 to -1.1 at BP-RP 2-2.5, so these counts are ~10% conservative): T<12: 11,610; T<13: 36,700; T<14: 74,331; T<15: 118,970. Within 30 pc: 5,025 total, 4,161 with T<~15.

### 3.1 Results (REALISTIC mode; STEP totals in brackets)

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

The known/model ratio of 0.36 is the empirical discount between "occurrence x idealised completeness" and "planets that actually survive detection, vetting and cataloguing" in the best-searched dataset in existence. I apply it below as a realism factor to the other datasets (range 0.36-1.0).

KBONUS eligible background M dwarfs (Kuber's 1,836; PSF photometry; 13 quarters assumed; noise 0.84x Kepler):

| G | n stars | transiting | detectable 1-4 Re | detectable 2-4 Re |
|---|---|---|---|---|
| 13-14 | 2 | 0.2 | 0.2 | 0.1 |
| 14-15 | 10 | 0.8 | 0.8 | 0.3 |
| 15-16 | 24 | 1.8 | 1.8 | 0.6 |
| 16-17 | 68 | 5.2 | 4.9 | 1.8 |
| >17 | 1,732 | 132.5 | 103.8 | 43.8 |
| total | 1,836 | 140 | 111.5 [121] | 46.6 |

With the 1.5x pessimistic noise the total drops to 80-94. Applying the 0.36 Kepler realism factor: ~40 planets (16 with Rp>2 Re), of which only ~3 lie among the 104 stars brighter than G 17 where vetting is tractable. Kuber searched 1,736 of these stars and reported 1 candidate plus unresolved leads (KIC 6629201 at 2.638 d, G 17.9; KIC 12164160 at 19.42 d; KIC 8743279 at 21.96 d; KIC 8803342 at 5.01 d; the latter three at G 18.4-18.8 with 0.3-1% depths). The gap between ~40 expected and ~1-5 harvested is explained by (a) a first screen that passes noise as often as data (Sect. 4), so SNR 9-15 signals cannot be separated from systematics, (b) the >=3-events-in-both-subsets holdout rule, which removes P > ~20-25 d and any star with few quarters, (c) uncertain M-dwarf status of colour-selected G 18-19 sources (many are likely reddened distant K dwarfs or giants for which DC15 rates do not apply), and (d) PSF-flux contamination (Sect. 2.2). The relaxed funnel (PSFFRAC>0.5, CROWDSAP>0.7, background, >=4 quarters) adds only ~500 stars (2,379 total; 133 detectable [147]), so KBONUS background M dwarfs are a small pool: realistic expectation of a few to ~15 vettable new candidates in total, essentially all fainter than G 17.

TESS FFIs, GCNS M dwarfs within 100 pc (dilution 0.85, 3-transit minimum; known = TOIs within 100 pc, 1-4 Re, P<50 d):

| Tmag | n stars | transiting | detectable (2 sec) | detectable (4 sec) | detectable (8 sec) | known TOIs | residual (4 sec) |
|---|---|---|---|---|---|---|---|
| <12 | 11,610 | 888 | 379 | 542 | 694 | 170 | 372 |
| 12-13 | 25,090 | 1,919 | 425 | 663 | 961 | 58 | 605 |
| 13-14 | 37,631 | 2,878 | 198 | 396 | 675 | 24 | 372 |
| 14-15 | 44,639 | 3,414 | 26 | 82 | 204 | 6 | 76 |
| total | 118,970 | 9,100 | 1,029 [1,297] | 1,684 [2,098] | 2,534 [3,058] | 258 | ~1,425 |

With the 0.36 realism factor the residual is ~370 (2 sectors) to ~820 (8 sectors) unreported detectable planets, ~70% at T 12-14, i.e. exactly the faint-FFI regime that SPOC/QLP (T<13.5) cover least well and that eumemic is now mining (9,979 K/M dwarfs and 107,311 faint M dwarfs searched; 645 counted candidates across all samples as of 2026-10-16, per https://github.com/eumemic/exoplanets README). A solo searcher entering this pool today competes with that effort; the differentiators are vetting depth and the fainter/longer-period corners (P 20-50 d needs >=6 sectors; eumemic's 20-100 d search covered only 3,661 K/M dwarfs).

K2 (M-dwarf proxy, 13,891 stars Kp<16.5):

| Kp | n stars | transiting | detectable | detectable 2-4 Re | known (approx. split of 132) | residual |
|---|---|---|---|---|---|---|
| <14 | 2,316 | 177 | 154 | 54 | 40 | 114 |
| 14-15 | 3,843 | 294 | 247 | 89 | 40 | 207 |
| 15-16 | 4,805 | 368 | 216 | 104 | 45 | 171 |
| 16-16.5 | 2,927 | 224 | 74 | 48 | 7 | 67 |
| total | 13,891 | 1,063 | 691 [754] | 295 | 132 | ~560 |

With the 0.36 realism factor: ~200 remaining, mostly P<10 d around Kp 14-16 stars, where the K2 campaigns' 80-d baselines give 8+ transits. K2 M dwarfs have been searched by several groups (K2SFF/EVEREST/Vanderburg, Dressing+2017/2019, Hardegree-Ullman+2020 "Scaling K2"), so the true residual is smaller; but the proxy count is itself uncertain by 2x (colour-only count 24,623), and K2 has no community candidate catalogue beyond the Exoplanet Archive, which makes "unreported" harder to establish than for TESS (CTOI list) or Kepler (KOI/TCE tables).

### 3.2 Assumptions list
1. Occurrence: DC15 Table 5 (Huber radii); density uniform in ln P inside each DC15 period interval and uniform in ln R inside each 0.5-Re radius sub-bin.
2. Host: 0.5 Msun, 0.5 Rsun for every star (GCNS BP-RP>1.84 sample is redder/smaller on average, which lowers occurrence of 2-4 Re planets but raises depth).
3. Transit probability 1.08 R*/a; SNR = depth / CDPP(T14) x sqrt(N_transits) (DC15 Eq. 2); CDPP(T) = CDPP(6 h) x (6/T)^0.16 for Kepler/KBONUS/K2 (DR25 3h/6h ratio), white-noise scaling for TESS.
4. Baselines/duty: Kepler 1,460 d x 0.92; KBONUS 13 quarters (1,170 d) x 0.90; K2 80 d x 0.95; TESS 27.4 d per sector x 0.90; >=3 transits required; P<50 d only.
5. Dilution 1.0 for Kepler/KBONUS/K2 (PSF/PDC deblended), 0.85 for TESS FFI apertures.
6. "Known" counts use Exoplanet Archive host Teff<4000 and exclude FP/FA (TESS) or FALSE POSITIVE (KOI); they are lower bounds on "already found" because TCEs/community candidates not promoted to KOI/TOI are not counted.
7. TESS sector count per star is a scenario parameter (2/4/8); the real distribution after 109 sectors is bimodal (most of the sky 2-5 sectors, CVZs 20-40).

## 4. The multiple-testing problem

Trial counts (`multiple_testing.py`; Ofir-style grid with frequency step dur_min x f / (OS x S), dur_min 1 h; 6 durations; ~P/2h independent phases per period):

| dataset | periods (OS=1 / OS=3) | effective independent trials | expected max SNR, white noise | SNR for <1 FA per 10k stars | P(SNR>=9) per star, white |
|---|---|---|---|---|---|
| Kepler 4 yr, 0.5-25 d | 137k / 411k | 6e7 | 5.99 | 6.97 | 7e-12 |
| Kepler 4 yr, 1-30 d | 119k / 358k | 7e7 | 6.02 | 6.99 | 8e-12 |
| TESS 2 sectors | 5k / 15k | 2e6 | 5.41 | 6.49 | 3e-13 |
| TESS 10 sectors over 5 yr | 177k / 532k | 8e7 | 6.03 | 7.00 | 9e-12 |
| K2 80 d | 7.5k / 22.5k | 3e6 | 5.48 | 6.55 | 4e-13 |

Measured in the sandbox (`bls_timing2.py`, astropy BoxLeastSquares, objective 'snr', pure Gaussian noise): max BLS SNR 5.38 (range 4.96-5.77, 8 realisations) for TESS 2 sectors; 5.43 (4.95-6.17) for one Kepler quarter; 5.82 (5.70-5.99, 3 realisations) for 4-yr Kepler at OS=1 and 6.13 at OS=2; the corresponding SDE ((max-mean)/std of the periodogram) was 3.9-5.5. TLS documentation from >10,000 noise-only runs: SDE 5.7 (FAP 10%), 6.1 (5%), 7.0 (1%), 8.3 (0.1%), 9.1 (0.01%), with the explicit warning that red noise makes these "too optimistic" ([TLS FAQ](https://transitleastsquares.readthedocs.io/en/latest/FAQ.html); Hippke & Heller 2019 [arXiv:1901.02015](https://arxiv.org/abs/1901.02015): SDE ~7 for 1% FPR for both BLS and TLS, TLS recovers 93% vs BLS 76% of injected Earths).

What this means for thresholds: in white noise, SDE >= 9.1 or BLS SNR >= 7 already gives <1 false alarm per 10,000 stars, and a single-star SNR>=9 event has probability ~1e-11. Real light curves are not white: Kuber's inverted (f -> 2-f) run passed his full first screen (BLS SNR>=9, >=3 events in both subsets, holdout SNR>=4, odd/even <3 sigma, no secondary, max single-event fraction <0.5) on 32 signals / 28 stars out of 1,736, versus 31 / 29 for the real data (binomial p = 1.00; no excess). So the first screen's empirical per-star false-alarm rate is 1.8% (90% CI 1.3-2.5%), i.e. ~180 false first-screen passes per 10,000 stars, 10 orders of magnitude above the white-noise expectation. Only an empirical null (inverted or scrambled light curves run through the identical pipeline, including detrending, holdout split and vetting) calibrates this; analytic FAPs are useless for red-noise-dominated Kepler/TESS photometry of faint stars.

Kuber's own numbers (editorial audit 2026-10-13, local file `analysis/editorial-audit-20261013/README.md`): at his candidate's strength (P>1.5 d, BLS SNR>=26.7, holdout SNR>=15.5) two inverted stars (KIC 8701308, KIC 5283912) still pass: per-star rate 1.2e-3 (90% CI 2.0e-4 to 3.6e-3). In a 1,736-star search the expected number of inverted-equivalent events at or above the candidate's scores is ~2 (90% CI 0.4-6.3), so the probability of seeing at least one such score by chance is ~86%. The detection statistics therefore carry essentially no global significance; the audit was right to retract the earlier "2/1,736 as an empirical false-alarm rate" framing. The case for KIC 5689121 rests entirely on the astrophysical tests (per-quarter depth tracking the target's aperture share with chi2 15.4/12 vs p~1e-36 for a contaminant host; raw-TPF centroid delta chi2 = 934 favouring the target; the neighbour's own KBONUS light curve flat at 24 +/- 41 ppm; odd/even and secondary null) plus TRICERATOPS FPP 0.069, which is a candidate, not a validation (<0.015 needed). A 1,736-star inverted null can only bound per-star FARs down to 3/1,736 = 1.7e-3 (95%); demonstrating <1 false alarm per 10,000 stars needs >=30,000 inverted stars with zero passes, which is why eumemic ran 798,769 inverted stars.

Chopra's numbers (https://github.com/paraschopra/tess-planet-candidates/blob/main/docs/GUIDE.md): 743 real vs 752 flipped signals pass the search cuts (SDE>=7, SNR>=6, >=3 transits; p = 0.84, noise-dominated); 31 vs 39 with no hard red flags; 10 vs 17 flag-free (no real excess); the strongest flag-free flipped signal is SDE 8.4 while the two real candidates are SDE 11.3 and 10.3; zero flag-free flipped signals reach SDE 9 in 526 stars, which caps the per-star FAR at 5.7e-3 (95%), i.e. "<3 fakes per survey" exactly as GUIDE.md now states after withdrawing the earlier 0.01 and 0.2 figures. The 38-scramble tests (best fakes SDE 8.5 and 9.3 vs 11.3/10.3) give p < 1/38 = 0.026 per candidate, uncorrected for the 526-star search. Both candidates are therefore ~2-3 sigma-equivalent events as detections and, with TRICERATOPS never completed, remain "weak candidates" by his own label. The flipped maximum of SDE 8.4 over 526 stars, versus the white-noise expectation of ~2.9-sigma-equivalent tails, again quantifies the red-noise inflation (~1.4x in SDE).

eumemic's calibration (README): 798,769 inverted stars scaled to 4,039,203 searched stars predict 96 false alarms among "clean" SNR>=9 signals (19 found before scaling) but 2,462 including signals with a neighbour within 15 arcsec or a centroid offset, against 9,938 real SNR>=9 signals (FA fraction ~25%); at SNR 7.5-9 the FA fraction is ~79% (2,756 vs 3,482). Lesson for the playbook: at SNR 7.5-9 most survivors are noise even after LEO-Vetter-style flags; SNR>=9 plus pixel-level tests is the floor; the inverted null must be large (>=10x the search sample) to say anything at the <1e-3 level.

## 5. Validation feasibility by magnitude (`followup_feasibility.py`)

| | A: Kp 16, M1.5 V, 300 pc (KIC 5689121-like) | B: Tmag 12, mid-M, 30 pc (TIC 4206066-like) |
|---|---|---|
| approximate V / Ic / J | 16.6 / 15.2 / 13.7 (+/-0.3) | 13.6 / 11.6 / 10.0 |
| depth, T14 | 2,378 ppm, 2.4 h (2.7 Re / 0.507 Rsun) | 1,509 ppm, 1.1 h (1.4 Re / 0.33 Rsun) |
| planet mass (Chen & Kipping 2017), RV K | 8.0 Me, 3.4 m/s | 2.6 Me, 2.4 m/s |
| 0.1 arcsec and 1 arcsec projected | 30 au, 300 au | 3 au, 30 au |
| photon noise per transit, 1-m / 0.4-m in Ic (my model: 40% throughput, 1,500 A, 3 arcsec aperture, dark sky) | 435 / 1,088 ppm | 99 / 248 ppm |
| single-transit SNR with assumed systematics floors (500 ppm 1-m, 800 ppm 0.4-m) | 3.6 (1-m), 1.8 (0.4-m) | 3.0 (1-m), 1.8 (0.4-m); 1 h is short, so 2-min cadence helps |
| ground photometry verdict | marginal: needs a 2-m-class telescope or 3-4 stacked 1-m transits; events every 13.8 d at a 2.4-h window and Dec +41 limit opportunities to ~5/yr from one site | routine TFOP SG1 work; 0.3-0.4 m observers quote 1.3-2.6 ppt per 10 min ([arXiv:2107.03480](https://arxiv.org/pdf/2107.03480)), a 1-m does better |
| speckle / AO | Gemini 'Alopeke/Zorro limiting magnitude "r ~18 in good observing conditions" but typical targets R 6-15 and "R > 16 needing 30-45 minutes" ([Scott et al. 2021](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2021.716560/pdf); [Howell & Furlan 2022](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.871163/pdf)); feasible in principle but competitive time, and even a Keck AO curve leaves FPP ~0.027 because companions inside 30 au survive (Kuber's own TRICERATOPS test); CFHT archival image saturates both cores | standard; 0.02-1.2 arcsec = 0.6-36 au excluded, which is what drives TRICERATOPS FPP below 0.015 for nearby hosts |
| RVs | impossible: K ~3 m/s on V 16.6 / J 13.7. MAROON-X's design goal is 1 m/s "for mid M dwarfs as faint as V=16.5 in 30 minutes or less" ([arXiv:1606.07140](https://www.arxiv-vanity.com/papers/1606.07140/)), a design goal on an 8-m telescope; CARMENES-class surveys target 324 bright M dwarfs with 3-4 m/s jitter ([arXiv:1711.06576](https://arxiv.org/abs/1711.06576)). Mass measurement would need ~30-50 h of 8-10 m time for a 3-sigma mass | hard but done routinely for V~13.6 M dwarfs with MAROON-X/ESPRESSO/CARMENES (K 2.4 m/s needs ~20-40 epochs at 1-1.5 m/s) |
| community catalogue route | none for Kepler-field stars (no KOI/TCE submission channel; CTOI requires a TIC and a TESS detection, which does not exist here: TESS folds give -483 to -2,075 ppm) | ExoFOP CTOI submission -> TFOP SG1/SG3 follow-up -> TOI promotion is the designed path ([tess.mit.edu/followup](https://tess.mit.edu/followup/) lists SG1 photometry, SG2 recon spectroscopy, SG3 high-resolution imaging, SG4 precise RVs) |
| realistic end state | candidate with Zenodo DOI + GitHub evidence pack; RNAAS/AJ note possible only with new pixel-level or photometric evidence (AAS82222 was declined); statistical validation is not reachable (FPP 0.069 vs 0.015 required; NFPP undefined without resolved neighbours) | CTOI -> TOI -> statistically validated or RV-confirmed planet with a TFOP co-author team; Rabtsevich's TIC 4206066 (FPP 0.03-0.04, flagged independently by LEO-Vetter and the TESS team) is already on this track |

Rule of thumb for the playbook: a candidate is "closable" by amateurs only if the host is T<~13 and within ~50 pc (SG1 photometry on 0.4-1 m telescopes, speckle in minutes, FPP<0.015 achievable with Gaia + contrast curves). For Kp>15 Kepler-field stars the ceiling is a carefully audited candidate.

## 6. Time, compute and data volume

Measured BLS runtimes (1 CPU, astropy 8.0.2 BoxLeastSquares with 'snr' objective, 6 durations 1-6 h, Ofir log-frequency grid; `bls_timing2.out`):

| light curve | points | periods | wall time |
|---|---|---|---|
| Kepler 1 quarter (90 d), 0.5-25 d OS=3 | 2,867 | 25,349 | 3.8 s |
| TESS 2 sectors 10-min (54 d), 0.5-25 d OS=3 | 6,997 | 15,208 | 2.8 s |
| TESS 10 sectors spread over 5.2 yr, OS=1 / OS=3 | 34,974 | 177k / 532k | 54 s / 178 s |
| Kepler 4 yr (46k points), 1-30 d OS=1 / OS=2 | 45,980 | 119k / 238k | 54 s / 99 s |
| Kepler 4 yr, 0.5-25 d OS=2 | 45,980 | 274k | 124 s |

Warning: `BoxLeastSquares.autoperiod` with its default frequency spacing (min_duration / baseline^2) produces ~5e7 trial periods for a 4-yr baseline and never finishes; always build your own Ofir/log-frequency grid (N ~ OS x S/dur_min x ln(f_max/f_min)).

Published/observed rates: Kuber's pilot median 18.7 s per star (mean 20.0 s; 64 stars; includes KBONUS download, cleaning, 0.5-25 d BLS with up to 2 peaks, holdout vetting; `evidence/codex/database_expansion/kbonus/pilot.log`); his 1,736-star search is therefore ~9-10 CPU-h. Chopra: 526 stars, ~1e6 trial periods each, 1.5 h on 12 cores = ~120 core-s per star (GUIDE.md). eumemic: stack-slide ~9 s per star (coherent 17 s, semi-coherent 4 s), 200 K/M dwarfs in 70 s on a 16-core c8g.4xlarge, Graviton4 2.6x faster than c6i at equal vCPU (README).

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

On the 4-CPU sandbox: 1,836 KBONUS stars at OS=1 = ~7 h wall (plus an inverted run = 14 h); 10k TESS 4-sector stars = ~4 h wall. On a 64-core cloud node (~$2-3/h spot) the 10k Kepler-length job is ~5 h; a 100k-star Kepler-length job is ~2 days and ~$150. Always run the inverted null (doubles cost) and leave budget for injection-recovery (x5-10 for 300 injections, as eumemic did).

Data volume (measured in the sandbox): KBONUS stitched light curve 5-11 MB per star (mean ~8.5 MB; the MAST M-dwarf bundle is 134 GB for 29,800 stars, i.e. 4.5 MB average); Kuber's 1,836 eligible stars = ~15 GB; Kepler PDC long-cadence file 0.47 MB per quarter (4,821 cool dwarfs x 17 quarters = 38 GB, or 2.3 GB for a 300-star pilot); TESS-SPOC FFI light curve 1.1 MB per sector (100k stars x 4 sectors = 440 GB; eumemic reports ~130 GB transferred / ~12 GB stored for 73,226 K/M-dwarf light curves after trimming); TESScut pixel cutouts are 10-50 MB per sector and only needed for the vetting stage. Download, not CPU, is the bottleneck: lightkurve `search_lightcurve` + `download` took ~30 s per KBONUS star here (MAST latency), so a 2,000-star pool is ~17 h of serial downloads unless parallelised or fetched from the bulk bundle.

## 7. Practitioner conclusions for the playbook

1. The arithmetic favours many-star, modest-noise pools: 76 transiting 1-4 Re planets per 1000 M dwarfs (26 with Rp>2 Re) means a search needs >=500 stars to expect a handful of real transiting planets, and >=2,000 to expect several detectable 2-4 Re ones. Kuber's pool was both small (1,836) and faint (median G 18.5), so his expected harvest was ~10-40 planets but his effective harvest after a null-dominated screen was 1 candidate plus a few unlocalised leads.
2. KBONUS background M dwarfs are nearly exhausted as a "bright, vettable" pool: only 104 eligible stars are brighter than G 17 and Kuber searched 1,736 of 1,836. New ground there means either looser catalogue cuts (PSFFRAC/CROWDSAP/quarter count: +500 stars), different detrending (his PSF vs SAP depths differ by 40%), or the non-M-dwarf background sources (406,548 total) at the price of a much lower occurrence rate and no TESS cross-check.
3. The largest residual of detectable, unreported planets is in TESS FFIs at T 12-14 (hundreds, even after a 0.36 realism discount), followed by K2 M dwarfs (~200) and Kepler primaries near threshold (~230 but already Robovetted). Only the TESS pool offers a path to validation (Sect. 5), and eumemic/Chopra are already in it, so differentiate on period range (20-50 d with >=6 sectors), on the 13.5<T<15 band, or on vetting quality.
4. Significance has to be measured, not computed: run the identical pipeline on inverted and on >=38 scrambled copies of every star; report the maximum null score and the per-star pass rate at the candidate's own thresholds (Kuber: 2/1,736; Chopra: 0/526 flag-free at SDE>=9); require the real candidate to exceed the entire null distribution, and size the null at >=10x the search sample if you want to quote a <1e-3 per-star FAR.
5. Sanity-check noise against magnitude before searching: a KBONUS light curve with 6-h CDPP below ~0.5x the DR25 cool-dwarf median at its G (e.g. 145 ppm at G 16.5, 270 ppm at G 18.5) is measuring the neighbour's flux; any transit there is a dilution puzzle, not a candidate.
6. Budget 1-2 CPU-min per Kepler-length star (BLS) and a few seconds per TESS star; the 4-CPU sandbox handles ~2,000 Kepler-length stars overnight including the inverted run; download latency, not CPU, dominates unless you use the bulk bundles.

## Files written (all under `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/feasibility-stats/`)
`yield_model.py` (+ `yield_model.out`, `yield_model_realistic.out`, `min_radius_table.csv`), `multiple_testing.py`, `bls_timing2.py` (+ `bls_timing2.out`), `kbonus_noise2.py` (+ `kbonus_target_noise.csv`), `kbonus_noise_from_cache.py` (+ `kbonus_noise_vs_G.csv`), `kbonus_funnel.py` (+ `kbonus_eligible_background_mdwarfs.csv`, `kbonus_noise_sample.csv`), `tap_dl.py`, `tap_cdpp2.py` (+ `dr25_cdpp_by_kp.csv`, `dl/dr25_cool_dwarfs.csv`), `known_bins.py`, `gaia_gcns3.py`, `tess_noise.py`, `tess_noise_measure.py`, `followup_feasibility.py`, downloaded catalogues/papers in `dl/` (`kbonus_mstars.fits`, `kbonus_source_catalog.fits/.parquet`, `koi_cumulative.csv`, `toi.csv`, `k2pandc.csv`, `christiansen2012.pdf/.txt`, `1501.01623.txt`, `2002.02573.txt`, `2209.04047.txt`, `1905.05900.txt`).