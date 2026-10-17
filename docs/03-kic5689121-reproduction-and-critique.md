# 03 · KIC 5689121 re-done from scratch: reproduction and red-team review

> Two agents worked on this on 2026-10-17 without access to Kuber's code: one rebuilt the detection end to end in this repo (`pipeline/run.py`, outputs in `runs/kic5689121/`), the other reviewed the candidate with a referee's eye using Kuber's public evidence pack plus the fresh reproduction. Their full reports follow verbatim; machine-readable results are in `docs/briefs/reproduction_and_critique.json`.

**Short version.** The 13.769-day signal is real, repeats 70 times, and the missing light sits on the M dwarf's pixel, not the brighter neighbour's. Kepler-169 e is excluded by CCD geometry, not only by epoch. The neighbour's *official* Kepler light curve even shows the predicted ~80 ppm diluted leak of the dip. What it is not: a validated planet. The dominant remaining scenario, a planet around an unresolved bound companion, cannot be excluded from any archive, and the host is too faint for the imaging that would settle it. The reviewer's "lessons for the developer" section is the design brief for this repo's own search.

---

## Part A · Independent reproduction

# KIC 5689121 end-to-end reproduction (2026-10-17, sandbox)

## 1. Environment
- venv `/home/user/ZyushgPlanet/.venv`, Python 3.13.16: lightkurve 2.6.0, astropy 8.0.2, astroquery 0.4.11, transitleastsquares 2.0.1, wotan 2.0, batman-package 2.5.3, numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2, pandas 2.3.3. `import batman` failed (`ModuleNotFoundError: distutils`, removed in Py3.12+) -> `uv pip install setuptools` (84.0.0) fixes it. Pinned in `/home/user/ZyushgPlanet/requirements.txt`.
- 4 CPUs shared with other agents (load average 5-7.5 during the run), 15 GB RAM.

## 2. Data access (commands actually run)
```python
import lightkurve as lk, astropy.units as u; from astropy.coordinates import SkyCoord
sr = lk.search_lightcurve(SkyCoord(285.907827,40.941988,unit='deg'), author='KBONUS-BKG', radius=20*u.arcsec)
# 17.8 s, 2 rows: target_name 'Gaia DR3 2103792461826460288' (distance 0.0") dataURL
# mast:HLSP/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits
# and 'Gaia DR3 2103792427466721792' -> .../005689116/hlsp_kbonus-bkg_kepler_kepler_kic-005689116_kepler_v1.0_lc.fits
from astroquery.mast import Observations
Observations.download_file(uri, local_path='data/kbonus/2103792461826460288/<name>')  # 8,847,360 bytes, 1.3 s
```
One multi-extension FITS per source holds all quarters (not one file per quarter). `lk.read(file)` works natively -> `KeplerLightCurve`, `time.format='bkjd'`, `scale='tdb'`, flux = FLUX (PSF).

## 3. FITS structure found (saved to `runs/kic5689121/fits_structure.json`)
28 HDUs: PRIMARY (header only) + LIGHTCURVE_STITCHED (50,182 rows) + 13 x (LIGHTCURVE_Q<n>, APERTURE_Q<n>) for Q1,2,3,4,6,7,8,10,11,12,14,15,16.
- PRIMARY: `TIMESYS='TDB'`, `XPOSURE=1765.464`, `KEPLERID=5689121`, `GAIAID='Gaia DR3 2103792461826460288'`, `TPFORG=5689116`, `KEPMAG=16.332`, `GMAG=16.282`, `BPMAG=17.423`, `RPMAG=15.229`, `PARALLAX=3.406`, `QDETECT='011110111011101110'`, `DOI='10.17909/7jbr-w430'`, `HLSPVER='V1.0'`, `ORIGIN='PSFMACHINE'`, `VERSION='1.1.2'`. No BJDREFI/BJDREFF/TIMEREF keywords.
- LIGHTCURVE_STITCHED columns: TIME (TUNIT `'BJD - 2454833'` = **BKJD**), CADENCENO, QUARTER, FLUX, FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR, PSF_FLAT_FLUX, PSF_FLAT_FLUX_ERR, SAP_QUALITY, FLATTEN_MASK. All 50,182 FLUX values finite; median FLUX 3991 e-/s, SAP 4604, PSF_FLAT 3987.
- LIGHTCURVE_Q<n> columns: CADENCENO, TIME, FLUX, FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR, CENTROID_COLUMN, CENTROID_ROW, SAP_QUALITY, SAP_BKG, RED_CHI2, PSF_FLUX_NOVA, PSF_FLUX_ERR_NOVA; header QUARTER, SEASON, MODULE, OUTPUT, CHANNEL, TPF_ROW/COL, SAP='optimal', FLFRCSAP, CROWDSAP, NPIXSAP, PSFFRAC, PERTRATI, PERTSTD, ITERNEG.
- **PSF flux = FLUX/FLUX_ERR** (corrected psfmachine LFD photometry; searched). **Aperture flux = SAP_FLUX** (diluted: 7-pixel quarters Q2/Q6/Q14 have SAP medians ~7100-7500 vs PSF ~4200-4300 e-/s). PSF_FLAT_FLUX = authors' flattened PSF flux; PSF_FLUX_NOVA = mean-shape-model flux.
- **Quality flags**: `SAP_QUALITY` = Kepler TPF QUALITY bitmask (values present: 0, 16, 128, 144, 1024, 2048, 2176, 2192, 8192, 8208, 8320, 8336, 262144, 262160, 262272, 262288, 270336, 270464, 393216, 393344); `FLATTEN_MASK` 0/1 (stitched only). 86-97 % of cadences per quarter have SAP_QUALITY==0.
- Per quarter: PSFFRAC 0.666-0.871, CROWDSAP 0.747-0.960, FLFRCSAP 0.688-0.955, NPIXSAP 2-7, time 131.51-1557.96 BKJD. Q5/9/13/17 missing (dead module 3, season 3).
- **Time conversion**: BJD_TDB = TIME + 2454833.0 (`fetch.BKJD_OFFSET`). Kuber's T0 = BJD_TDB 2455630.2454306 = BKJD 797.24543 = BTJD -1369.75457.

## 4. Pipeline written (`/home/user/ZyushgPlanet/pipeline/`, README there)
`fetch.py` (query/download/cache under `data/kbonus/<gaia_id>/`, `describe_fits`, `load_quarters`), `prep.py` (stitch with SAP_QUALITY==0, per-quarter median normalisation, wotan biweight, upward-only 4-sigma clip, noise metrics, transit mask for the second pass, invert), `search.py` (astropy BLS 1-60 d with P-dependent frequency grid, 210,253 periods, durations 1-6 h, objective='snr'; `bls_refine`; TLS 1-30 d with R=0.507, M=0.505, u=(0.2125,0.3984) from `catalog_info(TIC_ID=120686465)`, use_threads=3; `null_sde_excluding`), `vet.py` (box depth with robust errors, odd/even, secondary + phase scan, per-transit depths, random-phase empirical noise, fixed-ephemeris depth, ephemeris comparison), `report.py` (JSON + 9 PNGs), `run.py` (CLI). `.gitignore` ignores `.venv/ data/ __pycache__/ *.pyc runs/**/*.fits`.

**Detrending window**: for P=13.8 d, R=0.507 Rsun, M=0.505 Msun: a=0.0897 au, a/R*=38.0, T14(b=0,k=0.05)=2.91 h; x1.36 for e<=0.3 at apastron = 3.97 h; 3 x 3.97 h = 11.9 h = **0.496 d** (Hippke et al. 2019, AJ 158, 143, Sect. IV.6/VI: optimal biweight window ~3x transit duration; Sect. V.4: mask known transits).

**Two-pass design (discrepancy found and fixed)**: with the blind 0.496-d biweight, the window containing the ~2000-ppm, 2-h dip is pulled down ~250-300 ppm: measured 1-h-binned shoulders of +200 to +380 ppm within +/-5 h of transit, blind depth 1671 ppm, and the inverted-flux null's top peak in a 13-15 d grid was the shoulder itself (13.7698 d, SNR 7.98, box 2.2-3.7 h after mid-transit). Pass 2 re-detrends with `wotan.flatten(..., mask=transit_mask)` (mask width max(1.5 x box, 3 h)); depth becomes 2058 ppm and the null peak moves to 14.8173 d, SDE 4.84.

## 5. Commands run
```
. .venv/bin/activate
python -W ignore -m pipeline.run --gaia 2103792461826460288 --ra 285.907827 --dec 40.941988 \
  --neighbour 2103792427466721792 --out <scratch>/repro/smoke2 --no-tls --bls-pmin 13 --bls-pmax 15 \
  --compare-period 13.769448251 --compare-t0-bjd 2455630.2454306 --compare-duration-h 2.3     # 54 s, complete
nohup python -W ignore -m pipeline.run --gaia 2103792461826460288 --ra 285.907827 --dec 40.941988 \
  --neighbour 2103792427466721792 --out runs/kic5689121 --compare-period 13.769448251 \
  --compare-t0-bjd 2455630.2454306 --compare-duration-h 2.3 > runs/kic5689121/stdout.txt 2>&1 &   # started 11:59:34 UTC
python scratch/tls_time.py 12 16 3 1.1     # TLS timing on the pass-1 LC, 12-16 d
python scratch/bls_full.py                 # BLS 1-60 d timing/peaks on the pass-1 LC
```
At the time this output was forced (12:03:43 UTC) the full run had finished pass-1 prep and the 1-60 d BLS and was inside the target TLS; it will write `runs/kic5689121/results.json` and the nine PNGs when done (~20-30 min total). The smoke2 directory contains a complete BLS-only run with all figures and identical vetting code.

## 6. Numbers
**Light curve**: 46,091 points after SAP_QUALITY==0 + finite + err>0 (from 50,182 raw) in 13 quarters, baseline 1426.4 d; 432 points removed by the upward-only 4-sigma clip (964 ppm robust sigma). Noise: std 1031 ppm, robust sigma 952 ppm, **point-to-point 962 ppm**, median formal error 821 ppm; **6-h CDPP-like 259 ppm (binned robust) / 263 ppm (lightkurve estimate_cdpp, 12 cadences)**, white-noise expectation 278 ppm. Raw flux shows a stable ~1.7 % peak-to-peak rotational modulation at P_rot = 10.42 d (Lomb-Scargle; 10.2-11.2 d in every quarter), not harmonically related to 13.77 d.

**BLS (blind, pass 1, 1-60 d, 210,253 periods x 8 durations, 134.7 s)**: best **P = 13.769435 d, power (depth SNR) 35.25, SDE 33.21**; refined P = 13.769407 d, T0 = 136.31283 BKJD (BJD_TDB 2454969.31283), duration 2.0 h, depth 1685 +/- 48 ppm (formal). Scratch full-grid run (211,097 periods, 143.6 s): 13.769479 d, SDE 34.9; next distinct peaks 41.308 d (=3P, SDE 21.1), 4.5898 d (=P/3, 17.2), 55.078 d (=4P, 16.0), 1.967 d (13.4), 3.442 d (13.3) - i.e. the signal and its harmonics dominate.

**TLS (12-16 d, osf=3, dgs=1.1, 10,651 periods, 17.5 s, 1.64 ms/period on 3 threads)**: **P = 13.769272 d, SDE 48.70, SNR 25.78**, model depth 1438 ppm, duration 1.71 h, T0 = 136.32142 BKJD (BJD_TDB 2454969.32142), 104 epochs, 74 distinct transits with data. Projected 1-30 d grid: 180,843 periods -> ~300 s (running). TLS grid sizes: 1-30 d osf3 180,843; 1-60 d osf3 198,548; 0.5-25 d osf3 244,777.

**Pass 2 (masked re-detrend) BLS refine**: P = 13.769540 d, T0 = 136.30866 BKJD, 2.0 h, **depth 2071 +/- 48 ppm (SNR 43.0)**, odd 2077 / even 2064 ppm, phase-0.5 depth -0 +/- 48 ppm, 74 transits with data.

**Vetting (pass 2, robust errors)**: box depth **2058 +/- 67 ppm**; odd 2040 +/- 103 (n_in 130), even 2074 +/- 86 (n_in 147), **difference -0.26 sigma**; **secondary at phase 0.5: -62 +/- 60 ppm**, largest other dip 1.7 sigma at phase 0.175; **73 individual events, 72 with positive depth**, weighted mean 2076 +/- 55 ppm, one negative outlier (-700 ppm near BKJD 645); **empirical noise** from 400 random-phase folds: std 65 ppm -> **empirical SNR 31.6** (formal 31.1 from the robust error). Rp/R* = 0.045 -> 2.5 Re for 0.507 Rsun; Teq ~ 417 K.

**Inverted null (flux -> 2 - flux, pass-2 LC)**: 13-15 d BLS best P = 14.8173 d, power 7.92, **SDE 4.84** (vs 33-35 real). Full 1-60 d BLS + 1-30 d TLS nulls are in the running job (`inverted_null` key, plus `*_excluding_candidate_band`).

**Comparison to Kuber's ephemeris (P 13.769448251 d, T0 BJD_TDB 2455630.2454306)**: dP = +7.9 s (6.7 ppm) for the pass-2 refine, +0.8 s for the smoke refine, -1.1 s for the full-grid BLS, -15 s for TLS (grid steps ~5 s); T0 offset at mid-baseline epoch: 2.1 min; depth in a 2.3-h box at his exact ephemeris: 1887 +/- 69 ppm.

**Neighbour KIC 5689116 (Gaia DR3 2103792427466721792)**: 7.09" away, G 15.478 (2.10x brighter), parallax 0.302 mas, BP-RP 0.92 (background star). Own LC: 46,445 points, p2p 479 ppm, 6-h CDPP 105 ppm; BLS 13-15 d best 14.5957 d, SDE 5.51, 250 +/- 36 ppm at 1 h (noise). **Depth at the target ephemeris: -24 +/- 29 ppm**, versus ~981 ppm required to host the signal.

## 7. Comparison table (ours vs Kuber)
| quantity | this pipeline | Kuber |
|---|---|---|
| P (d) | 13.769435 (BLS) / 13.769272 (TLS) / 13.769540 (pass-2 refine) | 13.769448 +/- 0.000065 |
| T0 | BKJD 136.309-136.321 (offset 2.1 min at mid-baseline) | BJD_TDB 2455630.245431 |
| depth (ppm) | 1671 blind; 2058 +/- 67 masked; 2076 +/- 55 per-event; 1887 +/- 69 at his ephemeris (2.3 h) | 1438 +/- 105 (aperture) to 2165 +/- 98 / 2348 +/- 128 (PSF); 1830-2430 by box width |
| duration (h) | 2.0 (BLS box), 1.71 (TLS), ~2.3 visible | 2.1-2.4 |
| transits | 73-74 with data | 70 usable |
| SNR | BLS 35.3; empirical 31.6; TLS 25.8 | BLS 26.7; local 13.7 |
| odd/even | -0.26 sigma | agree |
| secondary | -62 +/- 60 ppm | none |
| neighbour | -24 +/- 29 ppm (needs 981) | 24 +/- 41 (needs ~1100) |

## 8. Figures
Final: `/home/user/ZyushgPlanet/runs/kic5689121/01_raw_stitched.png ... 09_per_transit.png` (written when the background job ends). Verified renderings from the identical code: `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/repro/smoke2/{01_raw_stitched,02_detrended,03_bls_periodogram,05_phase_fold,06_odd_even,07_secondary,08_neighbour_fold,09_per_transit}.png` (05 shows a clean ~2.3-h, ~2000-ppm dip with 30-min bins; 06 odd and even folds identical; 09 per-event depths flat at ~2100 ppm).

## 9. Discrepancies / caveats
1. Blind depth is biased low by ~20 % by the biweight (fixed by the masked second pass); any pipeline quoting a depth from a single wotan pass will under-report it and will see a spurious inverted-null peak at the true period.
2. astropy BLS `objective='likelihood'` is sign-agnostic; nulls must use `objective='snr'`.
3. TLS depth (1438 ppm) is a limb-darkened model minimum with a 1.71-h duration, not a box depth; do not compare it to Kuber's 1438 +/- 105 aperture depth (coincidence).
4. Formal +/-48 ppm errors are ~1.4x too small versus the random-phase estimate (65 ppm).
5. Full-grid TLS, inverted BLS/TLS over 1-60/1-30 d, and the neighbour's 1-60 d BLS are from the still-running job; values quoted are from measured narrower grids of the same code.
6. Nothing from the untrusted planethunt clone was executed or imported. Nothing was committed to git.

---

## Part B · Red-team review

# Red-team review: the KIC 5689121 13.77-day candidate and the method behind it

Reviewer stance: a referee's eye. Everything below is either (E) established by a measurement I ran in the sandbox or read from a primary source, (K) taken from Kuber Mehta's evidence repo (read-only, nothing executed), or (S) my inference/speculation, labelled.

Scratch analyses live in `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/critique/` (`q1.py` archive TAP, `q2.py` Gaia/TIC, `q3.py` Kepler LC headers, `a1.py` PSF/SAP/PDCSAP depths, `a2.py`/`a4.py` difference images + pixel geometry, `a5.py` per-event depths, `a7.py` TARS fold, `a8.py` rotation). The independent reproduction is in `/home/user/ZyushgPlanet/runs/kic5689121/results.json` (finished 12:13 UTC, 13.8 min).

## 0. Bottom line

The dip is real, periodic, transit-shaped and sits on (or within ~2") of the M dwarf. It is not Kepler-169 e leaking through the detector, and it is not the 7.1" neighbour. It is a **planet candidate** with a ~7 % false-positive probability dominated by a scenario (planet around an unresolved bound companion) that archival data cannot touch. The parts of Kuber's case that looked weakest to the editor (null-control parity, method-dependent depth, TESS non-detection, Barber flag) are either explained or are honest limitations rather than red flags; the parts that are genuinely weak are the host ambiguity, the host characterisation, and the fact that the *search* has no discriminating power at its thresholds.

## 1. Is the signal real and on target?

### 1.1 Reality of the signal (E)

* Blind recovery by the sandbox pipeline (`results.json`): BLS 1-60 d, 210,253 periods, best **P = 13.769435 d, power 35.2, SDE 33.2**; TLS 1-30 d, 180,849 periods, **P = 13.769643 d, SDE 127, SNR 25.5**, 74 distinct transits. Inverted-flux null on the same star: BLS best power 9.4 / SDE 5.15 (at 48.3 d), TLS SDE 8.2 — nothing at 13.77 d once the transits are masked during detrending.
* Pass-2 (transit-masked biweight) vetting: 2.3-h box depth **1916 +/- 68 ppm** (robust) / 1932 +/- 45 (formal); odd 1895 +/- 104 vs even 1936 +/- 84 (**-0.31 sigma**); phase-0.5 depth **-90 +/- 54 ppm**; largest other dip in a 195-phase scan 2.5 sigma; 400 random-phase folds give std 61 ppm (**empirical SNR 31**). Depth at Kuber's exact ephemeris: 1879 +/- 69 ppm.
* My per-event analysis (`a5.py`): 71 events with coverage, **70 positive**, weighted mean 1921 ppm, per-event scatter 867 ppm vs median formal error 512 ppm (red-noise factor **1.69**), empirical event-based SNR 18.3; first half 2023 +/- 145 vs second half 1749 +/- 147 ppm.
* **Not a perturbed-PSF artefact**: the KBONUS mean-PSF column `PSF_FLUX_NOVA` gives 2047 ppm vs 2000 ppm for the corrected `FLUX` (2.0-h box, pooled) — the check the Kepler Bonus paper explicitly recommends ("default to the photometry fitted with the mean PSF model" if the perturbed model misbehaves; arXiv:2310.17733).
* **Independent-pipeline confirmation in the pixels** (new): I downloaded the 13 official Kepler long-cadence light curves of the TPF host **KIC 5689116** (Kp 15.48, CDPP6 130-177 ppm) and folded `pdcsap_flux` at Kuber's ephemeris after a masked biweight: **80 +/- 26 ppm** (SAP: 84 +/- 25), random-phase std 23 ppm (3.4-3.6 sigma), phase 0.5: -13 ppm. That is the diluted leak of a ~2000-ppm dip on a star 1.8 px away that contributes a few percent of the flux in KIC 5689116's optimal aperture (CROWDSAP 0.93-0.98). It also explains why the Kepler pipeline never produced a TCE: predicted MES ~ 80 / (140 ppm x sqrt(6/2.3)) x sqrt(70) ~ 3, well below 7.1. NASA Exoplanet Archive TAP confirms **no DR24 or DR25 TCE on KIC 5689116** (`q1_q17_dr24_tce`, `q1_q17_dr25_tce`).

### 1.2 The 'method-dependent depth' is two different things (E)

Pooled PSF box depth vs box width (masked detrend, Kuber ephemeris): 1.5 h **2147 +/- 83**, 2.0 h 2000 +/- 77, 2.3 h 1846 +/- 77, 2.6 h 1707 +/- 79, 3.0 h 1501 +/- 76 ppm. Pooled SAP: 1638 -> 1119 ppm over the same widths; **SAP 2.3 h = 1394 +/- 71 ppm**, which is Kuber's "corrected aperture 1438 +/- 105". So:

* 2165 (local PSF fit), 2117 +/- 190 (MCMC k^2, `kfit_psf.json`), 1830-2430 (box-width/baseline range) are one quantity measured with different effective widths;
* 1438 is the **diluted** SAP depth ("corrected" refers to corrected *uncertainties*, see Kuber's `full-archive-review.md`), not a competing depth.
A referee should still object to quoting any of these as "the" depth without a limb-darkened model and a propagated R*; Rp is 2.4-3.1 Re depending on choices, and larger if the star is young (below).

### 1.3 Blend with the neighbour, and dilution (E)

Per-quarter table (`a1.py`, 2.3-h box):

| Q | season | npix | target share FLFRCSAP*PSF/SAP | PSF depth | SAP depth | SAP/PSF |
|---|---|---|---|---|---|---|
| 2,6,10,14 | 0 | 7,7,4,7 | 0.55-0.97 | 1341-2128 | 543-947 | 0.36-0.65 |
| 3,7,11,15 | 1 | 5,6,5,5 | 0.66-0.89 | 1491-2515 | 970-2234 | 0.65-0.91 |
| 4,8,12,16 | 2 | 4,4,3,2 | 0.96-1.03 | 1285-2464 | 1088-2198 | 0.85-1.02 |

Season-pooled: share 0.66 -> SAP/PSF 0.49; share 0.83 -> 0.85; share 1.01 -> 0.92. SAP depth tracks the target's share of aperture light (my crude fit: scale 0.85; Kuber's with per-event errors: alpha 0.91, chi2 15.4/12, p = 0.22; contaminant-hosted model p ~ 1e-36, `kdilution_pooled.json`). The 7-pixel quarters contain the most neighbour light and show the *shallowest* SAP dips — the signature of a target-hosted signal.

Neighbour KIC 5689116's own KBONUS photometry at the ephemeris: **PSF -27 +/- 31, NOVA -33 +/- 30, SAP +63 +/- 31 ppm** (the SAP leak is expected); hosting the dip would require >~950-1100 ppm. Reproduction: -31 +/- 28 ppm, needs ~913.

Difference images (`a4.py`), all 13 TPFs of KIC 5689116, per-event local out-minus-in images, flux-weighted centroid of the positive deficit vs Gaia DR3 positions through the TPF WCS: deficit centroid **0.35-0.75 px from the target (mean 0.56 px = 2.2") vs 1.04-1.46 px from the neighbour (mean 1.29 px = 5.1")**, closer to the target in **12/13 quarters** (Q2, season 0, is a tie: 1.18 vs 1.16 px); in every quarter the deficit peak pixel differs from the mean-image peak pixel (the neighbour). Total deficit, median **7.3 e/s**, equals 4000 e/s x 2000 ppm = 8 e/s expected for an on-target dip. Caveat (S): a positive-clipped centroid of a noisy 4x5 stamp is biased toward the stamp centre, which here is toward the neighbour, so the ~0.5-px offset from the target is not evidence for an intermediate source; a PRF fit (Kuber's official-PRF fits favour the target 13/13 quarters, sum dchi2 799; raw-TPF centroid dchi2 934 vs neighbour) is the proper version. None of these can separate the target from a companion inside ~1-2".

New fact a referee will notice: the neighbour has **Gaia RUWE 8.73 and astrometric_excess_noise_sig 976** (parallax 0.30 +/- 0.23 mas, BP-RP 0.92, Teff_gspphot 5484 K): a distant G star with a badly fitting astrometric solution, i.e. probably itself an unresolved multiple. Irrelevant to the candidate only because its own photometry is flat at the ephemeris and the deficit is on the target's pixel.

### 1.4 The Kepler-169 e coincidence (E)

* Archive values (TAP `cumulative`, `ps`): Kepler-169 e = KOI-505.01, **P = 13.767101 +/- 0.000030 d, T0 = BKJD 174.81089, depth 684.8 ppm, duration 2.91 h, DR25 MES 35.4**, host KIC 5689351 Kp 14.19, 4.39' away.
* Geometry: I read the Kepler LC headers for both stars in all 13 quarters: **same MODULE/OUTPUT/CHANNEL in every season** (Q1 mod 3/out 2/ch 6; season 0 mod 11/out 2/ch 34; season 1 mod 23/out 2/ch 78; season 2 mod 15/out 2/ch 50). From the Q6 TPFs: KIC 5689116 corner col 258/row 338, KIC 5689351 corner col 202/row 297, i.e. **~56 columns and ~41 rows apart (66 px = 263")**.
* Coughlin et al. 2014 (AJ 147, 119; arXiv:1401.1240) mechanisms and their geometry: direct PRF contamination requires separation < d_max = 50 sqrt(1e6 x 10^(-0.4 Kp) + 1)" = **88" for a Kp 14.19 parent** (263" here); **column anomaly** requires the same module and |dcol| <= 10 px (56 here); **CCD/video crosstalk** requires the same module, a *different* output, and the same row and column within 10 px (same output here, 41/56 px off); antipodal reflection requires opposite sides of the field. None applies. Furthermore Coughlin's match definition needs **both period and epoch** to agree (sigma_P >= 3.5 and sigma_T >= 2.0); here the periods differ by 0.00235 d (32 sigma) and the transit midpoints are never closer than 67 h (offset drifting 2.80 -> 3.04 d, `k169_overlap.json`). Contamination is simultaneous by construction, so a different epoch is physically impossible for any electronic or optical path. The box depths in the target's PSF light curve at all five Kepler-169 ephemerides are null (K). Verdict: a numerical coincidence (S: roughly a 1e-3 chance per candidate of landing within 0.0023 d of one of five nearby planets' periods in a 0.5-25 d search, i.e. unremarkable across tens of candidates).

### 1.5 Inverted controls and the search's discriminating power (E, from Kuber's tables)

I re-read `analysis/kbonus-extension-20261011/{,ranks577-1800/}signals_real_and_inverted.csv` (4,040 rows; 1,672 stars with both polarities). Top signal per star:

| cut | real stars | inverted stars |
|---|---|---|
| BLS SNR >= 9 | 342 | 355 |
| ... & P > 1.5 d & depth < 4000 ppm | 183 | 170 |
| ... & holdout SNR >= 4 | 16 | 18 |
| SNR >= 15 & P>1.5 & depth<4000 & holdout>=4 | 7 | 10 |
| SNR >= 26.7 & P>1.5 & depth<4000 & holdout>=4 | 3 | 3 |

The first screen has **zero** discriminating power: equal pass counts at every threshold. But the matching passes are not noise. The three real and three inverted stars at candidate strength are KIC 7428405 (1.857 d, SNR 32/34 in both polarities), KIC 11560248 (1.827 d, SNR 31/29 in both), KIC 7765766 (2.44 d, both), and KIC 8701308 (inverted 4.254 d SNR 29.7; real polarity 8.505 d = 2P): **periodic variables appearing in both polarities or at a harmonic**. Kuber's second "stronger-than-candidate" inverted star KIC 5283912 has a 7,883-ppm, 3.1-h "brightening" at 8.78 d (S: plausibly a neighbour's eclipse mis-assigned with the wrong sign by the simultaneous PSF fit, the "misattributed variability" failure mode the Kepler Bonus paper warns about).

A better-calibrated null: among the **1,330 stars whose real-polarity top SNR < 9** (quiet stars), the inverted SNR distribution has percentiles 50/90/99 = 6.8/9.2/14.5; with P > 1.5 d and depth < 4000 ppm only three exceed SNR 25 (KIC 8784301 34.9 but holdout 0.8 -> fails; KIC 8701308 29.7, caught by a harmonic veto; KIC 10287716 25.0, holdout 27.4, a brightening variable). Conversely among 1,317 stars whose inverted top SNR < 9, the strongest real transit-like signal is KIC 9177981 (SNR 16.4, odd/even 8.6 sigma = EB). The candidate's own inverted search peaks at SNR 10.7 (16.8 d) and fails holdout, so it is **not** a both-polarity variable. Reading: at the candidate's strength the pipeline's "noise plus unvetoed variability" rate is ~1-2 per 1,700 stars, and with a both-polarity/harmonic veto it is ~0-1 per 1,700. Kuber's conclusion ("the case rests on vetting, not detection strength") is correct; his earlier "2/1,736 empirical FAP" was over-interpreted and he has retracted it. The 31/29 vs 32/28 parity is an indictment of the first screen, not of the candidate.

Comparability note: Kuber's "BLS SNR 26.7" is computed on the 9 discovery quarters; scaling by sqrt(13/9) gives ~32, close to the reproduction's 35 on all 13 quarters.

### 1.6 Why TRICERATOPS says 6.9 % (E/K)

TRICERATOPS (Giacalone et al. 2021, arXiv:2002.00691, Table 1) defines **FPP = 1 - (P_TP + P_PTP + P_DTP)**; planets around a bound secondary (STP), background/foreground stars (BTP) and all EB flavours count as false positives; validation needs FPP < 0.015 and NFPP < 0.001. Kuber's runs (`fpp_summary.json`, four 1e6-draw runs with a Gaia G<21 field): TP 0.750, PTP 0.172, **STP 0.068**, DTP 0.009, BTP 6e-4, BEB 2e-4, all EB scenarios ~0; FPP 0.058-0.089, mean 0.069. Why STP: without a contrast curve TRICERATOPS assumes companions are only excluded beyond 2.2" (Gaia resolution), the bound-companion prior follows Moe & Di Stefano (2017), and a 2.1-2.4-h, 2000-ppm dip on a 0.5 Rsun primary can equally be a ~3.5-4 Re planet around a fainter bound M dwarf. EB scenarios die because odd/even match, there is no secondary and the Kepler-band depth/shape disfavours a diluted grazing EB. A hypothetical Keck AO curve only reaches 0.027 because companions inside ~0.1" (~30 au at 298 pc) remain. TRICERATOPS 1.1 does not use RUWE; MOLUSC would, but RUWE = 0.984 at 298 pc is a weak constraint (Castro-Ginard et al. 2024, arXiv:2404.14127: detection probability of unresolved binaries via RUWE falls steeply with distance; ~50 % at ~100 pc). NFPP = 0 is an artefact of not passing the resolved neighbours; the 7.1" neighbour is excluded empirically (1.3) and the other Gaia sources within 30" are G > 19 (would need > ~30 % eclipses) — acceptable, but a referee will want NFPP computed.

### 1.7 Gaia DR3 (E)

Target: G 16.28, BP-RP 2.19, parallax 3.406 +/- 0.040 mas (298 pc), pm (-5.40, +2.85) mas/yr, RUWE 0.984, ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.020, astrometric_excess_noise_sig 0.00, non_single_star 0, phot_bp_rp_excess_factor 1.364 (normal for BP-RP 2.2), Teff_gspphot 3753 K. Within 40": the 7.1" neighbour (above), six sources G 19.1-20.7 at 14.7-28", a G 18.8 at 30.5", two G 17.2-17.3 at 36-37". None shares the target's proper motion/parallax, so there is no resolved co-moving companion. TIC: Tmag 15.21, Teff 3633 K, R 0.507, M 0.505, rho 3.87 rho_sun, **contratio 0.33**, lumclass DWARF.

### 1.8 TESS cannot test it (E)

The only TESS products for TIC 120686465 are the ten **TARS** light curves (TESS All-Sky Rotation Survey, HLSP lead Andrew Boyle, DOI 10.17909/2sff-fn29, https://archive.stsci.edu/hlsp/tars; systematics-corrected FFI photometry for 7.5 M stars to Tmag 16, built for rotation, not transits). lightkurve 2.6.0 cannot parse them (`LightkurveError ... of type generic`); astropy reads the two-column TIME/FLUX table. Per-sector 2.3-h-binned scatter 2,600-5,100 ppm (S74 anomalously 790). Fold at Kuber's ephemeris: **+1181 ppm with random-phase std 994 ppm (1.2 sigma)**; per-event mean 642 +/- 603 ppm over 20 covered events. Expected diluted depth in the TESS band: target share of the pair plus faint neighbours ~0.41 (Tmag 15.21 vs 14.97 plus six G>18 stars) -> ~830 ppm -> expected SNR ~0.8 in TARS. Kuber's custom 1.5-px aperture (per-event scatter 1,294 ppm, predicted increment 987 ppm) would reach expected SNR ~3.4 over 20 events, and observed -1221 +/- 1236 ppm is in **1.8-sigma tension** with that prediction (1-px: 1.5 sigma; 2-px: 2.2 sigma; TGLC: 0.3 sigma). Ephemeris propagation uncertainty: 21 min at sector 14, 33 min at sector 81 (P uncertainty 6.5e-5 d over 221-353 cycles), i.e. ~1/5 of the duration — not the limiting factor. Conclusion: a mild, correlated, aperture-dependent tension; not a rejection, not a recovery. The 2023-2024 200-s sectors (74,75,80,81) are the only ones worth a pixel-level forward model.

### 1.9 Stellar density and duration (E)

a/R* = 39.7 (Kuber's MCMC median) -> rho* = 3pi/(G P^2) (a/R*)^3 = 6.24 g/cc = **4.43 rho_sun**; TIC gives 3.87; the posterior 4.4 (+2.3 -2.8) is consistent with the M dwarf but 1.2 sigma from a solar-density blend host, so it is weak evidence. The duration is more useful: T14(b=0) for P = 13.77 d around 0.505 Msun/0.507 Rsun is 2.9 h, matching 2.1-2.4 h at b ~ 0.5; a rho = 1 host would need b >~ 0.85. Kepler's 29.4-min cadence smears a ~15-min ingress, so b and rho are degenerate; do not over-sell this.

### 1.10 Seasonal / pixel-level consistency (E)

Per-quarter depths with the empirical per-event scatter as the error: chi2 = 15.6/12 (p = 0.21); seasons 0/1/2/3 = 1693 +/- 173, 1967 +/- 177, 2101 +/- 194, 1245 +/- 613 ppm. Season 0 (channel 34, 7-pixel SAP apertures) is where the KBONUS PSF model struggles for this pair: the neighbour's PSF flux has 11-14 % 5-95 % ranges in Q2/6/10/14 vs ~1 % elsewhere, and the target's depth and localisation are weakest there (Q2 centroid tie). The detection does not depend on season 0; a referee would ask for season-0 to be down-weighted or the PSF reliability discussed (`PERTRATI`/`PERTSTD` per quarter are in the headers).

### 1.11 Rotation, youth and the Barber et al. 2022 question (E/S)

Lomb-Scargle on the raw per-quarter PSF flux: **P_rot = 10.42 d globally, 10.3-11.2 d in 12/13 target quarters, absent from the neighbour** (whose best periods are 28-40 d, i.e. systematics). The CDS ReadMe for J/AJ/164/88 (fetched in-sandbox) defines the Table 5 flags as "Of sufficient quality for running the Notch transit-search pipeline" (TESS = N, Kepler = Y for this star), the identifier as Gaia EDR3, and the row columns as Voff = 0.795 km/s (tangential velocity offset), SpType M1.6, parallax 3.406 +/- 0.040 mas (the trailing 0.04 is e_plx, not a membership probability); no Prot, Li or RV. MELANGE-3 is a 105 +/- 10 Myr association (Barber et al. 2022, AJ 164, 88; arXiv:2206.08383). Because this Gaia source was never a Kepler target, the "Kepler light curve" Barber et al. had was most plausibly KIC 5689116's PDCSAP, assigned by positional crossmatch (S); in that light curve the signal is ~80 ppm (1.1) and Notch would never see it. So the flag does not imply prior coverage, but only the authors can confirm. A 10.4-d rotation period for an M1.6 dwarf is consistent with youth (S: field early-M dwarfs mostly rotate at tens of days; Pleiades-age early-M slow rotators are at ~5-10 d), which (i) supports candidate membership, (ii) implies a possibly inflated radius (Gaia GSP-Phot 0.555 Rsun vs TIC 0.507) and Rp closer to 3 Re, (iii) would make this a scientifically valuable young planet and raises the bar for host characterisation. The only other background survey of the field, Bienias & Szabo 2025 (A&A 694, A209; arXiv:2501.09152), ran a 3.5-sigma dip search on Q4 target-aperture pixels with a minimum detectable depth of 0.18-0.25 % and found only hot Jupiters at 2.5-4.8 d, so it had no sensitivity here.

### 1.12 What remains genuinely unexcluded

(E) A planet around an unresolved bound companion inside ~2" (STP, ~7 %); (E) a background star inside ~1" hosting a planet (BTP, <0.1 % by TRICERATOPS priors); (S) residual PSF-model mis-attribution is disfavoured but not quantified for season 0; (E) the host radius and age.

## 2. What would settle it

See the structured list. In priority order: AO/speckle imaging + MOLUSC/TRICERATOPS rerun (professional time); Gaia DR4 (scheduled December 2026 per ESA; verify) for bound companions; NIR reconnaissance spectroscopy for SB2 and host parameters; an independent PRF-fit difference-image localisation (amateur-feasible); a population-calibrated search FAP with injection-recovery (amateur-feasible); a TESS two-star forward model as a consistency check only (expected SNR <= 3). Multi-band transit photometry is not feasible for a 2000-ppm dip at G = 16.3.

## 3. What the developer should do differently

Summarised in the lessons list; the two decisive changes are **target selection for validatability** (brightness, isolation, imaging availability, TESS followability: Tmag <= 13, contratio < 0.3, predicted SNR >= 7) and **pre-registered, population-calibrated search statistics with a variability veto**, so the write-up can say "0 of N inverted/scrambled stars pass the full vetting" rather than "the first screen passes noise as often as data". Run the cheap on-target tests (neighbour LC, dilution model, difference image, official-pipeline leak) before any public statement, and compute NFPP by passing the Gaia neighbours. Publish as a candidate through a peer-reviewed venue or ExoFOP (TESS targets), not RNAAS, which describes itself as "moderated but not edited", non-peer-reviewed, <= 1,500 words (https://journals.aas.org/research-notes/).

## 4. Assessment of the reproduction

* **Recovers the same signal, blind**: P within 12 ppm of Kuber's (BLS 13.769435, pass-2 refine 13.769616, TLS 13.769643 vs 13.769448 +/- 0.000065 d); the pass-2 value is formally 2.6 sigma from Kuber's but the refinement grid, not the data, sets that difference; T0 offset at mid-baseline -1.4 min; 74 transits with data vs 70 usable; odd/even, secondary, neighbour null all agree; Rp/R* 0.044-0.046.
* **Disagreements/caveats**: (a) the structured summary quoted smoke-run numbers (2058 +/- 67 ppm, 2.0-h box) while the finished run gives 1916 +/- 68 ppm with a 2.3-h box — purely box-width dependence, but the narrative should be regenerated from `results.json`; (b) "BLS SNR" definitions differ (35 vs 26.7, see 1.5); (c) formal errors are 1.4-1.7x too small (the reproduction already flags this); (d) the single-star inverted null (SDE 5.2) is not a population FAP; (e) TLS depth (1600 ppm, 1.65 h) is a model minimum and must not be compared with box depths; (f) limb darkening came from the TESS-band TIC table; (g) no pixel-level, SAP-dilution or official-PDCSAP tests were in the pipeline — I added them here and they agree with Kuber (1.1, 1.3).
* **Useful discovery by the reproduction**: a blind 0.5-d biweight pulls the in-transit trend down by ~250-300 ppm, biasing the depth ~20 % low and manufacturing a +200-400-ppm shoulder that appears as an inverted-null peak at the true period (SNR ~8). A developer who sees "the signal is in both polarities" after a single-pass detrend could wrongly reject a real transit; the fix is the masked second pass (Hippke et al. 2019, AJ 158, 143).

## 5. Established vs speculation

Established: reality, periodicity, on-target-side localisation to ~2", non-contamination by Kepler-169, dilution behaviour, Gaia cleanliness, TESS impotence, 10.4-d rotation on the target, FPP mechanics. Speculation (labelled S above): the Barber light-curve identity, MELANGE-3 membership and youth-inflation, the nature of KIC 5283912's inverted signal, PSF mis-attribution in season 0, and the ~1e-3 chance of the Kepler-169 period coincidence.

## Sources

Martinez-Palomera, Hedges & Dotson 2023, AJ 166, 265, arXiv:2310.17733; MAST KBONUS-BKG page https://archive.stsci.edu/hlsp/kbonus-bkg (DOI 10.17909/7jbr-w430). Coughlin et al. 2014, AJ 147, 119, arXiv:1401.1240 (contamination geometry and match criteria). Giacalone et al. 2021, AJ 161, 24, arXiv:2002.00691 (TRICERATOPS scenarios, FPP = 1-(TP+PTP+DTP), thresholds). Barber et al. 2022, AJ 164, 88, arXiv:2206.08383; CDS ReadMe https://cdsarc.cds.unistra.fr/ftp/J/AJ/164/88/ReadMe. Bienias & Szabo 2025, A&A 694, A209, arXiv:2501.09152. TARS HLSP https://archive.stsci.edu/hlsp/tars (DOI 10.17909/2sff-fn29). Castro-Ginard et al. 2024, arXiv:2404.14127 (RUWE detectability). Hippke et al. 2019, AJ 158, 143 (wotan). NASA Exoplanet Archive TAP https://exoplanetarchive.ipac.caltech.edu/TAP/sync (tables `ps`, `cumulative`, `q1_q17_dr24_tce`, `q1_q17_dr25_tce`). Gaia DR3 via astroquery.gaia; TIC via astroquery.mast.Catalogs; Kepler LC/TPF headers via lightkurve 2.6.0. RNAAS scope https://journals.aas.org/research-notes/. Kuber's repo files cited inline (`analysis/validation-20261011/SCIENCE.md`, `results/kepler/*.json`, `results/triceratops/fpp_summary.json`, `analysis/editorial-audit-20261013/{README.md,verification.json}`, `analysis/kbonus-extension-20261011/**/signals_real_and_inverted.csv`, `evidence/codex/database_expansion/kbonus/pilot_summary.csv`).