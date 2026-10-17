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