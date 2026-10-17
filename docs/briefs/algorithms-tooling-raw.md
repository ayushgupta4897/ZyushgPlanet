# Algorithms and tooling for a transit search at scale on a 4-CPU box (checked 2026-10-17)

Everything below was either run in the sandbox venv (`/home/user/ZyushgPlanet/.venv`, Python 3.13.16) or read from a primary source; each claim carries its source. Scripts, logs and JSON results are in `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/algorithms-tooling/` (file names `s01_…` to `s09_…` below). Kuber's repo, Chopra's repo and eumemic's repo were read (never executed or imported) from `scratchpad/untrusted/planethunt` and `scratchpad/algorithms-tooling/upstream/`.

Two things changed days before this was written and matter for anyone copying 2024-25 tutorials: **TLS 2.0 (6 Oct 2026) / 2.0.1 (11 Oct 2026)** and **wotan 2.0 (11 Oct 2026)** were released ([TLS releases](https://github.com/hippke/tls/releases), [wotan releases](https://github.com/hippke/wotan/releases)). The venv already has both. TLS 2 has "18 bug fixes, 10 of which change results", a new default backend, a built-in limb-darkened model (batman no longer a runtime dependency) and is "30-60x faster than TLS 1.33"; wotan 2.0 fixed 31 bugs (robust clipping is now two-sided, slider windows include their last sample) and is 60x faster for the biweight.

---

## 0. Tool inventory: pip name, version, install status, docs

PyPI versions queried live on 2026-10-17 (`curl https://pypi.org/pypi/<pkg>/json`). "venv" = importable in the sandbox venv now.

| Tool | pip name | PyPI latest (upload date) | venv | Docs | Ran here? |
|---|---|---|---|---|---|
| lightkurve | `lightkurve` | 2.6.0 (2026-04-16) | 2.6.0 | https://lightkurve.github.io/lightkurve/reference/index.html | yes (search/download/TPF/TESScut/PRF/centroids) |
| astropy BLS | `astropy` | 8.0.2 (2026-10-15) | 8.0.2 | https://docs.astropy.org/en/stable/api/astropy.timeseries.BoxLeastSquares.html | yes |
| TLS | `transitleastsquares` | 2.0.1 (2026-10-11) | 2.0.1 | https://transitleastsquares.readthedocs.io/en/latest/Python%20interface.html | yes (benchmarks + real LC) |
| wotan | `wotan` | 2.0 (2026-10-11) | 2.0 | https://wotan.readthedocs.io / https://github.com/hippke/wotan | yes |
| batman | `batman-package` | 2.5.3 (2025-05-05) | installed (module reports `__version__` 2.5.1) | https://lkreidberg.github.io/batman/docs/html/tutorial.html | yes |
| emcee / corner | `emcee`, `corner` | 3.1.6 (2024-04-19) / 2.3.0 | installed by me | https://emcee.readthedocs.io | yes |
| celerite2 | `celerite2` | 0.3.3 (2026-07-12, py>=3.11) | installed by me | https://celerite2.readthedocs.io/en/latest/api/python/ | yes |
| dynesty / ultranest | `dynesty`, `ultranest` | 3.1.0 / 4.6.3 | dynesty installed | https://dynesty.readthedocs.io | no |
| juliet | `juliet` | 2.2.10 (2026-02-16) | installed by me (pulls celerite 0.4.3 and ~40 deps) | https://juliet.readthedocs.io/en/latest/user/priorsnparameters.html | import only |
| exoplanet | `exoplanet` (+`exoplanet-core`, `pymc`) | 0.6.0 (2024-04-16); exoplanet-core 0.4.0 (2026-07-25); pymc 6.3.2 (2026-09-08) | not installed (dry-run resolves) | https://docs.exoplanet.codes/en/latest/user/install/ | no |
| PyTransit | `pytransit` | 2.9.3 (2026-10-15) | not installed (dry-run: 33 pkgs incl. arviz/xarray) | https://github.com/hpparvi/PyTransit | no |
| allesfitter | `allesfitter` | 1.2.10 (2022-07-25, stale) | not installed | https://www.allesfitter.com | no |
| LDTk | `ldtk` | 1.9.3 (2026-10-15) | installed with `--no-deps` + `tenacity` (its metadata pins pandas 1.5.3, which will not build on py3.13) | https://github.com/hpparvi/ldtk | import only |
| LEO-Vetter | `leo-vetter` | 1.2.0 (2026-07-29) | installed by me (+lmfit, brokenaxes) | https://github.com/mkunimoto/LEO-vetter | yes (flux-level) |
| exovetter | `exovetter` | 0.0.16 (2026-04-27) | installed by me | https://exovetter.readthedocs.io | import only |
| TRICERATOPS | `triceratops` | 1.1.0 (2026-09-02) | NOT installed: bare `uv pip install triceratops` resolved to 1.0.21 and failed on pandas 1.5.3; `triceratops==1.1.0` dry-run resolves but pins `pytransit==2.2` ([setup.py](https://raw.githubusercontent.com/stevengiacalone/triceratops/master/setup.py)) -> use a separate venv | https://github.com/stevengiacalone/triceratops | no |
| psfmachine | `psfmachine` | 1.1.4 (2023-12-15, py<4,>=3.8) | FAILED: depends on pandas 1.5.3 which cannot build on Python 3.13 -> needs a py<=3.11 venv | https://ssdatalab.github.io/psfmachine/ | no |
| nuance | `nuance` | 0.8.1 (2024-07-31) | not installed (dry-run: jax 0.11.2, jaxlib, jaxopt, jaxtyping, tinygp 0.3.1, 12 pkgs) | https://nuance.readthedocs.io | no |
| cuvarbase (GPU BLS) | `cuvarbase` | 0.2.5 (2023-10-23) | not installable (no GPU; needs PyCUDA + CUDA toolkit; docs tested CUDA 8 and Python 2.7-3.6) | https://johnh2o2.github.io/cuvarbase/bls.html | no |
| transit-diffImage | GitHub only (`pip install .`) | n/a | not installed; requires `tess-point>=0.8.0, astropy, astroquery, scipy, scanf` (repo requirements.txt) | https://github.com/stevepur/transit-diffImage | no |
| tess-point | `tess-point` | 0.9.5.post1 (2026-09-10) | not installed | https://github.com/christopherburke/tess-point | no |
| SHERLOCK | `sherlockpipe` | 1.5.1 (2026-05-29) | dry-run FAILED (dependency `ellc` build needs numpy without build isolation; docs say install numpy and Cython first, Python 3.10) | https://sherlock-ph.readthedocs.io/en/latest/installation.html | no |
| MonoTools | `MonoTools` | 1.0.0 (2025-06-06) | not installed (PyMC5/pytensor backend) | https://github.com/hposborn/MonoTools | no |
| ExoMiner | container `ghcr.io/nasa/exominer:latest` (Podman) | repo has no release tag | not run (needs SPOC TCE/DV products) | https://github.com/nasa/ExoMiner | no |
| Astronet (Shallue & Vanderburg) | GitHub (archived 2023-08-10) | n/a | no weights found | https://github.com/google-research/exoplanet-ml | no |
| Astronet-Triage(-v2) | GitHub | n/a | README shows training your own model; no released checkpoints seen | https://github.com/mdanatg/Astronet-Triage ; data https://zenodo.org/record/7411579 | no |
| qusi (RAMjET) | `qusi` | 1.5.6 (2025-04-21) | not installed; "framework for producing neural networks"; no bundled pretrained models evident | https://qusi.readthedocs.io | no |
| QATS (Carter & Agol) | C++ tarball only | n/a | not installed. Warning: PyPI `qats` 5.5.0 is DNV's unrelated time-series library | https://faculty.washington.edu/agol/QATS/ (ascl:1712.014) | no |
| Kepler SOC pipeline / Robovetter | GitHub | n/a | pipeline is MATLAB/Java and "not expected that the reader will be able to build or run"; Robovetter is C++11 and needs DR25 metric input files | https://github.com/nasa/kepler-pipeline ; https://github.com/nasa/kepler-robovetter | no |

Other ecosystem packages checked on PyPI (not installed): `foldedleastsquares` 1.1.12 (SHERLOCK's TLS fork), `lcbuilder` 0.27.0, `dearwatson` 1.1.0 (WATSON vetting), `tglc` 0.8.1, `eleanor` 2.0.5 (2022), `stella` 0.1.0 (2021, flare CNN), `astrobase` 0.5.3 (2021), `jax` 0.11.2, `tinygp` 0.3.1, `numpyro` 0.22.0. Not on PyPI: `molusc`, `exominer`, `tesscut` (it is `astroquery.mast.Tesscut`).

---

## 1. Loading and cleaning

### 1.1 lightkurve search/download (ran: `s01_kbonus_load.py`, `s02_kbonus_download.py`)

```python
import lightkurve as lk, astropy.units as u
from astropy.coordinates import SkyCoord
c = SkyCoord(285.907827, 40.941988, unit='deg')
sr = lk.search_lightcurve(c, author='KBONUS-BKG', radius=20*u.arcsec)   # 26 s
# 2 products: target_name 'Gaia DR3 2103792461826460288' (dist 0.0") and 'Gaia DR3 2103792427466721792' (3.1")
sr = sr[sr.target_name == 'Gaia DR3 2103792461826460288']
lcc = sr.download_all(download_dir=DATA)   # 8.44 MB FITS: .../hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits
```

Note the `dataURL` names the file by KIC (`kic-005689121`) while `target_name` is the Gaia designation; the MAST product description says light curves are named by KIC ID when one exists, otherwise `gaia-dr3-<id>` ([MAST HLSP page](https://archive.stsci.edu/hlsp/kbonus-bkg)). Kuber's batch code bypassed lightkurve and hit the MAST `Mast.Caom.Filtered` service directly with `provenance_name='KBONUS-BKG'` and a list of `target_name` Gaia designations, then fetched `https://mast.stsci.edu/api/v0.1/Download/file?uri=<dataURL>` (his `pilot.py`, `batch2.py`); that is the right pattern for thousands of stars because one POST returns all URLs.

### 1.2 KBONUS-BKG file anatomy (verified with `astropy.io.fits`)

- `LIGHTCURVE_STITCHED` (HDU 1, 50,182 rows): `TIME` (BJD-2454833 = BKJD), `CADENCENO`, `QUARTER`, `FLUX`, `FLUX_ERR` (e-/s), `SAP_FLUX`, `SAP_FLUX_ERR`, `PSF_FLAT_FLUX`, `PSF_FLAT_FLUX_ERR`, `SAP_QUALITY`, `FLATTEN_MASK`. MAST defines `Flux` = "PSF flux from stitched quarters", `PSF_flat_flux` = "PSF flux ... after flattening", `SAP_quality` = "Quality flag from the TPF", `Flatten_mask` = "Quality flag from the flattening process" ([MAST](https://archive.stsci.edu/hlsp/kbonus-bkg)).
- `LIGHTCURVE_Q<n>` per quarter: `CADENCENO, TIME, FLUX` ("Corrected PSF flux"), `FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR, CENTROID_COLUMN, CENTROID_ROW, SAP_QUALITY, SAP_BKG, RED_CHI2, PSF_FLUX_NOVA, PSF_FLUX_ERR_NOVA`; and `APERTURE_Q<n>` masks. Per-quarter headers carry `PSFFRAC`, `CROWDSAP`, `FLFRCSAP` (Kuber's `pilot.py` reads `hd.header['PSFFRAC']` and drops quarters with PSFFRAC < 0.5). The primary header has `OBJECT='KIC 5689121'`, `KEPMAG 16.332`, `GMAG 16.282`, `TPFORG 5689116` (the TPF the source was extracted from), `DOI 10.17909/7jbr-w430`, `HLSPVER V1.0`.
- For this star: quarters 1,2,3,4,6,7,8,10,11,12,14,15,16 (13 quarters; the originating TPF exists only in those quarters), 46,523 good points after `SAP_QUALITY==0`, baseline 1,426 d, robust point-to-point scatter 1,033 ppm per 29.4-min cadence (`s04_e2e.log`).
- Kuber's choices (pilot.py): use `FLUX` (corrected PSF), `SAP_QUALITY == 0`, finite and positive errors, `PSFFRAC >= 0.5`, >= 300 good points per quarter, normalise each quarter by its median, convert BKJD to BTJD with `-2167`. Kepler Bonus itself (Martinez-Palomera et al. 2023, [arXiv:2310.17733](https://arxiv.org/abs/2310.17733)) extracted 606,900 light curves with the Linearized Field Deblending PSF method (psfmachine), 406,548 of them new background sources.

Practical: `lightkurve`'s generic reader hands you one table extension; for the per-quarter extensions and header metrics use `astropy.io.fits` directly (as in `s04_kbonus_e2e.py`).

### 1.3 Quality bitmasks (values printed from the venv, `lightkurve.utils`)

- Kepler/K2: `KeplerQualityFlags.DEFAULT_BITMASK = 1130799` = Attitude tweak, Safe mode, Coarse point, Earth point, Desaturation event, Manual exclude, Detector anomaly, No data, **Thruster firing** (the K2 6-h thruster-fire flag is in the default mask); `HARD_BITMASK = 1664431`; `HARDEST_BITMASK = 2096639`. Docs: https://lightkurve.github.io/lightkurve/reference/api/lightkurve.KeplerQualityFlags.html.
- TESS: `TessQualityFlags.DEFAULT_BITMASK = 17087` = Attitude tweak, Safe mode, Coarse point, Earth point, Argabrightening, **Desaturation event** (momentum dumps), Manual exclude, Impulsive outlier, Bad Calibration Exclude; `HARD = 24319`; `HARDEST = 65535`. Docs: https://lightkurve.github.io/lightkurve/reference/api/lightkurve.TessQualityFlags.html; bit definitions in the TESS SDPDD ([EXP-TESS-ARC-ICD-TM-0014-Rev-F](https://archive.stsci.edu/files/live/sites/mast/files/home/missions-and-data/active-missions/tess/_documents/EXP-TESS-ARC-ICD-TM-0014-Rev-F.pdf)).
- Kuber's `pixel_diagnostic.py` used `(QUALITY & (17087|8192)) == 0` on TESS TPFs; on Kepler TPFs `QUALITY == 0`. For a shallow-transit search the default mask is the right compromise; `QUALITY==0` (what KBONUS `SAP_QUALITY==0` gives you) is stricter and loses a few percent of cadences.

### 1.4 Quarter / campaign / sector boundaries

- Kepler: quarters are separated by ~3-day rolls and have ~1-day monthly downlink gaps; each quarter sits on a different CCD channel (different aperture, flux fraction, crowding), so **normalise and detrend per quarter and never across a boundary**. Reference: Kepler Data Characteristics Handbook KSCI-19040-005 (https://archive.stsci.edu/kepler/manuals/Data_Characteristics.pdf). For KBONUS the per-quarter `PSFFRAC`/`CROWDSAP` change, so SAP depths vary with quarter (Kuber turned that into his dilution test).
- K2: campaign-level files; 6-h thruster firings produce the dominant systematic; use `quality_bitmask='default'` (includes Thruster firing) and either the pre-corrected EVEREST/K2SFF HLSPs or `lightkurve.correctors.SFFCorrector` (https://lightkurve.github.io/lightkurve/reference/api/lightkurve.correctors.SFFCorrector.html); K2 Handbook: https://archive.stsci.edu/k2/manuals/k2_handbook.pdf.
- TESS: one file per sector; two orbits per sector with a ~1 d downlink gap; momentum dumps every few days (flag bit 32). eumemic found a QLP-specific systematic in Sectors 94-104: a ramp at each orbit's start/end that detrending turns into a dip repeating at ~13.7 d (43% of new signals landed at 12-14 d); their mask drops the first 1 d and last 0.5 d of each orbit's QLP data in those sectors ([validation/qlp_orbit_edges.md](https://github.com/eumemic/exoplanets/blob/main/validation/qlp_orbit_edges.md)). Chopra's triage also flags periods near the 13.7-d TESS orbit as a soft flag ([GUIDE.md](https://github.com/paraschopra/tess-planet-candidates/blob/main/docs/GUIDE.md)).
- Stitching: eumemic prefers SPOC 2-min PDCSAP, then TESS-SPOC FFI PDCSAP, then QLP, one product per sector; Kuber's `acquire.py` does the same ranking and stores one compressed `.npz` per star with `time, flux, flux_err, sector, exptime, cx, cy, bkg`. Per-sector noise rejection (Chopra `lcsearch.reject_bad_sectors`): drop a sector if its 1-h point-to-point noise exceeds 1% or 3x the star's median sector.

### 1.5 Outlier rejection and normalisation

- Asymmetric clipping (Chopra `lcsearch.clean`): after detrending, clip **upward** residuals at `FLARE_CLIP_SIGMA = 3.0` and **downward** only beyond `DROPOUT_CLIP_SIGMA = 10.0` ("downward points could be transits"). eumemic clips only upward at 4 sigma. `lightkurve.LightCurve.remove_outliers(sigma_upper=3, sigma_lower=10)` does the same ([docs](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.LightCurve.remove_outliers.html)). Note Chopra's inverted-null caveat: because the clipping is asymmetric, flipping the flux is not a perfect null.
- Normalisation: divide each quarter/sector by its own median before stitching (all three pipelines). Error bars: do not trust pipeline errors for faint stars; use the robust point-to-point scatter `1.4826*median(|diff(f)|)/sqrt(2)` per segment (Chopra) or a rolling +-0.5 d robust std floored at half the global scatter (eumemic `rolling_std`).
- Our run: after `SAP_QUALITY==0`, biweight (1.0 d per quarter) and the 3/10-sigma clip, 45,648 points remain of 46,523.

---

## 2. Detrending

### 2.1 wotan 2.0 (ran: `s04_kbonus_e2e.py`, in-line snippets)

```python
from wotan import flatten, transit_mask
flat, trend = flatten(t, f, method='biweight', window_length=1.0, break_tolerance=0.5,
                      edge_cutoff=0.0, return_trend=True)
mask = transit_mask(time=t, period=P, duration=2*dur, T0=T0)   # then flatten(..., mask=mask)
```

Signature in the venv: `flatten(time, flux, window_length=None, edge_cutoff=0, break_tolerance=None, cval=None, return_trend=False, method='biweight', kernel=None, kernel_size=None, kernel_period=None, proportiontocut=0.1, robust=False, max_splines=50, stdev_cut=2, return_nsplines=False, mask=None, verbose=False)`. Methods: biweight, huber, hampel, andrewsinewave, welsch, median, medfilt, mean, trim_mean, winsorize, rspline, hspline, pspline (needs `pygam`, not installed here), cofiam, cosine, savgol, lowess, supersmoother, ridge/lasso/elasticnet, gp (squared_exp, matern, periodic, periodic_auto) ([README](https://github.com/hippke/wotan)).

Rules from the wotan paper (Hippke et al. 2019, [arXiv:1906.00966](https://arxiv.org/html/1906.00966v3)): "the optimal width of the sliding window is about three times the transit duration that is searched for" (Sec. VI); for w/T14 >~ 2.2 "most (>~98%) of the flux integral is preserved" (Fig. 8); discarding half a window at segment edges is "slightly preferable" for automated searches (ΔSDE ~0.1, Sec. IV.4); iterative GPs and splines are "inferior to robust methods for the purpose of transit searches" (GPs missed 10-11 of 100 K2 planets vs 6 for the biweight). `break_tolerance` splits the series at gaps longer than the value so each segment is detrended alone (set it at 0.5 d so Kepler monthly gaps and TESS orbit gaps split).

Measured speed (wotan 2.0, venv): biweight on 67,000 points: 0.10 s (window 0.5 or 1.0 d), 0.16 s (2.0 d); median 0.10 s; rspline 0.34 s; hspline 2.95 s; cosine 3.5 s on 4,400 points; `gp` matern 55 s on 4,400 points (unusable in a batch search). First call pays JIT warm-up (4.4 s for the 13-quarter real LC loop).

### 2.2 Over-detrending, shown numerically

Synthetic 90-d, 30-min series (1,200 ppm white noise + 3,000 ppm/4.3 d and 1,000 ppm/1.7 d sinusoids) with 2,000 ppm, 2.3-h transits at 13.77 d. Recovered in-transit depth after detrending: biweight 1.0 d -> 1,970 ppm; biweight 0.3 d -> 1,443 ppm (28% eaten); cosine 1.0 d -> 1,515; wotan `gp` matern 1.0 -> 1,096; biweight 0.3 d **with the transit masked** -> 2,176; celerite2 Matern-3/2 GP fit with transits masked -> 2,232 (+-noise). The transit-mask re-detrend trick (find signal -> `transit_mask` -> re-flatten with `mask=` -> re-measure depth) recovers what the blind window removed; Chopra's `vet.masked_clean` and eumemic's "transit-masked re-detrending" in vetting do exactly this.

### 2.3 Concrete parameter choices

Expected central transit duration for a 0.505 Msun / 0.507 Rsun M dwarf (rho = 3.9 rho_sun), using `T14 ~ 13 h (P/365.25 d)^(1/3) rho^(-1/3)` (Chopra `duration_days`): P = 0.5 d -> 0.9 h, 1 d -> 1.2 h, 5 d -> 2.0 h, 13.77 d -> 2.8 h (observed 2.1-2.4 h at b ~ 0.5), 30 d -> 3.6 h.

- **Kepler 30-min, Kp 16 M dwarf, P 0.5-30 d:** 3 x T14_max = 0.45 d is the paper's rule, but all three practitioners went longer because the 30-min cadence gives only 4-5 points per transit and the window must contain enough points: Kuber 1.0 d (per quarter, `break_tolerance=0.5`), Chopra 0.75 d (chosen by injection-recovery: 16 planted 1.5 R_earth planets recovered 3/16, **10/16**, 8/16, 5/16 for windows 0.5, 0.75, 1.0, 1.5 d; his code comment notes in-transit points pull the trend by ~duration/window, ~18% for a 2.9-h dip in 0.75 d), eumemic `window = clip(3*T14(P_max), 0.5, 1.0)` d plus Lomb-Scargle prewhitening (4 harmonics, 3*P_rot chunks) when P_rot < 3 d and amplitude > 3x noise. Recommendation: 0.75-1.0 d biweight, `break_tolerance=0.5`, `edge_cutoff=0` (keep data; edges are handled by the per-event local-baseline vetting later), clip +3/-10 sigma, and verify with injection-recovery on your own sample as Chopra did.
- **TESS 200-s / 10-min FFIs (K/M dwarfs, P <= 30 d):** T14 <= ~4 h -> 0.5 d window (eumemic's floor); fit the biweight to 10-minute means and interpolate (eumemic: 14x faster trend, residuals differ by ~30 ppm, 30 vs 28 of 100 injections detected); bin to 10-20 min before the period search (Chopra `BIN_S = 600`; eumemic coarse 20 min, fine 10 min); clip flares upward at 3-4 sigma.
- GP alternative for active stars: celerite2 (`terms.SHOTerm(sigma, rho, Q|tau)`, `RotationTerm(sigma, period, Q0, dQ, f)`, `Matern32Term(sigma, rho)`, `GaussianProcess(kernel, mean=)`, `gp.compute(t, diag=yerr**2+jitter, quiet=True)`, `gp.log_likelihood(y)`, `gp.predict(y, t=, return_var=)`, O(N) — [API](https://celerite2.readthedocs.io/en/latest/api/python/)). Our Matern-3/2 maximum-likelihood fit on 4,409 masked points took 0.37 s. Use it for vetting/fitting (transits masked), not for the blind search. TLS 2 also offers analytic SDE background detrending: `power(SDE_detrend="hybrid")` (default `"median"`, kernel 30).

---

## 3. Period search

### 3.1 astropy `BoxLeastSquares` (ran: `s04_kbonus_e2e.py`, `s07_null_throughput.py`)

API ([docs](https://docs.astropy.org/en/stable/api/astropy.timeseries.BoxLeastSquares.html)): `BoxLeastSquares(t, y, dy=None)`; `power(period, duration, objective=None, method=None, oversample=10)` with `objective='likelihood'` (default, maximises the log-likelihood) or `'snr'` ("signal-to-noise of the transit depth", recommended when noise is correlated), `method='fast'` (Cython) or `'slow'`; `autoperiod(duration, minimum_period=None, maximum_period=None, minimum_n_transit=3, frequency_factor=1.0)` with `df = frequency_factor * min(duration) / (max(t)-min(t))**2`; `compute_stats(period, duration, transit_time)` returns `depth, depth_odd, depth_even, depth_half, depth_phased, harmonic_amplitude, harmonic_delta_log_likelihood, transit_times, per_transit_count, per_transit_log_likelihood` (verified keys).

**Do not use `autoperiod` on multi-year baselines.** Measured grid sizes (0.5-30 d): Kepler 4-yr with durations [0.05,0.1,0.2] d -> 83.8 million periods; TESS 10 sectors over 7 yr -> 307 million; TESS 2 consecutive sectors -> 143k. Use a stellar-density-aware grid instead (Chopra `period_grid`, eumemic `period_grid`): step frequency so the accumulated phase drift over the baseline stays below half the shortest plausible duration divided by an oversampling factor:

```python
def dur_days(P, rho): return 13/24*(P/365.25)**(1/3)*rho**(-1/3)
def grid(baseline, rho, pmin=0.5, pmax=30, oversample=3):
    f, fq = [], 1/pmin
    while fq > 1/pmax:
        f.append(fq); P = 1/fq
        fq -= 0.5*dur_days(P, rho)/(oversample*baseline*P)
    return np.sort(1/np.array(f))
```

Sizes: 512k periods (1,460 d baseline, rho 3.9), 325k (rho 1.0), 877k (2,500 d), 19k (54 d). TLS's Ofir-style grid for the same star and baseline is 256k (oversampling 3) and the actual TLS run used 180,849 periods for 1-30 d. Duration grid: 6-8 geometric durations between 0.35 x T14(P_min) and 1.3 x T14(P_max) (Chopra `duration_grid`), or per period chunk (eumemic `run_bls` uses 7 log-period chunks each with its own 6 durations 0.4-1.6 x T14).

Measured runtime, real LC (45,648 points), 499,259 periods, 8 durations, `oversample=5`, 1 core: `objective='snr'` 186 s, `'likelihood'` 150 s. With 6 durations and `oversample=3`: 152 s, peak RSS 115 MB. Cost is close to linear in N_points x N_periods x N_durations x oversample. Chopra's `sde()` normalises BLS power by a running median in log-period (60 bins) before taking the peak; eumemic's `sde()` uses a running median/MAD in 200 log-period bins.

Our BLS on KIC 5689121 (`objective='snr'`): P = 13.76953 d, t0 = BKJD 136.3098, duration 2.45 h, box depth 1,649 +- 56 ppm, odd 1,634 +- 81 vs even 1,662 +- 78 ppm, half-period depth 801 +- 40 ppm (as expected for a real P), phase-0.5 "secondary" -71 +- 57 ppm, `harmonic_delta_log_likelihood = -430` (transit strongly preferred over a sinusoid), 74 transits with data.

### 3.2 transitleastsquares 2.0.1 (ran: `s03_tls_bench.py`, `s04_kbonus_e2e.py`, `s07`)

```python
from transitleastsquares import transitleastsquares, catalog_info
ab, M, Mmin, Mmax, R, Rmin, Rmax = catalog_info(TIC_ID=120686465)   # ((0.2125, 0.3984), 0.5054, 0.021, 0.021, 0.5075, 0.016, 0.016)
model = transitleastsquares(t, flux, dy)
r = model.power(period_min=1, period_max=30, R_star=0.507, M_star=0.505,
                R_star_min=0.3, R_star_max=0.8, M_star_min=0.3, M_star_max=0.8,
                u=list(ab), use_threads=4, oversampling_factor=3, duration_grid_step=1.1,
                show_progress_bar=False, verbose=False)
r.period, r.period_uncertainty, r.T0, r.SDE, r.snr, r.FAP, r.depth, r.duration, r.rp_rs,
r.odd_even_mismatch, r.transit_count, r.distinct_transit_count, r.transit_depths, r.snr_per_transit
```

Defaults from `tls_constants.py` in the venv: `OVERSAMPLING_FACTOR=3`, `DURATION_GRID_STEP=1.1`, `FRACTIONAL_TRANSIT_DURATION_MAX=0.12`, `N_TRANSITS_MIN=2`, `T0_FIT_MARGIN=0.01`, `TRANSIT_DEPTH_MIN=10 ppm`, `R_STAR_MIN/MAX=0.13/3.5`, `M_STAR_MIN/MAX=0.1/1.0`, `SDE_DETREND='median'` (kernel 30; `'hybrid'` available), `DEFAULT_BACKEND='fused-pl'` (piecewise-linear template approximation, "max deviation 2% of the depth", "no loss of sensitivity in injection-recovery tests"), other backends `'fused-binned'`, `'fused'` (exact), `'numba'` (TLS 1.33 algorithm); override per call `power(backend=...)` or env `TLS_BACKEND`. Valid kwargs: `R_star, R_star_min, R_star_max, M_star, M_star_min, M_star_max, period_min, period_max, n_transits_min, per, rp, a, inc, b, ecc, w, u, limb_dark, duration_grid_step, transit_depth_min, oversampling_factor, SDE_detrend, T0_fit_margin, T0_search_margin, use_threads, show_progress_bar, transit_template ('default'|'grazing'|'box'), verbose, backend`. `use_threads` defaults to all cores — set it explicitly in a multiprocessing pipeline. TLS 2 pre-bins dense light curves internally per trial period (`BIN_FRACTION=0.03`, `BIN_FACTOR_MAX=64`), which is why user-side binning bought nothing below.

SDE thresholds: the TLS FAQ gives, from >10,000 white-noise runs, SDE 5.7 (1-FAP 0.9), 6.1 (0.95), **7.0 (0.99)**, 8.3 (0.999), **9.1 (0.9999)**; "An SDE above 9 is commonly treated as the threshold"; `results.FAP` returns NaN above 0.1 and the table is optimistic for red noise ([FAQ](https://transitleastsquares.readthedocs.io/en/latest/FAQ.html)). The TLS paper reports 93% vs 76% recovery (TLS vs BLS) of Earth-sized planets around Sun-like stars at 110 ppm/30 min, both at a 1% false-positive rate at SDE ~ 7 ([arXiv:1901.02015](https://arxiv.org/abs/1901.02015)). These are **per-light-curve** numbers; see Section 5 for the survey-level correction.

`catalog_info` caveats: `KIC_ID=5689121` raised `IndexError` (the KBONUS background star is not in the Vizier KIC-derived table TLS queries) while `TIC_ID=120686465` and `KIC_ID=5689116` worked; TLS ships Claret limb-darkening tables (`ld_claret_tess.csv`, `JAA546A14limb1-4.csv` for Kepler). For KBONUS stars, query by TIC.

**Measured TLS runtimes** (synthetic 4-yr Kepler-like LC, 67,535 points, 1,200 ppm noise, injected 2,000 ppm / 2.3 h at 13.769 d; 4 shared vCPUs, other agents loading the machine):

| run | config | periods | seconds | result |
|---|---|---|---|---|
| A | 1-30 d, `use_threads=2` | 185,100 | 121.7 | P 13.76952, SDE 116, SNR 36.6, 106 transits |
| B | 1-10 d, 2 threads | 146,253 | 116.2 | found **P/2 = 6.8847 d** (alias: the true period was outside the range) |
| C | binned to 1 h (33,121 pts), 1-30 d, 2 threads | 185,099 | 121.6 | SDE 106 — no speed-up vs A (TLS 2 pre-bins internally) |
| D | 1 thread, `oversampling_factor=2`, `duration_grid_step=1.15` | 123,400 | 109.7 | SDE 117 |
| E | 2 threads, `backend='fused'` (exact) | 185,100 | 363.8 | SDE 116 |
| real LC | KIC 5689121, 45,648 pts, 1-30 d, 1 thread | 180,849 | 198 / 191 | see below |
| real LC | same, `use_threads=4` | 180,849 | 100 | identical result |

Lessons: restricting the period range saves little (most trial periods are short); the only cheap levers are `oversampling_factor` (linear) and threads (~2x from 1 to 4 under load); binning does not help; the exact backend costs 3x. Peak RSS of the TLS main process: 225 MB (`s09_rss.py`).

Our TLS on the real cleaned LC: **P = 13.769643 +- 0.001299 d, T0 = BKJD 136.30486, SDE 105.3, SNR 27.5, depth 1,999 ppm (TLS template bottom), duration 1.55 h, rp/rs 0.040, odd/even mismatch 0.21 sigma, 104 transits (74 distinct with data), FAP 8.0e-5 (table floor)**; the next four peaks are all aliases within 0.001 d. This reproduces Kuber's 13.769448 d detection from the public product in ~3 minutes on one core.

### 3.3 Semi-coherent and stack-slide BLS (eumemic `fastbls.py`, used by Kuber)

Problem: multi-year TESS coverage forces a coherent period grid ~ baseline/season times finer than one season needs. eumemic's `fastbls.py` (numba `@njit`) implements (a) **semi-coherent**: per-season BLS Δχ² on a grid set by the longest season, summed over seasons separately per duration, then the top 40 local maxima are refined coherently on all data (`search.semicoherent_peak`); (b) **stack-slide**: each season folded once per coarse frequency, then for each fine frequency the folded sums are shifted by the phase accumulated since the season midpoint and summed — phase-coherent across seasons at near-coarse cost (`_stackslide`). Kuber's `pilot.py` imported these (`seasons(gap=20 d, max_len=90 d)` to cut Kepler into ~quarter-length blocks, `semicoherent_peak(rho=5, 0.5-25 d, k_peaks=40)` then `refine`). eumemic's injection test at nominal SNR 10 on 300 TESS stars ([README](https://github.com/eumemic/exoplanets)): coherent astropy BLS recovered 60% (31% at SNR >= 9), median 17.5 s/star; stack-slide 62% (29%), 8.9 s; semi-coherent 52% (26%), 4.2 s; "6 were found only by the coherent search and 2 only by stack-slide". The MIT-licensed `fastbls.py` is 242 lines and is the single most valuable piece of code to port for a 4-core Kepler search (the two numba kernels `_season_scan`, `_coherent_best`, `_stackslide`).

### 3.4 GPU, correlated-noise and matched-filter alternatives

- **cuvarbase** (`eebls_gpu`, `eebls_transit_gpu` with Keplerian `qmin_fac/qmax_fac/rho` duration limits; the Keplerian assumption "reduces the number of frequencies by over 1,000" for a 10-yr baseline; docs tested CUDA 8 and Python 2.7-3.6) — effectively unmaintained (last PyPI 2023) and irrelevant without a GPU ([docs](https://johnh2o2.github.io/cuvarbase/bls.html), [install](https://johnh2o2.github.io/cuvarbase/install.html)).
- **nuance** (Garcia et al. 2024, [arXiv:2402.06835](https://arxiv.org/abs/2402.06835); [docs](https://nuance.readthedocs.io/en/latest/markdown/ls.html)): JAX/tinygp; `linear_search(time, flux, gp=None, X=None, model=None, positive=True, progress=True, backend=None, batch_size=None)` over epochs and durations gives single-event log-likelihoods under a GP noise model; `periodic_search(epochs, durations, ls, snr_function)(periods)` combines them; mask with `utils.phase` and re-run for a second planet. Paper numbers: best method in 93% of 438 rapidly rotating M-dwarf TESS cases vs five detrend+BLS pipelines; TOI-540 (5 sectors) took 5 min 35 s on 12 Apple M2 Max cores vs 1 min 49 s for 15 biweight windows + BLS; the authors "expect nuance to perform between one and two orders of magnitude slower" and "do not recommend using nuance in the general case" but for stars with specific correlated noise. CPU works (`backend="cpu"`). Not installed here (JAX stack).
- **SPOC/TPS-style wavelet matched filter** (Jenkins 2002, ApJ 575, 493 — reference from memory, not re-fetched): no maintained public Python implementation; the Kepler SOC pipeline source is MATLAB/Java and NASA says readers are not expected to be able to build or run it ([nasa/kepler-pipeline](https://github.com/nasa/kepler-pipeline)). The practical substitute is LEO-Vetter's SES/MES time series (`TCELightCurve.get_SES_MES`) plus the Robovetter-style metrics (Section 5.5).
- **ML triage**: **ExoMiner** ([nasa/ExoMiner](https://github.com/nasa/ExoMiner)) ships a Podman image (`podman pull ghcr.io/nasa/exominer:latest`) and `exominer_pipeline/run_pipeline.py` with `--tic_ids` / `--tic_ids_fp` (CSV with `tic_id, sector_run`), `--data_collection_mode 2min|FFI`, `--exominer_model single|cv_ensemble|full_cv_ensemble|<keras path>`, `--num_processes`, `--download_spoc...`; but it scores **TESS SPOC TCEs** (it downloads SPOC light curves and DV XMLs), so it cannot score a KBONUS or a home-made signal; amd64/arm64 environment files suggest CPU inference, not stated explicitly; NOSA licence. **Astronet** (archived 2023-08-10) lists no pretrained weights; **Astronet-Triage-v2** code is at mdanatg/Astronet-Triage with training data on Zenodo 7411579 and a README that trains your own model (recall 99.6% at 75.7% precision in Tey et al. 2023, [arXiv:2301.01371](https://arxiv.org/abs/2301.01371)); **qusi** 1.5.6 is a framework with no bundled models evident. Conclusion: no off-the-shelf CPU model will triage arbitrary Kepler-background or TESS-FFI signals today; rule-based LEO-Vetter + an inverted-trained classifier (eumemic) is the realistic route.

---

## 4. Multi-planet iteration, single transits, TTVs

- **Mask-and-repeat**: all three pipelines mask |phase| < 1-2 x duration (`BoxLeastSquares.transit_mask` or `wotan.transit_mask`, Chopra factor 2.0, eumemic/Kuber 1.0 x duration on the detrended series) and re-run up to 3 signals (Chopra `MAX_SIGNALS=3`; eumemic 3, 4 for TOI hosts; Kuber 2). Stop rules: SDE < 7 or SNR < 6 (Chopra), BLS SNR < 7 (eumemic), SNR < 9 (Kuber). eumemic's `repeat` test (same period within 0.5%, or half/double within 0.2%, re-found after masking) marks leftover variability, not a new planet; `known_alias` catches 1/3...3 harmonics of a masked TOI; TOI masks are widened by 3 sigma of the propagated ephemeris error `sqrt(eT0^2 + (n eP)^2)` (`toi_windows`). TLS: `transit_mask` then `cleaned_array` then a new `transitleastsquares(...)` object.
- **Single-transit search**: nuance's `linear_search` is literally a single-event scan; LEO-Vetter computes the SES series; MonoTools 1.0.0 fits mono/duo-transits with a PyMC5/pytensor + celerite backend ([repo](https://github.com/hposborn/MonoTools)); SHERLOCK's `MoriartySearch` uses the SANTO neural network ([docs](https://sherlock-ph.readthedocs.io/en/latest/_autosummary/sherlockpipe.single_transits.search.MoriartySearch.html)). The simplest robust version is Kuber's `independent_vet.local_events`: per event, fit intercept + slope + exposure-integrated box in a window of max(3 x duration, 0.15 d), requiring >= 4 out-of-transit points on each side, with noise floored by the local MAD; a single-event scan is the same fit slid along the time axis.
- **TTV tolerance**: QATS (Carter & Agol 2013) allows the inter-transit spacing to vary within [DeltaMin, DeltaMax]; the only code is C++ (`make`, `./call_qats data.dat 794 797 14`, `call_qats_indices` for transit indices; data must be detrended, median-subtracted, sigma-normalised with gaps zero-filled; cite Carter & Agol 2012 and Kel'Manov & Jeon 2004) ([readme](https://faculty.washington.edu/agol/QATS/readme)). Cheaper in practice: the semi-coherent BLS tolerates phase drift within a season by construction, and the coherent refine step with `oversample >= 3` plus a fine-period scan (eumemic `refine`: +-3 coarse steps in 241 sub-steps, 12 durations) absorbs TTVs of order the cadence.

---

## 5. Statistics of large searches

### 5.1 The look-elsewhere effect, measured on this star (ran: `s07_null_throughput.py`)

Identical cleaning + astropy BLS (`objective='snr'`, 499k periods, 6 durations, `oversample=3`) on the real LC, its inversion (f -> 2 - f) and seven quarter-scrambled copies (quarter order permuted, each quarter's series circularly shifted): the **real** signal gives SNR 30.3 (global SDE 29.3) at 13.76953 d with 74 transits; the **inverted** LC gives SNR 11.0 (SDE 6.1) at 29.04 d; the seven **scrambles** give SNR 9.6-11.9 (SDE 5.1-6.8) at 19.7-29.0 d, all with the longest trial duration (4.56 h) and 35-54 "transits". So on a 45k-point, 1,000-ppm Kepler light curve searched over half a million periods, a BLS depth-SNR of 10-12 is what noise produces, and the noise peaks pile up at the long-period / long-duration corner of the grid — exactly the 12-25 d pile-up Chopra reported and the reason Kuber's "BLS SNR >= 9" first screen passed inverted stars as often as real ones (31 signals/29 stars real vs 32/28 inverted over 1,736 stars; his [kbonus-extension README](https://github.com/Kuberwastaken/planethunt/blob/main/analysis/kbonus-extension-20261011/README.md): 10/9 vs 10/9 for ranks 65-576 and 19/18 vs 22/19 for ranks 577-1800). Depth/depth_err from BLS is a *local* statistic that ignores the number of trials; only an empirical null (below) converts it into a survey false-alarm rate.

### 5.2 Calibrating thresholds with inverted and scrambled light curves

- **Chopra** (526 stars, ~1,500 searches, [GUIDE.md](https://github.com/paraschopra/tess-planet-candidates/blob/main/docs/GUIDE.md), `06_inverted.py`): flip the raw flux *before* cleaning so the asymmetric clipping is applied identically; run the full search and the full triage. Real vs flipped: 743 vs 752 signals pass the search; 31 vs 39 have no hard flags; 10 vs 17 are flag-free; strongest flag-free flipped signal SDE 8.4 vs his candidates at 11.3 and 10.3; no flag-free flipped signal reached SDE 9. A second null scrambled each candidate's own light curve 38 times (mask transits, slide each sector randomly) and reached SDE 8.5 and 9.3. His review then **withdrew** the extrapolated false-alarm counts (0.01 and 0.2 per survey): "38 samples can't support a tail extrapolation"; zero fakes in 526 stars only caps the rate at about 3 per survey (95%).
- **Kuber**: matched inverted run on every star; at the candidate's own strength (P > 1.5 d, BLS SNR >= 26.7, holdout SNR >= 15.5) two inverted stars (KIC 8701308, KIC 5283912) still pass — detection strength alone did not single out KIC 5689121.
- **eumemic** ([validation/false_alarm_classifier.md](https://github.com/eumemic/exoplanets/blob/main/validation/false_alarm_classifier.md), README Validation): inverted light curves of **798,769** random full-sky stars (two samples of 399,408 and 399,361) through the identical pipeline (`EXO_INVERT=1` negates the detrended residuals); expected false alarms in the real search = inverted count x (4,039,203 / 798,769). Among signals passing LEO-Vetter and the cuts with clean pixels: SNR >= 9 real 726 vs 96 expected false alarms; SNR 7.5-9: 314 expected. False alarms concentrate at low MES: 46% of clean SNR >= 9 signals with MES < 10 are expected to be noise, 16% at MES 10-12, 2% above — so they require MES >= 10 for tier N and moved 78 candidates with MES < 10 to tier C. A `HistGradientBoostingClassifier` (features: BLS SNR, MES, SDE, duration/expected duration, n good events, max single-event fraction, odd/even, secondary, per-sector depth consistency, LEO metrics CHI/DMM/Fred/SHP/chases, Tmag, contamination, radius, Teff; centroid metrics deliberately excluded because inversion does not invert centroids) scores real vs inverted with inverted rows weighted by the star ratio; at a 0.80 cut 69 low-band signals are selected with ~6 expected false alarms (92% purity), validated with TIC-grouped folds after an outside contributor found star leakage between folds. On the earlier 9,979 K/M-dwarf search, 1,500 inverted stars gave 157 pre-filter passes, 1 passing all LEO-Vetter tests, 0 with consistent per-sector depths.
- Caveats everyone states: inversion measures noise and systematics, not astrophysical false positives (background EBs); asymmetric clipping and selection rules tuned after seeing real signals make the null imperfect; report "in one inverted survey through the identical pipeline no signal reached X" rather than a p-value.

### 5.3 Injection-recovery for completeness

Chopra `inject_test.py` and eumemic `inject_test.py` inject a trapezoid (15% ingress) into the **raw** flux at the physically expected duration (`expected_duration(P, rho) x U(0.6, 1.0)`), depth set for a nominal SNR, then run the untouched pipeline; a hit = period within 0.3% of P, P/2 or 2P. eumemic: nominal SNR 10 injections are recovered at SNR >= 9 only 22% of the time and at >= 7.5 48% ("red noise in full-frame light curves lowers the measured SNR"). Use the same machinery to pick the detrending window (Chopra's 0.75 d) and to calibrate vetting flags: two of Chopra's red-flag tests fired on 15-25% of planted planets until the noise model was fixed (next item).

### 5.4 Empirical noise with random off-transit boxes (Chopra's fix) and Kuber's phase controls

Chopra `vet.depth_emp`: measure the box depth at the fixed ephemeris, then at 40 random off-transit phases (`NULL_OFFSETS`), and take 1.4826 x MAD of those 40 depths as the error; the Codex audit found the controls had been detrended with only the real transits masked, so each control "helped fit its own trend" (noise ~10% low) — `depth_emp_matched` re-detrends per control phase with that phase masked, and the candidates' significance fell from ~8.2 sigma to 7.2 and 7.7 sigma. Kuber's `independent_vet.vet(controls=24|64)` fits the per-event local-baseline box at `controls` phases excluding +-max(0.06, 2 x duration/P) around phase 0 and 0.5, excluding the real transit cadences, and reports `depth_over_control_scatter` and an add-one rank `(1 + #controls >= depth)/(1 + N)` with the caveat "rank is not a survey false-alarm probability". Both also report `max_positive_depth_contribution` / `max_ses_frac` (a single event may not carry > 0.5 of the signal) and a leave-one-event-out SNR.

### 5.5 Model comparison: transit vs sinusoid vs EB

- BLS `compute_stats['harmonic_delta_log_likelihood']` (> 0 favours a sinusoid; -430 for our candidate) and `depth_odd/depth_even/depth_phased` give the first-pass EB tests.
- LEO-Vetter 1.2.0 (installed; Kunimoto et al. 2025, AJ 170, 280) computes 121 metrics in 8.2 s on the 45k-point LC (`s08_leo.py`): `TCELightCurve(ID, time, raw, flux, flux_err, period, epoch, duration)`, `compute_flux_metrics(star)` needs `star = dict(rad, e_rad, mass, e_mass, Teff, e_Teff, rho, u1, u2, id, tic, ra, dec, raDegrees, decDegrees, planetData=None)` (keys grepped from the source; `sector, cam, ccd` for pixel vetting, TESS only), then `check_thresholds(metrics, "FA")` / `"FP"`. Default thresholds (`thresholds.py`): MES 6.2, N_transit 3, SHP 0.6, MS1 0.2, MS2 0.8, MS3 0.8, chases 0.78, DMM 1.5, max_SES_to_MES 0.8, AIC1 -60, AIC2 -30, SWEET 15, ASYM 8, CHI 7.8, frac_gap 0.5, V_shape 1.5, size 22, offset 15". For KIC 5689121: MES 24.2 (new_MES 23.0 after bad-event removal), 74 transits, SHP 0.17, CHI 32.3, DMM 0.96, Fred 1.27, sig_pri 24.4 vs sig_sec 3.0 / sig_ter 3.0 / sig_pos 4.8, max_SES 5.6, SWEET sine 0.54, transit AIC 187 vs line AIC 752, Rp 2.80 +- 0.13 R_earth, Teq 381 K — **no FA and no FP flag fires**. (eumemic replaced LEO's per-cadence loop with a compiled equivalent for 2-5x faster vetting.)
- exovetter 0.0.16 (installed) exposes `ModShift, Sweet, Lpp, OddEven, Centroid, TransitPhaseCoverage, VizTransits, LeoTransitEvents` with `vetter.run(tce, lightcurve)` on a `Tce(period=P*u.day, epoch=..., epoch_offset=exo_const.bkjd, depth=..*exo_const.ppm, duration=..*u.hr)` (from the installed `tce.py` docstring).
- Full Bayesian evidence (transit vs sinusoid vs EB with secondary) via nested sampling: juliet (dynesty/multinest, `dataset.fit(n_live_points=300)`, priors `P_p1, t0_p1, r1_p1/r2_p1 or p_p1, b_p1, q1_/q2_<inst>, ecc_p1/omega_p1, rho, mdilution_<inst>, mflux_<inst>, sigma_w_<inst>, GP_sigma/GP_timescale/GP_rho` — [priors page](https://juliet.readthedocs.io/en/latest/user/priorsnparameters.html)) or dynesty/ultranest directly. Budget minutes per candidate, not per star; TLS's SDE and LEO's AIC are the batch-level proxies.

---

## 6. Transit fitting

### 6.1 batman + emcee with 30-min exposure integration (ran: `s06_fit_vet.py`, 109 s)

```python
import batman, emcee
params = batman.TransitParams(); params.t0=T0; params.per=P; params.rp=0.046; params.a=39.7
params.inc=89.3; params.ecc=0; params.w=90; params.limb_dark='quadratic'; params.u=[0.2125, 0.3984]
m = batman.TransitModel(params, t, supersample_factor=7, exp_time=1765.464/86400)  # Kepler LC exposure in days
# free: t0, per, rp/rs, a/rs, b (inc = arccos(b/aRs)), f0, ln jitter; 32 walkers x 2500 steps
```

Kipping 2010 ([arXiv:1004.3741](https://arxiv.org/abs/1004.3741)): finite integration (Kepler 30 min) distorts the shape and systematically underestimates the light-curve stellar density unless the model is integrated numerically; batman's `supersample_factor`/`exp_time` do that ([tutorial](https://lkreidberg.github.io/batman/docs/html/tutorial.html)). Result on 2,637 points within +-0.4 d of transit: **t0 = BKJD 797.24544 (+0.0012/-0.0011) = BJD 2455630.2454 (Kuber: 2455630.245431), P = 13.769444 +- 0.000038 d (Kuber: 13.769448 +- 0.000065), rp/rs = 0.0443 (+0.0017/-0.0010), a/rs = 43.7 (+2.9/-7.5), b = 0.35 (+0.28/-0.24), rho_star = 5.9 (+1.3/-2.6) rho_sun vs 3.9 expected for 0.505 Msun/0.507 Rsun (consistent, as Kuber found: 4.4 +2.3/-2.8), Rp = 2.45 (+0.09/-0.06) R_earth for R* = 0.507 Rsun without dilution correction, (rp/rs)^2 = 1,959 (+152/-91) ppm.** Acceptance 0.38. Add a `rho_star` prior from TIC (`catalog_info(TIC_ID=...)` or `astroquery.mast.Catalogs.query_object('TIC 120686465', catalog='TIC')`) and a dilution term `mdilution` when the aperture contains the brighter neighbour.

### 6.2 Limb darkening

- Fast path: `transitleastsquares.catalog_info(TIC_ID=...)` returns Claret quadratic coefficients (Kepler: Claret & Bloemen 2011 table shipped as `JAA546A14limb1-4.csv`; TESS: `ld_claret_tess.csv`). For TIC 120686465 (Teff 3633 K): u = (0.2125, 0.3984). The KIC path failed for this background star.
- LDTk 1.9.3: `from ldtk import LDPSetCreator, BoxcarFilter; sc = LDPSetCreator(teff=(3633,160), logg=(4.73,0.1), z=(0,0.1), filters=[...]); ps = sc.create_profiles(); cq, eq = ps.coeffs_qd(do_mc=True)`; it "downloads the uncached stellar spectra from the Husser et al. server automatically" (PHOENIX specific intensities; expect hundreds of MB on first use) ([README](https://github.com/hpparvi/ldtk)). Installed here with `--no-deps` because of its pandas 1.5.3 pin; not run (download).

### 6.3 Other fitters (status)

- **exoplanet 0.6.0**: `pip install "exoplanet[pymc]"`, "updated to work with PyMC version 5" ([install docs](https://docs.exoplanet.codes/en/latest/user/install/)); current PyMC is 6.3.2 and exoplanet-core 0.4.0 was released 2026-07-25 — compatibility with PyMC 6 is **unverified** here.
- **juliet 2.2.10** installed (dynesty backend; celerite GP kernels); good for evidence comparison; heavy dependency tree.
- **PyTransit 2.9.3** (`RoadRunnerModel('quadratic').set_data(times).evaluate(k, ldc, t0, p, a, i)`, Numba) — active (2026-10-15) but TRICERATOPS 1.1.0 pins `pytransit==2.2`, so keep TRICERATOPS in its own venv (Chopra hit the same pins).
- **allesfitter 1.2.10**: last release 2022-07-25; use only if you need its EB/RV joint models.

---

## 7. Pixel-level tools (ran: `s05_pixels.py`)

```python
tpf = lk.search_targetpixelfile('KIC 5689116', author='Kepler', cadence='long', quarter=10).download(quality_bitmask='default')
# KeplerTargetPixelFile 4x5 px, 4,447 good cadences (623 flagged), channel 34, 1.70 MB, 3.2 s
xy_t = tpf.wcs.world_to_pixel(SkyCoord(285.907827, 40.941988, unit='deg'))   # (1.15, 0.08)
xy_n = tpf.wcs.world_to_pixel(SkyCoord(285.905233, 40.942079, unit='deg'))   # (2.26, 1.47); 7.06" apart = 1.8 px
tpf.plot_pixels(aperture_mask=tpf.pipeline_mask, show_flux=True)            # per-pixel light curves
col, row = tpf.estimate_centroids(aperture_mask='all')                       # flux-weighted centroids
from lightkurve.prf import KeplerPRF
prf = KeplerPRF(channel=tpf.channel, shape=tpf.shape[1:], column=tpf.column, row=tpf.row)
model = prf(center_col=tpf.column+xy_t[0], center_row=tpf.row+xy_t[1], flux=1, scale_col=1, scale_row=1, rotation_angle=0)
tc = lk.search_tesscut('TIC 120686465', sector=14).download(cutout_size=11)   # 1,237 x 11 x 11, 30-min, 2.98 MB, 38 s
```

Docs: [plot_pixels](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.KeplerTargetPixelFile.plot_pixels.html), [estimate_centroids](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.KeplerTargetPixelFile.estimate_centroids.html), [search_tesscut](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.search_tesscut.html); `lightkurve.prf.KeplerPRF` works but the `tpfmodel` submodule needs `oktopus` (warning printed). astroquery's `Tesscut.get_sectors / get_cutouts / download_cutouts(coordinates|object_name, size, sector)` recommends cutouts <= 30 px and notes TICA cutouts were discontinued in August 2025 — all cutouts are now SPOC FFIs ([docs](https://astroquery.readthedocs.io/en/latest/mast/mast_cut.html)).

What a naive difference image shows: averaging the 33 in-transit vs 128 near-transit out-of-transit cadences of the 7 Q10 transits gives a per-pixel max SNR of **0.4** (target pixel 0.2, neighbour pixel 0.1) — a single quarter of a 1,000-ppm, Kp 16 star cannot localise a 0.2% dip. Kuber's `pixel_localize.py` is the template for what works: for each event fit intercept + slope + exposure-integrated box **per pixel** with noise floored by the local MAD, inverse-variance combine events (cluster-robust error), stack all quarters, then fit Kepler PRF templates (`lightkurve.prf.KeplerPRF`, calibration files cached locally) at the two Gaia positions with a common WCS offset calibrated on the mean image; compare χ² (his Δχ² = 934 in favour of the target, with the stated caveat that pixel covariances and PRF mismatch are not calibrated). Chopra's TESS version (`10_pixels.py`): 15x15 TESScut per sector, difference = (before + after)/2 - in-transit per transit, per-pixel noise from 40 off-transit phases, joint all-sector Gaussian-spot fit at a common sky position with per-sector amplitude, Δχ² against each Gaia neighbour moved to the sector epoch — after his audit found the first version "could locate a brightening as easily as a dimming". eumemic's per-sector centroid F-test (transit-diffImage PRF fits in up to 12 sectors, scatter floored at 1" per axis) never rejects the true host at p < 0.001 on 565 known planets and rejects it 73% of the time for 175 neighbour-hosted signals ([centroid_sectors.md](https://github.com/eumemic/exoplanets/blob/main/validation/centroid_sectors.md)).

Tools: **transit-diffImage** (`tessDiffImage` class with `make_ffi_difference_image`, `make_planet_difference_image`, `get_tess_cut`, `find_transits`, `make_difference_image`, `make_stellar_scene`, `draw_difference_image`; Kepler via `keplerDiffImage.py` which "Requires" Kepler-RaDex2Pix; install by cloning and `pip install .`); **LEO-Vetter `pixel_vetting(tlc, star, sectors)`** wraps it (TESS only); **psfmachine** (`TPFMachine.from_TPFs(tpfs, n_r_knots=10, n_phi_knots=12).fit_lightcurves()`, Gaia-seeded PRF photometry — the engine behind Kepler Bonus) needs a Python <= 3.11 venv here. For TESS sectors with 200-s FFIs an 11x11 cutout will be ~9x the 3 MB measured for 30-min data.

---

## 8. Orchestration at scale on 4 CPUs

### 8.1 Pipeline design (what all three working pipelines converge on)

1. **Targets**: catalogue query frozen before any photometry is seen, written to CSV with a rank (Kuber sorted by Gaia G; Chopra by expected SNR of a 1.5 R_earth planet).
2. **Download -> cache**: one MAST POST for all product URLs; workers fetch (processes, not threads: eumemic saw a 3x gain because FITS parsing is GIL-bound); verify `SIMPLE` header and SHA-256; write one compressed `.npz` per star (`time, flux, flux_err, sector/quarter, exptime, cx, cy`) plus a manifest JSON. Sizes measured: KBONUS LC FITS 8.44 MB, cleaned npz 1.39 MB, Kepler TPF 1.70 MB/quarter, TESScut 11x11x1,237 2.98 MB.
3. **Detrend + clean**: per segment biweight (0.5-1.0 d), asymmetric clip, empirical errors (<0.2 s per star).
4. **Search**: coarse semi-coherent/stack-slide BLS on 20-min bins -> coherent refine of top peaks on 10-min bins -> `event_stats` (per-event depths, odd/even, secondary, max single-event fraction) -> up to 3 signals by masking; TLS only on survivors.
5. **Triage**: hard flags (odd/even, secondary, repeat period, duration vs expected, single-event dominance, gap edges, per-season depth consistency, broad depression, same period in >= 3 other stars, too deep) and soft flags (period in the 12-25 d and 0.5-1 d pile-ups) — Chopra's list; LEO-Vetter FA/FP; MES >= 10 (eumemic).
6. **Vet**: holdout seasons reserved **before** the search (Kuber: last max(2, ceil(0.25 N)) quarters; Chopra: hide-a-season after the fact, labelled retrospective), fixed-ephemeris local-baseline depths, phase controls, SAP/PSF consistency, pixel localisation, catalogue/literature matching, inverted null of the whole sample.
7. **Provenance**: per-star JSON results keyed by `input_sha256` and `pipeline_version` so reruns skip unchanged work (Kuber `run_search.py`), `skip if <out>/<tic>.json exists` (eumemic, Chopra), and for anything you publish a DAG that re-runs downstream stages when code or inputs change (eumemic `pipeline/dag/`, `dag.lock.json`).

Parallelism: `ProcessPoolExecutor`/`multiprocessing.Pool` with `imap_unordered(chunksize=1)`, **one process per core and `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`** set before importing numpy (all three pipelines do this), `use_threads=1` for TLS inside workers; print progress as JSON lines.

### 8.2 Measured per-star costs (this box, 4 shared vCPUs, load average 5-13 from other agents, so treat as upper bounds)

| stage | Kepler 4-yr LC (46-67k pts) | TESS 5 sectors (19k pts, binned to 20 min = 9.4k) |
|---|---|---|
| download | 8.4 MB FITS, seconds (MAST-bound) | 5 x ~2 MB FITS |
| wotan + clean | < 5 s incl. JIT | < 1 s |
| coherent astropy BLS, density grid, 6 durations, oversample 3 | 443 s/star in 4-worker test (152 s solo) -> **0.5 stars/min on 4 workers** | 50 s/star -> **4.7 stars/min on 4 workers** (265k periods) |
| semi-coherent numba BLS (eumemic/Kuber code, not re-run here) | eumemic: 4.2 s median per TESS star (12 vCPU box); Kuber's Codex run: 5 workers | 4.2 s (semi), 8.9 s (stack-slide), 17.5 s (coherent) |
| TLS 1-30 d (survivors only) | 100 s at 4 threads / 190 s at 1 thread | ~20-40 s (fewer points) |
| batman+emcee fit (32 x 2,500) | 109 s | similar |
| LEO-Vetter flux metrics | 8 s | 3-8 s |
| RAM | BLS worker 115 MB; TLS 225 MB; 4 workers < 2 GB | smaller |

### 8.3 Run-time estimates (search stage only; add ~10% for vetting of the few-percent survivors)

- **Fast tier** (semi-coherent numba BLS at ~10 s per Kepler star or ~5 s per TESS star, per core): 1k stars ~45 min / ~20 min; 10k ~7 h / ~3.5 h; 100k ~3 days / ~1.5 days on 4 cores -> overnight is fine up to ~10k, 100k wants a cloud box.
- **Full coherent astropy-BLS tier** (~150-450 s per Kepler star): 1k stars 10-30 h; 10k 4-13 days; 100k ~8,000-12,000 core-hours -> cloud only. TESS 5-sector stars at 50 s: 1k 3.5 h, 10k 35 h, 100k ~15 days on 4 cores.
- **TLS on every star** (~190 s/core for Kepler): 1k 13 h, 10k 5.5 days — reserve TLS for signals that pass the BLS screen.
- Kuber's actual 1,800-star KBONUS run used the semi-coherent code with 2-5 workers over two days of wall time (per his READMEs); eumemic reports 200 TESS stars in 70 s on a 16-core AWS c8g.4xlarge after optimisations (trend on 10-min means 14x faster, division-free stack-slide kernel 1.42x, Graviton4 2.6x) and downloads 20x faster inside AWS us-east-1 than from MAST.
- Disk: 10k Kepler stars = 84 GB of KBONUS FITS or 14 GB of npz (delete the FITS after parsing, keep SHA-256 and URL in the manifest); TPFs only for candidates (13 quarters x 1.7 MB = 22 MB per star).

What to run overnight locally: target selection, downloads (network-bound, 4 processes), detrending, the fast-tier search of up to ~10k stars, the inverted null of the same stars (doubles the cost), triage, and TLS + fits for the handful of survivors. What to push to a cloud box (or Modal): coherent full-grid BLS on > 1k Kepler stars, TLS on > 1k stars, the 100k-star tiers, and TRICERATOPS/MOLUSC runs (eumemic ran those on AWS spot instances with a separate pinned venv).

---

## 9. Minimal reproducible end-to-end (what `s04_kbonus_e2e.py` does, 1 core, ~7 min total)

load KBONUS per-quarter FLUX with `SAP_QUALITY==0` -> normalise per quarter -> wotan biweight 1.0 d (`break_tolerance=0.5`) -> clip +3/-10 sigma -> density-aware grid (499k periods, 0.5-30 d) -> astropy BLS `objective='snr'` (186 s) -> `compute_stats` -> TLS 1-30 d with TIC limb darkening (198 s) -> P = 13.76964 d, SDE 105, 74 transits, odd/even 0.2 sigma. Then `s06_fit_vet.py` (batman+emcee, 109 s) and `s08_leo.py` (LEO-Vetter, 8 s) give the physical parameters and a clean rule-based vetting pass, and `s07_null_throughput.py` shows the same light curve inverted/scrambled never exceeds BLS SNR 12 versus 30 for the real signal. None of this is a validation (no pixel localisation at depth, no FPP), but it is the complete, free, 4-CPU toolchain, with every tool's version pinned above.