# Algorithms and tooling for a transit search at scale on a 4-CPU box (verified 2026-10-17)

Scope: the free, CPU-only toolchain for (1) loading and cleaning Kepler / K2 / TESS light curves, (2) detrending, (3) period search, (4) multi-planet / single-transit / TTV handling, (5) survey statistics (look-elsewhere, inverted and scrambled nulls, injection-recovery), (6) transit fitting, (7) pixel-level tools, (8) orchestration and run-time budgets. Every number below was either run in the sandbox venv (`/home/user/ZyushgPlanet/.venv`, Python 3.13.16) or read from a primary source; `[verified]` means re-run or re-fetched during this fact-check, `[reported]` means taken from the previous agent's logs in `scratchpad/algorithms-tooling/` and not re-run, `[unverified]` means neither. Scratch files of this check are in `/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/verify-algorithms-tooling/` (`verification_log.txt`, `tls_in_pool.py`). Kuber's, Chopra's and eumemic's repos were **read, never executed or imported**, from `scratchpad/untrusted/planethunt` and `scratchpad/algorithms-tooling/upstream/`.

**Two releases days old change every 2024-25 tutorial.** TLS 2.0 (6 Oct 2026) and 2.0.1 (11 Oct 2026) and wotan 2.0 (11 Oct 2026) are major rewrites ([TLS releases](https://github.com/hippke/tls/releases), [wotan releases](https://github.com/hippke/wotan/releases)) `[verified]`. TLS 2: "18 bug fixes, 10 of which change results", new default numba backend `fused-pl` (piecewise-linear template, validated against the exact statistic on hundreds of injection-recovery cases), `backend="fused"` reproduces TLS 1.x up to rounding, batman is no longer a runtime dependency, opt-in `power(SDE_detrend="hybrid")`, "30-60x faster than TLS 1.33 on warm benchmarks". wotan 2.0: 31 bugs fixed (BUG-001..031 in [BUGLOG.md](https://raw.githubusercontent.com/hippke/wotan/master/BUGLOG.md)), 563-test suite; speed-ups measured on Kepler-10 (52k points, 0.5-d window, 2 cores): biweight 60x, median 21x, mean 19x, hspline 14x, lowess 9x, rspline 2.3x, cosine 1.1x. Result-changing wotan fixes: slider windows now include their last sample (BUG-002), non-positive flux is kept (BUG-003), biweight tolerance is unit-independent (BUG-009), `slide_clip` defaults really use MAD/median (BUG-006), robust clipping is two-sided (BUG-028). The venv has both. Kuber's search ran in "Python 3.10 with unpinned recent numpy/scipy/astropy/lightkurve" (his kbonus-extension README), i.e. TLS/wotan 1.x — his numbers are not bit-reproducible with TLS 2.0.1.

---

## 0. Tool inventory: pip name, version, install status, docs

PyPI queried live 2026-10-17 (`curl https://pypi.org/pypi/<pkg>/json`) `[verified]`. "venv" = what `uv pip list` shows in the sandbox venv now `[verified]`.

| Tool | pip name | PyPI latest (upload) | requires_python | venv | Docs | Ran here? |
|---|---|---|---|---|---|---|
| lightkurve | `lightkurve` | 2.6.0 (2026-04-16) | >=3.9,<4 | 2.6.0 | https://lightkurve.github.io/lightkurve/reference/index.html | yes (search/download/TPF/TESScut/PRF/centroids) |
| astropy BLS | `astropy` | 8.0.2 (2026-10-15) | >=3.11 | 8.0.2 | https://docs.astropy.org/en/stable/api/astropy.timeseries.BoxLeastSquares.html | yes |
| astroquery | `astroquery` | 0.4.11 | | 0.4.11 | https://astroquery.readthedocs.io/en/latest/mast/mast.html | yes |
| TLS | `transitleastsquares` | 2.0.1 (2026-10-11) | >=3.9 | 2.0.1 | https://transitleastsquares.readthedocs.io/en/latest/Python%20interface.html | yes |
| wotan | `wotan` | 2.0 (2026-10-11) | — | 2.0 | https://github.com/hippke/wotan | yes |
| batman | `batman-package` | 2.5.3 (2025-05-05) | — | 2.5.3 (module reports `__version__` 2.5.1) | https://lkreidberg.github.io/batman/docs/html/tutorial.html | yes |
| emcee / corner | `emcee`, `corner` | 3.1.6 (2024-04-19) / 2.3.0 (2026-07-05) | — / >=3.9 | 3.1.6 / 2.3.0 | https://emcee.readthedocs.io | yes |
| celerite2 | `celerite2` | 0.3.3 (2026-07-12) | >=3.11 | 0.3.3 (+ celerite 0.4.3 via juliet) | https://celerite2.readthedocs.io/en/latest/api/python/ | yes |
| dynesty / ultranest | `dynesty`, `ultranest` | 3.1.0 (2026-07-17) / 4.6.3 (2026-10-13) | >=3.8 / any | **both installed** | https://dynesty.readthedocs.io ; https://johannesbuchner.github.io/UltraNest/ | import only |
| juliet | `juliet` | 2.2.10 (2026-02-16) | >=2.7 | 2.2.10 | https://juliet.readthedocs.io/en/latest/user/priorsnparameters.html | import only |
| exoplanet | `exoplanet` (+`exoplanet-core`, `pymc`) | 0.6.0 (2024-04-16); core 0.4.0 (2026-07-25, py>=3.12); pymc 6.3.2 (2026-09-08, py>=3.12) | | not installed | https://docs.exoplanet.codes/en/latest/user/install/ (says "updated to work with PyMC version 5"; PyMC 6 untested) | no |
| PyTransit | `pytransit` | 2.9.3 (2026-10-15) | >=3.10 | not installed | https://github.com/hpparvi/PyTransit | no |
| allesfitter | `allesfitter` | 1.2.10 (2022-07-25, stale) | | not installed | https://www.allesfitter.com | no |
| LDTk | `ldtk` | 1.9.3 (2026-10-15) | >=3.8 | 1.9.3 (PyPI metadata lists unpinned `pandas`; the earlier "pins pandas 1.5.3" claim is wrong) | https://github.com/hpparvi/ldtk | import only |
| LEO-Vetter | `leo-vetter` | 1.2.0 (2026-07-29) | — | 1.2.0 (+ lmfit 1.3.4, brokenaxes 0.6.2) | https://github.com/mkunimoto/LEO-vetter ; paper arXiv:2509.10619 ("accepted for publication in AJ"; software DOI 10.5281/zenodo.16945146) | yes (flux-level) |
| exovetter | `exovetter` | 0.0.16 (2026-04-27) | >=3.7 | 0.0.16 | https://exovetter.readthedocs.io | import only |
| TRICERATOPS | `triceratops` | 1.1.0 (2026-09-02) | >=3.8 | NOT installed: `setup.py` pins `pytransit==2.2`, `setuptools<81`, and for py<3.13 `numpy<2`, `scipy<1.14` ([setup.py](https://raw.githubusercontent.com/stevengiacalone/triceratops/master/setup.py)); dry-run today resolves 1.1.0 + pytransit 2.2.0 -> separate venv | https://github.com/stevengiacalone/triceratops | no |
| psfmachine | `psfmachine` | 1.1.4 (2023-12-15) | >=3.8,<4 | NOT installed: metadata pins `numpy<2`, `pandas<2`, `astropy<6`; dry-run would downgrade numpy 2.5.4 -> 1.26.4 and pandas 2.3.3 -> 1.5.3 (no cp313 wheel; wheels cp38-cp311) -> Python 3.10/3.11 venv | https://ssdatalab.github.io/psfmachine/ | no |
| nuance | `nuance` | 0.8.1 (2024-07-31) | >=3.9 | not installed (JAX 0.11.2, tinygp 0.3.1, jaxopt...) | https://nuance.readthedocs.io | no |
| cuvarbase (GPU BLS) | `cuvarbase` | 0.2.5 (2023-10-23) | | not installable (no GPU; docs: CUDA 8.0, Python 2.7/3.4-3.6 tested, needs PyCUDA + scikit-cuda) | https://johnh2o2.github.io/cuvarbase/bls.html | no |
| transit-diffImage | GitHub only (`pip install .`) | n/a | | not installed; needed by `leo_vetter.pixel` (`import leo_vetter.pixel` -> ModuleNotFoundError: transitDiffImage) | https://github.com/stevepur/transit-diffImage | no |
| tess-point | `tess-point` | 0.9.5.post1 (2026-09-10) | | not installed | https://github.com/christopherburke/tess-point | no |
| SHERLOCK | `sherlockpipe` | 1.5.1 (2026-05-29) | >=3.11 (docs say Python 3.10!) | not installed (earlier dry-run reportedly failed on the `ellc` build `[reported]`; docs: install numpy and Cython first) | https://sherlock-ph.readthedocs.io/en/latest/installation.html | no |
| MonoTools | `MonoTools` | 1.0.0 (2025-06-06) | | not installed; pins `lightkurve==1.11.0`, `scipy==1.12`, needs `exoplanet[pymc]`, `pymc`, `celerite2` -> own venv | https://github.com/hposborn/MonoTools | no |
| ExoMiner | Podman image via repo manifest (image name not stated in README) | no release tag | | not run (scores TESS SPOC TCEs) | https://github.com/nasa/ExoMiner | no |
| Astronet (Shallue & Vanderburg) | GitHub (archived 2023-08-10) | n/a | | no weights | https://github.com/google-research/exoplanet-ml | no |
| Astronet-Triage(-v2) | GitHub | n/a | | README trains your own model; Releases empty; data https://zenodo.org/record/7411579 | https://github.com/mdanatg/Astronet-Triage | no |
| qusi (RAMjET) | `qusi` | 1.5.6 (2025-04-21) | >=3.9 | not installed; a framework, no bundled models evident `[unverified]` | https://qusi.readthedocs.io | no |
| QATS (Carter & Agol) | C++ tarball only | n/a | | not installed. PyPI `qats` 5.5.0 (2026-10-14) is DNV's unrelated library | https://faculty.washington.edu/agol/QATS/ | no |
| Kepler SOC pipeline / Robovetter | GitHub | n/a | | MATLAB/Java; Robovetter C++11 needs DR25 metric inputs | https://github.com/nasa/kepler-pipeline ; https://github.com/nasa/kepler-robovetter | no |

Other ecosystem packages on PyPI (not installed): `foldedleastsquares` 1.1.12 (SHERLOCK's TLS fork), `lcbuilder` 0.27.0 (py>=3.11), `dearwatson` 1.1.0 (WATSON vetting, py>=3.11), `tglc` 0.8.1 (**py<3.13, pins `wotan~=1.9`, `numpy<2`, `astroquery==0.4.7`; a dry-run install into the search venv would downgrade wotan to 1.10 and numpy to 1.26.4** `[verified]`), `eleanor` 2.0.5 (2022), `stella` 0.1.0 (2021), `astrobase` 0.5.3 (2021), `jax` 0.11.2, `tinygp` 0.3.1, `numpyro` 0.22.0. Not on PyPI: `molusc`, `exominer`, `tesscut` (use `astroquery.mast.Tesscut`).

**Rule:** the search venv holds exactly lightkurve 2.6.0, astropy 8.0.2, astroquery 0.4.11, transitleastsquares 2.0.1, wotan 2.0, numba 0.68.0, numpy 2.5.4, scipy 1.18.1, pandas 2.3.3, batman-package 2.5.3, emcee 3.1.6, corner 2.3.0, celerite2 0.3.3, dynesty 3.1.0, ultranest 4.6.3, juliet 2.2.10, leo-vetter 1.2.0, exovetter 0.0.16, ldtk 1.9.3. TRICERATOPS, psfmachine, tglc, MonoTools, SHERLOCK each get their own venv; any of them installed in the search venv changes wotan/numpy and therefore changes search results.

---

## 1. Loading and cleaning

### 1.1 lightkurve search/download (`s01_kbonus_load.py`, `s02_kbonus_download.py` `[reported]`; the API facts `[verified]`)

```python
import lightkurve as lk, astropy.units as u
from astropy.coordinates import SkyCoord
c = SkyCoord(285.907827, 40.941988, unit='deg')
sr = lk.search_lightcurve(c, author='KBONUS-BKG', radius=20*u.arcsec)   # ~26 s
# 2 products: target_name 'Gaia DR3 2103792461826460288' (0.0") and 'Gaia DR3 2103792427466721792' (3.1")
sr = sr[sr.target_name == 'Gaia DR3 2103792461826460288']
lcc = sr.download_all(download_dir=DATA)   # 8.44 MB: .../hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits
```

The file is named by KIC when the source is in the KIC, otherwise `gaia-dr3-<id>` (MAST: `hlsp_kbonus-bkg_kepler_kepler_<target>_kepler_v1.0_<prodtype>.<ext>`, [HLSP page](https://archive.stsci.edu/hlsp/kbonus-bkg) `[verified]`), while `target_name` is the Gaia designation. For thousands of stars do what Kuber's `batch2.py`/`acquire.py` do `[verified]`: one POST to `Mast.Caom.Filtered` with `filters=[{'paramName':'provenance_name','values':['KBONUS-BKG']},{'paramName':'target_name','values':[<Gaia designations>]}]`, `pagesize 50000`, then fetch `https://mast.stsci.edu/api/v0.1/Download/file?uri=<dataURL>`; for TESS, `provenance_name in ['SPOC','TESS-SPOC','QLP']`, `dataproduct_type='timeseries'`, drop SPOC rows with `t_exptime<100` (20-s products), rank SPOC < TESS-SPOC < QLP, one product per sector.

### 1.2 KBONUS-BKG file anatomy (`astropy.io.fits`, `[verified]`)

- 28 HDUs. `LIGHTCURVE_STITCHED` (50,182 rows): `TIME` (BKJD = BJD_TDB - 2454833), `CADENCENO`, `QUARTER`, `FLUX`, `FLUX_ERR` (e-/s), `SAP_FLUX`, `SAP_FLUX_ERR`, `PSF_FLAT_FLUX`, `PSF_FLAT_FLUX_ERR`, `SAP_QUALITY`, `FLATTEN_MASK`. MAST definitions: `Flux` = "PSF flux from stitched quarters", `PSF_flat_flux` = "PSF flux from stitched quarters after flattening", `SAP_quality` = "Quality flag from the TPF", `Flatten_mask` = "Quality flag from the flattening process" `[verified]`.
- `LIGHTCURVE_Q<n>`: `CADENCENO, TIME, FLUX` (corrected PSF flux), `FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR, CENTROID_COLUMN, CENTROID_ROW, SAP_QUALITY, SAP_BKG, RED_CHI2, PSF_FLUX_NOVA, PSF_FLUX_ERR_NOVA`; `APERTURE_Q<n>` masks. Per-quarter headers carry `PSFFRAC`, `CROWDSAP`, `FLFRCSAP`, `CHANNEL/MODULE/OUTPUT` (Q1: channel 6, module 3, output 2) but **no `EXPOSURE`, `TIMEDEL` or `BJDREFI`** — take those from Kepler documentation or the originating TPF. (The MAST page documents psffrac/crowdsap/flfrcsap only as catalogue columns, "Minimum detected PSF fraction / SAP crowding / SAP flux fraction"; the header keywords exist in the files regardless.)
- Primary header: `OBJECT 'KIC 5689121'`, `KEPMAG 16.332`, `GMAG 16.282`, `TPFORG 5689116`, `GAIAID 'Gaia DR3 2103792461826460288'`, `RA_OBJ 285.90784`, `DEC_OBJ 40.94198`, `TIMESYS TDB`, `DOI 10.17909/7jbr-w430`, `HLSPVER V1.0`; no TEFF/RADIUS/TICID.
- This star: quarters 1,2,3,4,6,7,8,10,11,12,14,15,16 (Q5, Q9, Q13 missing: the originating TPF's channel sat on failed Module 3 in those quarters), 46,523 good points after `SAP_QUALITY==0`, baseline 1,426 d, robust point-to-point scatter 1,033 ppm per 29.42-min cadence; per-quarter PSFFRAC 0.67-0.87, CROWDSAP 0.75-0.96, FLFRCSAP 0.69-0.95 (`s04_e2e.log`).
- Kuber's choices (`pilot.py` `[verified]`): `FLUX` (corrected PSF), `SAP_QUALITY == 0`, finite positive errors, `PSFFRAC >= 0.5`, >= 300 good points per quarter, per-quarter median normalisation, BKJD -> BTJD by `-2167` d, manifest JSON with `fits_sha256`. Kepler Bonus ([arXiv:2310.17733](https://arxiv.org/abs/2310.17733) `[verified]`) extracted 606,900 light curves with the Linearized Field Deblending (LFD) PSF method (psfmachine), 406,548 of them new background sources; data licence CC BY 4.0, cite the paper and the DOI.
- **Time systems:** BKJD = BJD_TDB - 2454833 (Kepler TPF header `BJDREFI 2454833`, `TIMESYS TDB`); BTJD = BJD_TDB - 2457000 (TESScut header `BJDREFI 2457000`); Kuber's T0 = BJD 2455630.245431 = BKJD 797.245431 = BTJD -1369.754569 `[verified]`.
- **Exposure vs cadence (needed for every model fit):** Kepler LC TPF header `[verified]`: `INT_TIME 6.0198 s`, `READTIME 0.5189 s`, `NUM_FRM 270`, `TIMEDEL 0.0204336 d` (1765.46 s) -> **integration 1625.35 s** inside a 1765.46-s cadence. TESScut Sector 14 header `[verified]`: `EXPOSURE 0.0165 d = 1425.6 s`, `TIMEDEL 1800 s`, `NUM_FRM 900 x 1.98 s` (cosmic-ray mitigation drops 20%); 10-min and 200-s FFIs scale the same way (475.2 s, 158.4 s). Always read `EXPOSURE`/`INT_TIME*NUM_FRM` from the product you fit.

Practical: lightkurve's generic reader hands you one table extension; for per-quarter extensions and header metrics use `astropy.io.fits` directly.

### 1.3 Quality bitmasks (values printed from the venv `[verified]`)

- Kepler/K2 `KeplerQualityFlags.DEFAULT_BITMASK = 1130799` = Attitude tweak (1), Safe mode (2), Coarse point (4), Earth point (8), Desaturation event (32), Manual exclude (256), Detector anomaly (16384), No data (65536), Thruster firing (1048576); `HARD = 1664431`, `HARDEST = 2096639`. **Not in default:** Zero crossing (16), Argabrightening (64), Cosmic ray in optimal aperture (128), Sudden sensitivity dropout (1024), Impulsive outlier (2048), Argabrightening on CCD (4096), Cosmic ray in collateral (8192), No fine point (32768), **Rolling band in optimal aperture (131072), Rolling band in full mask (262144)**, Possible thruster firing (524288). Docs: https://lightkurve.github.io/lightkurve/reference/api/lightkurve.KeplerQualityFlags.html.
- TESS `TessQualityFlags.DEFAULT_BITMASK = 17087` = Attitude tweak (1), Safe mode (2), Coarse point (4), Earth point (8), Argabrightening (16), Desaturation event (32, momentum dumps), Manual exclude (128), Impulsive outlier (512), Bad Calibration Exclude (16384); `HARD = 24319`, `HARDEST = 65535`. **Not in default:** Cosmic ray (64), Discontinuity corrected (256), Cosmic ray in collateral (1024), Straylight (2048), Straylight2 (4096), **Planet Search Exclude (8192)**, Insufficient Targets for Error Correction Exclude (32768). Docs: https://lightkurve.github.io/lightkurve/reference/api/lightkurve.TessQualityFlags.html; bit definitions in the TESS SDPDD [EXP-TESS-ARC-ICD-TM-0014-Rev-F](https://archive.stsci.edu/files/live/sites/mast/files/home/missions-and-data/active-missions/tess/_documents/EXP-TESS-ARC-ICD-TM-0014-Rev-F.pdf).
- Kuber's `pixel_diagnostic.py` uses `(QUALITY & (17087|8192)) == 0` on TESS TPFs `[verified]` and `QUALITY == 0` on Kepler TPFs. KBONUS `SAP_QUALITY==0` is stricter than lightkurve's default and loses a few percent of cadences; for rolling-band-prone Kepler channels add bits 131072/262144 to your mask.

### 1.4 Quarter / campaign / sector boundaries

- **Kepler** (measured on this star `[verified]`): quarter-boundary gaps 0.88-4.78 d (Q1->Q2 4.78, Q2->Q3 1.76, Q3->Q4 2.88, Q6->Q7 0.88, Q10->Q11 0.94, Q11->Q12 1.08, Q14->Q15 2.17, Q15->Q16 0.98), 15.8 d for Q7->Q8 (safe mode), 92-105 d where Module-3 quarters are missing; 57 intra-quarter gaps > 0.5 d (monthly downlinks, mostly 0.5-0.9 d). Each quarter sits on a different CCD channel (different aperture, PSFFRAC, CROWDSAP), so **normalise and detrend per quarter, never across a boundary**; Kuber turned the per-quarter flux-fraction change into his dilution test. Reference: Kepler Data Characteristics Handbook KSCI-19040-005 (https://archive.stsci.edu/kepler/manuals/Data_Characteristics.pdf, reachable `[verified]`). Suspect-period list for Kepler: ~93 d (quarter cycle) and harmonics, 372.5 d (orbit), periods commensurate with the 29.424-min cadence, and anything that disappears when rolling-band-flagged cadences are removed.
- **K2:** campaign files; 6-h thruster firings dominate; use `quality_bitmask='default'` (includes Thruster firing) and the EVEREST/K2SFF HLSPs or `lightkurve.correctors.SFFCorrector` (https://lightkurve.github.io/lightkurve/reference/api/lightkurve.correctors.SFFCorrector.html); K2 Handbook https://archive.stsci.edu/k2/manuals/k2_handbook.pdf.
- **TESS:** one file per sector, two orbits with a ~1-d downlink gap, momentum dumps (bit 32). eumemic found a QLP-specific ramp at orbit edges in Sectors 94-104 that detrending turns into a dip repeating at ~13.7 d (43% of new signals in a 9,112-star pilot landed at 12-14 d); their `load_star` drops the first 1.0 d and last 0.5 d of each orbit's QLP data in those sectors (`validation/qlp_orbit_edges.md`, read from the clone `[verified]`). Chopra flags periods in the 12-25 d pile-up near the 13.7-d TESS orbit and at 0.5-1 d as soft flags (`docs/GUIDE.md` `[verified]`).
- **Stitching:** eumemic prefers SPOC 2-min PDCSAP, then TESS-SPOC FFI PDCSAP, then QLP, one product per sector; Kuber's `acquire.py` does the same and stores one compressed `.npz` per star with `time, flux, flux_err, sector, exptime, cx, cy, bkg` `[verified]`. Chopra's `reject_bad_sectors` drops a sector if its 1-h point-to-point noise exceeds 1% (always) or, when the star has >= 3 sectors, 3x the star's median sector `[verified]`.

### 1.5 Outlier rejection and normalisation

- Asymmetric clipping (Chopra `lcsearch.clean` `[verified]`): after detrending, clip **upward** residuals at `FLARE_CLIP_SIGMA = 3.0` and **downward** only beyond `DROPOUT_CLIP_SIGMA = 10.0` ("downward points could be transits"). eumemic clips only upward at 4 sigma and first drops normalised flux outside (0, 5) `[verified]`. `lightkurve.LightCurve.remove_outliers(sigma_upper=3, sigma_lower=10)` does the same ([docs](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.LightCurve.remove_outliers.html)). Chopra's caveat: because the clipping is asymmetric, flipping the flux is not a perfect null.
- Normalisation: divide each quarter/sector by its own median before stitching (all three pipelines). Errors: do not trust pipeline errors for faint stars; use `1.4826*median(|diff(f)|)/sqrt(2)` per segment (Chopra) or eumemic's `rolling_std` (+-0.5 d robust std floored at half the global scatter).
- Our run: after `SAP_QUALITY==0`, biweight 1.0 d per quarter and the +3/-10 sigma clip, 45,648 of 46,523 points remain (`s04_e2e.log` `[reported]`).

---

## 2. Detrending

### 2.1 wotan 2.0

```python
from wotan import flatten, transit_mask
flat, trend = flatten(t, f, method='biweight', window_length=1.0, break_tolerance=0.5,
                      edge_cutoff=0.0, return_trend=True)
mask = transit_mask(time=t, period=P, duration=2*dur, T0=T0)   # then flatten(..., mask=mask)
```

Signatures in the venv `[verified]`: `flatten(time, flux, window_length=None, edge_cutoff=0, break_tolerance=None, cval=None, return_trend=False, method='biweight', kernel=None, kernel_size=None, kernel_period=None, proportiontocut=0.1, robust=False, max_splines=50, stdev_cut=2, return_nsplines=False, mask=None, verbose=False)`; `transit_mask(time, period, duration, T0)`. Methods: biweight, huber, hampel, andrewsinewave, welsch, median, medfilt, mean, trim_mean, winsorize, rspline, hspline, pspline (needs `pygam`, not installed), cofiam, cosine, savgol, lowess, supersmoother, ridge/lasso/elasticnet, gp (squared_exp, matern, periodic, periodic_auto) ([README](https://github.com/hippke/wotan)). wotan 2.0 BUG-024 fixed `mask` being silently ignored by pspline/cofiam/supersmoother/savgol/medfilt — in 1.x the transit-mask trick did nothing for those methods.

Rules from the wotan paper (Hippke et al. 2019, [arXiv:1906.00966](https://arxiv.org/html/1906.00966v3) `[verified]`): "the optimal width of the sliding window is about three times the transit duration" (Sec. VI; Kepler test Sec. IV.6); the window "can be slightly shorter (w/T14 >~ 2.2) to preserve almost all (>~98%) of the flux integral" (white-noise injections); discarding half a window at segment edges is best "only by a slight margin (dSDE ~0.1)" (Sec. IV.4); on the 100-planet K2 sample GPs "are worse than all robust methods", missing 10-11 planets vs 6 for biweight/Welsch. `break_tolerance` splits at gaps longer than the value so each segment is detrended alone; 0.5 d splits every Kepler quarter boundary and most monthly gaps (measured gap list above).

Measured speed (wotan 2.0, venv `[reported]`): biweight on 67,000 points 0.10 s (window 0.5 or 1.0 d), 0.16 s (2.0 d); median 0.10 s; rspline 0.34 s; hspline 2.95 s; cosine 3.5 s on 4,400 points; `gp` matern 55 s on 4,400 points (unusable in a batch search). First call pays JIT warm-up (4.4 s for the 13-quarter loop).

### 2.2 Over-detrending, shown numerically (`[reported]`)

Synthetic 90-d, 30-min series (1,200 ppm white noise + 3,000 ppm/4.3 d and 1,000 ppm/1.7 d sinusoids), 2,000-ppm 2.3-h transits at 13.77 d. Recovered in-transit depth: biweight 1.0 d -> 1,970 ppm; biweight 0.3 d -> 1,443 ppm (28% eaten); cosine 1.0 d -> 1,515; wotan `gp` matern 1.0 -> 1,096; biweight 0.3 d **with transits masked** -> 2,176; celerite2 Matern-3/2 fit with transits masked -> 2,232 (0.37 s on 4,409 points). The transit-mask re-detrend (detect -> `transit_mask` -> `flatten(mask=)` -> re-measure) recovers what the blind window removed; Chopra's `vet.masked_clean` (mask_factor 2.0) and eumemic's transit-masked re-detrending do exactly this.

### 2.3 Concrete parameter choices

Expected central duration for a 0.505 Msun / 0.507 Rsun M dwarf (rho = 3.87 rho_sun), `T14 ~ 13 h (P/365.25 d)^(1/3) rho^(-1/3)` (Chopra `duration_days`) `[verified]`: P = 0.5 d -> 0.92 h, 1 d -> 1.16 h, 5 d -> 1.98 h, 13.77 d -> 2.78 h (observed 2.1-2.4 h at b ~ 0.5), 30 d -> 3.6 h.

- **Kepler 30-min, Kp 16 M dwarf, P 0.5-30 d:** 3 x T14_max = 0.45 d is the paper's rule, but all three practitioners went longer because 30-min cadence gives only 4-5 points per transit: Kuber 1.0 d per quarter (`break_tolerance` 0.5), Chopra 0.75 d chosen by injection-recovery (8 stars x 2 planted 1.5 R_earth planets: 3/16, **10/16**, 8/16, 5/16 recovered for 0.5, 0.75, 1.0, 1.5 d; his comment: in-transit points pull the trend by ~duration/window, ~18% for a 2.9-h dip in 0.75 d) `[verified]`, eumemic `window = clip(3*T14(P_max), 0.5, 1.0)` d with Lomb-Scargle prewhitening (0.1-10 d on 30-min bins, 4 harmonics, chunked) only when `P_rot < 3.0 d and var_snr > 3.0` `[verified]`. Recommendation: 0.75-1.0 d biweight, `break_tolerance=0.5`, `edge_cutoff=0` (keep data; edges are handled by per-event local-baseline vetting), clip +3/-10 sigma, and confirm the window by injection-recovery on your own sample.
- **TESS 200-s / 10-min FFIs (K/M dwarfs, P <= 30 d):** T14 <= ~4 h -> 0.5 d window (eumemic's floor); fit the biweight to 10-min means and interpolate (eumemic: 14x faster trend, residuals differ by ~30 ppm, 30 vs 28 of 100 injections detected `[verified]`); bin to 10-20 min before the period search (Chopra `BIN_S = 600`; eumemic coarse 20 min, fine 10 min); clip flares upward at 3-4 sigma.
- GP alternative for active stars: celerite2 (`terms.SHOTerm(sigma=, rho=, Q=|tau=)`, `terms.RotationTerm(sigma=, period=, Q0=, dQ=, f=)`, `terms.Matern32Term(sigma=, rho=)`, `GaussianProcess(kernel, t=None, mean=0.0)`, `gp.compute(t, yerr=|diag=, quiet=True)`, `gp.log_likelihood(y)`, `gp.predict(y, t=, return_var=)`; O(N); [API](https://celerite2.readthedocs.io/en/latest/api/python/) `[verified]`). Use it for vetting/fitting with transits masked, not for the blind search. TLS 2 also offers `power(SDE_detrend="hybrid")` (default `"median"`, `SDE_MEDIAN_KERNEL_SIZE 30`, hybrid edge fraction 0.05, kernel factor 1.5 `[verified]`).

---

## 3. Period search

### 3.1 astropy `BoxLeastSquares` (`[verified]` API; runs `[reported]`)

API ([docs](https://docs.astropy.org/en/stable/api/astropy.timeseries.BoxLeastSquares.html)): `BoxLeastSquares(t, y, dy=None)`; `power(period, duration, objective=None, method=None, oversample=10)` with `objective='likelihood'` (default) or `'snr'` ("signal-to-noise of the transit depth"), `method='fast'` (Cython) or `'slow'`; `autoperiod(duration, minimum_period=None, maximum_period=None, minimum_n_transit=3, frequency_factor=1.0)` with `df = frequency_factor * min(duration) / (max(t)-min(t))**2`; `compute_stats(period, duration, transit_time)` returns `depth, depth_odd, depth_even, depth_half, depth_phased, harmonic_amplitude, harmonic_delta_log_likelihood (>0 favours the sinusoid), transit_times, per_transit_count, per_transit_log_likelihood`; `transit_mask(t, period, duration, transit_time)`. Durations must be strictly shorter than the smallest trial period; `dy` must be the empirical scatter because `'snr'` power scales with it (the median-normalised SDE does not).

**Do not call `autoperiod` on multi-year baselines** `[verified]`: 0.5-30 d with durations [0.05,0.1,0.2] d -> 83.84 M periods for a 1,460-d Kepler baseline (79.98 M for this star's 1,426 d), 256.8 M for 7 yr of TESS sectors, ~143k for two consecutive sectors. Use a stellar-density-aware grid (Chopra `period_grid`, eumemic `period_grid`): step frequency so the phase drift over the baseline stays below half the shortest plausible duration divided by an oversampling factor:

```python
def dur_days(P, rho): return 13/24*(P/365.25)**(1/3)*rho**(-1/3)
def grid(baseline, rho, pmin=0.5, pmax=30, oversample=3):
    f, fq = [], 1/pmin
    while fq > 1/pmax:
        f.append(fq); P = 1/fq
        fq -= 0.5*dur_days(P, rho)/(oversample*baseline*P)
    return np.sort(1/np.array(f))
```

Sizes `[verified]`: 512,101 periods (1,460 d, rho 3.9), 325,337 (rho 1.0), 876,886 (2,500 d), 18,940 (54 d), 500,176 for this star's 1,426 d. TLS's Ofir-style grid for the same star was 180,849 periods for 1-30 d (oversampling 3). Duration grid: Chopra `duration_grid(n=7)` geometric between 0.35 T14(P_min) and min(1.3 T14(P_max), 0.4 P_min) `[verified]`; eumemic uses 7 log-period chunks each with 6 durations 0.4-1.6 x T14 `[reported]`.

Runtime, real LC (45,648 points), 499,259 periods, 8 durations, `oversample=5`, 1 core `[reported]`: `'snr'` 186 s, `'likelihood'` 150 s; 6 durations + `oversample=3`: 152 s, peak RSS 115 MB. Cost ~ N_points x N_periods x N_durations x oversample. Chopra's `sde()` divides power by a running median in 60 log-period bins; eumemic's uses median/MAD in 200 bins `[verified]`.

Our BLS on KIC 5689121 (`objective='snr'`) `[reported]`: P = 13.76953 d, t0 = BKJD 136.3098, duration 2.45 h, box depth 1,649 +- 56 ppm, odd 1,634 +- 81 vs even 1,662 +- 78 ppm, half-period depth 801 +- 40 ppm, phase-0.5 "secondary" -71 +- 57 ppm, `harmonic_delta_log_likelihood = -430`, 74 transits with data. Note the box depth (1,649 ppm) sits well below the physical depth (TLS 1,999 ppm; batman (rp/rs)^2 1,959 ppm; Kuber 2,165-2,348 ppm PSF): box depths are for consistency tests, not for the headline number.

### 3.2 transitleastsquares 2.0.1 (`[verified]` API/constants; runs `[reported]`)

```python
from transitleastsquares import transitleastsquares, catalog_info
ab, M, Mmin, Mmax, R, Rmin, Rmax = catalog_info(TIC_ID=120686465)   # ((0.2125, 0.3984), 0.5054, 0.0212, 0.0212, 0.5075, 0.0163, 0.0163), 1.3 s
model = transitleastsquares(t, flux, dy)
r = model.power(period_min=1, period_max=30, R_star=0.507, M_star=0.505,
                R_star_min=0.3, R_star_max=0.8, M_star_min=0.3, M_star_max=0.8,
                u=list(ab), use_threads=1, oversampling_factor=3, duration_grid_step=1.1,
                show_progress_bar=False, verbose=False)
r.period, r.period_uncertainty, r.T0, r.SDE, r.snr, r.FAP, r.depth, r.duration, r.rp_rs,
r.odd_even_mismatch, r.transit_count, r.distinct_transit_count, r.transit_depths, r.snr_per_transit
```

Defaults from `tls_constants.py` `[verified]`: `OVERSAMPLING_FACTOR 3`, `DURATION_GRID_STEP 1.1`, `FRACTIONAL_TRANSIT_DURATION_MAX 0.12`, `N_TRANSITS_MIN 2`, `T0_FIT_MARGIN 0.01`, `TRANSIT_DEPTH_MIN 1e-5`, `R_STAR_MIN/MAX 0.13/3.5`, `M_STAR_MIN/MAX 0.1/1.0`, `SDE_DETREND 'median'` (kernel 30), `DEFAULT_BACKEND 'fused-pl'` (others `'fused-binned'`, `'fused'` exact, `'numba'` = TLS 1.33 algorithm; per call `power(backend=...)` or env `TLS_BACKEND`), `BIN_FRACTION 0.03`, `BIN_FACTOR_MAX 64` (internal pre-binning), `WORKER_IDLE_TIMEOUT 300 s`, `TLS_VERSION 'Transit Least Squares TLS 2.0.1 (11 Oct 2026)'`. `VALID_PARAMETERS`: `R_star, R_star_min, R_star_max, M_star, M_star_min, M_star_max, period_min, period_max, n_transits_min, per, rp, a, inc, b, ecc, w, u, limb_dark, duration_grid_step, transit_depth_min, oversampling_factor, SDE_detrend, T0_fit_margin, T0_search_margin, use_threads, show_progress_bar, transit_template ('default'|'grazing'|'box'), verbose, backend`.

**Threads and pools `[verified]`:** `use_threads` defaults to `multiprocessing.cpu_count()`; TLS 2 runs its search on a persistent pool of forked worker processes kept alive 300 s between `power()` calls (~15 ms to start with `fork`, ~3 s with `spawn`/`forkserver`, the Linux default from Python 3.14). Inside a `multiprocessing.Pool` worker, `use_threads=2` raises `AssertionError('daemonic processes are not allowed to have children')`; `use_threads=1` works there, and `use_threads>1` works inside `concurrent.futures.ProcessPoolExecutor` (non-daemonic workers). So: either one TLS per `ProcessPoolExecutor` worker with `use_threads=1`, or one TLS at a time with `use_threads=4`.

SDE thresholds (TLS [FAQ](https://transitleastsquares.readthedocs.io/en/latest/FAQ.html) `[verified]`): from >10,000 white-noise runs, SDE 5.7 (1-FAP 0.9), 6.1 (0.95), **7.0 (0.99)**, 8.3 (0.999), **9.1 (0.9999)**; "often required to be >9"; in red noise "high SDE values occur more often" so the table is optimistic. `results.FAP` is a table lookup `[verified]`: SDE 5 -> NaN, 5.7 -> 0.092, 6 -> 0.057, 7 -> 0.0094, 8 -> 0.0014, and **every SDE >= 9.1 returns the floor 8.0032e-05** — FAP carries no information for strong signals. TLS paper ([arXiv:1901.02015](https://arxiv.org/abs/1901.02015), A&A 623, A39 `[verified]`): recovery ~93% (TLS) vs ~76% (BLS) for Earth-sized planets around Sun-like stars at 1% false-positive rate, SDE~7 in both cases — per-light-curve numbers; Section 5 gives the survey-level correction.

`catalog_info` `[verified]`: `KIC_ID=5689121` raises `IndexError` (the KBONUS background star is absent from the Vizier KIC table TLS queries); `TIC_ID=120686465` works; `KIC_ID=5689116` also works but returns the **neighbour's** parameters (u = (0.4703, 0.1893), M 1.047, R 0.96) — never fall back to `TPFORG`. TLS ships Claret tables (`ld_claret_tess.csv`; `JAA546A14limb1-4.csv` for Kepler).

**Measured TLS runtimes** (`[reported]`; synthetic 4-yr Kepler-like LC, 67,535 points, 1,200 ppm, injected 2,000 ppm / 2.3 h at 13.769 d; 4 shared vCPUs under load):

| run | config | periods | seconds | result |
|---|---|---|---|---|
| A | 1-30 d, `use_threads=2` | 185,100 | 121.7 | P 13.76952, SDE 116, SNR 36.6, 106 transits |
| B | 1-10 d, 2 threads | 146,253 | 116.2 | **P/2 = 6.8847 d** (alias: true period outside the range) |
| C | binned to 1 h (33,121 pts), 1-30 d | 185,099 | 121.6 | SDE 106 — no speed-up (TLS 2 pre-bins internally) |
| D | 1 thread, `oversampling_factor=2`, `duration_grid_step=1.15` | 123,400 | 109.7 | SDE 117 |
| E | 2 threads, `backend='fused'` (exact) | 185,100 | 363.8 | SDE 116 |
| real LC | KIC 5689121, 45,648 pts, 1-30 d, 1 thread | 180,849 | 198 / 191 | below |
| real LC | same, `use_threads=4` | 180,849 | 100 | identical |

Lessons: restricting the period range saves little (most trial periods are short); cheap levers are `oversampling_factor` (linear) and threads; binning does not help; the exact backend costs 3x; peak RSS 225 MB. Always test P/2, 2P, 3P/2, 2P/3 after a detection (eumemic `known_alias` uses 1/3..3 harmonics).

Our TLS on the real cleaned LC `[reported]`: **P = 13.769643 +- 0.001299 d, T0 = BKJD 136.30486, SDE 105.3, SNR 27.5, depth 1,999 ppm (template bottom), duration 1.55 h, rp/rs 0.040, odd/even mismatch 0.21 sigma, 104 transits (74 distinct with data), FAP 8.0e-5 (floor)**; the next four peaks are aliases within 0.001 d. This reproduces Kuber's 13.769448 d from the public product in ~3 min on one core.

### 3.3 Semi-coherent and stack-slide BLS (eumemic `fastbls.py`, used by Kuber) `[verified]` structure

`pipeline/fastbls.py` is 242 lines, MIT, with numba kernels `_season_scan`, `_coherent_best`, `_stackslide` and wrappers `seasons(t, gap=20, max_len=None)`, `semicoherent`, `coherent`, `stackslide`. (a) **semi-coherent**: per-season BLS dchi2 on a grid set by the longest season, summed over seasons per duration; `search.semicoherent_peak(t, y, dy, rho, pmin, pmax, k_peaks=40)` takes the top 40 local maxima (plus top 20 by raw statistic) and refines coherently; (b) **stack-slide**: each season folded once per coarse frequency, then fine frequencies shift the folded sums by the phase accumulated since the season midpoint. Kuber's `pilot.py` imports these: `short_blocks = seasons(gap=20 d, max_len=90 d)` to cut Kepler into quarter-length blocks, `semicoherent_peak(tt, yy, ee, 5., .5, 25., k_peaks=40)`, then `refine(tt, yy, ee, p, 5., np.ptp(tt))`. eumemic `refine` `[verified]`: step = 0.5 T14/P/span, 241 frequencies in +-3 coarse steps, 12 durations between max(0.3 T14, 0.4 h) and min(2 T14, 0.3 P, 10 h), astropy BLS `objective='likelihood'`, `oversample=10`, SNR = sqrt(2 x power). eumemic injection test at nominal SNR 10 on 300 TESS stars (README `[verified]`): coherent astropy BLS 60% recovered (31% at SNR >= 9), 17.5 s/star; stack-slide 62% (29%), 8.9 s; semi-coherent 52% (26%), 4.2 s; "6 were found only by the coherent search and 2 only by stack-slide; 18 only by coherent and 4 only by semi-coherent". Porting `fastbls.py` (numba 0.68.0 is in the venv) is the single most valuable piece of engineering for a 4-core Kepler search.

### 3.4 GPU, correlated-noise and matched-filter alternatives

- **cuvarbase** (`eebls_gpu`, `eebls_transit_gpu` with `qmin_fac/qmax_fac/rho`; Keplerian assumption cuts a 10-yr grid from 5.35e8 to 2.7e5 frequencies; docs tested CUDA 8.0 and Python 2.7/3.4-3.6; needs PyCUDA + scikit-cuda) `[verified]` — unmaintained (last PyPI 2023), irrelevant without a GPU ([bls](https://johnh2o2.github.io/cuvarbase/bls.html), [install](https://johnh2o2.github.io/cuvarbase/install.html)).
- **nuance** (Garcia et al., [arXiv:2402.06835](https://arxiv.org/abs/2402.06835) `[verified]`): JAX/tinygp; `linear_search(time, flux, gp=None, X=None, model=None, positive=True, progress=True, backend=None|'cpu'|'gpu', batch_size=None)` gives single-event log-likelihoods under a GP; `periodic_search` combines them ([docs](https://nuance.readthedocs.io/en/latest/markdown/ls.html)). Best method in 93% of 438 rapidly rotating M-dwarf TESS cases vs five detrend+BLS pipelines; TOI-540 (sectors 4,5,6,31,32) took 5 min 35 s on 12 Apple M2 Max cores vs 1 min 49 s for 15 biweight windows + BLS; authors "expect nuance to perform between one and two orders of magnitude slower" and "do not recommend using nuance in the general case". Not installed (JAX stack).
- **SPOC/TPS wavelet matched filter** (Jenkins 2002, ApJ 575, 493; confirmed via the Kepler pipeline overview [arXiv:1001.0258](https://arxiv.org/pdf/1001.0258) which also gives the 7.1-sigma MES threshold): no maintained public Python implementation; the SOC pipeline is MATLAB/Java ([nasa/kepler-pipeline](https://github.com/nasa/kepler-pipeline)). Practical substitute: LEO-Vetter's SES/MES series plus Robovetter-style metrics (Section 5.5).
- **ML triage** `[verified]`: ExoMiner's `exominer_pipeline/run_pipeline.py` takes `--tic_ids` (format `{tic}_{sector_run},...`) or `--tic_ids_fp`, `--data_collection_mode 2min|FFI`, `--num_processes`, `--num_jobs`, `--dv_xml_data_repository`, `--lc_data_repository`, `--stellar_parameters_source ticv8|tess-spoc|<csv>`, `--ruwe_source gaiadr2|unavailable|<csv>`, `--task phot-vetting|planet-validation`, `--exominer_model`, `--max_model_workers`; it scores **TESS SPOC TCEs** from SPOC light curves + DV XML, NOSA licence; the README does not name the container image (get it "using the manifest" in the getting-started docs) and states no hardware requirement. Astronet (archived 2023-08-10) ships no weights; Astronet-Triage's README trains your own model, Releases is empty, data on Zenodo 7411579, Tey et al. 2023 report 99.6% recall at 75.7% precision ([arXiv:2301.01371](https://arxiv.org/abs/2301.01371), DOI 10.3847/1538-3881/acad85); qusi is a framework. Conclusion unchanged: no off-the-shelf CPU model scores arbitrary Kepler-background or home-made FFI signals; rule-based LEO-Vetter plus an inverted-trained classifier (eumemic) is the realistic route.

---

## 4. Multi-planet iteration, single transits, TTVs

- **Mask-and-repeat**: mask |phase| < 1-2 x duration (`BoxLeastSquares.transit_mask` or `wotan.transit_mask`; Chopra factor 2.0, eumemic/Kuber 1.0 x duration on the detrended series) and re-run up to 3 signals (Chopra `MAX_SIGNALS=3` `[verified]`; eumemic 3, 4 for TOI hosts; Kuber 2). Stop rules: SDE < 7 or SNR < 6 (Chopra), BLS SNR < 7 (eumemic), SNR < 9 (Kuber). eumemic's `repeat` test (same period within 0.5%, or half/double within 0.2%, re-found after masking) marks leftover variability; `known_alias` catches 1/3..3 harmonics of a masked TOI; TOI masks widen by 3 sigma of `sqrt(eT0^2 + (n eP)^2)`. TLS: `transit_mask` -> `cleaned_array` -> new `transitleastsquares(...)`.
- **Single-transit search**: nuance `linear_search`; LEO-Vetter's SES series; MonoTools 1.0.0 (PyMC/exoplanet + celerite2; pins lightkurve 1.11.0 -> own venv); SHERLOCK's `MoriartySearch` ([docs](https://sherlock-ph.readthedocs.io/en/latest/_autosummary/sherlockpipe.single_transits.search.MoriartySearch.html)). Simplest robust version is Kuber's `independent_vet.local_events` `[verified]`: per event, fit intercept + slope + exposure-integrated box in a window `min(max(3*duration, 0.15 d), 0.23*P)`, requiring >= 12 points in the window, >= 2 in transit and >= 4 beyond +-duration on each side, split by sector/quarter; noise floored by the local MAD; a single-event scan is the same fit slid along the time axis.
- **TTV tolerance**: QATS (Carter & Agol) allows inter-transit spacing in [DeltaMin, DeltaMax]; C++ only (`make`; `./call_qats <file> <DeltaMin> <DeltaMax> <q>`, e.g. `./call_qats kep36.dat 794 794 14` -> S_Best 82.6976, M_Best 55; `./call_qats_indices <file> <M> <DeltaMin> <DeltaMax> <q>`; single-column data, detrended, median-subtracted, sigma-normalised, gaps zero-filled; cite Carter & Agol 2012 and Kel'manov & Jeon 2004) ([readme](https://faculty.washington.edu/agol/QATS/readme) `[verified]`). Cheaper: the semi-coherent BLS tolerates within-season phase drift by construction, and the coherent refine (241 sub-steps, `oversample=10`) absorbs TTVs of order the cadence.

---

## 5. Statistics of large searches

### 5.1 The look-elsewhere effect on this star (`s07_null_throughput.py` `[reported]`)

Identical cleaning + astropy BLS (`objective='snr'`, 499k periods, 6 durations, `oversample=3`) on the real LC, its inversion (f -> 2 - f) and seven quarter-scrambled copies: **real** SNR 30.3 (global SDE 29.3) at 13.76953 d, 74 transits; **inverted** SNR 11.0 (SDE 6.1) at 29.04 d; **scrambles** SNR 9.6-11.9 (SDE 5.1-6.8) at 19.7-29.0 d, all but one at the longest trial duration (4.56 h) with 35-54 "transits". On a 45k-point, 1,000-ppm Kepler light curve searched over half a million periods, a BLS depth-SNR of 10-12 is what noise produces, piling up in the long-period / long-duration corner — the 12-25 d pile-up Chopra reported and the reason Kuber's "BLS SNR >= 9" first screen passed inverted stars as often as real ones (31 signals/29 stars real vs 32/28 inverted over 1,736 stars; stage 1: 10/9 vs 10/9 for ranks 65-576 over 497 stars; stage 2: 19/18 vs 22/19 for ranks 577-1800 over 1,175 stars; [kbonus-extension README](https://github.com/Kuberwastaken/planethunt/blob/main/analysis/kbonus-extension-20261011/README.md) `[verified]` from the clone). Depth/depth_err is a *local* statistic; only an empirical null converts it into a survey false-alarm rate. Do not mix statistics: BLS SNR (depth/err), TLS SNR (depth/std x sqrt(N_in)), and MES (matched-filter, Kepler threshold 7.1) have different null distributions.

### 5.2 Calibrating thresholds with inverted and scrambled light curves (`[verified]` from the clones)

- **Chopra** (526 stars, ~1,500 searches, `docs/GUIDE.md`, `06_inverted.py`): flip raw flux *before* cleaning. Real vs flipped: 743 vs 752 signals pass the search (SDE >= 7, SNR >= 6, >= 3 transits); 31 vs 39 have no hard flags; 10 vs 17 are flag-free; strongest flag-free flipped SDE 8.4 vs candidates 11.3 and 10.3. A second null scrambled each candidate's light curve 38 times (mask transits, slide each sector) and reached SDE 8.5 and 9.3. Review **withdrew** the extrapolated false-alarm counts (0.01 and 0.2 per survey): "38 samples can't support a tail extrapolation"; zero fakes in 526 flipped stars only caps the rate at ~3 per survey (95%).
- **Kuber**: matched inverted run on every star; at the candidate's strength (P > 1.5 d, BLS SNR >= 26.7, holdout SNR >= 15.5) two inverted stars (KIC 8701308, KIC 5283912) still pass — strength alone did not single out KIC 5689121 `[reported from his audit]`.
- **eumemic** (`validation/false_alarm_classifier.md`, README): inverted light curves of **798,769** random full-sky stars (two samples) through the identical pipeline (`EXO_INVERT=1` negates the detrended residuals); expected false alarms = inverted count x (4,039,203 / 798,769). Clean-pixel SNR >= 9: 726 real vs 96 expected false alarms; SNR 7.5-9: 313.5 expected (62 x 4,039,203/798,769). 46% of clean SNR >= 9 signals with MES < 10 are expected to be noise, 16% at MES 10-12, 2% above -> MES >= 10 for tier N; 78 candidates moved to tier C. `HistGradientBoostingClassifier` (BLS SNR, MES, SDE, duration ratio, n good events, max single-event fraction, odd/even, secondary, per-sector depth consistency, LEO CHI/DMM/Fred/SHP/chases, Tmag, contamination, radius, Teff; centroids excluded because inversion does not invert centroids), inverted rows weighted by the star ratio; cut 0.80 -> 69 signals with 5.85 expected false alarms (~92% purity), TIC-grouped folds after a contributor found star leakage. Earlier 9,979 K/M-dwarf search: 1,500 inverted stars gave 157 pre-filter passes, 1 passing all LEO tests, 0 with consistent per-sector depths.
- Caveats all three state: inversion measures noise and systematics, not astrophysical false positives; asymmetric clipping and rules tuned after seeing real signals make the null imperfect; report "in one inverted survey through the identical pipeline no signal reached X", not a p-value.

### 5.3 Injection-recovery for completeness

Chopra `inject_test.py` and eumemic `inject_test.py` inject a trapezoid (15% ingress) into the **raw** flux at `expected_duration(P, rho) x U(0.6, 1.0)`, depth set for a nominal SNR, then run the untouched pipeline; a hit = period within 0.3% of P, P/2 or 2P `[reported]`. eumemic: nominal-SNR-10 injections recovered at SNR >= 9 only 22% of the time, 48% at >= 7.5 ("red noise in full-frame light curves lowers the measured SNR") `[verified]`. Use the same machinery to pick the detrending window (Chopra's 0.75 d) and to calibrate vetting flags: two of Chopra's red-flag tests fired on 15-25% of planted planets until the noise model was fixed, dropping to 5% `[verified]`.

### 5.4 Empirical noise with off-transit boxes (Chopra) and Kuber's phase controls (`[verified]`)

Chopra `vet.depth_emp`: box depth at the fixed ephemeris, then at 40 off-transit phases (`NULL_OFFSETS` = 20 in 0.1-0.4 plus 20 in 0.6-0.9); error = 1.4826 x MAD of the 40 depths. The Codex audit found controls detrended with only the real transits masked, so each control "helped fit its own trend" — `depth_emp_matched` re-detrends per control phase with that phase masked; significance fell from ~8.2 sigma to 7.2 and 7.7 sigma. Kuber's `independent_vet.vet(controls=24|64)` fits per-event local-baseline boxes at control phases excluding |phase| and |phase-0.5| < max(0.06, 2 x duration/P) and the real transit cadences, reporting `depth_over_control_scatter` and an add-one rank `(1 + #controls >= depth)/(1 + N)` with "rank is not a survey false-alarm probability"; both report `max_positive_depth_contribution` / `max_ses_frac` (< 0.5) and leave-one-event-out SNR.

### 5.5 Model comparison: transit vs sinusoid vs EB

- BLS `compute_stats['harmonic_delta_log_likelihood']` (> 0 favours a sinusoid; -430 here) and `depth_odd/depth_even/depth_phased/depth_half` are the first-pass EB tests.
- **LEO-Vetter 1.2.0** `[verified]`: `TCELightCurve(tic, time, raw, flux, flux_err, per, epo, dur_days, planetno=1)`; `compute_flux_metrics(star)` needs `star = dict(rad, e_rad, mass, e_mass, Teff, e_Teff, rho, u1, u2, id, tic, ra, dec, raDegrees, decDegrees, planetData=None)` (KeyError `'e_rad'` is the first failure if omitted); `check_thresholds(metrics, "FA"|"FP")`. Defaults (`thresholds.py`): MES 6.2, N_transit 3, SHP 0.6, MS1 0.2, MS2 0.8, MS3 0.8, MS4 0, MS5 -1, MS6 -1, chases 0.78 (only if N_transit <= 5), DMM 1.5, max_SES_to_MES 0.8 (only if N_transit <= 10), AIC1 -60, AIC2 -30, SWEET 15, ASYM 8, CHI 7.8, frac_gap 0.5, V_shape 1.5, size 22, offset 15". KIC 5689121 (re-run, 4.7 s, 121 metrics): MES 24.21 (new_MES 22.98 after bad-event removal; 74 -> 66 transits), SHP 0.17, CHI 32.3, DMM 0.96, Fred 1.27, sig_pri 24.4 vs sig_sec 3.0 / sig_ter 3.0 / sig_pos 4.8, max_SES 5.6, SWEET 0.54, transit AIC 187 vs line AIC 752, Rp 2.80 +- 0.13 R_earth, Teq 381 K — **no FA, no FP flag**. (LEO's quick fit gives b 0.90, a/Rs 20.4: a reminder it is not a substitute for the MCMC fit.) `pixel_vetting(tlc, star, sectors, quality_flags=None, tdi_dir='.', n_bad=0, max_sectors=inf)` is TESS-only and needs `transitDiffImage` (not installed: `import leo_vetter.pixel` fails). eumemic replaced LEO's per-cadence loop with a compiled equivalent (2-5x) `[reported]`.
- **exovetter 0.0.16** `[verified]`: `ModShift, Sweet, Lpp, OddEven, Centroid, TransitPhaseCoverage, VizTransits, LeoTransitEvents`; `Tce(period=P*u.day, epoch=..., epoch_offset=exo_const.bkjd|btjd|bjd|mbjd, depth=d*exo_const.ppm, duration=h*u.hr)`.
- Bayesian evidence via nested sampling: juliet (`dataset.fit(sampler='dynesty'|'multinest'|'dynamic_dynesty'|'ultranest'|'slicesampler_ultranest'|'emcee'|'zeus', n_live_points=...)`; priors `P_p1, t0_p1, p_p1` or `r1_p1/r2_p1`, `b_p1` or `a_p1`, `rho` (kg/m^3), `ecc_p1/omega_p1` or `(s)esinomega/(s)ecosomega`, `q1_<inst>/q2_<inst>`, `mdilution_<inst>`, `mflux_<inst>`, `sigma_w_<inst>`, GP `GP_sigma/GP_rho/GP_timescale/GP_B/GP_C/GP_L/GP_Prot/GP_S0/GP_omega0/GP_Q`; results in `results.posteriors['posterior_samples'|'lnZ'|'lnZerr']` — [priors page](https://juliet.readthedocs.io/en/latest/user/priorsnparameters.html) + installed `fit.py` `[verified]`) or dynesty/ultranest directly (both installed). Budget minutes per candidate, not per star.

---

## 6. Transit fitting

### 6.1 batman + emcee with exposure integration (`s06_fit_vet.py` `[reported]`; exposure facts `[verified]`)

```python
import batman, emcee
params = batman.TransitParams(); params.t0=T0; params.per=P; params.rp=0.046; params.a=39.7
params.inc=89.3; params.ecc=0; params.w=90; params.limb_dark='quadratic'; params.u=[0.2125, 0.3984]
m = batman.TransitModel(params, t, supersample_factor=7, exp_time=1625.35/86400)  # Kepler LC INTEGRATION (270 x 6.0198 s), not the 1765.46-s cadence
# free: t0, per, rp/rs, a/rs, b (inc = arccos(b/aRs)), f0, ln jitter; 32 walkers x 2500 steps
```

`TransitModel(params, t, max_err=1.0, nthreads=1, fac=None, transittype='primary', supersample_factor=1, exp_time=0.0)` `[verified]`. Kipping 2010 ("Binning is sinning", MNRAS 408, 1758, [arXiv:1004.3741](https://arxiv.org/abs/1004.3741) `[verified]`): finite integration distorts the shape and systematically underestimates the light-curve stellar density unless integrated numerically. The previous agent's fit used `exp_time=1765.464 s` (the cadence, 8% too long) and should be re-run; its result on 2,637 points within +-0.4 d of transit: **t0 = BKJD 797.24544 (+0.0012/-0.0011) = BJD 2455630.2454 (Kuber: 2455630.245431), P = 13.769444 +- 0.000038 d (Kuber: 13.769448 +- 0.000065), rp/rs = 0.0443 (+0.0017/-0.0010), a/rs = 43.7 (+2.9/-7.5), b = 0.35 (+0.28/-0.24), rho_star = 5.9 (+1.3/-2.6) rho_sun vs 3.87 expected (Kuber: 4.4 +2.3/-2.8), Rp = 2.45 (+0.09/-0.06) R_earth for R* = 0.507 Rsun without dilution, (rp/rs)^2 = 1,959 (+152/-91) ppm**, acceptance 0.38, 109 s. Add a `rho_star` prior from TIC (`catalog_info(TIC_ID=...)` or `astroquery.mast.Catalogs.query_object('TIC 120686465', catalog='TIC')`) and a dilution term when the aperture contains the brighter neighbour; record the emcee seed.

### 6.2 Limb darkening

- Fast path: `catalog_info(TIC_ID=...)` -> Claret quadratic coefficients; for TIC 120686465 (Teff 3633 K) u = (0.2125, 0.3984) `[verified]`. KIC path fails for this background star; the neighbour's KIC returns (0.4703, 0.1893).
- LDTk 1.9.3: `from ldtk import LDPSetCreator, BoxcarFilter; sc = LDPSetCreator(teff=(3633,160), logg=(4.73,0.1), z=(0,0.1), filters=[...]); ps = sc.create_profiles(); cq, eq = ps.coeffs_qd(do_mc=True)`; downloads PHOENIX specific intensities from the Husser et al. server on first use (hundreds of MB) ([README](https://github.com/hpparvi/ldtk)). Installed; not run (download) `[unverified]`.

### 6.3 Other fitters (status `[verified]` metadata)

- **exoplanet 0.6.0**: `pip install "exoplanet[pymc]"`; docs say "updated to work with PyMC version 5"; current PyMC is 6.3.2 and exoplanet-core 0.4.0 (2026-07-25, py>=3.12) — PyMC 6 compatibility untested.
- **juliet 2.2.10**: installed (dynesty/ultranest/multinest/emcee/zeus backends; celerite GP kernels); heavy dependency tree.
- **PyTransit 2.9.3**: active (2026-10-15, py>=3.10) but TRICERATOPS pins `pytransit==2.2` -> keep TRICERATOPS in its own venv.
- **allesfitter 1.2.10**: last release 2022-07-25; only for EB/RV joint models.
- **MonoTools 1.0.0**: pins `lightkurve==1.11.0`, `scipy==1.12` -> own venv.

---

## 7. Pixel-level tools (`s05_pixels.py` `[reported]`; API `[verified]`)

```python
tpf = lk.search_targetpixelfile('KIC 5689116', author='Kepler', cadence='long', quarter=10).download(quality_bitmask='default')
# KeplerTargetPixelFile 4x5 px, 4,447 good cadences (623 flagged), channel 34, 1.70 MB; 10.5 s today
xy_t = tpf.wcs.world_to_pixel(SkyCoord(285.907827, 40.941988, unit='deg'))   # (1.15, 0.08)
xy_n = tpf.wcs.world_to_pixel(SkyCoord(285.905233, 40.942079, unit='deg'))   # (2.26, 1.47); 7.06" = 1.8 px
tpf.plot_pixels(aperture_mask=tpf.pipeline_mask, show_flux=True)
col, row = tpf.estimate_centroids(aperture_mask='all', method='moments')
from lightkurve.prf import KeplerPRF
prf = KeplerPRF(channel=tpf.channel, shape=tpf.shape[1:], column=tpf.column, row=tpf.row)
model = prf(center_col=tpf.column+xy_t[0], center_row=tpf.row+xy_t[1], flux=1, scale_col=1, scale_row=1, rotation_angle=0)
tc = lk.search_tesscut('TIC 120686465', sector=14).download(cutout_size=11)   # 1,237 x 11 x 11, 30-min, 2.98 MB
```

Docs: [plot_pixels](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.KeplerTargetPixelFile.plot_pixels.html), [estimate_centroids](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.KeplerTargetPixelFile.estimate_centroids.html), [search_tesscut](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.search_tesscut.html). `lightkurve.prf.KeplerPRF` works; the `tpfmodel` submodule needs `oktopus` (warning on import). astroquery `Tesscut.get_sectors / get_cutouts / download_cutouts(coordinates|object_name, size (default 5 px), sector)` recommends <= 30 px, has a 600-s timeout, and TICA cutouts were discontinued August 2025 (SPOC FFIs only) ([docs](https://astroquery.readthedocs.io/en/latest/mast/mast_cut.html) `[verified]`).

A naive single-quarter difference image (33 in-transit vs 128 near-transit out-of-transit cadences of the 7 Q10 transits) gives a per-pixel max SNR of **0.4** (target pixel 0.2, neighbour 0.1) `[reported]` — one quarter of a 1,000-ppm, Kp 16 star cannot localise a 0.2% dip. Kuber's `pixel_localize.py` `[verified]`: per event fit intercept + slope + exposure-integrated box **per pixel** with local-MAD noise floors, inverse-variance combine events, stack all quarters, fit `lightkurve.prf.KeplerPRF` templates (`LocalPRF(KeplerPRF)`) at the two Gaia positions through `WCS(h[2].header)` with a common offset calibrated on the quiescent image; compare chi2 (his sum dchi2 = 934 in favour of the target; his own caveat: "Pixel covariances and PRF-model mismatch not calibrated; chi-square differences are not a planet false-positive probability. Gaia2016 coordinates not proper-motion propagated"). Chopra's TESS version (`10_pixels.py`): 15x15 TESScut per sector, (before+after)/2 - in-transit per transit, per-pixel noise from 40 off-transit phases, joint all-sector Gaussian-spot fit with per-sector amplitude, dchi2 against each Gaia neighbour — rewritten after the audit found the first version "could locate a brightening as easily as a dimming" `[verified]`. eumemic's per-sector centroid F-test (transit-diffImage PRF fits, scatter floored at 1" per axis) rejects the true host at p < 0.001 for 0 of 565 known planets and for 127 of 175 (73%) neighbour-hosted signals; weak (38%) for sources within 15" (`validation/centroid_sectors.md` `[verified]`).

Tools: **transit-diffImage** (`tessDiffImage`: `make_ffi_difference_image`, `make_planet_difference_image`, `get_tess_cut`, `find_transits`, `make_difference_image`, `make_stellar_scene`, `draw_difference_image`; `keplerDiffImage.py` needs Kepler-RaDex2Pix; clone + `pip install .`; requirements `tess-point>=0.8.0, astropy, astroquery, scipy, scanf` `[reported]`); **LEO-Vetter `pixel_vetting`** wraps it (TESS only); **psfmachine** (`TPFMachine.from_TPFs(tpfs, n_r_knots=10, n_phi_knots=12).fit_lightcurves()`, Gaia-seeded PRF photometry, the Kepler Bonus engine) needs a Python <= 3.11 venv. For 200-s FFI sectors an 11x11 cutout is ~9x the 3 MB measured at 30 min.

---

## 8. Orchestration at scale on 4 CPUs

### 8.1 Pipeline design (what all three working pipelines converge on)

1. **Targets**: catalogue query frozen before any photometry is seen, CSV with a rank and SHA-256 (Kuber by Gaia G; Chopra by expected SNR of a 1.5 R_earth planet).
2. **Download -> cache**: one MAST POST for all product URLs; workers fetch (processes, not threads: eumemic saw 3x because FITS parsing is GIL-bound `[verified]`); verify `SIMPLE` and SHA-256; one compressed `.npz` per star (`time, flux, flux_err, sector/quarter, exptime, cx, cy`) plus a manifest JSON. Sizes: KBONUS LC FITS 8.44 MB, cleaned npz 1.39 MB, Kepler TPF 1.70 MB/quarter, TESScut 11x11x1,237 2.98 MB `[reported]`.
3. **Detrend + clean**: per segment biweight (0.5-1.0 d), asymmetric clip, empirical errors (< 0.2 s per star).
4. **Search**: coarse semi-coherent/stack-slide BLS on 20-min bins -> coherent refine of top peaks on 10-min bins -> `event_stats` (per-event depths, odd/even, secondary, max single-event fraction) -> up to 3 signals by masking; TLS only on survivors; test P/2, 2P harmonics.
5. **Triage**: hard flags (odd/even, secondary, repeat period, duration vs expected, single-event dominance, gap edges, per-season depth consistency, broad depression, same period in >= 3 other stars, too deep) and soft flags (12-25 d and 0.5-1 d pile-ups; for Kepler also ~93 d/372.5 d harmonics and cadence-commensurate periods); LEO-Vetter FA/FP; MES >= 10 (eumemic).
6. **Vet**: holdout seasons reserved **before** the search and written to `<id>_split.json` (Kuber: last max(2, ceil(0.25 N)) quarters; Chopra: hide-a-season after the fact, labelled retrospective), fixed-ephemeris local-baseline depths, phase controls, SAP/PSF/second-reduction consistency, pixel localisation, catalogue/literature matching, inverted null of the whole sample.
7. **Provenance**: per-star JSON keyed by `input_sha256` and `pipeline_version` (git SHA) so reruns skip unchanged work (Kuber `run_search.py`), `skip if <out>/<tic>.json exists` (eumemic, Chopra), `uv pip freeze` with versions (TLS 2 changed results in 10 fixes; numba 0.68.0), recorded seeds, and for anything you publish a DAG that re-runs downstream stages when code or inputs change (eumemic `pipeline/dag/`, `dag.lock.json`).

Parallelism: `ProcessPoolExecutor` (not `multiprocessing.Pool` if TLS runs inside workers with `use_threads>1`), `imap_unordered(chunksize=1)` or `map`, one process per core, `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` set before importing numpy, `use_threads=1` for TLS inside workers, `NUMBA_CACHE_DIR` on a writable path shared by workers, progress as JSON lines.

### 8.2 Measured per-star costs (`[reported]`; 4 shared vCPUs, load average 5-13 from other agents, so treat as upper bounds)

| stage | Kepler 4-yr LC (46-67k pts) | TESS 5 sectors (19k pts, binned to 20 min = 9.4k) |
|---|---|---|
| download | 8.4 MB FITS, seconds (MAST-bound) | 5 x ~2 MB FITS |
| wotan + clean | < 5 s incl. JIT | < 1 s |
| coherent astropy BLS, density grid, 6 durations, oversample 3 | 443 s/star in a 4-worker test (152 s solo) -> **0.5 stars/min on 4 workers** | 50 s/star -> **4.7 stars/min on 4 workers** (265k periods) |
| semi-coherent numba BLS (eumemic/Kuber code, not re-run here) | eumemic: 4.2 s median per TESS star (12-16 vCPU boxes); Kuber: 2-5 workers | 4.2 s (semi), 8.9 s (stack-slide), 17.5 s (coherent) |
| TLS 1-30 d (survivors only) | 100 s at 4 threads / 190 s at 1 thread | ~20-40 s |
| batman+emcee fit (32 x 2,500) | 109 s | similar |
| LEO-Vetter flux metrics | 4.7-8 s `[verified]` | 3-8 s |
| RAM | BLS worker 115 MB; TLS 225 MB; 4 workers < 2 GB | smaller |

### 8.3 Run-time estimates (search stage; add ~10% for vetting the few-percent survivors)

- **Fast tier** (semi-coherent numba BLS at ~10 s per Kepler star or ~5 s per TESS star per core): 1k stars ~45 min / ~20 min; 10k ~7 h / ~3.5 h; 100k ~3 days / ~1.5 days on 4 cores -> overnight is fine to ~10k; 100k wants a cloud box.
- **Full coherent astropy-BLS tier** (~150-450 s per Kepler star): 1k 10-30 h; 10k 4-13 days; 100k ~8,000-12,000 core-hours -> cloud. TESS 5-sector stars at 50 s: 1k 3.5 h, 10k 35 h, 100k ~15 days on 4 cores.
- **TLS on every star** (~190 s/core for Kepler): 1k 13 h, 10k 5.5 days — reserve TLS for BLS survivors.
- Kuber's 1,800-star KBONUS run used the semi-coherent code with 2-5 workers over two days of wall time (his READMEs `[verified]`); eumemic reports 200 TESS stars in 70 s on a 16-core AWS c8g.4xlarge after optimisations (trend on 10-min means 14x faster for the trend itself / 1.16x overall, division-free stack-slide kernel 1.42x, Graviton4 2.6x; batman built from source on Graviton) and downloads 20x faster inside AWS us-east-1 `[verified]`.
- Disk: 10k Kepler stars = 84 GB of KBONUS FITS or 14 GB of npz (delete FITS after parsing; keep URL + SHA-256); TPFs only for candidates (13 quarters x 1.7 MB = 22 MB per star).

Run overnight locally: target selection, downloads (4 processes), detrending, fast-tier search of up to ~10k stars, the inverted null of the same stars (doubles the cost), triage, TLS + fits for survivors. Push to a cloud box (or Modal): coherent full-grid BLS on > 1k Kepler stars, TLS on > 1k stars, 100k-star tiers, TRICERATOPS/MOLUSC (eumemic used AWS spot instances with a separate pinned venv).

---

## 9. Minimal reproducible end-to-end (`s04_kbonus_e2e.py`, 1 core, ~7 min `[reported]`)

load KBONUS per-quarter FLUX with `SAP_QUALITY==0` -> normalise per quarter -> wotan biweight 1.0 d (`break_tolerance=0.5`) -> clip +3/-10 sigma -> density-aware grid (499k periods, 0.5-30 d) -> astropy BLS `objective='snr'` (186 s) -> `compute_stats` -> TLS 1-30 d with TIC limb darkening (198 s) -> P = 13.76964 d, SDE 105, 74 transits, odd/even 0.2 sigma. Then `s06_fit_vet.py` (batman+emcee, 109 s; re-run with `exp_time=1625.35/86400`) and `s08_leo.py` (LEO-Vetter, 4.7 s `[verified]`) give physical parameters and a clean rule-based vetting pass, and `s07_null_throughput.py` shows the same light curve inverted/scrambled never exceeds BLS SNR 12 versus 30 for the real signal. None of this is a validation (no pixel localisation at depth, no FPP), but it is the complete, free, 4-CPU toolchain with every version pinned above.

---

## Sources (primary)

- TLS releases https://github.com/hippke/tls/releases ; TLS FAQ https://transitleastsquares.readthedocs.io/en/latest/FAQ.html ; TLS paper https://arxiv.org/abs/1901.02015 (A&A 623, A39)
- wotan releases https://github.com/hippke/wotan/releases ; BUGLOG https://raw.githubusercontent.com/hippke/wotan/master/BUGLOG.md ; wotan paper https://arxiv.org/html/1906.00966v3
- astropy BLS https://docs.astropy.org/en/stable/api/astropy.timeseries.BoxLeastSquares.html
- lightkurve quality flags https://lightkurve.github.io/lightkurve/reference/api/lightkurve.KeplerQualityFlags.html , https://lightkurve.github.io/lightkurve/reference/api/lightkurve.TessQualityFlags.html ; TESS SDPDD https://archive.stsci.edu/files/live/sites/mast/files/home/missions-and-data/active-missions/tess/_documents/EXP-TESS-ARC-ICD-TM-0014-Rev-F.pdf
- Kepler Bonus HLSP https://archive.stsci.edu/hlsp/kbonus-bkg (DOI 10.17909/7jbr-w430); Martinez-Palomera et al. https://arxiv.org/abs/2310.17733 ; Kepler Data Characteristics Handbook https://archive.stsci.edu/kepler/manuals/Data_Characteristics.pdf ; K2 Handbook https://archive.stsci.edu/k2/manuals/k2_handbook.pdf
- astroquery TESScut https://astroquery.readthedocs.io/en/latest/mast/mast_cut.html
- batman https://lkreidberg.github.io/batman/docs/html/tutorial.html ; Kipping 2010 https://arxiv.org/abs/1004.3741
- celerite2 https://celerite2.readthedocs.io/en/latest/api/python/ ; juliet https://juliet.readthedocs.io/en/latest/user/priorsnparameters.html ; exoplanet https://docs.exoplanet.codes/en/latest/user/install/
- LEO-Vetter https://github.com/mkunimoto/LEO-vetter , https://arxiv.org/abs/2509.10619 ; exovetter https://exovetter.readthedocs.io ; TRICERATOPS setup.py https://raw.githubusercontent.com/stevengiacalone/triceratops/master/setup.py
- nuance https://arxiv.org/abs/2402.06835 , https://nuance.readthedocs.io/en/latest/markdown/ls.html ; cuvarbase https://johnh2o2.github.io/cuvarbase/bls.html , https://johnh2o2.github.io/cuvarbase/install.html ; QATS https://faculty.washington.edu/agol/QATS/readme
- ExoMiner https://github.com/nasa/ExoMiner , https://raw.githubusercontent.com/nasa/ExoMiner/main/exominer_pipeline/run_pipeline.py ; Astronet-Triage https://github.com/mdanatg/Astronet-Triage , https://arxiv.org/abs/2301.01371 ; Kepler pipeline overview https://arxiv.org/pdf/1001.0258
- psfmachine https://ssdatalab.github.io/psfmachine/ ; transit-diffImage https://github.com/stevepur/transit-diffImage ; SHERLOCK https://sherlock-ph.readthedocs.io/en/latest/installation.html ; MonoTools https://github.com/hposborn/MonoTools
- Practitioner repos (read from local clones): https://github.com/Kuberwastaken/planethunt (pilot.py, batch2.py, independent_vet.py, pixel_localize.py, acquire.py, analysis/kbonus-extension-20261011/README.md), https://github.com/paraschopra/tess-planet-candidates (docs/GUIDE.md, pipeline/lcsearch.py, pipeline/vet.py), https://github.com/eumemic/exoplanets (README.md, pipeline/search.py, pipeline/fastbls.py, validation/*.md)
- PyPI JSON for every package in Section 0 (queried 2026-10-17)