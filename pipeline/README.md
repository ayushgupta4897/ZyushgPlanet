# `pipeline/` — Kepler Bonus transit search & vetting

A small, reusable Python package that reproduces, end to end, the kind of detection Kuber Mehta
reported for **KIC 5689121** (Gaia DR3 2103792461826460288): download the Kepler Bonus background-source
PSF light curve from MAST, clean and detrend it, search it blind with BLS and TLS, and vet the best
signal (odd/even, secondary, per-transit depths, empirical noise, inverted-flux null, neighbour control).

Everything runs from the repo's virtualenv:

```bash
cd /home/user/ZyushgPlanet
. .venv/bin/activate                 # Python 3.13; packages pinned in ../requirements.txt
python -m pipeline.run --gaia 2103792461826460288 --ra 285.907827 --dec 40.941988 \
    --neighbour 2103792427466721792 --out runs/kic5689121 \
    --compare-period 13.769448251 --compare-t0-bjd 2455630.2454306 --compare-duration-h 2.3
```

Wall time for that command in the sandbox (4 shared CPUs, `--threads 3`): see `runs/kic5689121/run.log`
(about 20-30 minutes; BLS 1-60 d ~2.5 min per light curve, TLS 1-30 d ~5 min per light curve, three
BLS runs + two TLS runs in total).

---

## Data product (what `fetch.py` finds)

* HLSP **KBONUS-BKG** v1.0 (Martinez-Palomera, Hedges & Dotson 2023, AJ 166, 265; arXiv:2310.17733;
  DOI 10.17909/7jbr-w430; format page https://archive.stsci.edu/hlsp/kbonus-bkg).
* MAST serves **one multi-extension FITS per source for all quarters**, e.g.
  `mast:HLSP/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits`
  (8.8 MB). Sources without a KIC id are named `gaia-dr3-<id>`. The CAOM `target_name` is the Gaia DR3
  designation, `provenance_name`/`author` is `KBONUS-BKG`.
* HDU layout (KIC 5689121: 28 HDUs):

  | HDU | content |
  |---|---|
  | `PRIMARY` | header only: `KEPLERID`, `GAIAID`, `TPFORG` (= KIC of the target pixel file the source was extracted from; here 5689116), `RA_OBJ/DEC_OBJ`, `GMAG/BPMAG/RPMAG/KEPMAG`, `PARALLAX`, `XPOSURE=1765.464` s, `TIMESYS='TDB'`, `QDETECT` (18-char Q0..Q17 detection string), `DOI`, `HLSPVER` |
  | `LIGHTCURVE_STITCHED` | all quarters, 50,182 rows: `TIME, CADENCENO, QUARTER, FLUX, FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR, PSF_FLAT_FLUX, PSF_FLAT_FLUX_ERR, SAP_QUALITY, FLATTEN_MASK` |
  | `LIGHTCURVE_Q<n>` (x13) | per quarter: `CADENCENO, TIME, FLUX, FLUX_ERR, SAP_FLUX, SAP_FLUX_ERR, CENTROID_COLUMN, CENTROID_ROW, SAP_QUALITY, SAP_BKG, RED_CHI2, PSF_FLUX_NOVA, PSF_FLUX_ERR_NOVA`; header `QUARTER, SEASON, MODULE, OUTPUT, CHANNEL, FLFRCSAP, CROWDSAP, NPIXSAP, PSFFRAC, PERTRATI, PERTSTD, ITERNEG, SAP='optimal'` |
  | `APERTURE_Q<n>` (x13) | uint8 SAP aperture mask in the originating TPF frame (e.g. 4x5) |

* **Which flux is which:** `FLUX`/`FLUX_ERR` = corrected **PSF** photometry (psfmachine LFD) — this is
  what we search. `SAP_FLUX` = simple aperture photometry inside `APERTURE_Q<n>` (diluted by neighbours;
  `CROWDSAP`, `FLFRCSAP` quantify that). `PSF_FLAT_FLUX` = the authors' flattened PSF flux (stitched
  extension only). `PSF_FLUX_NOVA` = PSF flux from the mean, un-perturbed shape model.
* **Quality flags:** `SAP_QUALITY` is the Kepler TPF `QUALITY` bitmask copied over (values present for this
  star: 0, 16, 128, 144, 1024, 2048, 2176, 2192, 8192, 8208, 8320, 8336, 262144, 262160, 262272, 262288,
  270336, 270464, 393216, 393344). `FLATTEN_MASK` (0/1) is the authors' flattening mask. There is no
  separate PSF-fit quality flag; `RED_CHI2` per cadence and `PSFFRAC`/`ITERNEG` per quarter are the
  extraction diagnostics.
* **Time:** `TIME` has `TUNIT='BJD - 2454833'` on `TIMESYS='TDB'`, i.e. **BKJD**. There is no
  `BJDREFI/BJDREFF/TIMEREF` keyword; `BJD_TDB = TIME + 2454833.0` (`fetch.BKJD_OFFSET`). `lightkurve.read()`
  parses the file natively as a `KeplerLightCurve` with `time.format='bkjd'`, `scale='tdb'`.
* KIC 5689121 has 13 quarters: Q1,2,3,4,6,7,8,10,11,12,14,15,16 (Q5/9/13/17 fall on the dead module 3 in
  season 3; Q0 and Q17 are not in the product), 131.5-1558.0 BKJD, PSFFRAC 0.67-0.87, CROWDSAP 0.75-0.96.

## Modules

### `fetch.py`
* `find_products(gaia_id, ra, dec)` — `astroquery.mast.Observations.query_criteria(provenance_name='KBONUS-BKG',
  target_name='Gaia DR3 <id>')`, falling back to `lightkurve.search_lightcurve(SkyCoord, author='KBONUS-BKG',
  radius=20")`.
* `download(...)` — `Observations.download_file(dataURL)` into `data/kbonus/<gaia_source_id>/<MAST filename>`
  (cached; `data/` is git-ignored) plus a `manifest.json`.
* `describe_fits(path)` — prints/returns the HDU list, primary header, columns, per-quarter metrics and the
  quality values seen (saved as `fits_structure.json` by `run.py`).
* `load_quarters(path, min_psffrac=0)` — list of `Quarter` objects (PSF flux, SAP flux, quality, centroids, header).
* CLI: `python -m pipeline.fetch --gaia 2103792461826460288`.

### `prep.py`
1. `stitch(quarters, quality_mode='hard')` — keep `SAP_QUALITY == 0` (`'default'` = lightkurve
   `KeplerQualityFlags.DEFAULT_BITMASK`, `'none'`), finite values, `err > 0`; normalise each quarter by its
   median PSF flux (levels differ by quarter because the PSF model/pixel footprint changes every season).
2. `detrend(lc, window_length, break_tolerance=0.5, method='biweight', mask=None)` — wotan biweight;
   segments split at gaps > 0.5 d; `edge_cutoff=0`.
   **Window = 3 x longest plausible transit.** For P = 13.8 d around a 0.507 Rsun / 0.505 Msun star:
   a = 0.0897 au, a/R* = 38.0, T14(b=0) = P/pi * asin(sqrt((1+k)^2-b^2)/(a/R*)) = **2.91 h**; allowing e <= 0.3
   at apastron (x sqrt((1+e)/(1-e)) = 1.36) gives 3.97 h; 3 x 3.97 h = 11.9 h = **0.496 d**
   (`default_window_days()`). (Hippke et al. 2019, AJ 158, 143 recommend ~3 x T14 for the biweight.)
   A blind search to 60 d would want ~0.8 d.
3. `clip_upward(lc, n_sigma=4)` — remove only points **above** the trend by > 4 robust sigma (flares,
   cosmic rays). Nothing below the trend is ever clipped, so no dip of any depth can be removed.
4. `noise_metrics(lc)` — std, robust (MAD) sigma, point-to-point scatter (std(diff)/sqrt2), median formal
   error, and a **6-h CDPP-like** number two ways: robust std of 6-h (12-cadence) bin means, and
   `lightkurve.LightCurve.estimate_cdpp(transit_duration=12)`.
5. `transit_mask(time, P, T0, dur)` + `prepare(..., mask_ephemeris=(P, T0, dur))` — **second ("masked")
   pass**: the trend is estimated from out-of-transit points only (wotan `mask=`). Without it, a window
   containing a 2-h, ~2000-ppm transit is pulled down by ~(n_in/n_window) x depth ~ 250-300 ppm, which
   (i) biases the depth low (1671 vs 2058 ppm for KIC 5689121) and (ii) creates +200-400 ppm "shoulders"
   within +/-5 h of each transit that an inverted-flux null picks up as a spurious 13.77-d dip. The blind
   searches use pass 1; all depths, odd/even, secondary and the null use pass 2.
6. `invert(lc)` — null-test copy with `flux -> 2 - flux`.

### `search.py`
* `bls_search(t, f, e, pmin=1, pmax=60, durations_h=(1,1.5,2,2.5,3,4,5,6), oversample=3, q_ref_h=2)` —
  `astropy.timeseries.BoxLeastSquares.power(periods, durations, objective='snr')`. Period grid: frequency
  step `df = q_ref/(oversample * baseline * P)` so a 2-h transit drifts by <= 1/3 of its duration over the
  1426-d baseline -> **210,246 periods between 1 and 60 d**, ~2.5 min on one core. `power` is the depth
  SNR of the best box at each period ("BLS SNR"); `sde = (power - median)/std`. The best peak is refined on
  a fine grid (+/-0.3 %, 1500 periods, 1-4 h in 0.25-h steps) and `compute_stats` is stored. Top distinct
  peaks (harmonics merged) are listed.
* `bls_refine(t, f, e, period)` — same fine grid + `compute_stats` on the masked light curve.
* `tls_search(t, f, e, pmin=1, pmax=30, oversampling_factor=3, duration_grid_step=1.1, use_threads=3,
  star=STAR)` — `transitleastsquares` with `R_star=0.507`, `M_star=0.505` (+/- ranges 0.45-0.56) and
  quadratic limb darkening `u=(0.2125, 0.3984)` (what `transitleastsquares.catalog_info(TIC_ID=120686465)`
  returns for this star; the TIC mass/radius it returns, 0.5054/0.5075, are the values used above).
  Grid: 180,843 periods for 1-30 d; measured **1.64 ms/period with 3 threads on 46k points => ~5 min**.
  The grid is restricted to 1-30 d (BLS covers 1-60 d) to keep a full run under ~30 min; no binning is
  needed because Kepler long cadence is already 29.4 min.
* `null_sde_excluding(periods, power, P)` — best SDE of a periodogram outside +/-2 % of P and its
  harmonics (P/4, P/3, P/2, 2P, 3P, 4P), for quoting an inverted null without the candidate band.

### `vet.py` (fixed ephemeris, pass-2 light curve)
* `box_depth` — `1 - mean(in-transit)` relative to the median of out-of-transit flux; error from the robust
  scatter of in- and out-of-transit points (not the formal white-noise error).
* `odd_even` — depths of odd and even epochs and their difference in sigma.
* `secondary_eclipse` — box depth at phase 0.5 plus a 200-step scan of all phases >= 3 half-durations from
  the primary (most significant other dip).
* `per_transit_depths` — local depth of every event (baseline = median of points 0.75-4 durations away);
  weighted mean, scatter, number of positive events, max single-event share, chi2 of constant depth.
* `empirical_noise` — fold at 400 random T0 offsets (>= 2 durations from the real transit) with the same
  box; the std of those depths is the empirical (red-noise-inclusive) error; `empirical_snr = depth / std`.
* `ephemeris_depth` — box depth at an *external* ephemeris (neighbour control; Kuber's ephemeris).
* `compare_ephemeris` — dP and T0 offset (minutes) at the epoch nearest mid-baseline.

### `report.py`
`write_json` (numpy-safe) and figures (`01_raw_stitched`, `02_detrended` with predicted-transit ticks,
`03_bls_periodogram` with the inverted null overlaid, `04_tls_periodogram`, `05_phase_fold` with 30-min bins
and the TLS model, `06_odd_even`, `07_secondary`, `08_neighbour_fold`, `09_per_transit`).

### `run.py`
CLI driver (see top). Options: `--quality {hard,default,none}`, `--window`, `--clip-sigma`,
`--bls-pmin/--bls-pmax/--bls-oversample`, `--tls-pmin/--tls-pmax/--tls-osf/--tls-dgs`, `--threads`,
`--no-tls`, `--no-tls-null`, `--tls-neighbour` (TLS on the neighbour too; default BLS only),
`--compare-period/--compare-t0-bjd/--compare-duration-h`, `--min-psffrac`.

## Outputs (`runs/<label>/`)
* `results.json` — everything: `primary_header`, `prep`, `noise`, `bls` (incl. `top_peaks`, `stats`),
  `tls`, `pass1_ephemeris`, `pass1_box_depth`, `bls_refine_pass2`, `adopted_ephemeris`, `vetting`
  (`box_depth`, `odd_even`, `secondary`, `per_transit`, `empirical_noise`), `inverted_null` (`bls`, `tls`,
  `*_excluding_candidate_band`), `comparison`, `neighbour` (own `bls`, `depth_at_target_ephemeris_pass2`,
  `depth_required_if_hosted_by_neighbour_ppm`), `software`, `runtime_total_s`.
  Depths inside `vetting`/`pass1_box_depth` are fractional (x1e6 for ppm); `*_ppm` keys are already ppm.
  Times: `*_bkjd` (BJD_TDB - 2454833) and `*_bjd_tdb`.
* `fits_structure.json`, `fits_structure_neighbour.json` — the FITS description.
* `lightcurve_pass1.csv`, `lightcurve_pass2_masked.csv`, `*_neighbour.csv` — time (BKJD and BJD_TDB),
  normalised flux, detrended flux, error, quarter, normalised SAP flux.
* `run.log` — the full console log with runtimes.

## Interpreting the numbers
* **SDE** (signal detection efficiency) is the peak height in units of the periodogram's own scatter. For
  white noise the TLS FAQ (https://transitleastsquares.readthedocs.io/en/latest/FAQ.html, >10,000 noise-only
  searches) maps SDE 7.0 -> FAP 1 %, SDE 8.3 -> 0.1 %, SDE 9.1 -> 0.01 %; with red noise these are
  optimistic, which is why the pipeline also reports the inverted-null SDE from the *same* light curve: a
  real transit should stand far above it.
* **BLS SNR / depth_snr** use the formal per-point errors and are optimistic; `empirical_snr` from
  random-phase folds includes correlated noise and is the number to quote.
* **Odd/even** difference > 3 sigma or a **secondary** at phase 0.5 > 3 sigma indicates an eclipsing binary
  at twice the period or a secondary eclipse.
* **Neighbour**: if the dip were on the brighter neighbour (KIC 5689116, 2.1x brighter in G), its own PSF
  light curve would have to show `depth_required_if_hosted_by_neighbour_ppm` at the same ephemeris.
* This pipeline establishes a **candidate**, not a planet: it does not do pixel-level centroid/PRF
  localisation, Gaia-neighbour dilution modelling, TRICERATOPS/LEO-Vetter false-positive probabilities, or
  stellar-density consistency fits. See the repo playbook for those steps.
