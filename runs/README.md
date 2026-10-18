# `runs/` — outputs behind every reported number

Each directory is one run id. Configs are copied in as `config_used.json`; the pre-registration they follow is `pipeline/prereg_poolA.json`.

| Run | What it is |
|---|---|
| `controls/` | Positive controls (KIC 5689121, Kepler-138 and the pixel-test controls) that must pass before pool stars are run |
| `controls_v1.0_failed/` | The first control run that failed; kept as the audit trail |
| `poolA_pilot/` | Pilot slice of the target list |
| `poolA_main/` | The full search of 6,465 Kepler Bonus background dwarfs (`targets/poolA_kbonus_bg_dwarfs_ranked.csv`) |
| `kic5689121/` | Reproduction of the KIC 5689121 detection with `pipeline.run` |
| `kic*_deep/` | Deep light-curve runs (`pipeline.run`) on the signals that survived stage 2 |

Inside `poolA_main/`:

| Path | Produced by | Holds |
|---|---|---|
| `results.jsonl` | `pipeline.run_batch` (on Modal) | First-screen BLS results for every star, real and inverted flux (`2 - f`) |
| `triage_all.csv`, `triage.md`, `survivors.csv` | `pipeline.triage` | First-screen pass/fail, real vs inverted counts, ephemeris matches |
| `vet/<star>/vet.json`, `vet_summary.csv`, `survivors_stage2.csv` | `pipeline.vet_batch` | Stage-2 vetting per signal; fold plots kept for the 22 stage-2 passes |
| `pixel_v15/`, `pixel_v16/`, `pixel_season/` | `pipeline.pixel_vet`, `pipeline.pixel_season_stack` | Pixel-level localisation (per quarter; per season, v1.8) |
| `adversarial/<star>/` | independent re-derivation and red-team review | `rederive_result.json`, `redteam.md`, `redteam_verdict.json`, the re-derivation scripts |
| `fpp*/`, `molusc/`, `gaia_contrast/` | `pipeline.fpp_triceratops`, `scripts/molusc_run.py`, `scripts/gaia_contrast_curve.py` | False-positive probabilities and their companion constraints |
| `transit_fit/` | `scripts/transit_fit.py` | Transit fits and the stellar-density check |
| `ephemeris_fieldwide/` | `scripts/ephemeris_fieldwide.py` | Field-wide ephemeris matching |
| `human_review.json` | — | Record of the human inspection of every plot |

Not committed, because the scripts above regenerate them: downloaded light curves and target-pixel files (`data/`, from MAST), periodogram arrays (`*.npz`, `*.npy`), light-curve dumps, MOLUSC companion tables (hashes in `molusc_summary.json`) and fold plots for signals that failed stage 2.
