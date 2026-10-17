# ZyushgPlanet

An archival exoplanet hunt run with Claude Code, built to do what Kuber Mehta did on
2026-10-11 (find an unreported transit candidate in old Kepler data with an AI agent),
but with the discipline his own audit showed was missing: a frozen target list, thresholds
registered before the search, inverted-light-curve null controls on every star, positive
controls that must pass, and ephemeris matching against every known object.

**Status (2026-10-17):** the first pool, 6,465 Kepler Bonus background dwarfs that were never
Kepler targets and have no published transit search, has been screened (6,365 on Modal, 100
locally; 31 had no usable cadences). First screen: 448 real vs 97 inverted passes. Stage 2:
11 survivors under the registered v1.4 rules, 23 under the post-hoc v1.5 rules, 22 under v1.7 (all reported).
Pixel localisation puts 7 on target. The independent re-derivation + red-team pass has cleared four so far,
two of them candidate-grade pending a TRICERATOPS run (`docs/10`); the remaining batches are running.
Nothing has been called a candidate.

## Read in this order

| Doc | What it is |
|---|---|
| `docs/00-executive-summary.md` | One page: what Kuber found, what we built, what we found, how close |
| `docs/01-anatomy-of-kubers-discovery.md` | What actually happened with KIC 5689121, from his evidence repo and our own archive queries |
| `docs/01b-kuber-method-verified.md` | His pipeline line by line, every count re-derived |
| `docs/02-landscape-and-open-ground.md` | The October-2026 wave (Rabtsevich, Chopra, eumemic, Kuber and about twenty more), what each covered, verified open-ground ranking |
| `docs/03-kic5689121-reproduction-and-critique.md` | Independent reproduction in this repo and a referee-style red-team review |
| `docs/04-kepler-bonus-dataset.md`, `docs/04b-other-datasets.md` | The dataset we search, and every other archive with what is left in each |
| `docs/05-vetting-and-validation.md` | The test battery, thresholds, false-positive zoo, and what "validated" means |
| `docs/06-publishing-and-community.md` | Venues, policies, precedents, how to approach professionals |
| `docs/07-feasibility-and-statistics.md` | Occurrence rates, noise, multiple testing, what can be validated at which magnitude |
| `docs/08-agentic-workflow-with-claude-code.md` | Repo architecture, guardrails, scale, templates |
| `docs/09-algorithms-and-tooling.md` | Detrending, BLS/TLS with measured runtimes, fitting, pixel tools |
| `docs/10-roadmap.md` | Where the search stands, the survivors, what each still needs, next milestones |
| `docs/briefs/` | Raw and verified research briefs with sources, JSON included |

## The search

| Path | Purpose |
|---|---|
| `targets/poolA_kbonus_bg_dwarfs_ranked.csv` | Frozen Pool A list (6,465 stars; SHA-256 in the prereg file) |
| `pipeline/prereg_poolA.json` | Pre-registered search parameters and first-screen thresholds, with a revision log |
| `pipeline/batch_search.py` | One star in, one JSON out: download, detrend, BLS, diagnostics, inverted null |
| `pipeline/run_batch.py` | Local multiprocessing runner (resumable) |
| `pipeline/modal_app.py`, `pipeline/modal_client.py` | Modal deployment driven over HTTPS; code shipped per run |
| `pipeline/triage.py`, `pipeline/ephemeris_match.py` | Apply the screen, compare real vs inverted pass rates per period bin, match survivors against KOIs, TCEs, EBs, confirmed planets |
| `pipeline/vet_batch.py`, `scripts/reverdict.py` | Stage-2 vetting (masked re-detrend, per-event depths, gap proximity, variability, TLS, known objects) and re-verdicting under a revised rule set |
| `pipeline/pixel_vet.py`, `scripts/pixel_vet_survivors.py` | Pixel-level localisation on the originating target-pixel files, with registered verdict rules |
| `pipeline/fetch_known.py` | Downloads the known-object catalogues into `data/known/` with a dated SHA-256 manifest (required before triage) |
| `scripts/incremental_vet.sh` | Fetch-triage-vet loop used while a Modal batch runs |
| `pipeline/run.py` (+ `fetch/prep/search/vet/report.py`) | Deep single-target reproduction (used on KIC 5689121) |
| `runs/<run-id>/` | Results, logs, triage tables; `runs/controls/` are the positive controls |

```sh
uv venv .venv && . .venv/bin/activate && uv pip install -r requirements.txt
python -m pipeline.fetch_known                                           # known-object catalogues + manifest
python -m pipeline.run_batch --config pipeline/prereg_poolA.json \
    --fnames kic-005689121 kic-007603200 --run-id controls --workers 2   # must pass
python -m pipeline.run_batch --config pipeline/prereg_poolA.json \
    --targets targets/poolA_kbonus_bg_dwarfs_ranked.csv --start 1 --stop 100 --run-id pilot --workers 3
python -m pipeline.triage --run-id pilot --config pipeline/prereg_poolA.json
python -m pipeline.vet_batch --run-id pilot --skip-done
python scripts/pixel_vet_survivors.py pilot
```

### Running on Modal (optional, for thousands of stars)

The Claude Code cloud sandbox cannot speak gRPC, so the Modal app is deployed once from a laptop
and driven over HTTPS; the search code and config are shipped with every run, so pipeline changes
need no redeploy.

```sh
pip install -r requirements-modal.txt            # on the laptop
export ZP_ENDPOINT_TOKEN='<shared secret>'       # any random string; the sandbox uses the same one
modal deploy pipeline/modal_app.py               # prints the start/status/results/stop URLs
# in the sandbox: .modal_endpoints.json = {"base": "https://<ws>--zyushg-planet-screen", "suffix": ".us-east.modal.run", "token": "<shared secret>"}
python -m pipeline.modal_client start --run-id poolA_main --start 101 --stop 6465
python -m pipeline.modal_client status --run-id poolA_main      # also writes runs/<id>/modal_status.json
python -m pipeline.modal_client fetch  --run-id poolA_main      # appends results.jsonl; writes config_used.json
modal app stop zyushg-planet-screen              # when no further batches are planned
```

Data caches (`data/`) are not committed; everything is re-downloadable from MAST
(Kepler Bonus HLSP, DOI 10.17909/7jbr-w430, CC BY 4.0) and the NASA Exoplanet Archive.

## Rules this repo follows

- Nothing is a "planet". Survivors of the screen are "first-screen signals"; after pixel
  localisation and ephemeris clearance they are "candidates"; validation needs FPP < 1.5 %.
- Thresholds are never changed after looking at results. Any later cut is labelled post-hoc.
- Every number in a write-up must trace to a file produced by code in this repo.
- Third-party evidence repos are read, never executed.
