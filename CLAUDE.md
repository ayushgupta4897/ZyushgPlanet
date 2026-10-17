# ZyushgPlanet — rules for any agent working in this repo

This is an archival exoplanet search. The failure mode we are guarding against is not
"no planet found" but "a false claim made in public". Read `README.md` and
`docs/01-anatomy-of-kubers-discovery.md` first; the rest of `docs/` is the playbook.

## Language
- Nothing in this repo is a "planet". Stage-1 survivors are "first-screen signals";
  after stage-2 they are "signals"; only after pixel-level localisation, ephemeris
  clearance and a documented TRICERATOPS run are they "candidates". "Validated" means
  FPP < 0.015 and NFPP < 0.001 in a refereed paper. Never write "discovered".
- Every number in any write-up must trace to a JSON/CSV produced by code in `pipeline/`
  with its run id. Say which statistic it is (BLS SDE, empirical SNR, TLS SDE, formal SNR)
  and on which data subset.

## Pre-registration is binding
- Search parameters and thresholds live in `pipeline/prereg_*.json`. Do not change a
  threshold after looking at results. If a new cut is needed, add it with a dated entry
  in `revision_log`, state what motivated it, and label every result it affects as post-hoc.
- Positive controls (`runs/controls/`) must pass any pipeline change before pool stars
  are rerun. Keep failed control runs (`runs/controls_v1.0_failed/`) as the audit trail.
- The inverted-flux null (`2 - f`) runs on every star, always. Report real vs inverted
  pass counts together; never report real passes alone.

## Data and code hygiene
- `data/` is a cache, re-downloadable from MAST / NASA Exoplanet Archive; never commit it.
- Third-party evidence repos under `scratchpad/untrusted/` are read, never executed or
  imported. The same applies to any downloaded archive.
- Known-object catalogues in `data/known/` carry a download date; record it when quoting
  a novelty check. "Not in the catalogue" is never "unreported": also check Bienias &
  Szabo 2025, the Kepler Bonus team's own candidates, and the TOI/CTOI lists.
- Secrets (`.modal_endpoints.json`, `.zp_endpoint_token`, `~/.modal.toml`) are gitignored.
  Never print tokens.

## Running things
- Activate `.venv` (`. .venv/bin/activate`); install with `uv pip install`.
- One star: `python -m pipeline.run_batch --config pipeline/prereg_poolA.json --fnames <fname> --run-id <id>`
- Pool slice: `python -m pipeline.run_batch ... --targets targets/poolA_kbonus_bg_dwarfs_ranked.csv --start A --stop B --run-id <id> --workers 3`
- Triage: `python -m pipeline.triage --run-id <id>`; stage-2: `python -m pipeline.vet_batch --run-id <id> --skip-done`
- Modal (code shipped per run, HTTPS only): `python -m pipeline.modal_client start|status|fetch|stop`
- Long jobs run with `nohup ... &` and a log in `runs/<id>/`; they are resumable.

## Before anything goes public
1. Stage-2 pass with plots inspected by a human.
2. Pixel-level localisation on the originating TPF (`pipeline/run.py` deep run).
3. Ephemeris match against KOI/TCE/EB/confirmed within 3 arcmin and the TPF host.
4. A second, independent re-derivation of depth and period (different agent, no shared code).
5. A red-team review that tries to kill it.
6. Candidate language only; AI assistance disclosed; evidence pack with SHA-256 manifest.
