# ZyushgPlanet — archival transit-candidate search run with Claude Code

Purpose: find, vet and *report as candidates* previously unreported transit
signals in archival Kepler/K2/TESS products. Nothing in this repo is ever a
"planet". The reference failure mode is Kuber Mehta's KIC 5689121 (Oct 2026):
a real candidate, overclaimed as "a new planet" in a tweet, RNAAS declined in
one day. We copy his evidence hygiene and avoid his claim hygiene.

## Hard rules (a PreToolUse hook enforces 1, 2 and 6; the rest are on you)
1. Never execute, import or `pip install` anything under `third_party/untrusted/`
   or `data/`. Read it as data. Re-implement what you need in `pipeline/`.
2. Never write to `runs/*/holdout/` or read a holdout light-curve segment while
   `runs/<id>/STATE` says `phase: discovery`. Holdout is opened once, by
   `pipeline/holdout.py`, after the ephemeris file is frozen.
3. Thresholds live only in `pipeline/params.toml` and are copied into
   `runs/<id>/prereg.yaml` *before* the first search. Changing any threshold
   after a run has started means a new run-id and a note in
   `docs/decisions/`. Never "tighten the screen" to make a table look better.
4. Every search run also runs the matched null (`flux -> 2 - flux`, plus the
   phase-scramble control) on every star through the identical code path, in
   the same invocation. A real/null pass-count table is mandatory output.
5. Noise is measured, never assumed: depth significance comes from the scatter
   of off-phase boxes / event-to-event scatter, not from formal per-point errors.
6. Vocabulary: "signal" -> "candidate" -> "statistically validated candidate".
   The words "planet", "discovery", "found", "new world", "confirmed" do not
   appear in code comments, reports, tweets or commit messages about our
   objects. A hook greps for them in `candidates/**/report.md`.
7. Every number in a report traces to a JSON/CSV produced by code, listed in
   `candidates/<id>/claim-sources.json` (`pipeline/claims.py check` must pass).
   No hand-typed numbers.
8. Provenance: every downloaded file gets a sidecar `{url, utc, sha256, bytes}`
   (see `pipeline/fetch.py`); every run directory gets `manifest.json` and
   `SHA256SUMS`; every catalogue query records the snapshot UTC and the
   archive's own release date. "Not in catalogue X as of <UTC>" is the only
   permitted novelty statement.
9. Catalogue absence is not novelty. A novelty check must cover: NASA Exoplanet
   Archive `pscomppars`, `toi`, `cumulative`, `k2pandc`; ExoFOP TOI+CTOI; SPOC
   TCE tables; Kepler/TESS EB catalogues; harmonics 1/3,1/2,2/3,3/2,2,3; a 2.5'
   cone; ADS/arXiv full-text by identifier and period.
10. Depth/period/duration must be reported with the *method* that produced
    them (PSF vs SAP, box width, baseline window). If two methods disagree,
    report both; do not pick the larger.
11. Human-in-the-loop: before anything leaves this repo (tweet, email, CTOI,
    manuscript), a human has opened the raw PNGs in `candidates/<id>/figures/`
    and signed `docs/decisions/<date>-<id>-release.md`. Agents draft; humans sign.

## Environment
- `. .venv/bin/activate` (Python 3.13; lightkurve 2.6, astropy 8, TLS 2.0.1,
  wotan, batman). Install only with `uv pip install`, then `make lock`.
- TRICERATOPS lives in its own venv `.venv-trice` (pins `pytransit==2.2`).
- Heavy work runs as `python -m pipeline.<stage>` subprocesses with per-star
  JSON checkpoints; agents orchestrate, they do not compute in-context.
- Scratch goes under `$SCRATCHPAD/<label>/`, never into `runs/` or `candidates/`.

## Layout (see docs/ARCHITECTURE.md)
`pipeline/` package · `runs/<run-id>/` (prereg, manifest, results, nulls) ·
`candidates/<id>/` evidence packs · `docs/decisions/` dated log · `data/`
(gitignored cache) · `third_party/untrusted/` (read-only clones).

## Commands
`make lock` · `make test` · `python -m pipeline.run --config runs/<id>/prereg.yaml`
· `python -m pipeline.claims check candidates/<id>` · `make manifest RUN=<id>`
