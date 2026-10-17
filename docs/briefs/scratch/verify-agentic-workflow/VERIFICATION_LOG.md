# verify-agentic-workflow — what was checked (2026-10-17)
Sandbox: Python 3.13.16, uv 0.11.32 (has --generate-hashes), Claude Code 2.1.302, node 22.22.0, jq 1.7, 4 CPUs.
- Kuber clone @09f4aa0 (2026-10-13): read run_search.py, pilot.py, acquire.py, watch_search.py, known_control.py, bootstrap.sh, requirements-*.txt,
  SEARCH_AUDIT.md, kbonus-extension README, editorial-audit README, independent_check README, TASKS.md, SCIENCE.md (FAP correction),
  publication/2026-10-13/{claim-sources,figure_sources,submission-record,word-count,bulk-data-manifest}.json, FINAL_STATUS.json, planethunt.md,
  catalogs/provenance.json (5 entries; pscomppars + exofop_toi recorded as errors), novelty_report.json keys, FOLLOWUP_PLAN.md (5 emails).
  FINDING: pilot.py one() caches with `if dest.exists(): return json.loads(...)` -- no status/version/sha check (run_search.py has it).
  Runtimes: kbonus 129 files median 16.0 s p90 27.0 max 35.4 (0.60 CPU-h); TESS regular 175 files median 6.1 p90 11.6 max 19.7; USP 176 files median 3.0.
- guard.sh tests: denies as claimed; BYPASS: `cd third_party/untrusted/... && python run_search.py` -> no decision. Read tool not matched by matcher.
- node --check of workflow body with stubbed hooks: OK.
- lightkurve: search_tesscut -> [14,26,40,41,53,54,74,75,80,81] in 28.4 s; KBONUS-BKG 2 products (0.0", 3.07") in 213.7 s; TESS LC 14 products TARS/TASOC in 30.3 s.
- PyPI: triceratops 1.1.0 (2026-09-02, >=3.8, pytransit==2.2, numpy<2 for py<3.13, numpy>=2.1 for py>=3.13); leo-vetter 1.2.0 (2026-07-29).
- NASA TAP tables present: pscomppars ps toi cumulative q1_q17_dr25_koi q1_q17_dr25_sup_koi q1_q17_dr25_tce k2pandc keplerstellar q1_q17_dr25_ks keplernames k2names. `fpwg`, `koi` -> ERROR.
- ExoFOP CSV endpoints download_toi.php / download_ctoi.php (output=csv) OK.
- KBONUS mstars catalogue Content-Length 2,390,400 B.
- Docs fetched: memory, hooks (3 chunks), workflows, tools-reference, cloud-environments, routines, sub-agents, skills, scheduled-tasks, headless,
  claude-code-on-the-web, permissions; ExoFOP news; Zenodo developers + github/describe-software; AAS LLM + RNAAS; KBONUS HLSP; MAST acks; Gaia credit;
  Chopra GUIDE.md/README.md/RESPONSE_TO_CODEX.md; eumemic README/cloud.py/params.toml; modal pricing; MAST TESS holdings; vantage c8g/c6i.16xlarge; CFF schema guide; uv pip compile; LEO-Vetter README.
