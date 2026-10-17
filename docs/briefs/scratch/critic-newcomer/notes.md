# Newcomer-lens critique notes (2026-10-17)

## Numbers checked against runs/
- poolA_main/triage.md: 6365 stars (ok 6365, errors 0, skipped 31 = too_few, n_after_quality 0); real 438 / inverted 95; mean 125 s/star. survivors.csv 438 rows.
- poolA_pilot/summary.json: 100 stars, real 10 / inverted 2; vet_summary: stage2_pass 0/10.
- poolA_main/vet_summary.csv: stage2_pass True = 11; names match docs/10 table; P, depth, dur, SNR, n_events, frac_positive, rp_earth_box all match (rounded).
- pixel.json summaries for all 11 match docs/10 col "Pixel localisation" (1.1"/12/13/7.1; 1.7"/13/17/5.2; 2.3"/8/13/5.7; 2.9"/9/12/3.1; 3.1"/10/12/3.7; 3.6"/8/13/3.1; 2/14 (other closer 12/14); 6 q, 4.6"; 7.9"/2.5; 0/13; 28.6"/1510 ppm).
- controls/pixel/kic-005689121: 13/13 quarters -> matches docs/10 reading guide.
- controls: kic-007603200 stage2_pass False (fails known_object); no pixel dir -> docs/10 row "Kepler-138 c + d pass stage 2 and pixel localisation" NOT backed. s2 signal 23.089 d exists in results.jsonl.
- target CSV sha256 matches prereg; 6465 rows; pilot ranks 1-100, main 101-6465.
- Wall time "2 h 17 m": modal_call.json 12:34:33Z; results complete between 14:50 and 15:00 (incremental_vet.log). Plausible, no file states it.
- requirements.txt lacks modal (modal_app.py imports it); modal not installed in venv.
- data/ gitignored incl. data/known/*.csv; ephemeris_match._load reads them with no downloader; triage imports ephemeris_match -> fresh clone cannot triage.
- triage --run-id poolA_main fails without --config (no config_used.json for Modal runs). CLAUDE.md gives the bare form.
- README --help: run_batch, triage, vet_batch, modal_client, run all parse. README status para stale ("is being screened").
