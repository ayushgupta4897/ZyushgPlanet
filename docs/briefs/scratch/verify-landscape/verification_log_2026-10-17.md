# Verification log (adversarial pass), 2026-10-17 ~12:00-12:30 UTC
Raw downloads: downloads/ (ctoi_20261017.csv, toi_20261017.csv, eumemic/ sparse clone @f1d17f0, paraschopra/ sparse clone @74d5624, zenodo/*.json, arxiv/*.xml, fx_*.json, others/*.md, newrepos/*.md, moquin_poster.pdf, gh_search2.json, gaia_m60.csv, gaia_m60100.csv, tess_obs.html)

## Confirmed exactly
- ExoFOP CTOI 5,138 rows / 4,935 TICs / 4,116 CTOI-only; TOI 8,151 rows / 7,826 TICs; 1,397 FP/FA-only TOI hosts; CTOI lastmod max 2026-10-08; Aug-2026 1,173 (lafargamagro); none of the AI-agent TICs is a CTOI; TIC 300381700 CTOI by salinas at 9.072 d.
- eumemic: 2,990 allsky rows (tiers as stated), 68 earlier rows, counted 587+58=645, 630 unique TICs, Tmag median 11.91, dist 10/43 within 60/100 pc, Teff<4000 16; pipeline/search.py P_MIN,P_MAX=0.5,30; README cuts for all six samples; 'Built with Claude Code' line 1034; ~17/4/9 s per star; issues #1,#2,#3,#5,#9,#10,#11 as described.
- paraschopra: 01_targets.py constants (60 pc, Tmag<13, R<1.5, Rcont<0.2, >=3 sectors, >=2 yr, 1.5 Re @10 d ranking); 526 downloads (342 ok+184 cached); 15,332 ranked; 0.5-25 d; CTOI_DRAFT 2026-08-19 statement; Substack 2026-10-09.
- planethunt: chain 29,800->9,204->6,477->2,495->1,836->1,800 (recomputed from MAST catalogues: 9204/6477/2495/1836 exact); 1,736 searched; 31 real vs 32 inverted; RNAAS AAS82222 declined 2026-10-13 (submission-record.json); KIC 5689121: G=16.28, PSFFRAC 0.666, CROWDSAP 0.747, 13 PSF quarters, tpforg 5689116; no KOI on kepid 5689121; Kepler-169 = KIC 5689351.
- Zenodo stats/dates for 22967456, 22999133, 23175179, 23257784, 23325261, 23249265, 23118499, 23134068, 19661443.
- X: Rabtsevich 1,187,670 views/8,742 likes/830 RT/304 replies/354 quotes; Mehta 9,650/239/8/15 ("Found a new planet with Codex this weekend").
- HN: 50002665 301 pts/132 comments; Chopra 2 pts; Mehta Show HN 1 pt.
- TESS: "Orbit 229, Sector 109. Sector 110 begins on October 31, 2026"; MAST: S107 non-FFI ingest COMPLETED (start 01 Oct 2026); TESScut through S108 orbit 2b; S109 orbit 1a early release IN PROGRESS (16 Oct 2026); TESS-SPOC newest bundle Year 7 (S84-96); QLP HLSP lists sectors up to 104.
- Gaia DR3 M-dwarf counts: 103,730 at G 14-16.5 within 100 pc; 83,586 at 60-100 pc (TAP, parallax_over_error>10, bp_rp>2, M_G>7.5).
- k2pandc: 2,390 confirmed / 1,372 candidate / 293 FP / 22 refuted.
- arXiv abstracts: T16, RAVEN, ExoMiner++ (v4), ExoMiner++ 2.0, LEO-Vetter, Eschen 2024, Tschudi 2026 (v3), DIAmante, Tey 2023, DTARPS II, Olmschenk 2021, Kunimoto 2022, 2610.04567 = math.FA paper.
- SciAm 2026-10-14 (Kovac), Gizmodo 2026-10-09 (Wright) quotes; Reddit quotes (ErrorlessQuaak, admiralrewd) from RSS.
- Comdex4: 23 merged PRs, 7 titled "Claude/tess transit hunter repo alu4hv" (#1,2,5,6,8,14,21).

## Corrected / refuted
- KBonus background planet searches exist: Bienias & Szabo 2025 (A&A 694, A209; arXiv:2501.09152) 7 hot-Jupiter candidates; Moquin et al. 2026 NExScI workshop poster (Poster_JulieMoquin_41.pdf): BOI-1.01 (P=3.37 d, 2.94 Re, M dwarf G=17.8) from a matched-filter TCE pipeline (TraSHDUMP/RECYCLEBin/INCInERATOr), "expected yield of several hundred new planet candidates".
- 27,378 = ALL background quality-cut KBonus sources (incl. 1,836 M-catalogue stars); non-M = 25,542 (F 733, G 17,087, K 7,449, bp_rp>2 273); G<14 1,257; G 14-16 4,894.
- DDT Program 100 not on tess.mit.edu/public/DDTs/DDT_list.html (highest entry 099).
- Missed projects: avijit-thawani (1,098 stars 60-100 pc, retracted), alexvaldex bot (5,000+ stars), brahle (445 dossiers), darshanahirrao (1,031 stars, none), moralomonte/mega-rhexas, alplix/planodia, ashtonwallen/ExoplanetWorkbench, emaadhasan (1,066 stars), coracle-us/duo-transits, klucilla/refute.
- Mehta's agents: Codex (Kepler) + "ChatGPT Work" (follow-up) per kuber.studio; Claude Code not credited there; 1,736 searches = 1,648 unique KICs.
- CTOI-only hosts were inside eumemic's full-sky sample (searched unmasked with up to 3 iteratively-masked signals), so "never searched" overstates the gap.

## Unverifiable
- ExoFOP policy date 2026-08-19 (only Chopra's doc; no Wayback snapshots; page undated).
- Bluesky (public API connection reset / 403), Reddit search endpoints blocked (0-byte), Mastodon search returned empty on mastodon.social/fosstodon; astrodon.social unreachable.
- Zenodo full-text search is noisy; no new candidate preprints seen beyond those listed, plus F. O. Krone eclipsing-companion preprints (2026-10-16).
