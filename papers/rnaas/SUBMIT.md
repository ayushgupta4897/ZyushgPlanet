# Submitting the RNAAS note

Files: `note.tex` (AASTeX 7.0.1, RNAAS option), `figure1.pdf`, `build/note.pdf`, `cover_letter.txt`.
Every number in the note is filled by `build_note.py` from `report/data.json` and the transit fits; `numbers.json` records them.
Current count (`word_count.json`): see file; about 930 words including title, headers, caption and references (limit 1,500); abstract 94 (limit 150); one figure; no line numbers; first person singular.

Submitting is public: RNAAS publishes within about 72 hours of acceptance, indexed in ADS, under CC BY 4.0. A decline is also part of the record. Do the steps in order.

## 1. Before anything goes public (this repo's rules, `CLAUDE.md`)
- [x] Inspect every plot yourself (done by Aayush Gupta, 2026-10-17; `runs/poolA_main/human_review.json`): `candidates/*/stage2/fold_events.png`, `candidates/*/pixels/diffimages_v16.png`, `candidates/*/pixels/season_stack/season_maps.png`, `candidates/*/deep_run/0[5-9]_*.png`, `candidates/*/transit_fit/fit.png`, and `figure1.png`.
- [ ] Read the note and the cover letter and change anything you would not say yourself. The AAS policy allows AI drafting only when the author supervises the text.
- [ ] Recommended: ask one professional astronomer to read the note first (`docs/06-publishing-and-community.md` §6 lists who and how). Kuber's note was declined at moderation within 11 hours.

## 2. Archive the data
- [x] Bundle built for Aayush Gupta (ORCID 0009-0006-1165-0828).
- [x] On zenodo.org: New upload, choose Dataset, upload `papers/rnaas/zenodo/*.zip`, copy title, description, license and keywords from `papers/rnaas/zenodo/.zenodo.json`, then **Reserve DOI** (reserved: 10.5281/zenodo.23442225).
- [x] Zenodo record published 2026-10-17; DOI resolves; file checksum matches the local bundle (md5 7c7dcf9b60badf009445308388cc043a).

## 3. Finalise the manuscript
```sh
. .venv/bin/activate
python papers/rnaas/build_note.py   # author, email, ORCID and DOI already set
```
- [x] DOI in `note.tex` and `cover_letter.txt`.

## 4. Submit
- [x] AAS Journals portal (https://aas.msubmit.net), journal **Research Notes of the AAS**, article type Research Note.
- [x] Upload `note.tex`, `figure1.pdf` and, if asked, the compiled PDF.
- [x] File types: `note.tex` = Manuscript file; `figure1.pdf` = Figure file (keep this exact name, `\plotone` loads it); `build/note.pdf` = Author PDF.
- [x] Data questions: MAST = **Yes** (Kepler Bonus HLSP 10.17909/7jbr-w430; Kepler long-cadence Q0-Q17 10.17909/T9488N). IPAC = **Yes** (NASA Exoplanet Archive: KOI cumulative 10.26133/NEA4, DR25 TCE 10.26133/NEA15, Kepler stellar 10.26133/NEA6, Planetary Systems 10.26133/NEA12; accessed 2026-10-17, `data/known/manifest.json`). All six appear in `note.tex` as `\dataset` links, with the Zenodo DOI 10.5281/zenodo.23442225.
- [x] When the portal asks about generative AI, state: "I used Claude (Anthropic) to help write analysis code and draft the text; I reviewed all results and the text and take responsibility for them."
- [x] Paste `cover_letter.txt` into the cover-letter field.

- [x] **Submitted 2026-10-18; AAS manuscript # AAS82486** (receipt acknowledged by the AAS Editorial Office). Intended upload: `note.tex` and `build/note.pdf` from commit 80544de (with the archive DOI sentence) and `figure1.pdf`; confirm in the portal PDF under "Check Manuscript Status". Digital assets: none. UAT: Exoplanets (498), Transit photometry (1709), Exoplanet detection methods (489), Mini Neptunes (1063), Light curves (918), Astrostatistics (1882). Portal abstract entered as plain text (same wording as `note.tex`).

## 5. After acceptance
- Until a decision: keep the report, screens and poster private; README stays silent.
- The note is the citable record; you can then post the screens and the report link with the RNAAS DOI.
- ExoFOP community candidates: ExoFOP accepts an RNAAS as the publication only if the search methodology was previously published in a refereed journal. This pipeline is not, so expect the upload to need a methods paper first (`docs/06` §0, item 5).
- Next scientific step: speckle or adaptive-optics imaging of the three hosts (G 16–17.4), ideally with a professional collaborator.
