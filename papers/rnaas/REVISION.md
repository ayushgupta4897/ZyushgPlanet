# Prepared revision of AAS82486 (not submitted)

The submitted files are frozen in `submitted_AAS82486/` (SHA-256 in its `SHA256SUMS`). `note.tex` and `build/note.pdf`
in this folder are the prepared revision. Use it only if the editor asks for changes or allows minor revisions before
publication; do not replace files unprompted.

## Changes (all numbers from run outputs; see `numbers.json`)

1. **Inverted-flux null carried through the later stages.** New sentence in Vetting: of the 97 inverted first-screen
   signals, 6 pass the same second-stage rules and none localizes on target in the pixel test (2 offset, 4 inconclusive),
   against 22 and 5 for the real signals. Sources: `runs/poolA_main/vet_inverted_summary.json`,
   `runs/poolA_main/pixel_v16_inverted_summary.json`; prereg `revision_log` v1.9 (post hoc, dated 2026-10-18).
2. **Detection statistics for each candidate.** Empirical S/N (depth over the scatter of same-duration boxes at random
   phases), number of transits, BLS SDE and the strongest inverted-flux SDE: 20.0 / 86 / 89 (14), 23.8 / 86 / 46 (13),
   15.2 / 253 / 33 (15). Source: `report/data.json` (stage-2 `vet.json` and first-screen results).
3. **TESS.** Each host was observed in 9–10 TESS sectors, but at T = 15.5–16.5 an undiluted, ideal-aperture estimate
   (Stassun et al. 2018 noise model) gives a combined transit S/N of at most 1.0, so TESS cannot test the signals.
   Source: `runs/poolA_main/tess_sensitivity.json` (`scripts/tess_sensitivity.py`). Adds one reference.

4. **Sample selection with counts.** 425,219 background sources → 118,648 (PSF fraction > 0.5) → 64,075 (crowding > 0.7)
   → 27,378 (≥ 8 quarters) → 13,696 (G < 17.5) → 9,513 (never Kepler targets) → 9,335 (excluding the 1,800 sources of
   Mehta 2026, now cited) → 6,465 Gaia DR3 dwarfs. Source: `docs/data/target_selection_counts.json`
   (`scripts/target_selection_counts.py`, reproduces the frozen list exactly).
5. **Re-derivation agreement.** The independent re-derivation reproduces each candidate's period within 10 s and depth
   within 7 % (largest differences 9.2 s and 6.2 %). Source: `report/data.json` (`rederive`).
6. **Scope.** One sentence: the search is not calibrated for completeness and no occurrence rate is inferred.

Length: 1,097 words by our count (was 934); 846 by the AAS count (everything after the abstract except references; limit
1,350); abstract unchanged (94 words); still one figure.

## The six points the same editorial office raised on AAS82207 (another AI-assisted note, 12–15 Oct 2026)

| AAS82207 point | Submitted AAS82486 | Prepared revision |
|---|---|---|
| 1. First person singular, including the AI statement | Yes | Yes |
| 2. Terms used consistently, in their established sense | Yes (first-screen signal / signal / candidate, defined) | Yes |
| 3. Explicit counts for every sample cut | Cuts listed, counts only at the end | Every cut counted (item 4) |
| 4. Validation of the photometry | Positive controls; independent re-derivation mentioned | Adds re-derivation agreement (item 5) |
| 5. What the search could have missed | Not stated | Scope sentence (item 6) |
| 6. Audit of failures | Rejection classes listed | Adds the inverted null through every stage (item 1) |

## Text for the editor, if a revision is requested

> I have added the following: the inverted-flux null carried through the second stage and the pixel test (6 of 97
> inverted signals pass the second stage and none localizes on target, against 22 and 5 of the real signals); the
> empirical S/N, number of transits and BLS SDE (with the inverted-flux value) for each candidate; and a TESS
> sensitivity estimate showing the hosts are too faint for TESS to test the signals (combined S/N at most 1.0). No
> result or conclusion has changed. I have also given the count after each sample cut (including the exclusion of an earlier
> public search of the same data), the agreement of the independent re-derivation, and a statement that no occurrence
> rate is inferred.
