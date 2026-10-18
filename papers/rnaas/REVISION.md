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

Length: 1,050 words (was 934; limit 1,500); abstract unchanged (94 words); still one figure.

## Text for the editor, if a revision is requested

> I have added three items: the inverted-flux null carried through the second stage and the pixel test (6 of 97
> inverted signals pass the second stage and none localizes on target, against 22 and 5 of the real signals); the
> empirical S/N, number of transits and BLS SDE (with the inverted-flux value) for each candidate; and a TESS
> sensitivity estimate showing the hosts are too faint for TESS to test the signals (combined S/N at most 1.0). No
> result or conclusion has changed.
