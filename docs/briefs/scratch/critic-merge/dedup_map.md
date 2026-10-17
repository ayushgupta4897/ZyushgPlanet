# Merge of audit A (newcomer/reproducibility lens) and audit B (referee lens), 2026-10-17

A: 15 gaps, 8 contradictions. B: 17 gaps, 10 contradictions.

## Gap merges (keep most specific fix)
- A3 (Kepler-138 control claim) == B9 -> one item, B9 fix (adds "c passes every stage-2 criterion except known-object veto", negative-control run), A3 detail "13/13 quarters", commit d04cfac kept.
- A4 (438/95 vs 6,465) == B8 -> one item; B8 percentages + triage.py multi-run fix; A4 detail "docs/10 later uses 6,365 denominator", "6,334 with usable data" kept.
- A12 (prereg name v1.1) folded into B3 (criterion changed post-registration without revision_log); B3 fix already sets the name.
- A15 (kic-001577049 SNR 1.9 unexplained) folded into B4 (tls_sde narrow window + no empirical SNR floor); B4 fix is the specific one, A15 footnote kept as the minimum text fix.
- B17 split: config_used.json half folded into A2 (blocking: triage command fails); pixel_vet docstring 12/13 vs 13/13 half kept as its own minor item.
- B13 (overclaims) kept separate from B1/B3 although two phrases overlap; it is the wording fix list.

Result: 32 -> 28 gaps (6 blocking, 12 important, 10 minor).

## Contradiction merges (none dropped)
- Kepler-138: A#1 + B#2 -> one line with both file citations.
- 438/95 pool: A#2 + B#3 -> one line.
- prereg v1.1 vs v1.2: A#7 + B#4 -> one line incl. B's registered_utc-vs-commit-time detail.
Result: 18 -> 15 contradictions.

## Ranking rule
blocking (science-claim blockers first, then reproducibility blockers) > important > minor.
