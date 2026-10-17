# 00 · Executive summary

**The milestone.** Do what Kuber Mehta did on 2026-10-11 (find an unreported transiting-planet signal in archival space-telescope data with an AI coding agent), but in a way that survives a referee.

**What his find really is.** A real, repeating, on-target transit signal (P = 13.769 d, ~2 ppt, 70 events) on a faint red dwarf that was never a Kepler target, pulled out of the 2023 Kepler Bonus background-star photometry. It is a *candidate* with a 6.9 % false-positive probability, it cannot be validated with any telescope an amateur can reach, his RNAAS note was declined in a day, and his own inverted-light-curve controls showed his first screen passed noise as often as data. The tweet said "planet"; the repo says "candidate". (`docs/01`, `01b`, `03`)

**The wave he is part of.** Between 2026-09-25 and 2026-10-16 at least twelve AI-agent or citizen archival searches went public; one (eumemic) swept ~4.17 million TESS dwarfs. The community-TOI portal has accepted only published candidates since August. The verified open-ground ranking put the non-M Kepler Bonus background dwarfs first for novelty and faint nearby TESS M dwarfs first for validatability. (`docs/02`, `04b`)

**What we built.** A frozen, hashed target list of 6,465 never-searched Kepler Bonus background dwarfs; a pre-registered two-stage screen with an inverted-flux null on every star; positive controls that caught three real bugs before any pool star was judged; a pixel-level localisation test validated on Kuber's star; a Modal deployment driven over HTTPS (the sandbox cannot speak gRPC); ephemeris matching against every known object; and a repo-level rule set (`CLAUDE.md`). (`pipeline/`, `docs/05`, `08`, `09`)

**What we found (2026-10-17).** 6,365 stars screened on Modal in 2 h 17 m plus 100 locally: 438 real vs 95 inverted first-screen passes; 11 signals survive stage-2 vetting; 3 of them put the missing light firmly on the target star in the Kepler pixels (kic-006699887 at 15.1 d, kic-009520080 at 2.4 d, kic-010213563 at 36.8 d), 3 are plausible, 5 are weak or off-target. Independent re-derivation and red-team review are in progress. Nothing is called a planet. (`docs/10`)

**How close.** To a defensible *candidate*: close, pending the adversarial pass and a false-positive-probability run. To a *validated planet*: far, for the same reason as Kuber (hosts at ~1 kpc, G 16–17). To a *validatable* find: that is Pool B (faint nearby TESS M dwarfs), already scoped. (`docs/07`, `10`)

**Honest caveats.** The pixel test uses an empirical point-source model, not the official Kepler PRF; per-event errors understate correlated noise (calibrated on controls); the stage-2 criteria were registered after a 208-star partial look at the batch (logged as such); host radii are colour-based estimates; a Goddard/KBonus team expects "several hundred" candidates from this same dataset and may publish first.
