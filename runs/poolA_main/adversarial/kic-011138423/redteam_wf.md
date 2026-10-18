# kic-011138423 red team (2026-10-17): verdict SURVIVES

Notes: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-011138423/redteam.md`; scripts, JSON, logs and `sha256_manifest.txt` in `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-011138423/redteam/`. Own code only; no `pipeline/` imports.

## Killed
1. **Bright TPF host KIC 11138418** (G 12.5, 11.4"): its own KBONUS light curve folded at the ephemeris gives -8.1+-8.9 ppm (PSF), -5.0+-3.7 ppm (SAP); with a measured 0.1-0.2 %/pixel PRF wing at the target it cannot supply the 1.7-1.9 e/s lost by the target.
2. **Resolved Gaia neighbours**: required eclipse depths >= 150 % (all G > 18.5) or 17-41 % (G = 16.79 star 17" away, whose position shows no deficit).
3. **Ephemeris contamination**: 203 KOI/TCE/EB harmonic matches within 0.2 %; the only period-coherent, phase-aligned one is Kepler-20 f (4P, 98 ppm, 426', different module), a chance hit that could populate only every fourth event. Only this star of 6,365 has a real detection near 4.894 d.
4. **Already known?** No KOI/TCE/EB/PS within 5', no TOI/PS in a 3' box, SIMBAD empty at the position, not in Bienias & Szabo 2025 (7 candidates) nor the Kepler Bonus public candidate.
5. **Variability/systematics**: no LS power at P or harmonics (amp ~100 ppm); per-event chi2/dof 1.38 is the 73rd percentile of random ephemerides on the same data; near-gap vs other events 1114+-215 vs 1239+-67 ppm; per-year t0 within +-7 min; flagged cadences 7.6 % in transit vs 8.0 % overall; PSF-fit chi2 unchanged in transit.
6. **Shape/phase**: odd/even 0.35 sigma; secondary +51+-106 ppm (< 27 % of primary at 3 sigma); no dip at any other phase; duration 1.8-2.0 h matches a central transit of the 0.53 Rsun dwarf.

## Finding that must enter the evidence pack (not a kill)
KBONUS fractional depth is **season-dependent**: S0 900+-105, S1 1243+-175, S2 1548+-213, S3 1698+-180 ppm (S0 vs S3: 3.8 sigma; per-quarter chi2 38.6/16 vs random-ephemeris max 26.5). Cause: the un-rescaled per-quarter PSF flux of the star varies 1017-1789 e/s with the stamp geometry (host-wing contamination; the stitched FLUX hides this because every quarter is rescaled to 1357 e/s). The **absolute lost light is constant**: raw PSF 1.73 e/s (chi2 14.4/16), raw SAP 1.58 e/s (15.1/16), NOVA same, pixel-level 1.80 e/s. A leaking neighbour would vary with geometry; a source fully inside the photometry does not. Radius must be quoted as **2.0-2.5 Re** (depth 1200-1770 ppm depending on the star's true flux), not 2.1 Re.

## Independent pixel localisation (season stacks, empirical PRF)
Target pixel is the strongest deficit pixel in all four seasons (SNR 7.7/7.4/4.8/5.7); PRF-at-target beats null by delta chi2 48/43/27/21 and beats PRF-at-host by 41/34/25/12; best positions 1.7-2.9" from the Gaia position in different sky directions each season; moving the source beyond the stamp edge is disfavoured in the three clean seasons (S0 is systematics-dominated, chi2_null 157/34 px).

## Left for validation
Sub-pixel blend (< 2-3"), twin EB at 2P, on-target grazing low-mass companion. Settle with TRICERATOPS (2-3" exclusion), speckle/AO, a resolved-shape transit, and a contamination-immune depth definition (absolute e/s over an independent flux).