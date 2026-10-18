# kic-006448637 — red team (2026-10-17)

**Verdict: survives** (signal, not candidate). Notes: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006448637/redteam.md`; pixel scripts/logs in `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006448637/redteam/`.

## Why the pipeline pixel test failed and why that is not a kill
`pixel.json` reports the deficit 70-115 arcsec from the target in 0/13 quarters. `diffimages.png` shows why: the TPF belongs to KIC 6448574, a Kp 8.75 saturated star with a 40-80-row bleed column; the target sits in the bottom 1-4 rows of a 52x13 stamp, 1.7 px from the bleed column (psffrac 0.50). The whole-stamp significance-weighted centroid is set by bleed pixels with event-differenced noise of hundreds of e/s, while the target's entire eclipse deficit is ~5 e/s. "Depth in aperture" (87 ppm) is the saturated host's aperture.

My localised re-analysis of the same TPF pixels (`redteam/tpf_local.py`, `tpf_radius.py`): target-box deficit positive in all 12 usable quarters (511-1228 ppm = 3-7 e/s), neighbour box (G = 16.0 star at 2.6 px) consistent with zero, deficit complete inside r <= 1.6 px with a null 2.5-4.5 px annulus every quarter, mean centroid offset dx = 0.00 +- 0.06, dy = +0.21 +- 0.07 px (< 0.9 arcsec). `pipeline/pixel_vet.py` needs a saturated-host / edge-of-stamp mode; I did not modify it.

## Blend tests
- KIC 6448650 (G 16.03, 10 arcsec) own KBONUS PSF flux at the ephemeris: -75 +- 121 ppm vs >= 1120 ppm required (~10 sigma). KIC 6448629 (G 18.4): 971 +- 694 vs 1.0 % required (~13 sigma). KIC 6448666 (Kp 11, 0.554-d TCE, 7.6 px): not a harmonic, and the deficit is compact at the target. Fainter Gaia sources at 1.9-5.8 px would need 2.4-12 % eclipses where the centroid excludes them.
- Shape: trapezoid chi2 grid (cadence-integrated) prefers a box; pure V disfavoured by delta chi2 = 24; T23/T14 >= ~0.55 at 2 sigma. Secondary < ~130 ppm (3 sigma) at any phase. Odd/even 0.7-2.0 sigma (0.4 without epoch 28). Depth by season 857/857/728/853 +- 65-101 ppm. SAP/PSF 0.96 +- 0.07.
- Not excluded by Kepler data: unresolved companion inside ~1-2 arcsec (needs >= 2.5 % eclipse), EB at 2P with equal eclipses.

## Systematics
2 of 24 events within 0.5 d of a gap (0.7 expected); depth 812 -> 794 ppm without them. Epoch 28 (Q13 start, 2268 +- 349 ppm, 4.2-sigma pull) drives the per-year depth drift (chi2 13/4 -> 4.5/4) and should be flagged; both depth sets must be quoted. No depth trend with row position; excluding bleed-column pixels changes nothing. LS variability 32 ppm at 0.74 d. Inverted null SDE 9.9 vs 50.6. The -8 h bump in `fold_events.png` is a single brightening (epoch 39), not a dip.

## Ephemeris contamination
Host TPF channels: ch 23 / 39 / 63 / 47 in seasons 0-3. **Kepler-1715 b (KOI-308.01, KIC 6291837, Kp 12.35, P = 35.5970 d, 684 ppm) is on the same channel** (its TPF: Q2 ch 23, Q4 ch 63) 18.6 arcmin away. But P differs by 0.0097 d = 20 sigma of the 35.58777 +- 0.00045 d timing-fit period, and its transits are 0.84-1.20 d from ours at every epoch (min 0.835 d = 4.5 durations) — the event trains never overlap, so crosstalk/column-anomaly contamination (which preserves timing) is excluded. No other KOI/TCE/EB at P, P/2, P/3, P/4, 2P, 3P, 4P is on-channel-and-on-time. Nearest KOI/EB: KIC 6448768 (16.49 d, 127 arcsec), unrelated.

## Known?
Absent from KOI, DR25 TCE, Kirk+2016 EB, confirmed planets (`data/known`, 2026-10-17), Bienias & Szabo 2025 (their 7 source apertures: KIC 11033410, 4543171, 8175131, 8293539, 9040849, 9674230, 4459924), Moquin 2026 BOI-1; SIMBAD empty within 40 arcsec except the Kp 11 star at 30 arcsec. The KBonus team's own unpublished search could still hold it.

## Host consistency
T14 4.1-4.8 h at 35.59 d -> rho*(b=0) 2.1-3.3 rhosun vs host 1.7-1.9 rhosun -> b ~ 0.3-0.6; Rp = 2.4-2.5 Re undiluted. Consistent; not discriminating against a background dwarf EB.

## What would settle it
Pipeline pixel_vet rerun in local mode; TRICERATOPS; high-resolution imaging inside 1 arcsec; explicit 2P test; both depth sets and SNR 10.0 quoted; live ExoFOP/VSX/Gaia NSS cross-match.