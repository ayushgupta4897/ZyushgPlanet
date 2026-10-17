# Red-team review: kic-009520080 (KIC 9520080, Gaia DR3 2127366472162960768)

Full notes: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-009520080/redteam.md`; scripts and outputs in `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-009520080/redteam_out/`.

## Verdict: SURVIVES, with two open items that must travel with it

No non-planet explanation or instrumental origin was found. Open items: (1) a marginal 2.2-2.8 sigma brightening at phase 0.5 (no planet produces one; the signal is a formal stage-2 FAIL under the post-hoc v1.5 rule (c), which `prereg_poolA.json` records was added because of this star; `stage2_pass_v14 = true`); (2) the generic unresolved-blend risk (G > 21 background EB inside the 4-arcsec pixel, or a bound companion EB), which only TRICERATOPS-style statistics and follow-up can bound.

## Key results
- **Pixel localisation (own code, 17 TPFs):** deficit in the target's pixel every quarter, combined 1.02 +- 0.04 e/s (24.5 sigma), matching the expected 1.5 e/s in-stamp missing light for 1570 ppm on the target; host-dominated pixels -0.33 +- 0.51 e/s (-11 ppm). Host blend excluded by > 50x. KBONUS light curves: KIC 9520088 (G 14.8) -6 +- 12 ppm, KIC 9520089 (G 13.9) -190 +- 210 ppm at the ephemeris. The G 20.9 star at 3.5 arcsec shows only PSF spill (0.24 vs 1.01 e/s) where its pixel differs from the target's.
- **Ephemeris contamination:** no KOI/TCE/EB/confirmed planet within 3 arcmin; closest field-wide period match of consequence is Kepler-113 c at 2:1 (748 ppm, 4.8 deg away, epoch phase 0.48, no odd/even alternation): not a contaminant. Channels per season 62/46/22/38; same depth on all four.
- **Systematics:** gap_frac_1d 0.09, no gap-distance dependence; the per-event chi2 inflation and the low third-quartile depth are one event (n=342, BKJD 945.09, isolated +6000-8000 ppm spikes, error under-estimated); without it reduced chi2 ~ 1.5 as in the controls.
- **Variability:** LS peak 0.286 d (~100 ppm), unrelated; power at P, P/2, 2P negligible.
- **Phase 0.5:** mean +192 ppm (z 2.84, rank 1/177 boxes), clipped mean +126 (z 2.24), median +51 (z 1.32), per-epoch-median +237 +- 69 vs control +46 +- 69 (z 2.7); 13 >3-sigma cadences at 13 different epochs (not one flare); consistent across seasons at 1.1-2.4 sigma each. Unexplained, marginal, not a false-transit mechanism.
- **Known?** Not a KOI/TCE/EB/confirmed planet, not in Bienias & Szabo 2025, no TOI/CTOI within 5 arcmin, SIMBAD clean; KIC 9520080 was never a Kepler target. Gaia DR3: RUWE 1.15, ipd_frac_multi_peak 0, non_single_star 0; neighbour list complete to 15 arcsec.
- **Host consistency:** a/Rs 9.5, T14(b=0) 1.99 h, observed 1.125-1.15 h -> b ~ 0.85, Rp ~ 3.0 Re (box) to ~3.3 Re (smearing-corrected); not discriminating against a blend on a denser star but not in tension with the K dwarf.

## What would settle it
Ground-based transit photometry at predicted times; TRICERATOPS FPP with a contrast curve; AO/speckle imaging; an independent look at the phase-0.5 brightening; and the v1.5 stage-2 FAIL carried in any report.