# Red-team review: kic-005609180 (KBONUS-BKG PSF light curve, 21.106 d first-screen signal)

Date 2026-10-17. Candidate language only; nothing here is a planet. Full notes with file references: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005609180/redteam.md`; data products under `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005609180/redteam/`.

## Verdict: KILL

The 21.106 d signal is a coincidence of three post-gap PSF-model recovery ramps in the KBONUS stitched PSF flux (epochs 4, 36, 39), all in season-0 quarters on channel 17 (module 7, output 1). With the first day after every gap masked, the period is undetectable (15-30 d BLS: SDE 1.66 at 21.1 d, rank > 40) while the identical inverted-flux null reaches SDE 8.31. The remainder, beyond the ~2 d ramp timescale, is a 3.1-3.3 sigma PSF-only residual (201-220 +- 65 ppm) that is 1.5-1.8 sigma in SAP (97-113 +- 63 ppm) despite CROWDSAP 0.95-1.0, has chi2 56-59 for 43-46 dof, is absent in season 2 (-110 +- 145 ppm), is carried by three events (161 +- 70 ppm without them) and is matched by a -265 +- 67 ppm brightening 3 h earlier in a phase scan.

## Key evidence

1. **Events 4, 36, 39 are ramps** (`redteam/stitched_vs_quarter.txt`, `redteam/quarter_diag.txt`): they lie 0.113, 0.061, 0.138 d after a Q2 gap, the Q9/Q10 roll and a Q10 downlink gap. At those cadences the stitched FLUX (searched) reads -43 / -14 / -47 ppt while per-quarter FLUX reads -32 / -4 / +29 ppt, PSF_FLUX_NOVA +17 / -28 / +16 ppt and SAP +7 / -7 / +8 ppt. RED_CHI2 of the PSF fit changes by -17000 and -6400; SAP centroid moves 0.02 px. The stitched PSF flux differs from the per-quarter PSF flux by 85-96 ppt MAD (max 180 ppt) in Q6, Q10, Q14 (all channel 17) vs 0.3-5 ppt in most other quarters.
2. **Ramps are generic** (`redteam/ramp_census.json`): 11 of 85 gaps give stitched-PSF ramps < -5 ppt (7 < -10 ppt); 8 of 11 on channel 17, 3 on channel 57, none on channels 65/25. The three largest (-38, -38, -25 ppt) have opposite-sign SAP ramps (+6.7, +7.2, +6.0 ppt).
3. **The period is the gap comb** (`redteam/residual_robustness.txt`): only 10.552 and 21.103 d (of 70 aliases in 1-50 d) phase all three big ramps within 0.15 d.
4. **Masked BLS** (`redteam/masked_bls.txt`): PSF all data 21.1069 d SDE 8.56 (rank 1); PSF with 1.0 d post-gap masked SDE 1.66 at 21.15 d (rank > 40, top 28.5 d at 5.86); inverted-flux masked top 23.16 d SDE 8.31; SAP masked 2.34 at 21.13 d.
5. **Residual vs time since gap** (`redteam/strict_threshold.txt`, 3 h box, own local-baseline depths): > 1 d: PSF 260 +- 63 / SAP 208 +- 60; > 2 d: PSF 201 +- 65 / SAP 97 +- 63; > 2.5-4 d: PSF 213-220 +- 65-67 / SAP 112-113 +- 63-64. Step at 2 d is epochs 5, 20, 67 (1.3-2.0 d after gaps with measured ramps).
6. **Other channels**: no LS power at 21.1 d or harmonics (PSF, SAP, PSF_FLAT_FLUX); odd/even and secondary null (meaningless on a ramp); 8 h box-SNR-preferred duration exceeds the 4.7-6.5 h b=0 maximum for a 0.93-1.3 Rsun host; pipeline's own gap_frac_1d 0.21, event_weighted_depth 318 vs 656 ppm, event_chi2_const 87/60 already pointed here.
7. **Known objects**: never a Kepler target; TPF host KIC 5609183 (Kp 13.74, 12.6 arcsec) has no KOI/TCE/EB; nothing within 180 arcsec except an unrelated 485 d TCE at 141 arcsec; no field-wide KOI/TCE/EB with matching period and epoch; live NEA and SIMBAD null; not in Bienias & Szabo 2025 or Moquin+2026. Novelty is not the issue.

## Disposition
Do not advance to pixel-level work or candidate status. If a pixel run is done anyway, it must exclude the twelve gap-edge epochs or it will localise the thermal ramp.