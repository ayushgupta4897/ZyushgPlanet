## Red team: kic-009655224 — verdict **survives** (signal, not a candidate)

Notes: `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-009655224/redteam.md`; scripts, logs, JSON and `MANIFEST.sha256` in `/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-009655224/redteam/`. No kill found; no test exceeded ~2.6 sigma.

**Kill attempts that failed**
- *Known/contaminating object*: no KOI, DR25 TCE or Kepler EB row for KIC 9655224/9655236/9655226; no ephemeris match at any harmonic field-wide (nearest: KIC 8757619 P 6.70842 at 89', epoch off 8.6 h; EB KIC 9719636 P 3.35157 at 52', off 4.5 h and no phase-0.5 ghost); SIMBAD empty at target; Gaia DR3 target clean (RUWE 1.04, no multi-peak/NSS), nothing within 5".
- *Pixel origin (own code, `redteam/pixel_fit.log`)*: 3-source PRF fit on 10 host TPFs: target 1.05 +/- 0.15 e/s, G=15.9 neighbour 0.39 +/- 0.18, host -0.06 +/- 0.14; target-only beats neighbour-only by dchi2 = 46. Whole-stamp deficit 1.32 +/- 0.34 e/s vs 1.07 lost in the PSF LC and 1.04 in SAP. Pipeline's own pixel.json summed over 10 quarters agrees (target chi2 299.5 vs 328.6 neighbour, 340.9 host). Its weak centroid statistic is edge-biased; `deficit_total_e_s` is a positive-clipped sum, not light lost.
- *Neighbours' own LCs at the ephemeris, all quarters*: G=15.9 star 0.43 +/- 0.23 e/s (needs >= 2 e/s to be the source), host -0.15 +/- 0.12 e/s.
- *EB tests*: no secondary at any phase (max 1.7 sigma), odd/even 0.3 sigma, flat bottom over V by dchi2 9, duration consistent with a 0.83 Rsun dwarf at b ~ 0.65, Rp ~ 2.25 Re.
- *Systematics*: no gap/trend/flag correlation (desats z -1.1; max flag z 2.6 of 11), present in all four channels, not carried by outlier events (429 +/- 75 ppm without the top 5; median depth = mean depth).
- *Variability*: nothing at P or harmonics in Lomb-Scargle.

**Residual weaknesses (carry forward)**: season-1 (ch 23) depth 368 +/- 112 vs ~700 ppm (p ~ 0.1); 0.4 +/- 0.2 e/s echo in the G=15.9 neighbour consistent with leakage but unproven; per-event errors ~20% underestimated; blend within ~4" not excludable by Kepler pixels; host G 17.1 at 2.2 kpc, so validation can only be statistical.

**To settle**: Kepler-PRF (not Gaussian) difference-image fit per season; TRICERATOPS FPP/NFPP with the Gaia neighbours; PRF-based per-season re-extraction; KIC lookup in Bienias & Szabo 2025 / Moquin 2026 lists.