## kic-003963519 red team (2026-10-17): verdict **survives**

Own code only (`runs/poolA_main/adversarial/kic-003963519/redteam/{lc_checks.py,pixel_checks.py}` -> `lc_checks.json`, `pixel_checks.json`, `ephem_matches.log`, `gaia_45arcsec.csv`); full notes in `runs/poolA_main/adversarial/kic-003963519/redteam.md`. Candidate language only.

**Blend / background EB.** Independent pixel localisation (season-stacked event-differenced deficit maps on the 14 KIC 3963530 TPFs, host-PRF template, grid incl. 3 px beyond the stamp): target dchi2 vs null **70.1** over 4 seasons; best alternatives: TPF host (10.4") 11.4, G 20.3 star (8.8") 11.3, nearest neighbour G 18.67 at 7.0" **2.0**. Margin 59 (pipeline's own per-quarter fits: 140 vs 57, margin 84); v16 rule met. Target inside the 1-sigma position region in 3/4 seasons, 3-sigma in 4/4. Pixel-implied depth on the target **1053 +- 135 ppm** = light-curve depth. Shape flat-bottomed (T23/T14 0.97, V worse by dchi2 54); odd/even 0.24 sigma; secondary 59 ppm (err 124); strongest dip at any phase 3.0 sigma post-hoc. SAP/PSF depth 0.85 +- 0.14 (PSF 1194, SAP 1009, PSF_FLAT 1200, NOVA 1201 ppm), as expected for on-target light. Caveat: the target sits on the stamp edge (60-71% of PRF on silicon); localisation across the edge is weaker; Gaia shows no off-stamp source within 3 px along the edge normal; an unresolved companion within 1-2" is untested.

**Variability.** LS on detrended flux: 0.62027 d, 100 ppm half-amplitude (the periodicity the inverted null nearly passed on); P/P_LS = 10.99, cannot build the box; depth after notching 1177 +- 146 ppm. Raw LS top peak 52.5 d (systematics); power at P 9e-6.

**Systematics.** Per-event chi2/dof 1.56 vs random-ephemeris calibration median 1.50 (p84 1.69; 41% of random ephemerides scatter more). Far from gaps 1296 +- 101 ppm (136 events) vs near gaps 793 +- 370 (12). Without the 10 strongest events 1025 +- 102. Quarter thirds 1432/1232/1131 +- 170. Per-quarter chi2 10.1/13; per-channel 1428 +- 462 (ch8), 1003 +- 274 (ch36), 1393 +- 207 (ch80), 1171 +- 203 (ch52). Per-year T0 shifts within +-12 min. RED_CHI2 in/out +0.05%.

**Ephemeris contamination.** Channels 52/8/36/80 by season. No KOI/TCE/EB/confirmed object anywhere in the field with |P/(mP)-1| < 5e-4 (m = 1/4..4) and aligned epoch; nearest period KOI-1765.01 = EB KIC 8043714 (6.814989 d: 20 sigma off, epoch 3.2 d off, 314' away). Nothing within 3'. Catalogues fetched 2026-10-17T15:46:52Z.

**Known?** Not KOI/TCE/EB/PS; never a Kepler target; SIMBAD 60" empty; outside Bienias & Szabo 2025 scope (HJs 2.5-4.8 d); not Moquin 2026 BOI-1.01.

**Host / duration.** Gaia DR3 plx 0.229 +- 0.061 (SNR 3.8), RUWE 1.02, GSP-Phot 5110 K, logg 4.69, R 0.60 (0.55-0.74) Rsun; T14_max 2.5-3.0 h vs observed 2.37-2.44 h (b 0.2-0.6); implied radius 2.1-2.8 Re (up to ~7 Re if the star is evolved).

**Residual risks:** unresolved companion/hierarchical EB (needs TRICERATOPS + imaging), stamp-edge truncation, modest SNR (8.5-9.6), host radius uncertain, year-3/ch52 depth 1.5 sigma low.