# critic-referee notes (2026-10-17) — numbers re-derived from runs/poolA_main/*.csv and vet/pixel JSON

## Inverted null, per period bin (triage_all.csv, ok==True, pass_real / pass_inv by best-peak period)
P 0.5-2 d: real 24 / inv 25 ; P 2-10 d: 33 / 36 ; P 10-20 d: 12 / 4 ; P 20-40 d: 369 / 30
-> at P<10 d the screen has no discriminating power; the 12x excess at 20-40 d is the gap-systematics regime inversion cannot mimic.
Inverted best-peak medians: SDE 11.5 (>9), snr_emp 7.09 (~7.1), P 0.62 d, dur 7.05 h; 5841/6365 inverted vs 2902/6365 real best peaks fail duty cycle.
Pass SDE&SNR alone: real 1732, inv 145.

## Pilot (ranks 1-100): 10 real / 2 inverted, 0 stage-2 passes. poolA_main triage_all rank range 101-6465 (n=6365).
-> pool totals are 448 real / 97 inverted over 6,465, not 438/95.

## Prereg history (git, +0000): a5f8fa6 12:02 v1.0 (json says 12:05Z); c0b354d 12:12 v1.1 (12:20Z); a30c2fc 12:18 v1.2 (12:40Z);
90ae39a 12:44 stage2 with event_chi2_const_pvalue_min 0.001 (json 13:05Z); d04cfac 12:47 -> event_reduced_chi2_max 3.0 (json 13:20Z), not in revision_log. name field still v1.1.
Control KIC 5689121 stage-2: chi2 140.5 / 72 dof -> p ~ 3e-6, would have failed the p>=0.001 version.

## Stage-2 fail counts (438): tls_sde 349, event_chi2 347, events_positive 253, odd_even 136, depth_pass_ratio 136, ...
tls_sde is TLS SDE in a +/-2% period window (vet_batch._tls frac=0.02): SDE undefined there.
kic-001577049: snr_emp_pass2 1.88 (box_noise 703 ppm), tls_snr 18.5 -> passes (no stage-2 empirical SNR floor).

## Quartile depths chi2 (3 dof, formal errors): 009520080 31.6 (1957/1601/1091/1566); 011285493 14.7; 003963519 11.9; 001577049 16.2; 006367225 9.0; 006448637 8.8; 003831969 8.6
No criterion in prereg stage2_vetting; docs/10 L13 lists "per-quartile consistency" as a stage-2 test.

## Pixel
Target within 0.5 px of stamp edge or outside: 006699887 7/13 quarters; 009520080 12/17; 010213563 13/13.
Sum chi2 Gaussian at target vs at brightest in-stamp neighbour: 006699887 536 vs 1047 (G15.1); 009520080 839 vs 1568 (G14.8); 010213563 1143 vs 1556 (G14.3) — computed per quarter but NOT in summary or docs.
006448637: stamps 79x12 / 51x13 / 40x13; target_px y = -10.3, -3.7, -3.3 in Q3/Q6/Q7 (outside); 2-3 events/quarter; amp_gauss_at_target = -1.9e28 in Q3; median target->deficit 95" -> "not testable", not "off target".
Noise fallback e7dc7a8 14:52 added after control validation ad3cc1f 14:15 and first survivor runs a8236aa 14:18; no field in pixel.json records which branch was used.
Controls: runs/controls/pixel has only kic-005689121 (13/13). controls/vet_summary: kic-007603200 stage2_fails = known_object.
