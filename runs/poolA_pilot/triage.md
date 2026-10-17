# Triage: run `poolA_pilot`

- stars with results: 100 (ok 100, errors 0, skipped 0)
- first-screen passes: **real 10**, **inverted 2** (pre-registered thresholds: {"sde_min": 9.0, "snr_emp_min": 7.1, "n_transits_min": 3, "max_single_event_frac": 0.5, "odd_even_sigma_max": 3.0, "secondary_sigma_max": 3.0, "depth_ppm_max": 50000, "period_min_d": 0.5, "max_duty_cycle": 0.15})
- mean runtime 167 s/star; median point-to-point scatter 477 ppm

- sde: real p50/p90/p99/max = [10.3, 18.6, 85.9, 86.0]; inverted = [11.5, 30.5, 86.7, 96.4]
- snr_emp: real p50/p90/p99/max = [11.8, 27.6, 1141.6, 3565.6]; inverted = [5.9, 10.5, 14.8, 15.2]

## Real first-screen survivors (10)

| rank | fname | G | type | d (pc) | P (d) | depth (ppm) | dur (h) | SDE | SNR_emp | Ntr | o/e σ | sec σ | maxSE | inv SDE | both-pol | known |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 45 | gaia-dr3-2079295441538666496 | 13.8 | K | 529 | 37.6555 | 473±7 | 3.0 | 11.9 | 16.2 | 33 | 0.1 | 0.1 | 0.47 | 11.2 |  | none:  |
| 41 | kic-010001229 | 14.9 | K | 411 | 32.7565 | 525±33 | 1.5 | 13.5 | 15.2 | 40 | 2.2 | 0.5 | 0.40 | 14.6 |  | none:  |
| 51 | gaia-dr3-2128830403176666880 | 13.7 | G | 584 | 34.4885 | 233±7 | 3.0 | 11.2 | 14.3 | 37 | 0.7 | -0.6 | 0.50 | 9.7 |  | none:  |
| 9 | kic-006192017 | 16.5 | mid-late-M | 250 | 32.1054 | 726±62 | 7.1 | 10.5 | 11.3 | 26 | 2.9 | 0.1 | 0.38 | 16.6 |  | none:  |
| 70 | kic-007660076 | 16.1 | early-M | 441 | 31.1973 | 1073±56 | 1.5 | 10.5 | 10.8 | 38 | 0.5 | 0.8 | 0.45 | 9.3 |  | none:  |
| 55 | kic-008903463 | 14.9 | K | 436 | 38.2427 | 518±27 | 1.5 | 11.4 | 10.6 | 29 | 1.5 | -0.8 | 0.50 | 10.9 |  | none:  |
| 60 | gaia-dr3-2077399234963970560 | 13.8 | G | 411 | 27.5015 | 76±6 | 7.1 | 10.4 | 10.4 | 49 | 0.6 | 0.2 | 0.32 | 11.3 |  | none:  |
| 40 | kic-010263791 | 15.4 | early-M | 533 | 29.6989 | 256±30 | 3.0 | 10.1 | 10.0 | 40 | 0.8 | -0.4 | 0.42 | 11.1 |  | none:  |
| 100 | kic-009835496 | 14.8 | early-M | 474 | 23.9351 | 200±18 | 7.1 | 10.9 | 7.7 | 26 | 0.3 | 1.7 | 0.36 | 8.0 |  | none:  |
| 22 | kic-008004947 | 17.2 | mid-late-M | 310 | 31.1750 | 1332±83 | 1.5 | 10.7 | 7.2 | 39 | 1.7 | 0.2 | 0.48 | 10.5 |  | NEARBY: TCE:0 1:1 sep=2.90' TCE dphase=0.493 |

## Inverted (null) passes (2) — the false-alarm floor

| fname | P (d) | depth (ppm) | SDE | SNR_emp | Ntr |
|---|---|---|---|---|---|
| kic-007703185 | 0.5149 | 310 | 34.4 | 15.2 | 2443 |
| kic-008903463 | 27.7897 | 177 | 10.9 | 10.0 | 50 |
