#!/usr/bin/env python
"""
Cross-check of the kic-004474607 signal against the eclipsing binary that hosts the TPF.

The originating TPF star KIC 4474637 is KOI-3829.01 (FALSE POSITIVE, koi_fpflag_ss=1) and a
Kepler EB (Kirk+2016, J/AJ/151/68): P = 3.8867 d, depth ~37%, duration ~16 h, Kp 16.81, i.e. as
bright as our target and ~20 arcsec away. Catalogue values are read from data/known/*.csv
(downloaded 2026-10-17 per data/known/manifest.json). This script asks:
  (a) does the EB period (or its harmonics) leak into my detrended PSF flux of the target?
  (b) does each claimed 32.9 d event fall at an EB eclipse phase?
  (c) is 32.94 d a rational alias of 3.8867 d?
Reads only my own npz outputs and results.json (no pipeline code).
"""
import json, csv
import numpy as np
from astropy.timeseries import BoxLeastSquares

OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-004474607/rederive"
KNOWN = "/home/user/ZyushgPlanet/data/known"

# EB ephemeris from the KOI table (KIC 4474637 = K03829.01)
with open(f"{KNOWN}/koi_cumulative.csv") as fh:
    for row in csv.DictReader(fh):
        if row["kepid"] == "4474637":
            P_EB = float(row["koi_period"]); T0_EB = float(row["koi_time0bk"])
            DUR_EB_h = float(row["koi_duration"]); DEP_EB = float(row["koi_depth"])
            disp = row["koi_disposition"]; flags = {k: row[k] for k in row if k.startswith("koi_fpflag")}
            break

res = json.load(open(f"{OUT}/results.json"))
P = res["ephemeris"]["period_d"]; t0 = res["ephemeris"]["t0_bkjd"]; dur = res["ephemeris"]["duration_box_d"]
S = np.load(f"{OUT}/bls_search.npz"); periods = S["periods"]; pw = S["real_power"]
D = np.load(f"{OUT}/detrended_psf_masked.npz"); t, y, ye = D["t"], D["y"], D["ye"]

out = dict(eb=dict(kepid=4474637, koi="K03829.01", disposition=disp, fpflags=flags, period_d=P_EB,
                   t0_bkjd=T0_EB, duration_h=DUR_EB_h, depth_ppm=DEP_EB,
                   catalogue="data/known/koi_cumulative.csv fetched 2026-10-17"))

# (c) rational alias check
r = P / P_EB
best = min(((abs(r - n / d), n, d) for n in range(1, 40) for d in range(1, 8)), key=lambda x: x[0])
out["alias"] = dict(ratio=r, nearest_rational=f"{best[1]}/{best[2]}", frac_mismatch=best[0] / r,
                    drift_over_baseline_d=float(best[0] / r * P * (t.max() - t.min()) / P))

# (a) BLS power at EB period and harmonics within my own 1-50 d periodogram
med = np.median(pw); sd = np.std(pw)
leak = {}
for lab, Pl in [("P_EB", P_EB), ("2P_EB", 2 * P_EB), ("3P_EB", 3 * P_EB), ("4P_EB", 4 * P_EB), ("P_EB/2", P_EB / 2)]:
    if periods.min() <= Pl <= periods.max():
        m = np.abs(periods / Pl - 1) < 0.002
        leak[lab] = dict(period_d=Pl, max_power=float(pw[m].max()), sde_local=float((pw[m].max() - pw.mean()) / sd))
out["bls_power_at_eb_harmonics"] = leak
out["bls_power_at_best"] = dict(max_power=float(pw.max()), sde=float((pw.max() - pw.mean()) / sd))

# direct fold of my detrended target flux on the EB ephemeris: any eclipse-shaped residual?
ph = ((t - T0_EB) / P_EB + 0.5) % 1.0 - 0.5
inb = np.abs(ph * P_EB) < DUR_EB_h / 24 / 2
oot = np.abs(ph * P_EB) > DUR_EB_h / 24
d_eb = float(np.mean(y[oot]) - np.mean(y[inb])); e_eb = float(np.std(y[oot]) / np.sqrt(inb.sum()))
sec = np.abs((ph - 0.5 + 0.5) % 1.0 - 0.5) * P_EB < DUR_EB_h / 24 / 2
d_eb2 = float(np.mean(y[oot]) - np.mean(y[sec])); e_eb2 = float(np.std(y[oot]) / np.sqrt(sec.sum()))
out["target_flux_folded_on_eb"] = dict(primary_depth_ppm=d_eb * 1e6, primary_err_ppm=e_eb * 1e6,
                                       primary_sigma=d_eb / e_eb, secondary_depth_ppm=d_eb2 * 1e6,
                                       secondary_sigma=d_eb2 / e_eb2, n_in_primary=int(inb.sum()),
                                       note="fraction of the 37% EB eclipse that leaks into the target's PSF flux: "
                                            f"{d_eb / (DEP_EB * 1e-6) * 100:.3f}% of the EB depth")

# (b) EB phase of each covered 32.9 d event
events = []
for ev in res["transits"]["events"]:
    if not ev["covered"]:
        continue
    phe = ((ev["tc"] - T0_EB) / P_EB) % 1.0
    in_ecl = min(phe, 1 - phe) * P_EB * 24 < DUR_EB_h / 2 or abs(phe - 0.5) * P_EB * 24 < DUR_EB_h / 2
    events.append(dict(epoch=ev["epoch"], tc=ev["tc"], depth_ppm=ev["depth_ppm"], err_ppm=ev["err_ppm"],
                       eb_phase=phe, inside_eb_eclipse_window=bool(in_ecl)))
n_in = sum(e["inside_eb_eclipse_window"] for e in events)
frac_window = 2 * DUR_EB_h / 24 / P_EB   # fraction of EB cycle inside primary+secondary windows
out["events_vs_eb_phase"] = dict(events=events, n_events=len(events), n_inside_eb_eclipse=n_in,
                                 expected_if_random=len(events) * frac_window,
                                 depth_inside_ppm=float(np.mean([e["depth_ppm"] for e in events if e["inside_eb_eclipse_window"]])) if n_in else None,
                                 depth_outside_ppm=float(np.mean([e["depth_ppm"] for e in events if not e["inside_eb_eclipse_window"]])) if n_in < len(events) else None)

# does removing EB-eclipse-coincident cadences change the stacked depth?
ph_eb_all = ((t - T0_EB) / P_EB) % 1.0
in_eb_any = (np.minimum(ph_eb_all, 1 - ph_eb_all) * P_EB * 24 < DUR_EB_h / 2) | (np.abs(ph_eb_all - 0.5) * P_EB * 24 < DUR_EB_h / 2)
phs = ((t - t0) / P + 0.5) % 1.0 - 0.5
inbox = np.abs(phs * P) < dur / 2; oot2 = np.abs(phs * P) > 1.5 * dur
keep = ~in_eb_any
d_clean = float(np.mean(y[oot2 & keep]) - np.mean(y[inbox & keep]))
e_clean = float(np.std(y[oot2 & keep]) / np.sqrt((inbox & keep).sum()))
out["depth_excluding_eb_eclipse_cadences"] = dict(depth_ppm=d_clean * 1e6, err_analytic_ppm=e_clean * 1e6,
                                                   n_in=int((inbox & keep).sum()), n_in_full=int(inbox.sum()))

json.dump(out, open(f"{OUT}/eb_crosscheck.json", "w"), indent=2)
print(json.dumps(out, indent=1)[:5000])
