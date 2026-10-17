#!/usr/bin/env python
"""
Supplement to rederive_kic011285493.py (same author, same loader, nothing from pipeline/).

Purpose: the main re-derivation measured the box depth with a box equal to the
trapezoid T14 (5.9 h), which includes ingress/egress and therefore reads lower
than the pipeline's 4.5 h box.  To compare like with like, measure the box depth
(with the same random-phase empirical error) as a function of box width, on both
PSF FLUX and SAP_FLUX, at the re-derived ephemeris.  Also convert each depth to an
absolute deficit in e/s so the two apertures can be compared independent of dilution.

Writes: depth_vs_box_width.csv and merges a 'depth_vs_box_width' + 'comparison_with_pipeline'
block into results.json.
"""
import json
import os
import sys
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rederive_kic011285493 as R   # own code

res = json.load(open(os.path.join(HERE, "results.json")))
P = res["ephemeris"]["period_d"]
t0 = res["ephemeris"]["t0_bkjd"]
dur = res["ephemeris"]["duration_bls_h"] / 24.0      # mask half-width used in the main run
T14 = res["ephemeris"]["duration_trapezoid_T14_h"]
T23 = res["ephemeris"]["duration_trapezoid_T23_h"]

widths_h = sorted(set([round(T23, 2), 3.0, round(dur * 24, 2), 4.0, 4.5, 5.0, round(T14, 2)]))

with fits.open(R.FITS) as h:
    d = h["LIGHTCURVE_STITCHED"].data
    sq = np.array(d["SAP_QUALITY"], int)
    good = sq == 0
    med_psf = float(np.nanmedian(np.array(d["FLUX"], float)[good]))
    med_sap = float(np.nanmedian(np.array(d["SAP_FLUX"], float)[good]))

table = []
for col, med in (("FLUX", med_psf), ("SAP_FLUX", med_sap)):
    t, f_raw, e_raw, q = R.load(col)
    tmask = np.abs(R.phase_days(t, P, t0)) < 1.0 * dur
    f2, keep, mad = R.detrend(t, f_raw, window=1.0, mask=tmask)
    tm, fm = t[keep], f2[keep]
    for w in widths_h:
        wd = w / 24.0
        dep, n_in = R.box_depth(tm, fm, P, t0, wd, exclude=[(P / 2, wd)])
        err, _ = R.empirical_depth_error(tm, fm, P, t0, wd, ntrial=1500, avoid=[(0.0, 1.5 * wd)])
        row = dict(column=col, box_h=w, depth_ppm=dep * 1e6, err_ppm=err * 1e6,
                   snr=dep / err, n_in=n_in, median_flux_e_s=med,
                   deficit_e_s=dep * med, deficit_err_e_s=err * med)
        table.append(row)
        print(f"{col:9s} box={w:4.2f} h  depth={dep*1e6:6.0f} +- {err*1e6:4.0f} ppm  "
              f"SNR={dep/err:4.1f}  n_in={n_in:4d}  deficit={dep*med:5.2f} +- {err*med:4.2f} e/s", flush=True)

import csv
with open(os.path.join(HERE, "depth_vs_box_width.csv"), "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(table[0].keys()))
    wr.writeheader(); wr.writerows(table)

def pick(col, w):
    return [r for r in table if r["column"] == col and abs(r["box_h"] - w) < 0.01][0]

psf45, sap45 = pick("FLUX", 4.5), pick("SAP_FLUX", 4.5)
ratio45 = sap45["depth_ppm"] / psf45["depth_ppm"]
ratio45_err = ratio45 * np.hypot(sap45["err_ppm"] / sap45["depth_ppm"], psf45["err_ppm"] / psf45["depth_ppm"])
def_ratio = sap45["deficit_e_s"] / psf45["deficit_e_s"]
def_ratio_err = def_ratio * np.hypot(sap45["deficit_err_e_s"] / sap45["deficit_e_s"],
                                     psf45["deficit_err_e_s"] / psf45["deficit_e_s"])

res["depth_vs_box_width"] = {
    "note": "box depth (mean out - mean in, out = beyond 1 box width either side, secondary window excluded) "
            "vs box width at the re-derived ephemeris; error = robust std of 1500 random-phase boxes",
    "rows": table,
    "psf_4p5h_box": psf45, "sap_4p5h_box": sap45,
    "sap_over_psf_depth_ratio_4p5h": ratio45, "sap_over_psf_depth_ratio_4p5h_err": ratio45_err,
    "psf_over_sap_median_flux": med_psf / med_sap,
    "sap_over_psf_absolute_deficit_ratio_4p5h": def_ratio,
    "sap_over_psf_absolute_deficit_ratio_4p5h_err": def_ratio_err,
}
json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=float)
print("merged into results.json")
