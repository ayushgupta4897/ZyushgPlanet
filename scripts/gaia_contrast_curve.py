"""Empirical Gaia DR3 contrast sensitivity around a target, usable as a TRICERATOPS contrast curve.

    python scripts/gaia_contrast_curve.py <fname> <ra> <dec> [--radius-deg 1.0] [--out runs/poolA_main/gaia_contrast/<fname>]

Method (after Brandeker & Cataldi 2019, A&A 621, A86, who did this for DR2 over the whole sky): take every Gaia DR3
source within `radius` of the target, find all pairs closer than 6 arcsec, and bin them in (separation, delta G).
At separations of 4-6 arcsec Gaia resolves everything, so the pair density there (per unit sky area, per delta-G bin)
is the reference. The completeness at a smaller separation is the ratio of the observed pair density to that
reference; the contrast curve is, for each separation, the largest delta G at which completeness is still >= 0.5.
This is the magnitude difference down to which Gaia would have *catalogued* a companion at that separation; a
companion fainter than the curve is not excluded.

Caveats: this is a statistical property of the catalogue in this field, not an observation of the target; it says
nothing inside ~0.4 arcsec; it is in Gaia G (passed to TRICERATOPS as filt="Vis"). Output: contrast_curve.csv
(sep_arcsec, delta_mag), completeness.json, and a PNG.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("fname")
ap.add_argument("ra", type=float)
ap.add_argument("dec", type=float)
ap.add_argument("--radius-deg", type=float, default=1.0)
ap.add_argument("--gmax", type=float, default=21.0)
ap.add_argument("--out", default=None)
a = ap.parse_args()
out = a.out or f"runs/poolA_main/gaia_contrast/{a.fname}"
os.makedirs(out, exist_ok=True)

from astroquery.gaia import Gaia  # noqa: E402

Gaia.ROW_LIMIT = -1
q = f"""SELECT source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {a.ra}, {a.dec}, {a.radius_deg}))
AND phot_g_mean_mag < {a.gmax}"""
job = Gaia.launch_job_async(q)
t = job.get_results().to_pandas()
t.columns = [c.lower() for c in t.columns]
print(f"{len(t)} Gaia DR3 sources with G < {a.gmax} within {a.radius_deg} deg")

# pairs within 6 arcsec via a KD-tree on tangent-plane coordinates
from scipy.spatial import cKDTree  # noqa: E402

cosd = np.cos(np.radians(a.dec))
x = (t.ra.values - a.ra) * cosd * 3600.0
y = (t.dec.values - a.dec) * 3600.0
tree = cKDTree(np.c_[x, y])
pairs = tree.query_pairs(6.0, output_type="ndarray")
i, j = pairs[:, 0], pairs[:, 1]
sep = np.hypot(x[i] - x[j], y[i] - y[j])
g = t.phot_g_mean_mag.values
dg = np.abs(g[i] - g[j])
bright = np.minimum(g[i], g[j])
# restrict primaries to roughly the target's brightness class so the faint-companion limit is representative
keep = (bright > 13) & (bright < 18.5)
sep, dg = sep[keep], dg[keep]
print(f"{len(sep)} pairs within 6 arcsec with a 13 < G < 18.5 primary")

sep_edges = np.arange(0.0, 6.01, 0.25)
dg_edges = np.arange(0.0, 8.01, 0.5)
H, _, _ = np.histogram2d(sep, dg, bins=[sep_edges, dg_edges])
area = np.pi * (sep_edges[1:] ** 2 - sep_edges[:-1] ** 2)  # annulus area per bin
dens = H / area[:, None]
ref = dens[(sep_edges[:-1] >= 4.0)].mean(axis=0)  # reference density per delta-G bin at 4-6 arcsec
with np.errstate(divide="ignore", invalid="ignore"):
    comp = np.where(ref > 0, dens / ref, np.nan)
comp = np.clip(comp, 0, 1.5)
curve = []
for k in range(len(sep_edges) - 1):
    smid = 0.5 * (sep_edges[k] + sep_edges[k + 1])
    if smid > 4.0:
        break
    ok = np.where(comp[k] >= 0.5)[0]
    # contiguous from delta G = 0: the largest bin index such that all bins up to it are >= 0.5 completeness
    lim = 0.0
    for b in range(len(dg_edges) - 1):
        if np.isfinite(comp[k, b]) and comp[k, b] >= 0.5:
            lim = dg_edges[b + 1]
        else:
            break
    curve.append((smid, lim, int(H[k].sum())))
cc = pd.DataFrame(curve, columns=["sep_arcsec", "delta_mag", "n_pairs_in_bin"])
# enforce monotonic non-decreasing contrast with separation (physical for a resolution-limited catalogue)
cc["delta_mag"] = np.maximum.accumulate(cc.delta_mag.values)
cc = cc[cc.delta_mag > 0]
cc[["sep_arcsec", "delta_mag"]].to_csv(os.path.join(out, "contrast_curve.csv"), index=False, header=False, sep=",")  # TRICERATOPS reads comma-delimited
cc.to_csv(os.path.join(out, "contrast_curve_with_counts.csv"), index=False)
json.dump({"fname": a.fname, "ra": a.ra, "dec": a.dec, "radius_deg": a.radius_deg, "gmax": a.gmax, "n_sources": int(len(t)),
           "n_pairs": int(len(sep)), "primary_G_range": [13, 18.5], "sep_edges": sep_edges.tolist(), "dg_edges": dg_edges.tolist(),
           "completeness": np.nan_to_num(comp, nan=-1).tolist(), "reference_density_per_arcsec2": ref.tolist(),
           "curve": cc.to_dict(orient="records"), "method": "Brandeker & Cataldi 2019 style: resolved-pair density vs 4-6 arcsec reference, 50 % completeness",
           "band": "Gaia G (TRICERATOPS filt='Vis')"}, open(os.path.join(out, "completeness.json"), "w"), indent=1)
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    im = ax[0].imshow(comp.T, origin="lower", aspect="auto", extent=[sep_edges[0], sep_edges[-1], dg_edges[0], dg_edges[-1]], vmin=0, vmax=1.2, cmap="viridis")
    ax[0].set_xlabel("separation (arcsec)"); ax[0].set_ylabel("delta G"); ax[0].set_title(f"Gaia DR3 pair completeness, {a.fname} field"); plt.colorbar(im, ax=ax[0])
    ax[1].step(cc.sep_arcsec, cc.delta_mag, where="mid"); ax[1].set_xlabel("separation (arcsec)"); ax[1].set_ylabel("delta G at 50 % completeness"); ax[1].set_title("contrast curve"); ax[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(out, "contrast_curve.png"), dpi=120)
except Exception as ex:  # noqa: BLE001
    print("plot failed:", ex)
print(cc.to_string(index=False))
