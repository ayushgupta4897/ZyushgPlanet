"""Field-wide (whole Kepler field, any channel) period + epoch match of a signal against KOI, DR25 TCE and the
Kepler EB catalogue, with harmonics 1:1, 2:1, 1:2, 3:1, 1:3, 4:1, 1:4. Motivated by kic-011285493, where a 12.5 %-deep
EB 5.6 arcmin away (outside the registered 3-arcmin radius) was imprinted channel-wide on KBONUS PSF photometry.

    python scripts/ephemeris_fieldwide.py poolA_main kic-006699887

A row is a *coincidence* when |dP/P| < 1e-3; it is flagged `epoch_match` when the signal's mid-transit falls within
two transit durations of the catalogue object's eclipse at the catalogue period. Output: runs/<run>/ephemeris_fieldwide/<fname>.json
"""
import json, os, sys
import numpy as np, pandas as pd
BK = 2454833.0
run, fname = sys.argv[1], sys.argv[2]
v = json.load(open(f"runs/{run}/vet/{fname}/vet.json"))
P, t0, D = float(v["P"]), float(v["t0_bjd"]) - BK, float(v["dur_h"]) / 24
man = json.load(open("data/known/manifest.json"))
cats = [("KOI", pd.read_csv("data/known/koi_cumulative.csv"), "koi_period", "koi_time0bk", "kepoi_name", "koi_depth"),
        ("DR25_TCE", pd.read_csv("data/known/dr25_tce.csv"), "tce_period", "tce_time0bk", "kepid", "tce_depth"),
        ("Kepler_EB_Kirk2016", pd.read_csv("data/known/J_AJ_151_68_catalog.csv"), "Per", "BJD0", "KIC", None)]
rows = []
for name, df, pc, tc, idc, dc in cats:
    per = pd.to_numeric(df[pc], errors="coerce").values
    for ratio, lab in [(1, "1:1"), (2, "2:1"), (0.5, "1:2"), (3, "3:1"), (1 / 3, "1:3"), (4, "4:1"), (0.25, "1:4")]:
        frac = np.abs(per / (P * ratio) - 1)
        for i in np.where(frac < 1e-3)[0]:
            r = df.iloc[i]
            their_t0 = pd.to_numeric(r[tc], errors="coerce")
            off = None
            if np.isfinite(their_t0):
                their_t0 = float(their_t0) - (BK if their_t0 > 2400000 else 0)
                off = float(((t0 - their_t0 + per[i] / 2) % per[i]) - per[i] / 2)
            rows.append({"catalogue": name, "id": str(r[idc]), "ratio": lab, "period_d": float(per[i]), "dP_over_P": float(frac[i]),
                         "t0_offset_d": off, "epoch_match": bool(off is not None and abs(off) < 2 * D),
                         "depth_ppm": None if dc is None else float(pd.to_numeric(r[dc], errors="coerce")),
                         "ra": None if "ra" not in df.columns else float(r["ra"]), "dec": None if "dec" not in df.columns else float(r["dec"])})
# physical viability: a 1:1 twin contaminates directly; a 2:1 twin (their P = 2 x ours) would alternate odd/even depths;
# a 1:2 twin (their P = ours / 2) would put an equal dip at phase 0.5; 3:1 / 1:3 / 4:1 / 1:4 likewise leave phase-scan signatures.
oe = abs(v.get("odd_even_sigma") or 0); sec = abs(v.get("secondary_sigma") or 0)
ra0, dec0 = float(v["ra"]), float(v["dec"])
for r in rows:
    if r["ra"] is not None:
        r["sep_arcmin"] = 60 * float(np.hypot((r["ra"] - ra0) * np.cos(np.radians(dec0)), r["dec"] - dec0))
    viable = r["epoch_match"] and (r["ratio"] == "1:1" or (r["ratio"] in ("2:1", "3:1", "4:1") and oe > 3) or (r["ratio"] in ("1:2", "1:3", "1:4") and sec > 3))
    r["viable_contaminant"] = bool(viable)
out = {"fname": fname, "odd_even_sigma": v.get("odd_even_sigma"), "secondary_sigma": v.get("secondary_sigma"),
       "n_viable_contaminants": int(sum(r["viable_contaminant"] for r in rows)), "viable_contaminants": [r for r in rows if r["viable_contaminant"]], "P_d": P, "t0_bkjd": t0, "dur_d": D, "catalogues_fetched_utc": man.get("fetched_utc"),
       "n_period_coincidences": len(rows), "n_epoch_matches": int(sum(r["epoch_match"] for r in rows)),
       "epoch_matches": [r for r in rows if r["epoch_match"]], "all_period_coincidences": rows,
       "note": "Whole Kepler field, any channel. A period coincidence without an epoch match cannot contaminate a folded light curve over a 4-year baseline."}
os.makedirs(f"runs/{run}/ephemeris_fieldwide", exist_ok=True)
json.dump(out, open(f"runs/{run}/ephemeris_fieldwide/{fname}.json", "w"), indent=1)
print(f"{fname}: {len(rows)} period coincidences, {out['n_epoch_matches']} with matching epoch, {out['n_viable_contaminants']} viable given odd/even {oe:.1f} sigma and secondary {sec:.1f} sigma:", [(r['catalogue'], r['id'], r['ratio'], round(r['t0_offset_d'], 3), round(r.get('sep_arcmin') or -1)) for r in out['viable_contaminants']])
print("   epoch matches by harmonic:", [(r['ratio'], r['id'], round(r.get('sep_arcmin') or -1, 1)) for r in out['epoch_matches']])
