"""Stage-2 vetting of the inverted-flux (2 - f) first-screen passes: the null for the later stages.

    python scripts/vet_inverted.py [--workers 3]

Every star whose inverted light curve passed the first screen (runs/poolA_main and runs/poolA_pilot triage_all.csv)
is vetted with pipeline.vet_batch.vet_one(invert=True) at the inverted signal's period, duration and epoch, under the
current registered stage-2 rules. Writes runs/poolA_main/vet_inverted/<fname>/vet.json, vet_inverted_summary.csv and
vet_inverted_summary.json (counts compared with the real-flux stage-2 passes).
"""
import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.batch_search import BKJD_OFFSET, load_config  # noqa: E402
from pipeline.vet_batch import vet_one  # noqa: E402

OUT = "runs/poolA_main/vet_inverted"
CFG = "pipeline/prereg_poolA.json"


def one(row):
    cfg = load_config(CFG)
    r = vet_one(row["fname"], row["inv_period_d"], row["inv_duration_h"], row["inv_t0_bjd"] - BKJD_OFFSET, row["star"], cfg,
                os.path.join(OUT, row["fname"]), make_plots=False, invert=True)
    return {"fname": row["fname"], "run": row["run"], "P": r.get("P"), "depth_pass2_ppm": r.get("depth_pass2_ppm"),
            "snr_emp_pass2": r.get("snr_emp_pass2"), "n_events": r.get("n_events"), "ok": r.get("ok"), "error": r.get("error"),
            "stage2_pass": r.get("stage2_pass"), "stage2_fails": ",".join(r.get("stage2_fails", []))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    tg = pd.read_csv("targets/poolA_kbonus_bg_dwarfs_ranked.csv").set_index("fname")
    rows = []
    for run in ("poolA_main", "poolA_pilot"):
        t = pd.read_csv(f"runs/{run}/triage_all.csv")
        t = t[t.pass_inv.fillna(False).astype(bool)]
        for _, s in t.iterrows():
            star = tg.loc[s.fname].to_dict() if s.fname in tg.index else {}
            star = {k: (None if pd.isna(v) else v) for k, v in star.items() if isinstance(v, (int, float, str))}
            both = bool(s.get("pass_real")) and abs(float(s.real_period_d) / float(s.inv_period_d) - 1) < 0.01
            star["_both_polarity"] = both
            rows.append({"fname": s.fname, "run": run, "inv_period_d": float(s.inv_period_d), "inv_duration_h": float(s.inv_duration_h),
                         "inv_t0_bjd": float(s.inv_t0_bjd), "star": star})
    os.makedirs(OUT, exist_ok=True)
    print(f"{len(rows)} inverted first-screen passes", flush=True)
    res = []
    with ProcessPoolExecutor(a.workers) as ex:
        for r in ex.map(one, rows):
            res.append(r)
            print(f"{r['fname']}: {'PASS' if r['stage2_pass'] else 'fail ' + r['stage2_fails']}", flush=True)
    df = pd.DataFrame(res)
    df.to_csv("runs/poolA_main/vet_inverted_summary.csv", index=False)
    fails = df[~df.stage2_pass.astype(bool)].stage2_fails.str.split(",").explode().value_counts().to_dict()
    summ = {"inverted_first_screen": len(df), "inverted_stage2_pass": int(df.stage2_pass.astype(bool).sum()),
            "errors": int((~df.ok.astype(bool)).sum()), "fail_reasons": fails, "rules": "prereg stage2_vetting (current version)",
            "passing": df[df.stage2_pass.astype(bool)].fname.tolist()}
    json.dump(summ, open("runs/poolA_main/vet_inverted_summary.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
