"""Recompute stage-2 verdicts for every existing vet.json under the CURRENT prereg (no re-run), rewrite vet_summary.csv,
and write survivors_stage2.csv. Keeps the previously stored verdict as stage2_pass_prev / stage2_fails_prev."""
import glob, json, os, sys
sys.path.insert(0, os.getcwd())
import pandas as pd
from pipeline.batch_search import load_config
from pipeline.vet_batch import stage2_verdict
cfg = load_config(sys.argv[2] if len(sys.argv) > 2 else "pipeline/prereg_poolA.json")
run = sys.argv[1]
rows = []
sv = pd.read_csv(f"runs/{run}/survivors.csv").set_index("fname") if os.path.exists(f"runs/{run}/survivors.csv") else None
for vf in sorted(glob.glob(f"runs/{run}/vet/*/vet.json")):
    r = json.load(open(vf))
    bp = bool(sv.loc[r["fname"], "both_polarity"]) if sv is not None and r["fname"] in sv.index and "both_polarity" in sv.columns else False
    if "stage2_pass_v14" not in r:
        r["stage2_pass_v14"], r["stage2_fails_v14"] = r.get("stage2_pass"), r.get("stage2_fails")
    ok, fails = stage2_verdict(r, cfg, both_polarity=bp)
    r["stage2_pass"], r["stage2_fails"], r["stage2_rules"] = ok, fails, cfg["name"]
    json.dump(r, open(vf, "w"), indent=1, default=str)
    flat = {k: v for k, v in r.items() if not isinstance(v, (list, dict))}
    if isinstance(r.get("tls"), dict):
        flat.update({f"tls_{k}": v for k, v in r["tls"].items()})
    flat["stage2_fails"] = ",".join(fails); flat["stage2_fails_v14"] = ",".join(r.get("stage2_fails_v14") or [])
    rows.append(flat)
df = pd.DataFrame(rows); df.to_csv(f"runs/{run}/vet_summary.csv", index=False)
surv = df[df.stage2_pass == True].sort_values("snr_emp_pass2", ascending=False)
surv.to_csv(f"runs/{run}/survivors_stage2.csv", index=False)
print(f"{run}: vetted {len(df)}; pass under {cfg['name']}: {len(surv)}; pass under v1.4 rules: {int(df.stage2_pass_v14.fillna(False).sum())}")
print(surv[["fname", "P", "depth_pass2_ppm", "snr_emp_pass2", "secondary_sigma", "gap_frac_1d", "duration_ratio", "stage2_fails_v14"]].to_string(index=False))
