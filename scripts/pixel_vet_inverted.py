"""Pixel localisation (v1.6 per-quarter test) of the inverted-flux signals that pass stage 2: the null for that stage.

    python scripts/pixel_vet_inverted.py     # -> runs/poolA_main/pixel_v16_inverted/<fname>/pixel.json + summary json

Same test and rules as the real signals (pipeline.pixel_vet, as run by scripts/pixel_vet_survivors.py ... pixel_v16),
with the difference image mirrored (--invert) so that excess light during the inverted event plays the role of
missing light.
"""
import glob
import json
import os
import subprocess
import sys

OUT = "runs/poolA_main/pixel_v16_inverted"
done = {}
for vf in sorted(glob.glob("runs/poolA_main/vet_inverted/*/vet.json")):
    r = json.load(open(vf))
    if not r.get("stage2_pass"):
        continue
    fname, out = r["fname"], os.path.join(OUT, r["fname"])
    if not os.path.exists(os.path.join(out, "pixel.json")):
        print(f"pixel_vet --invert {fname} P={r['P']:.4f}", flush=True)
        cmd = [sys.executable, "-m", "pipeline.pixel_vet", "--fname", fname, "--period", str(r["P"]), "--t0-bjd", str(r["t0_bjd"]),
               "--dur-h", str(r["dur_h"]), "--ra", str(r["ra"]), "--dec", str(r["dec"]), "--tpforg", str(int(r["tpforg"])), "--out", out, "--invert"]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=2400)
        if p.returncode != 0:
            print("  ERROR:", p.stderr[-400:], flush=True)
            done[fname] = "error"
            continue
    done[fname] = json.load(open(os.path.join(out, "pixel.json")))["summary"]["verdict"]
    print(f"  {fname}: {done[fname]}", flush=True)
summ = {"inverted_stage2_pass": len(done), "verdicts": done, "on_target": sorted(k for k, v in done.items() if v == "on-target")}
json.dump(summ, open("runs/poolA_main/pixel_v16_inverted_summary.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
