"""Run pixel-level localisation on every stage-2 survivor of a run (skips ones already done)."""
import glob, json, os, subprocess, sys
run = sys.argv[1] if len(sys.argv) > 1 else "poolA_main"
outroot = sys.argv[2] if len(sys.argv) > 2 else "pixel"
for vf in sorted(glob.glob(f"runs/{run}/vet/*/vet.json")):
    r = json.load(open(vf))
    if not (r.get("stage2_pass") or r.get("stage2_pass_v14")):
        continue
    fname = r["fname"]
    out = f"runs/{run}/{outroot}/{fname}"
    if os.path.exists(os.path.join(out, "pixel.json")):
        continue
    print(f"pixel_vet {fname} P={r['P']:.4f}", flush=True)
    cmd = [sys.executable, "-m", "pipeline.pixel_vet", "--fname", fname, "--period", str(r["P"]), "--t0-bjd", str(r["t0_bjd"]),
           "--dur-h", str(r["dur_h"]), "--ra", str(r["ra"]), "--dec", str(r["dec"]), "--tpforg", str(int(r["tpforg"])), "--out", out]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=1500)
        print("\n".join(l for l in p.stdout.splitlines() if l.startswith((" ", "{", "}", '"'))), flush=True)
        if p.returncode != 0:
            print("  ERROR:", p.stderr[-400:], flush=True)
    except subprocess.TimeoutExpired:
        print("  TIMEOUT", flush=True)
