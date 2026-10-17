"""
Modal fan-out for the KBONUS background-star transit screen.

Needs MODAL_TOKEN_ID / MODAL_TOKEN_SECRET in the environment (set them in the cloud
environment settings, not in chat).  Then:

    modal run pipeline/modal_app.py --start 1 --stop 6465 --run-id poolA_full

Each container downloads one star's FITS from MAST (~9 MB), runs the identical
search_star() as the local runner (same prereg config, same thresholds), and
returns the dict.  Results are written locally to runs/<run-id>/results.jsonl.
Cost estimate: ~90 CPU-s/star -> 6,465 stars ~ 160 CPU-h ~ US$30-40 on Modal at
~US$0.19/CPU-h; wall-clock ~1-2 h with 100 concurrent containers.
"""
from __future__ import annotations

import json
import os

import modal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy", "scipy", "astropy", "wotan", "requests", "pandas")
    .add_local_dir(os.path.join(ROOT, "pipeline"), remote_path="/root/pipeline")
    .add_local_file(os.path.join(ROOT, "pipeline", "prereg_poolA.json"), remote_path="/root/prereg.json")
)

app = modal.App("zyushg-planet-screen", image=image)


@app.function(cpu=2.0, memory=2048, timeout=1800, max_containers=100, retries=1)
def screen_one(fname: str) -> dict:
    import sys
    sys.path.insert(0, "/root")
    from pipeline.batch_search import load_config, passes_first_screen, search_star
    cfg = load_config("/root/prereg.json")
    res = search_star(fname, cfg, cache_dir="/tmp/kb", keep_fits=False)
    if res.get("ok") and "real" in res:
        res["pass_real_s1"] = bool(passes_first_screen(res["real"]["s1"], cfg)[0])
        res["pass_inv_s1"] = bool(passes_first_screen(res["inverted"]["s1"], cfg)[0]) if "inverted" in res else False
    return res


@app.local_entrypoint()
def main(start: int = 1, stop: int = 100, run_id: str = "poolA_modal",
         targets: str = "targets/poolA_kbonus_bg_dwarfs_ranked.csv"):
    import pandas as pd
    df = pd.read_csv(os.path.join(ROOT, targets))
    df = df[(df["rank"] >= start) & (df["rank"] <= stop)]
    run_dir = os.path.join(ROOT, "runs", run_id)
    os.makedirs(run_dir, exist_ok=True)
    out = os.path.join(run_dir, "results.jsonl")
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["fname"] for l in open(out)}
    todo = [f for f in df["fname"].tolist() if f not in done]
    print(f"{len(todo)} stars to screen on Modal (run {run_id})")
    n = npr = npi = 0
    with open(out, "a") as fh:
        for res in screen_one.map(todo, order_outputs=False, return_exceptions=True):
            if isinstance(res, Exception):
                res = {"fname": "?", "ok": False, "error": repr(res)}
            fh.write(json.dumps(res, default=str) + "\n")
            fh.flush()
            n += 1
            npr += bool(res.get("pass_real_s1"))
            npi += bool(res.get("pass_inv_s1"))
            if n % 25 == 0:
                print(f"  {n}/{len(todo)} done; first-screen pass real={npr} inverted={npi}")
    print(f"done: {n} stars; real passes {npr}, inverted passes {npi} -> {out}")
