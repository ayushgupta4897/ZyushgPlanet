"""
Local batch runner for the KBONUS background-star transit screen.

    . .venv/bin/activate
    python -m pipeline.run_batch --config pipeline/prereg_poolA.json \
        --targets targets/poolA_kbonus_bg_dwarfs_ranked.csv --start 1 --stop 200 \
        --run-id poolA_pilot --workers 3

Writes runs/<run-id>/results.jsonl (one line per star, append-only; safe to resume)
and runs/<run-id>/summary.json.  Also accepts --fnames for ad-hoc lists (controls).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.batch_search import load_config, passes_first_screen, search_star  # noqa: E402


def _worker(args):
    fname, cfg, cache_dir, keep = args
    return search_star(fname, cfg, cache_dir=cache_dir, keep_fits=keep)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--targets", default=None)
    ap.add_argument("--fnames", nargs="*", default=None)
    ap.add_argument("--start", type=int, default=1, help="first rank (1-based, inclusive)")
    ap.add_argument("--stop", type=int, default=None, help="last rank (inclusive)")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--keep-fits", action="store_true")
    ap.add_argument("--cache-dir", default="data/kbonus_batch")
    a = ap.parse_args()

    cfg = load_config(a.config)
    run_dir = os.path.join("runs", a.run_id)
    os.makedirs(run_dir, exist_ok=True)
    out_path = os.path.join(run_dir, "results.jsonl")
    with open(os.path.join(run_dir, "config_used.json"), "w") as fh:
        json.dump(cfg, fh, indent=1)

    if a.fnames:
        fnames = list(a.fnames)
    else:
        df = pd.read_csv(a.targets)
        df = df[df["rank"] >= a.start]
        if a.stop is not None:
            df = df[df["rank"] <= a.stop]
        fnames = df["fname"].tolist()
    done = set()
    if os.path.exists(out_path):
        with open(out_path) as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["fname"])
                except Exception:  # noqa: BLE001
                    pass
    todo = [f for f in fnames if f not in done]
    print(f"[{a.run_id}] {len(fnames)} targets, {len(done)} done, {len(todo)} to run, workers={a.workers}", flush=True)

    t0 = time.time()
    n_pass_real = n_pass_inv = 0
    with Pool(a.workers) as pool, open(out_path, "a") as fh:
        for i, res in enumerate(pool.imap_unordered(_worker, [(f, cfg, a.cache_dir, a.keep_fits) for f in todo]), 1):
            if res.get("ok") and "real" in res:
                pr, _ = passes_first_screen(res["real"]["s1"], cfg)
                pi, _ = passes_first_screen(res["inverted"]["s1"], cfg) if "inverted" in res else (False, [])
                res["pass_real_s1"], res["pass_inv_s1"] = bool(pr), bool(pi)
                n_pass_real += pr
                n_pass_inv += pi
            fh.write(json.dumps(res, default=_json_default) + "\n")
            fh.flush()
            if i % 5 == 0 or i == len(todo):
                el = time.time() - t0
                s1 = res.get("real", {}).get("s1", {}) if res.get("ok") else {}
                print(f"  {i}/{len(todo)}  {el/60:.1f} min  {el/i:.0f} s/star  pass(real)={n_pass_real} pass(inv)={n_pass_inv}"
                      f"  last={res['fname']} P={s1.get('period_d', float('nan')):.4f} sde={s1.get('sde', float('nan')):.1f}"
                      f" snr_emp={s1.get('snr_emp', float('nan')):.1f} rt={res.get('runtime_s')}s"
                      + (f" ERR={res['error'][:80]}" if res.get('error') else ""), flush=True)
    summarize(run_dir, cfg)


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def summarize(run_dir: str, cfg: dict):
    rows = [json.loads(l) for l in open(os.path.join(run_dir, "results.jsonl"))]
    ok = [r for r in rows if r.get("ok") and "real" in r]
    real = [r for r in ok if r.get("pass_real_s1")]
    inv = [r for r in ok if r.get("pass_inv_s1")]
    summ = {
        "n_rows": len(rows), "n_ok": len(ok), "n_error": sum(1 for r in rows if r.get("error")),
        "n_skipped": sum(1 for r in rows if r.get("skipped")),
        "first_screen_pass_real": len(real), "first_screen_pass_inverted": len(inv),
        "mean_runtime_s": (sum(r.get("runtime_s", 0) for r in rows) / max(1, len(rows))),
        "real_passes": [{"fname": r["fname"], **{k: r["real"]["s1"].get(k) for k in ("period_d", "depth_ppm", "duration_h", "sde", "snr_emp", "n_transits", "odd_even_sigma", "secondary_sigma")}} for r in real],
        "inverted_passes": [{"fname": r["fname"], **{k: r["inverted"]["s1"].get(k) for k in ("period_d", "depth_ppm", "duration_h", "sde", "snr_emp", "n_transits")}} for r in inv],
    }
    with open(os.path.join(run_dir, "summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, default=_json_default)
    print(json.dumps({k: v for k, v in summ.items() if not k.endswith("passes")}, indent=1))
    print("real passes:", json.dumps(summ["real_passes"], indent=1, default=_json_default))
    print("inverted passes:", json.dumps(summ["inverted_passes"], indent=1, default=_json_default))


if __name__ == "__main__":
    main()
