"""
Sandbox-side driver for the deployed Modal app (plain HTTPS, no gRPC).

Config file .modal_endpoints.json (gitignored) in the repo root:
    {"base": "https://<workspace>--zyushg-planet-screen", "suffix": ".us-east.modal.run", "token": "<ZP_ENDPOINT_TOKEN>"}
Endpoint URLs are   <base>-start<suffix>   <base>-status<suffix>   etc.

    python -m pipeline.modal_client start  --run-id poolA_pilot --start 1 --stop 300
    python -m pipeline.modal_client status --run-id poolA_pilot
    python -m pipeline.modal_client fetch  --run-id poolA_pilot      # -> runs/<run-id>/results.jsonl
    python -m pipeline.modal_client stop   --call-id fc-...
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _cfg():
    p = os.path.join(ROOT, ".modal_endpoints.json")
    if not os.path.exists(p):
        sys.exit("missing .modal_endpoints.json (base url + token)")
    return json.load(open(p))


def _url(c, name):
    return f"{c['base']}-{name}{c.get('suffix', '.modal.run')}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["start", "status", "fetch", "stop"])
    ap.add_argument("--run-id")
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--stop", type=int, default=100)
    ap.add_argument("--fnames", nargs="*")
    ap.add_argument("--targets", default="targets/poolA_kbonus_bg_dwarfs_ranked.csv")
    ap.add_argument("--config", default="pipeline/prereg_poolA.json")
    ap.add_argument("--call-id")
    a = ap.parse_args()
    c = _cfg()
    h = {"X-ZP-Token": c["token"]}

    if a.cmd == "start":
        import pandas as pd
        if a.fnames:
            fnames = a.fnames
        else:
            df = pd.read_csv(os.path.join(ROOT, a.targets))
            fnames = df[(df["rank"] >= a.start) & (df["rank"] <= a.stop)]["fname"].tolist()
        code = {"batch_search.py": open(os.path.join(ROOT, "pipeline", "batch_search.py")).read()}
        cfg = json.load(open(os.path.join(ROOT, a.config)))
        r = requests.post(_url(c, "start"), headers=h, json={"run_id": a.run_id, "fnames": fnames, "code": code, "config": cfg}, timeout=120)
        r.raise_for_status()
        print(r.json())
        run_dir = os.path.join(ROOT, "runs", a.run_id)
        os.makedirs(run_dir, exist_ok=True)
        with open(os.path.join(run_dir, "modal_call.json"), "w") as fh:
            json.dump(r.json(), fh)
        with open(os.path.join(run_dir, "config_used.json"), "w") as fh:
            json.dump(cfg, fh, indent=1)
    elif a.cmd == "status":
        r = requests.get(_url(c, "status"), headers=h, params={"run_id": a.run_id}, timeout=60)
        r.raise_for_status()
        print(json.dumps(r.json(), indent=1))
        run_dir = os.path.join(ROOT, "runs", a.run_id)
        if os.path.isdir(run_dir):
            with open(os.path.join(run_dir, "modal_status.json"), "w") as fh:
                json.dump(r.json(), fh, indent=1)
    elif a.cmd == "fetch":
        run_dir = os.path.join(ROOT, "runs", a.run_id)
        os.makedirs(run_dir, exist_ok=True)
        out = os.path.join(run_dir, "results.jsonl")
        if not os.path.exists(os.path.join(run_dir, "config_used.json")):
            with open(os.path.join(run_dir, "config_used.json"), "w") as fh:
                json.dump(json.load(open(os.path.join(ROOT, a.config))), fh, indent=1)
        off_p = os.path.join(run_dir, ".fetch_offset")
        offset = int(open(off_p).read()) if os.path.exists(off_p) else 0
        r = requests.get(_url(c, "results"), headers=h, params={"run_id": a.run_id, "offset": offset}, timeout=300)
        r.raise_for_status()
        with open(out, "a") as fh:
            fh.write(r.text)
        nxt = r.headers.get("X-Next-Offset", str(offset + len(r.content)))
        open(off_p, "w").write(nxt)
        n = sum(1 for _ in open(out))
        print(f"fetched {len(r.content)} bytes; {n} result lines now in {out}")
    elif a.cmd == "stop":
        r = requests.post(_url(c, "stop"), headers=h, json={"call_id": a.call_id}, timeout=60)
        r.raise_for_status()
        print(r.json())


if __name__ == "__main__":
    main()
