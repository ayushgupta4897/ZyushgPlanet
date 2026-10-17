"""
Modal fan-out for the KBONUS background-star transit screen, driven over HTTPS.

Why HTTPS endpoints: the Claude Code cloud sandbox cannot speak gRPC (Modal's
client protocol), so the app is deployed ONCE from a laptop and then controlled
through three authenticated web endpoints.  Containers only exist while a batch
is running; the deployment itself costs nothing when idle.

One-time deploy (laptop, inside the repo):

    pip install modal
    modal token set --token-id ... --token-secret ...      # already done
    export ZP_ENDPOINT_TOKEN='<shared secret from the session>'
    modal deploy pipeline/modal_app.py

    # prints three URLs like https://<workspace>--zyushg-planet-screen-start.modal.run

Tear-down when the experiments are over:

    modal app stop zyushg-planet-screen

Endpoints (all require header  X-ZP-Token: $ZP_ENDPOINT_TOKEN):
    POST /start   {"run_id": "poolA_full", "start": 1, "stop": 6465}   -> {"call_id": ...}
    GET  /status?run_id=poolA_full                                      -> progress json
    GET  /results?run_id=poolA_full&offset=0                            -> JSONL text (append-only)
    POST /stop    {"call_id": "..."}                                    -> cancels a running batch

Cost: ~90 CPU-s per star (real + inverted BLS) -> 6,465 stars ~ 160 CPU-h,
roughly US$30-40 at Modal's ~US$0.19/CPU-h; ~1-2 h wall-clock at 100 containers.
"""
from __future__ import annotations

import json
import os

import modal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APP_NAME = "zyushg-planet-screen"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy", "scipy", "astropy", "wotan", "requests", "pandas", "fastapi[standard]")
    .add_local_dir(os.path.join(ROOT, "pipeline"), remote_path="/root/pipeline")
    .add_local_file(os.path.join(ROOT, "pipeline", "prereg_poolA.json"), remote_path="/root/prereg.json")
    .add_local_file(os.path.join(ROOT, "targets", "poolA_kbonus_bg_dwarfs_ranked.csv"), remote_path="/root/targets.csv")
)

app = modal.App(APP_NAME, image=image)
vol = modal.Volume.from_name("zyushg-planet-runs", create_if_missing=True)
secret = modal.Secret.from_local_environ(["ZP_ENDPOINT_TOKEN"])
RUNS = "/runs"

with image.imports():
    from fastapi import HTTPException, Request
    from fastapi.responses import PlainTextResponse


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


@app.function(cpu=1.0, memory=1024, timeout=24 * 3600, volumes={RUNS: vol})
def run_batch(run_id: str, fnames: list[str]) -> dict:
    """Map screen_one over fnames and stream results into the Volume (resumable)."""
    import time
    d = os.path.join(RUNS, run_id)
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, "results.jsonl")
    done = set()
    if os.path.exists(out):
        with open(out) as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["fname"])
                except Exception:  # noqa: BLE001
                    pass
    todo = [f for f in fnames if f not in done]
    status = {"run_id": run_id, "total": len(fnames), "already_done": len(done), "todo": len(todo),
              "done": 0, "pass_real": 0, "pass_inv": 0, "errors": 0, "state": "running", "started": time.time()}

    def flush():
        with open(os.path.join(d, "status.json"), "w") as fh:
            json.dump(status, fh)
        vol.commit()

    flush()
    last = time.time()
    with open(out, "a") as fh:
        for res in screen_one.map(todo, order_outputs=False, return_exceptions=True):
            if isinstance(res, Exception):
                res = {"fname": "?", "ok": False, "error": repr(res)}
            fh.write(json.dumps(res, default=str) + "\n")
            fh.flush()
            status["done"] += 1
            status["pass_real"] += bool(res.get("pass_real_s1"))
            status["pass_inv"] += bool(res.get("pass_inv_s1"))
            status["errors"] += bool(res.get("error"))
            if time.time() - last > 60:
                flush()
                last = time.time()
    status["state"] = "finished"
    status["finished"] = time.time()
    flush()
    return status


def _auth(request) -> None:
    if request.headers.get("x-zp-token", "") != os.environ.get("ZP_ENDPOINT_TOKEN", "__unset__"):
        raise HTTPException(status_code=401, detail="bad token")


@app.function(secrets=[secret], volumes={RUNS: vol}, timeout=120)
@modal.fastapi_endpoint(method="POST")
def start(body: dict, request: Request):
    _auth(request)
    import pandas as pd
    run_id = body["run_id"]
    if body.get("fnames"):
        fnames = list(body["fnames"])
    else:
        df = pd.read_csv("/root/targets.csv")
        df = df[(df["rank"] >= int(body.get("start", 1))) & (df["rank"] <= int(body.get("stop", 100)))]
        fnames = df["fname"].tolist()
    call = run_batch.spawn(run_id, fnames)
    return {"call_id": call.object_id, "run_id": run_id, "n": len(fnames)}


@app.function(secrets=[secret], volumes={RUNS: vol}, timeout=60)
@modal.fastapi_endpoint(method="GET")
def status(run_id: str, request: Request):
    _auth(request)
    vol.reload()
    p = os.path.join(RUNS, run_id, "status.json")
    if not os.path.exists(p):
        return {"run_id": run_id, "state": "unknown"}
    return json.load(open(p))


@app.function(secrets=[secret], volumes={RUNS: vol}, timeout=120)
@modal.fastapi_endpoint(method="GET")
def results(run_id: str, request: Request, offset: int = 0):
    _auth(request)
    vol.reload()
    p = os.path.join(RUNS, run_id, "results.jsonl")
    if not os.path.exists(p):
        return PlainTextResponse("", headers={"X-Next-Offset": "0"})
    with open(p, "rb") as fh:
        fh.seek(int(offset))
        data = fh.read()
    return PlainTextResponse(data.decode("utf-8", "replace"), headers={"X-Next-Offset": str(int(offset) + len(data))})


@app.function(secrets=[secret], timeout=60)
@modal.fastapi_endpoint(method="POST")
def stop(body: dict, request: Request):
    _auth(request)
    modal.FunctionCall.from_id(body["call_id"]).cancel()
    return {"cancelled": body["call_id"]}


@app.local_entrypoint()
def main(start: int = 1, stop: int = 100, run_id: str = "poolA_modal"):
    """Laptop-only alternative: `modal run pipeline/modal_app.py --start 1 --stop 300 --run-id poolA_pilot`
    streams results straight into runs/<run-id>/results.jsonl on the laptop."""
    import pandas as pd
    df = pd.read_csv(os.path.join(ROOT, "targets", "poolA_kbonus_bg_dwarfs_ranked.csv"))
    df = df[(df["rank"] >= start) & (df["rank"] <= stop)]
    run_dir = os.path.join(ROOT, "runs", run_id)
    os.makedirs(run_dir, exist_ok=True)
    out = os.path.join(run_dir, "results.jsonl")
    done = {json.loads(l)["fname"] for l in open(out)} if os.path.exists(out) else set()
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
