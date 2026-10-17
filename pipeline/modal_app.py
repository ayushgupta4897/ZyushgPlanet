"""
Modal fan-out for the transit screen, driven over HTTPS, code shipped per run.

Why: the Claude Code cloud sandbox cannot speak gRPC (Modal's client protocol), so
this app is deployed ONCE from a laptop and then controlled through authenticated
web endpoints.  The search code itself (pipeline/batch_search.py) and the
pre-registered config are sent in the /start request and stored in a Modal Volume,
so pipeline changes never require a redeploy.  Containers exist only while a batch
is running; an idle deployment costs nothing.

One-time deploy (laptop, inside the repo, branch claude/planet-hunt-playbook):

    pip install modal
    export ZP_ENDPOINT_TOKEN='<shared secret from the session>'
    modal deploy pipeline/modal_app.py
    # copy the printed URLs (…--zyushg-planet-screen-start.modal.run etc.) back to the session

Tear-down when the experiments are over:

    modal app stop zyushg-planet-screen

Endpoints (all require header  X-ZP-Token: $ZP_ENDPOINT_TOKEN):
    POST /start   {"run_id": "...", "fnames": [...], "code": {"batch_search.py": "..."}, "config": {...}}
                  -> {"call_id": "...", "n": N}
    GET  /status?run_id=...                      -> progress json
    GET  /results?run_id=...&offset=0            -> JSONL text; header X-Next-Offset
    POST /stop    {"call_id": "..."}             -> cancels a running batch

Cost: ~90 CPU-s per star (real + inverted BLS) -> 6,465 stars ~ 160 CPU-h,
roughly US$30-40 at Modal's ~US$0.19/CPU-h; ~1-2 h wall-clock at 100 containers.
The sandbox-side driver is pipeline/modal_client.py.
"""
from __future__ import annotations

import json
import os

import modal

APP_NAME = "zyushg-planet-screen"
RUNS = "/runs"

image = modal.Image.debian_slim(python_version="3.12").pip_install(
    "numpy", "scipy", "astropy", "wotan", "requests", "pandas", "fastapi[standard]"
)
app = modal.App(APP_NAME, image=image)
vol = modal.Volume.from_name("zyushg-planet-runs", create_if_missing=True)
secret = modal.Secret.from_local_environ(["ZP_ENDPOINT_TOKEN"])

with image.imports():
    from fastapi import HTTPException, Request
    from fastapi.responses import PlainTextResponse


_LOADED: dict[str, object] = {}


def _load_pipeline(run_id: str):
    """Import batch_search.py + config shipped for this run (cached per container)."""
    if run_id in _LOADED:
        return _LOADED[run_id]
    import importlib.util
    import sys
    vol.reload()
    code_dir = os.path.join(RUNS, run_id, "code")
    spec = importlib.util.spec_from_file_location(f"batch_search_{run_id}", os.path.join(code_dir, "batch_search.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    cfg = mod.load_config(os.path.join(code_dir, "config.json"))
    _LOADED[run_id] = (mod, cfg)
    return _LOADED[run_id]


@app.function(cpu=2.0, memory=2048, timeout=1800, max_containers=100, retries=1, volumes={RUNS: vol})
def screen_one(fname: str, run_id: str) -> dict:
    mod, cfg = _load_pipeline(run_id)
    res = mod.search_star(fname, cfg, cache_dir="/tmp/kb", keep_fits=False)
    if res.get("ok") and "real" in res:
        res["pass_real_s1"] = bool(mod.passes_first_screen(res["real"]["s1"], cfg)[0])
        res["pass_inv_s1"] = bool(mod.passes_first_screen(res["inverted"]["s1"], cfg)[0]) if "inverted" in res else False
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
        for res in screen_one.map(todo, kwargs={"run_id": run_id}, order_outputs=False, return_exceptions=True):
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
    run_id = body["run_id"]
    fnames = list(body["fnames"])
    code_dir = os.path.join(RUNS, run_id, "code")
    os.makedirs(code_dir, exist_ok=True)
    for name, text in body["code"].items():
        with open(os.path.join(code_dir, os.path.basename(name)), "w") as fh:
            fh.write(text)
    with open(os.path.join(code_dir, "config.json"), "w") as fh:
        json.dump(body["config"], fh)
    vol.commit()
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
