"""
Download the known-object catalogues used by ephemeris_match.py / triage.py / vet_batch.py
into data/known/ and write data/known/manifest.json (URL, UTC time, SHA-256, row count).

    . .venv/bin/activate && python -m pipeline.fetch_known

Sources (all public):
  * NASA Exoplanet Archive TAP  https://exoplanetarchive.ipac.caltech.edu/TAP/sync
      - cumulative (KOI table, all dispositions)      -> koi_cumulative.csv
      - q1_q17_dr25_tce (Kepler DR25 TCEs)            -> dr25_tce.csv
      - ps (confirmed planets in the Kepler field)     -> ps_keplerfield.csv
      - keplerstellar (all Kepler targets; kepid only) -> kepler_target_kepids.csv
  * VizieR J/AJ/151/68 (Kirk et al. 2016 Kepler EB catalogue)   -> J_AJ_151_68_catalog.csv
  * VizieR J/ApJS/258/16 (Prsa et al. 2022 TESS EB catalogue)   -> J_ApJS_258_16_tess-ebs.csv
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import sys

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "known")
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

QUERIES = {
    "koi_cumulative.csv": "select kepid,kepoi_name,kepler_name,koi_disposition,koi_pdisposition,koi_period,koi_time0bk,koi_duration,koi_depth,koi_fpflag_nt,koi_fpflag_ss,koi_fpflag_co,koi_fpflag_ec,ra,dec,koi_kepmag from cumulative",
    "dr25_tce.csv": "select kepid,tce_plnt_num,tce_period,tce_time0bk,tce_duration,tce_depth,tce_max_mult_ev,ra,dec from q1_q17_dr25_tce",
    "ps_keplerfield.csv": "select pl_name,hostname,tic_id,pl_orbper,pl_tranmid,pl_trandur,pl_trandep,ra,dec,disc_facility from ps where default_flag=1 and ra between 279 and 302.5 and dec between 36 and 53",
    "keplerstellar.csv": "select kepid,st_quarters,kepmag,teff,radius,ra,dec from keplerstellar",
}
VIZIER = {
    "J_AJ_151_68_catalog.csv": ("J/AJ/151/68", "catalog"),
    "J_ApJS_258_16_tess-ebs.csv": ("J/ApJS/258/16", "tess-ebs"),
}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    os.makedirs(OUT, exist_ok=True)
    manifest = {"fetched_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "files": {}}
    for name, q in QUERIES.items():
        path = os.path.join(OUT, name)
        r = requests.get(TAP, params={"query": q, "format": "csv"}, timeout=900)
        r.raise_for_status()
        with open(path, "wb") as fh:
            fh.write(r.content)
        n = sum(1 for _ in open(path)) - 1
        manifest["files"][name] = {"source": TAP, "query": q, "rows": n, "sha256": _sha(path)}
        print(f"{name}: {n} rows", flush=True)
    # derived: list of KIC ids that were ever observed as Kepler targets
    import pandas as pd
    ks = pd.read_csv(os.path.join(OUT, "keplerstellar.csv"))
    obs = ks[ks.st_quarters.astype(str).str.contains("1")][["kepid"]].drop_duplicates()
    p = os.path.join(OUT, "kepler_target_kepids.csv")
    obs.to_csv(p, index=False)
    manifest["files"]["kepler_target_kepids.csv"] = {"derived_from": "keplerstellar.csv (any observed quarter)", "rows": int(len(obs)), "sha256": _sha(p)}
    print(f"kepler_target_kepids.csv: {len(obs)} rows", flush=True)
    from astroquery.vizier import Vizier
    v = Vizier(columns=["**"], row_limit=-1)
    for name, (cat, table) in VIZIER.items():
        res = v.get_catalogs(cat)
        t = [x for x in res if x.meta.get("name", "").endswith("/" + table)]
        if not t:
            print(f"{name}: table {table} not found in {cat}", file=sys.stderr)
            continue
        path = os.path.join(OUT, name)
        t[0].to_pandas().to_csv(path, index=False)
        manifest["files"][name] = {"source": f"VizieR {cat}/{table}", "rows": int(len(t[0])), "sha256": _sha(path)}
        print(f"{name}: {len(t[0])} rows", flush=True)
    with open(os.path.join(OUT, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    print("manifest written:", os.path.join(OUT, "manifest.json"))


if __name__ == "__main__":
    main()
