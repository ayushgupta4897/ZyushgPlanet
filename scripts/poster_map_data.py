"""Sky positions and search stage of every Pool A star, for the star-field poster -> report/map.json.

    python scripts/poster_map_data.py

stage: 0 searched, 1 first-screen pass (real flux), 2 stage-2 pass, 3 candidate.
Sources: targets/poolA_kbonus_bg_dwarfs_ranked.csv, runs/poolA_{main,pilot}/triage_all.csv, runs/poolA_main/vet_summary.csv,
report/data.json (candidates).
"""
import json

import pandas as pd

tg = pd.read_csv("targets/poolA_kbonus_bg_dwarfs_ranked.csv")
s1 = set()
for run in ("poolA_main", "poolA_pilot"):
    t = pd.read_csv(f"runs/{run}/triage_all.csv")
    s1 |= set(t[t.pass_real.fillna(False).astype(bool)].fname)
vs = pd.read_csv("runs/poolA_main/vet_summary.csv")
s2 = set(vs[vs.stage2_pass == True].fname)  # noqa: E712
cands = set(json.load(open("report/data.json"))["candidates"])
stage = [3 if f in cands else 2 if f in s2 else 1 if f in s1 else 0 for f in tg.fname]
out = {"n": len(tg), "counts": {str(k): stage.count(k) for k in range(4)},
       "stars": [[round(r, 4), round(d, 4), s] for r, d, s in zip(tg.ra, tg.dec, stage)],
       "candidates": {f: [round(float(tg.loc[tg.fname == f, "ra"].iloc[0]), 5), round(float(tg.loc[tg.fname == f, "dec"].iloc[0]), 5)] for f in sorted(cands)}}
json.dump(out, open("report/map.json", "w"), separators=(",", ":"))
print(out["n"], out["counts"])
