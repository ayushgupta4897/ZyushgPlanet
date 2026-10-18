"""RNAAS turnaround from Crossref (AAS deposits date_received / date_accepted / date_epub per note).

    python scripts/rnaas_turnaround.py [--since 2024-01-01] [--out docs/briefs/rnaas_turnaround.json]

Only published notes carry these dates, so declined submissions are invisible here (survivorship).
"""
import argparse
import datetime as dt
import json
import re
import statistics as st
import urllib.parse
import urllib.request

ap = argparse.ArgumentParser()
ap.add_argument("--since", default="2024-01-01")
ap.add_argument("--out", default="docs/briefs/rnaas_turnaround.json")
a = ap.parse_args()
items, off = [], 0
while True:
    q = urllib.parse.urlencode({"filter": f"issn:2515-5172,from-pub-date:{a.since}", "rows": 1000, "offset": off, "select": "DOI,title,assertion"})
    m = json.load(urllib.request.urlopen("https://api.crossref.org/works?" + q, timeout=180))["message"]
    items += m["items"]
    off += 1000
    if off >= m["total-results"]:
        break
rows = []
for it in items:
    d = {x["name"]: x["value"] for x in it.get("assertion", []) if "name" in x}
    try:
        r, acc, e = (dt.date.fromisoformat(d[k][:10]) for k in ("date_received", "date_accepted", "date_epub"))
    except (KeyError, ValueError):
        continue
    rows.append({"doi": it["DOI"], "title": re.sub(r"<[^>]+>|\s+", " ", (it.get("title") or [""])[0]).strip(), "received": str(r), "accept_days": (acc - r).days, "online_days": (e - r).days})


def summary(rs):
    out = {"n": len(rs)}
    for k in ("accept_days", "online_days"):
        v = sorted(x[k] for x in rs)
        out[k] = {"median": st.median(v), "p10": v[int(0.1 * (len(v) - 1))], "p25": v[int(0.25 * (len(v) - 1))], "p75": v[int(0.75 * (len(v) - 1))], "p90": v[int(0.9 * (len(v) - 1))], "max": v[-1], "within_3d": sum(x <= 3 for x in v)}
    return out


planet = re.compile(r"(?i)exoplanet|planet|transit|TOI|KOI|hot jupiter|neptune|super-earth")
res = {"source": "api.crossref.org works, issn 2515-5172", "fetched_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "since": a.since,
       "all": summary(rows), "y2026": summary([x for x in rows if x["received"] >= "2026"]),
       "y2026_planet_titles": summary([x for x in rows if x["received"] >= "2026" and planet.search(x["title"])]),
       "rows": rows}
json.dump(res, open(a.out, "w"), indent=1)
for k in ("all", "y2026", "y2026_planet_titles"):
    s = res[k]
    print(f"{k}: n={s['n']}  accept median {s['accept_days']['median']} d (p25-p75 {s['accept_days']['p25']}-{s['accept_days']['p75']}, p90 {s['accept_days']['p90']}); online median {s['online_days']['median']} d (p90 {s['online_days']['p90']}); online within 3 d of receipt {s['online_days']['within_3d']}")
