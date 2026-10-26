"""Research Notes of the AAS in our corridor (exoplanets, astrobiology), with review times, from Crossref.

    python scripts/rnaas_corridor.py [--since 2025-01-01] [--out docs/data/rnaas_corridor.json]

AAS does not deposit the topical corridor in Crossref, so notes are classed by keywords in title + abstract. Review
time is date_accepted - date_received as AAS deposits them; only published notes carry dates (declines are invisible).
"""
import argparse
import datetime as dt
import html
import json
import re
import statistics as st
import urllib.parse
import urllib.request

ap = argparse.ArgumentParser()
ap.add_argument("--since", default="2025-01-01")
ap.add_argument("--out", default="docs/data/rnaas_corridor.json")
a = ap.parse_args()

items, cur = [], "*"
while True:
    q = urllib.parse.urlencode({"filter": f"issn:2515-5172,from-pub-date:{a.since}", "rows": 1000, "cursor": cur,
                                "select": "DOI,title,abstract,assertion,author,published"})
    m = json.load(urllib.request.urlopen("https://api.crossref.org/works?" + q, timeout=180))["message"]
    items += m["items"]; cur = m.get("next-cursor")
    if not m["items"] or len(items) >= m["total-results"]:
        break

EXO = re.compile(r"(?i)exoplanet|exomoon|transit(s|ing)?\b|\bTOI[- ]?\d|\bKOI[- ]?\d|hot (jupiter|neptune|saturn)|super-earth|sub-neptune|"
                 r"mini-neptune|radial[- ]velocit|planet(ary)? (candidate|host|system|companion|mass|radius|radii|occurrence)|"
                 r"planet-hosting|circumbinary planet|protoplanetary|debris disk|direct imaging|microlensing planet|"
                 r"\b(Kepler|K2|TESS|WASP|HAT-P|KELT|TRAPPIST|LHS|GJ|LTT|Wolf|TOI|XO|Qatar|NGTS|HATS)-\s?\d+\s?[b-h]\b")
ASTROBIO = re.compile(r"(?i)astrobiolog|biosignature|technosignature|\bSETI\b|habitab|prebiotic|origin(s)? of life|"
                      r"extraterrestrial (life|intelligence)|life detection|biomarker|\bbiotic\b")
CANDIDATE = re.compile(r"(?i)candidate|validat|new (planet|transit)|discover|single[- ]transit|duo[- ]?transit|false[- ]positive|"
                       r"search for (transit|planet)|vetting|TOI|CTOI|planet hunters")
SURVEY = re.compile(r"(?i)\bKepler\b|\bK2\b|\bTESS\b|Gaia|archival")

rows = []
for it in items:
    d = {x["name"]: x["value"] for x in it.get("assertion", []) if "name" in x}
    title = re.sub(r"<[^>]+>|\s+", " ", (it.get("title") or [""])[0]).strip()
    abstract = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", it.get("abstract") or ""))).strip()
    text = f"{title} {abstract}"
    exo, bio = bool(EXO.search(text)), bool(ASTROBIO.search(text))
    if not (exo or bio):
        continue
    try:
        r, acc, e = (dt.date.fromisoformat(d[k][:10]) for k in ("date_received", "date_accepted", "date_epub"))
    except (KeyError, ValueError):
        r = acc = e = None
    rows.append({"doi": it["DOI"], "title": title, "n_authors": len(it.get("author", [])),
                 "first_author": " ".join(filter(None, [(it.get("author") or [{}])[0].get("given"), (it.get("author") or [{}])[0].get("family")])),
                 "exoplanets": exo, "astrobiology": bio, "candidate_report": bool(exo and CANDIDATE.search(text) and SURVEY.search(text)),
                 "received": str(r) if r else None, "accepted": str(acc) if acc else None, "online": str(e) if e else None,
                 "accept_days": (acc - r).days if r else None, "online_days": (e - r).days if r else None})


def summary(rs):
    v = sorted(x["accept_days"] for x in rs if x["accept_days"] is not None)
    if not v:
        return {"n": 0}
    q = lambda f: v[int(f * (len(v) - 1))]  # noqa: E731
    return {"n": len(v), "median": st.median(v), "p25": q(0.25), "p75": q(0.75), "p90": q(0.9), "max": v[-1],
            "share_over": {d: round(sum(x > d for x in v) / len(v), 3) for d in (7, 10, 14, 21, 28)}}


y26 = [x for x in rows if (x["received"] or "") >= "2026"]
res = {"source": "api.crossref.org works, issn 2515-5172", "fetched_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
       "since": a.since, "n_notes_fetched": len(items), "n_corridor": len(rows),
       "classification": "keywords in title + abstract (AAS does not deposit the corridor in Crossref)",
       "review_days": {"corridor_all": summary(rows), "corridor_2026": summary(y26),
                       "exoplanets_2026": summary([x for x in y26 if x["exoplanets"]]),
                       "astrobiology_2026": summary([x for x in y26 if x["astrobiology"]]),
                       "candidate_reports_all": summary([x for x in rows if x["candidate_report"]]),
                       "single_author_2026": summary([x for x in y26 if x["n_authors"] == 1])},
       "rows": sorted(rows, key=lambda x: x["received"] or "", reverse=True)}
json.dump(res, open(a.out, "w"), indent=1)
print(f"{len(items)} notes since {a.since}; {len(rows)} in corridor "
      f"({sum(x['exoplanets'] for x in rows)} exoplanet, {sum(x['astrobiology'] for x in rows)} astrobiology, "
      f"{sum(x['candidate_report'] for x in rows)} candidate reports)")
for k, s in res["review_days"].items():
    print(k, json.dumps(s))
