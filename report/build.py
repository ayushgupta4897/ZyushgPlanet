"""Inline report/data.json, report/map.json and report/src/viz.js into the HTML templates -> report/dist/*.html.

    python report/build.py
"""
import os

here = os.path.dirname(os.path.abspath(__file__))
data = open(os.path.join(here, "data.json")).read()
viz = open(os.path.join(here, "src", "viz.js")).read()
mp = os.path.join(here, "map.json")
star_map = open(mp).read() if os.path.exists(mp) else "null"  # scripts/poster_map_data.py
os.makedirs(os.path.join(here, "dist"), exist_ok=True)
for name in ("report.html", "screens.html", "poster.html", "poster-grand.html"):
    src = os.path.join(here, "src", name)
    if not os.path.exists(src):
        continue
    s = open(src).read().replace("/*__DATA__*/null", data).replace("/*__VIZ__*/", viz).replace("/*__MAP__*/null", star_map)
    open(os.path.join(here, "dist", name), "w").write(s)
    print("built", name, len(s) // 1024, "KB")
