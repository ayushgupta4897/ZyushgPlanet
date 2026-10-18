"""Label the artist's impressions in report/art/raw/ -> report/art/<fname>.png.

The raw images are AI-generated illustrations (OpenAI gpt-image-2.5-sunburst; prompts in report/art/prompts.json),
written from each candidate's measured radius, period and host temperature. They are not photographs: the candidates
have never been imaged and are known only from transits in Kepler photometry. Every labelled copy says so.
"""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, "data.json")))
for k in ["kic-006699887", "kic-003831969", "kic-011138423"]:
    c, h = D["candidates"][k], D["candidates"][k]["host"]
    img = Image.open(os.path.join(here, "art", "raw", f"{k}.png")).convert("RGB")
    w, ht = img.size
    fig = plt.figure(figsize=(w / 150, ht / 150 + 1.35), dpi=150, facecolor="black")
    ax = fig.add_axes([0, 1.35 / (ht / 150 + 1.35), 1, (ht / 150) / (ht / 150 + 1.35)])
    ax.imshow(img)
    ax.set_axis_off()
    ax.text(0.015, 0.975, "ARTIST'S IMPRESSION · NOT A PHOTOGRAPH", transform=ax.transAxes, color="white", fontsize=10,
            fontweight="bold", va="top", alpha=0.9)
    ly = round(h["dist_pc"] * 3.26156, -2)
    fig.text(0.015, 0.125, f"KIC {int(c['kic'])}  ·  transit candidate (unconfirmed)  ·  ~{ly:,.0f} light-years",
             color="white", fontsize=15, va="center")
    extra = "  ·  second on-target signal at 47.7 d shown as the distant dot" if k == "kic-006699887" else ""
    fig.text(0.015, 0.075, f"P = {c['period_d']:.3f} d  ·  R = {c['rp_earth']['median']:.1f} $R_\\oplus$  ·  host {h['teff_gspphot']:,.0f} K, "
             f"G = {h['G']:.1f}  ·  FPP {100 * c['fpp']['molusc']['mean']:.2f}%{extra}", color="#c8c8c8", fontsize=10, va="center")
    fig.text(0.015, 0.032, "AI-generated illustration (OpenAI gpt-image-2.5) drawn from the measured parameters. "
             "The candidate has never been imaged; it is known only from transits in Kepler data.", color="#8a8a8a", fontsize=8.5, va="center")
    fig.savefig(os.path.join(here, "art", f"{k}.png"), facecolor="black")
    plt.close(fig)
    print("wrote", k)
