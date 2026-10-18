"""Real sky images of the three candidate host stars, labelled with their Gaia DR3 distance.

    python report/sky_cutouts.py        # -> report/sky/<fname>.png, report/sky/hosts.png, report/sky/plates/<fname>.jpg

Images: Pan-STARRS1 DR1 colour (z, g+z, g) via the CDS hips2fits service; one panel per host, 1.5'
field, north up, east left. The transit candidates themselves are not visible in any image:
they were detected only as periodic dimming of the host star in Kepler photometry.
Distances: Gaia DR3 (report/data.json, host.dist_pc).
"""
import io
import json
import os
import urllib.parse
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, "data.json")))
K = ["kic-006699887", "kic-003831969", "kic-011138423"]
FOV_ARCMIN, NPIX = 1.5, 900
os.makedirs(os.path.join(here, "sky"), exist_ok=True)


def cutout(ra, dec, npix=NPIX, fov_arcmin=FOV_ARCMIN):
    q = urllib.parse.urlencode({"hips": "CDS/P/PanSTARRS/DR1/color-z-zg-g", "width": npix, "height": npix, "fov": fov_arcmin / 60,
                                "projection": "TAN", "coordsys": "icrs", "ra": ra, "dec": dec, "format": "png"})
    return Image.open(io.BytesIO(urllib.request.urlopen("https://alasky.cds.unistra.fr/hips-image-services/hips2fits?" + q, timeout=120).read())).convert("RGB")


def panel(ax, c, img):
    ly = c["host"]["dist_pc"] * 3.26156
    ax.imshow(img, origin="upper")
    ax.set_axis_off()
    m, g, r = NPIX / 2, NPIX * 0.035, NPIX * 0.09
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ax.plot([m + dx * g, m + dx * r], [m + dy * g, m + dy * r], color="white", lw=1.6)
    bar = NPIX / FOV_ARCMIN * (10 / 60)
    ax.plot([NPIX * 0.06, NPIX * 0.06 + bar], [NPIX * 0.93] * 2, color="white", lw=2)
    ax.text(NPIX * 0.06 + bar / 2, NPIX * 0.905, '10"', color="white", ha="center", fontsize=9)
    ax.text(NPIX * 0.94, NPIX * 0.93, "N up, E left", color="white", ha="right", fontsize=8, alpha=0.8)
    ax.set_title(f"KIC {int(c['kic'])}\n~{round(ly, -2):,.0f} light-years ({c['host']['dist_pc']:,.0f} pc)", color="white", fontsize=12, loc="left")
    ax.text(0, -0.04, f"Candidate: P = {c['period_d']:.3f} d, R = {c['rp_earth']['median']:.1f} $R_\\oplus$   ·   host G = {c['host']['G']:.1f}, "
            f"{c['host']['teff_gspphot']:,.0f} K", transform=ax.transAxes, color="#bbbbbb", fontsize=8.5, va="top")


imgs = {k: cutout(D["candidates"][k]["ra"], D["candidates"][k]["dec"]) for k in K}
for k in K:
    fig, ax = plt.subplots(figsize=(6, 6.6), facecolor="black")
    panel(ax, D["candidates"][k], imgs[k])
    fig.text(0.02, 0.015, "Pan-STARRS1 colour image of the host star (crosshair). The candidate itself is not resolved; it was detected as transits in Kepler data.",
             color="#888888", fontsize=6.5, wrap=True)
    fig.savefig(os.path.join(here, "sky", f"{k}.png"), dpi=200, facecolor="black", bbox_inches="tight")
    plt.close(fig)
fig, axs = plt.subplots(1, 3, figsize=(16, 6.4), facecolor="black")
for ax, k in zip(axs, K):
    panel(ax, D["candidates"][k], imgs[k])
fig.text(0.01, 0.02, "Host stars of the three transit candidates, Pan-STARRS1 DR1 colour (z, g+z, g) via CDS hips2fits; 1.5' fields, crosshair on the host. "
         "Distances from Gaia DR3. The candidates are not resolved in any image: they were found as periodic dimming in Kepler Bonus photometry.",
         color="#888888", fontsize=8)
fig.savefig(os.path.join(here, "sky", "hosts.png"), dpi=170, facecolor="black", bbox_inches="tight")
print("wrote", os.listdir(os.path.join(here, "sky")))
# unannotated square plates for the poster (1.2' field; the poster draws its own crosshair, scale bar and labels)
os.makedirs(os.path.join(here, "sky", "plates"), exist_ok=True)
for k in K:
    cutout(D["candidates"][k]["ra"], D["candidates"][k]["dec"], npix=1200, fov_arcmin=1.2).save(os.path.join(here, "sky", "plates", f"{k}.jpg"), quality=93)
print("wrote plates")
