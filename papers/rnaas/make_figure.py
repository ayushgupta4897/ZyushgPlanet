"""Figure 1 of the RNAAS note, from report/data.json (itself exported from runs/ by scripts/report_data.py).

    python papers/rnaas/make_figure.py   -> papers/rnaas/figure1.pdf, figure1.png
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, "..", "..", "report", "data.json")))
C = [D["candidates"][k] for k in ("kic-006699887", "kic-003831969", "kic-011138423")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
fig, ax = plt.subplots(2, 3, figsize=(7.1, 4.9), gridspec_kw={"height_ratios": [1.15, 1], "hspace": 0.36, "wspace": 0.42})
for j, c in enumerate(C):
    a = ax[0, j]
    p = np.array(c["lc_points"], float); b = np.array(c["lc_binned"], float); m = np.array(c["lc_model"], float)
    a.plot(p[:, 0], p[:, 1] / 1e3, ".", ms=1.2, color="0.72", rasterized=True)
    a.errorbar(b[:, 0], b[:, 1] / 1e3, b[:, 2] / 1e3, fmt="o", ms=2.2, color="k", elinewidth=0.6, capsize=0)
    a.plot(m[:, 0], m[:, 1] / 1e3, "-", color="#b2182b", lw=1.1)
    d = c["depth_ppm"] / 1e3
    a.set_ylim(-2.4 * d, 1.2 * d); a.set_xlim(b[0, 0], b[-1, 0])
    a.set_xlabel("hours from mid-transit")
    if j == 0:
        a.set_ylabel("relative flux $-$ 1 (ppt)")
    a.set_title(f"KIC {c['kic']}   P = {c['period_d']:.4f} d", fontsize=8)
    a.text(0.03, 0.04, f"{c['n_events']} transits\n$R_p$ = {c['rp_earth']['median']:.2f} $R_\\oplus$\nFPP = {100 * c['fpp']['molusc']['mean']:.2f}%",
           transform=a.transAxes, fontsize=6.8, va="bottom")
    mp = max(c["season_maps"], key=lambda s: np.max(s["snr"]))
    snr = np.array(mp["snr"], float); H, W = snr.shape
    a2 = ax[1, j]
    im = a2.imshow(snr, origin="lower", cmap="Greys", vmin=-3, vmax=max(6, np.nanmax(snr)), extent=(-0.5, W - 0.5, -0.5, H - 0.5))
    tx, ty = mp["target"]
    a2.plot(tx, ty, "o", mfc="none", mec="#b2182b", ms=9, mew=1.3)
    for nx, ny in mp["neighbours"]:
        if np.hypot(nx - tx, ny - ty) > 0.3 and -0.5 <= nx <= W - 0.5 and -0.5 <= ny <= H - 0.5:
            a2.plot(nx, ny, "x", color="#2166ac", ms=4.5, mew=1.0)
    a2.set_xlim(-0.5, W - 0.5); a2.set_ylim(-0.5, H - 0.5)
    a2.set_xticks([]); a2.set_yticks([])
    ps = c["pixel_season"]
    a2.set_title(f"season {mp['season']};  $\\Delta\\chi^2$ = {ps['sum_dchi2_null_minus_target']:.0f} / {ps['sum_dchi2_bestalt_minus_target']:.0f}", fontsize=7, loc="left")
    cb = fig.colorbar(im, ax=a2, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=6); cb.set_label("S/N of missing flux", fontsize=6.5)
fig.savefig(os.path.join(here, "figure1.pdf"), bbox_inches="tight", dpi=300)
fig.savefig(os.path.join(here, "figure1.png"), bbox_inches="tight", dpi=200)
print("wrote figure1.pdf/png")
