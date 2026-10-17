#!/usr/bin/env python
"""Fold figure for the independent re-derivation of kic-008620445 (reads results.json and the
scratch fold_data.npz written by rederive_kic008620445.py; no pipeline code)."""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-008620445/rederive"
SCRATCH = ("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/"
           "scratchpad/rederive-kic008620445")
R = json.load(open(os.path.join(OUT, "results.json")))
Z = np.load(os.path.join(SCRATCH, "fold_data.npz"))

# palette: fixed categorical order (blue, orange, aqua), neutral inks, recessive grid
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.facecolor": SURF,
                     "figure.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False})

t, f, Q = Z["t"], Z["f"], Z["Q"]
ts, fs = Z["ts"], Z["fs"]
Pb, Db, t0b = float(Z["Pb"]), float(Z["Db"]), float(Z["t0b"])


def phase_time(tt, P, t0):
    return ((tt - t0 + 0.5 * P) % P) - 0.5 * P


def binned(x, y, width, lo, hi):
    edges = np.arange(lo, hi + width, width)
    idx = np.digitize(x, edges) - 1
    xb, yb, eb = [], [], []
    for i in range(len(edges) - 1):
        m = idx == i
        if m.sum() >= 3:
            xb.append(0.5 * (edges[i] + edges[i + 1]))
            yb.append(np.mean(y[m]))
            eb.append(np.std(y[m]) / np.sqrt(m.sum()))
    return np.array(xb), np.array(yb), np.array(eb)


fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
ax = axes.ravel()

# (a) periodogram with inverted null
P, pw, Pi, pwi = Z["P"], Z["power"], Z["Pi"], Z["poweri"]
ax[0].plot(Pi, pwi, color=INK2, lw=0.5, alpha=0.45, label="inverted flux (2-f) null")
ax[0].plot(P, pw, color=C1, lw=0.6, label="PSF flux")
for k, pk in enumerate(R["bls_psf"]["top3_raw"]):
    ax[0].annotate(f"{k+1}: {pk['period_d']:.4f} d", (pk["period_d"], pk["power"]),
                   xytext=(4, 6 - 10 * k), textcoords="offset points", fontsize=8, color=INK)
ax[0].set_xscale("log")
ax[0].set_xlabel("period (d)")
ax[0].set_ylabel("BLS power (delta log-likelihood)")
ax[0].set_title(f"Blind BLS 1-50 d, {R['bls_psf']['n_periods']} periods; "
                f"SDE raw {R['bls_psf']['sde_raw']:.1f}, inverted {R['bls_inverted_null']['sde_raw']:.1f}",
                fontsize=9)
ax[0].legend(frameon=False, fontsize=8, loc="upper left")

# (b) primary fold
dt = phase_time(t, Pb, t0b)
m = np.abs(dt) < 0.5
ax[1].plot(dt[m] * 24, (f[m] - 1) * 1e6, ".", color=INK2, ms=2, alpha=0.25, rasterized=True)
xb, yb, eb = binned(dt[m], f[m], 0.02, -0.5, 0.5)
ax[1].errorbar(xb * 24, (yb - 1) * 1e6, eb * 1e6, fmt="o", color=C1, ms=4, lw=1, capsize=0, label="30-min bins")
dep = R["primary"]["depth_ppm"]
ax[1].plot([-Db / 2 * 24, -Db / 2 * 24, Db / 2 * 24, Db / 2 * 24], [0, -dep, -dep, 0], color=C2, lw=1.5,
           label=f"box {dep:.0f} +- {R['primary']['empirical_err_ppm']:.0f} ppm")
ax[1].set_ylim(-4000, 2500)
ax[1].set_xlabel("hours from mid-transit")
ax[1].set_ylabel("relative flux (ppm)")
ax[1].set_title(f"PSF fold  P={Pb:.5f} d  t0={t0b:.4f} BKJD  dur={Db*24:.2f} h  "
                f"empirical SNR {R['primary']['empirical_snr']:.1f}", fontsize=9)
ax[1].legend(frameon=False, fontsize=8, loc="lower right")

# (c) odd / even
ep = np.round((t - t0b) / Pb).astype(int)
for sel, col, lab in [(ep % 2 == 1, C1, "odd"), (ep % 2 == 0, C2, "even")]:
    mm = sel & (np.abs(dt) < 0.5)
    xb, yb, eb = binned(dt[mm], f[mm], 0.03, -0.5, 0.5)
    ax[2].errorbar(xb * 24, (yb - 1) * 1e6, eb * 1e6, fmt="o", color=col, ms=4, lw=1, capsize=0, label=lab)
oe = R["odd_even"]
ax[2].set_ylim(-4000, 2500)
ax[2].set_xlabel("hours from mid-transit")
ax[2].set_ylabel("relative flux (ppm)")
ax[2].set_title(f"odd {oe['odd_depth_ppm']:.0f}+-{oe['odd_err_ppm']:.0f}  even {oe['even_depth_ppm']:.0f}"
                f"+-{oe['even_err_ppm']:.0f} ppm  ({oe['sigma']:.1f} sigma)", fontsize=9)
ax[2].legend(frameon=False, fontsize=8, loc="lower right")

# (d) secondary at phase 0.5
dts = phase_time(t, Pb, t0b + Pb / 2)
m = np.abs(dts) < 0.5
ax[3].plot(dts[m] * 24, (f[m] - 1) * 1e6, ".", color=INK2, ms=2, alpha=0.25, rasterized=True)
xb, yb, eb = binned(dts[m], f[m], 0.02, -0.5, 0.5)
ax[3].errorbar(xb * 24, (yb - 1) * 1e6, eb * 1e6, fmt="o", color=C1, ms=4, lw=1, capsize=0)
ax[3].axvspan(-Db / 2 * 24, Db / 2 * 24, color=C2, alpha=0.15, lw=0)
ax[3].set_ylim(-4000, 2500)
ax[3].set_xlabel("hours from phase 0.5")
ax[3].set_ylabel("relative flux (ppm)")
ax[3].set_title(f"secondary at phase 0.5: {R['secondary']['depth_ppm']:.0f} ppm "
                f"({R['secondary']['sigma']:.1f} sigma)", fontsize=9)

# (e) per-event and per-year depths
ev = np.genfromtxt(os.path.join(OUT, "events.csv"), delimiter=",", names=True)
ok = ev["n_in"] >= 1
ax[4].errorbar(ev["tc_bkjd"][ok], ev["depth_ppm"][ok], ev["err_ppm"][ok], fmt="o", color=INK2, ms=3,
               lw=0.7, capsize=0, alpha=0.6, label="per-event (local baseline)")
yrs = R["per_year"]
ycent = {"Y1_Q1-4": 335, "Y2_Q5-8": 720, "Y3_Q9-12": 1090, "Y4_Q13-17": 1430}
for k, v in yrs.items():
    ax[4].errorbar(ycent[k], v["depth_ppm"], v["empirical_err_ppm"], fmt="s", color=C1, ms=7, lw=1.6,
                   capsize=3)
    ax[4].annotate(f"{v['depth_ppm']:.0f}", (ycent[k], v["depth_ppm"]), xytext=(6, 4),
                   textcoords="offset points", fontsize=8, color=INK)
ax[4].axhline(dep, color=C2, lw=1, label="global box depth")
ax[4].axhline(0, color=INK2, lw=0.6)
ax[4].set_ylim(-6000, 8000)
ax[4].set_xlabel("time (BKJD)")
ax[4].set_ylabel("depth (ppm)")
ax[4].set_title(f"per-event and per-year depths; {R['transits']['n_with_ge1_point']} of "
                f"{R['transits']['n_epochs_in_baseline']} epochs have in-box data", fontsize=9)
ax[4].legend(frameon=False, fontsize=8, loc="upper right")

# (f) SAP vs PSF fold
dtsap = phase_time(ts, Pb, t0b)
for x, y, col, lab in [(dt, f, C1, "PSF FLUX"), (dtsap, fs, C3, "SAP_FLUX")]:
    mm = np.abs(x) < 0.5
    xb, yb, eb = binned(x[mm], y[mm], 0.02, -0.5, 0.5)
    ax[5].errorbar(xb * 24, (yb - 1) * 1e6, eb * 1e6, fmt="o", color=col, ms=4, lw=1, capsize=0, label=lab)
sv = R["sap_vs_psf"]
ax[5].set_ylim(-3000, 2000)
ax[5].set_xlabel("hours from mid-transit")
ax[5].set_ylabel("relative flux (ppm)")
ax[5].set_title(f"SAP depth {R['sap']['depth_ppm']:.0f}+-{R['sap']['empirical_err_ppm']:.0f} ppm; "
                f"SAP/PSF = {sv['ratio']:.2f}+-{sv['ratio_err']:.2f} (CROWDSAP {sv['expected_ratio_if_on_target_crowdsap']:.2f})",
                fontsize=9)
ax[5].legend(frameon=False, fontsize=8, loc="lower right")

fig.suptitle("kic-008620445: independent re-derivation (own loader, 1 d running-median detrend, astropy BLS) "
             "- first-screen signal, not a candidate", fontsize=10, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(os.path.join(OUT, "fold_rederive.png"), dpi=130)
print("wrote", os.path.join(OUT, "fold_rederive.png"))
