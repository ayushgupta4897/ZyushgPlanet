"""Results JSON + PNG figures.

Figures (all PNG, 150 dpi, light surface):
  01_raw_stitched.png      per-quarter normalised PSF flux with the wotan trend
  02_detrended.png         detrended flux with predicted transit times marked
  03_bls_periodogram.png   BLS SNR periodogram, best period and harmonics marked
  04_tls_periodogram.png   TLS SDE periodogram
  05_phase_fold.png        fold at the best period with 30-min bins and the TLS model
  06_odd_even.png          odd and even folds
  07_secondary.png         fold centred on phase 0.5
  08_neighbour_fold.png    neighbour's own light curve folded at the TARGET ephemeris
  09_per_transit.png       individual event depths vs time
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# fixed categorical palette (dataviz skill reference palette), assigned by role, never cycled
C = {"raw": "#9a9995", "trend": "#eb6834", "data": "#2a78d6", "binned": "#0b0b0b", "model": "#eb6834",
     "marker": "#1baf7a", "odd": "#2a78d6", "even": "#eb6834", "null": "#6250d6", "neigh": "#e87ba4"}
plt.rcParams.update({"figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "axes.edgecolor": "#c3c2b7",
                     "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 9, "axes.titlesize": 10, "legend.frameon": False})


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return str(o)


def write_json(obj, path: Path):
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default))
    return str(path)


def _bin(x, y, width):
    edges = np.arange(x.min(), x.max() + width, width)
    idx = np.digitize(x, edges) - 1
    xb, yb, eb = [], [], []
    for i in range(len(edges) - 1):
        m = idx == i
        if m.sum() >= 3:
            xb.append(x[m].mean()); yb.append(np.mean(y[m])); eb.append(np.std(y[m]) / np.sqrt(m.sum()))
    return np.array(xb), np.array(yb), np.array(eb)


def fig_raw(lc, out, title):
    fig, ax = plt.subplots(figsize=(13, 3.6))
    for q in np.unique(lc.quarter):
        m = lc.quarter == q
        ax.plot(lc.time[m], lc.flux[m], ".", ms=1.5, color=C["raw"], alpha=0.7)
        ax.text(lc.time[m].mean(), np.nanpercentile(lc.flux[m], 99.7) + 0.004, f"Q{q}", ha="center", fontsize=7, color="#52514e")
    if lc.trend is not None:
        ax.plot(lc.time, lc.trend, "-", lw=0.8, color=C["trend"], label="wotan biweight trend")
        ax.legend(loc="lower left")
    ax.set_xlabel("Time (BKJD = BJD_TDB - 2454833)"); ax.set_ylabel("PSF flux / quarter median")
    ax.set_title(title); fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def fig_detrended(lc, out, title, period, t0, duration):
    fig, ax = plt.subplots(figsize=(13, 3.6))
    ax.plot(lc.time, lc.flat, ".", ms=1.5, color=C["data"], alpha=0.6, label="detrended PSF flux")
    n = np.arange(np.floor((lc.time.min() - t0) / period), np.ceil((lc.time.max() - t0) / period) + 1)
    tc = t0 + n * period
    has = np.array([np.any(np.abs(lc.time - x) < 0.5 * duration) for x in tc])
    ax.vlines(tc[has], 0.985, 0.99, color=C["marker"], lw=1, label=f"predicted transits with data ({has.sum()})")
    ax.vlines(tc[~has], 0.985, 0.99, color=C["raw"], lw=0.8, alpha=0.6, label=f"in gaps ({(~has).sum()})")
    ax.set_ylim(0.98, 1.012); ax.set_xlabel("Time (BKJD)"); ax.set_ylabel("Relative flux")
    ax.legend(loc="upper left", ncol=3); ax.set_title(title); fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def fig_periodogram(periods, power, best, out, title, ylabel, harmonics=True, null=None):
    fig, ax = plt.subplots(figsize=(13, 3.4))
    ax.plot(periods, power, "-", lw=0.6, color=C["data"], label="search")
    if null is not None:
        ax.plot(null[0], null[1], "-", lw=0.5, color=C["null"], alpha=0.7, label="inverted-flux null")
    ax.axvline(best, color=C["model"], lw=1, ls="--", label=f"best P = {best:.5f} d")
    if harmonics:
        for k in (2, 3, 4, 0.5, 1 / 3):
            if periods.min() < best * k < periods.max():
                ax.axvline(best * k, color=C["marker"], lw=0.6, ls=":", alpha=0.8)
    ax.set_xscale("log"); ax.set_xlabel("Period (d)"); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.legend(loc="upper right"); fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def fig_fold(t, f, period, t0, duration, out, title, model=None, depth_ppm=None, window_h=12, color=None,
             bin_min=30.0, center=0.0, label="detrended flux"):
    from .vet import phase
    ph = (phase(t, period, t0) - center + 0.5) % 1.0 - 0.5
    x = ph * period * 24
    m = np.abs(x) < window_h
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot(x[m], f[m], ".", ms=2, color=color or C["data"], alpha=0.35, label=label)
    xb, yb, eb = _bin(x[m], f[m], bin_min / 60)
    ax.errorbar(xb, yb, eb, fmt="o", ms=4, color=C["binned"], lw=1, label=f"{bin_min:.0f}-min bins")
    if model is not None:
        mp, my = model
        mx = mp * period * 24
        o = np.argsort(mx); ax.plot(mx[o], my[o], "-", lw=1.5, color=C["model"], label="TLS model")
    ax.axvspan(-duration * 12, duration * 12, color=C["marker"], alpha=0.08, lw=0)
    ax.set_xlim(-window_h, window_h)
    lo = np.nanpercentile(f[m], 0.3); hi = np.nanpercentile(f[m], 99.7)
    ax.set_ylim(min(lo, 1 - 2.5 * (depth_ppm or 1500) * 1e-6), hi)
    ax.set_xlabel(f"Hours from {'mid-transit' if center == 0 else 'phase %.2f' % center}"); ax.set_ylabel("Relative flux")
    ax.set_title(title); ax.legend(loc="lower left"); fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def fig_odd_even(t, f, period, t0, duration, oe, out, title, window_h=8):
    from .vet import phase, epoch
    ep = epoch(t, period, t0); x = phase(t, period, t0) * period * 24
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, name, sel in ((axes[0], "odd", ep % 2 == 1), (axes[1], "even", ep % 2 == 0)):
        m = sel & (np.abs(x) < window_h)
        ax.plot(x[m], f[m], ".", ms=2, color=C[name], alpha=0.35)
        xb, yb, eb = _bin(x[m], f[m], 0.5)
        ax.errorbar(xb, yb, eb, fmt="o", ms=4, color=C["binned"], lw=1)
        d = oe[name]
        ax.axhline(1 - d["depth"], color=C[name], lw=1, ls="--")
        ax.set_title(f"{name}: depth {d['depth']*1e6:.0f} +/- {d['depth_err']*1e6:.0f} ppm (n_in={d['n_in']})")
        ax.set_xlabel("Hours from mid-transit"); ax.set_xlim(-window_h, window_h)
    axes[0].set_ylabel("Relative flux")
    lo = np.nanpercentile(f[np.abs(x) < window_h], 0.3); hi = np.nanpercentile(f[np.abs(x) < window_h], 99.7)
    axes[0].set_ylim(lo, hi)
    fig.suptitle(f"{title}  (odd-even difference {oe['difference_sigma']:.2f} sigma)")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def fig_per_transit(pt, out, title):
    ev = pt["events"]
    if not ev:
        return
    tc = np.array([e["tc_bkjd"] for e in ev]); d = np.array([e["depth_ppm"] for e in ev]); er = np.array([e["depth_err_ppm"] for e in ev])
    fig, ax = plt.subplots(figsize=(13, 3.4))
    ax.errorbar(tc, d, er, fmt="o", ms=4, color=C["data"], lw=1, label="single-event box depth")
    s = pt["summary"]
    ax.axhline(s["weighted_mean_depth_ppm"], color=C["model"], lw=1, ls="--",
               label=f"weighted mean {s['weighted_mean_depth_ppm']:.0f} +/- {s['weighted_mean_err_ppm']:.0f} ppm")
    ax.axhline(0, color="#52514e", lw=0.6)
    ax.set_xlabel("Mid-transit time (BKJD)"); ax.set_ylabel("Depth (ppm)"); ax.set_title(title); ax.legend(loc="upper right")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
