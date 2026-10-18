#!/usr/bin/env python
"""
Follow-up to rederive_kic005609180.py: robustness of the 21.1 d signal to
(i) gap-edge epochs, (ii) edge trimming before the blind BLS, (iii) choice of box
duration, and (iv) the SAP-vs-PSF dilution check per year and per event.
Imports only the functions defined in my own rederive script (no pipeline code).
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rederive_kic005609180 import (load, detrend, segments, run_bls, sde_of, top_peaks, phase_fold,
                                   box_depth, empirical_depth_error, transit_epochs, per_event_depth,
                                   refine_ephemeris, OUT_DIR, PIPELINE_P, PIPELINE_DEPTH_PPM, PIPELINE_DUR_H)

res0 = json.load(open(os.path.join(OUT_DIR, "rederive_results.json")))
P = res0["ephemeris"]["period_d"]; T0 = res0["ephemeris"]["t0_bkjd"]
out = {"based_on": "rederive_results.json", "P": P, "T0_bkjd": T0}

t, f, fe, q = load("FLUX")
ts, fs, fes, qs = load("SAP_FLUX")


def trim_edges(t, pad=0.5):
    keep = np.ones_like(t, bool)
    for (a, b) in segments(t):
        keep[a:b] &= (t[a:b] - t[a] > pad) & (t[b - 1] - t[a:b] > pad)
    return keep


def prep(t, f, fe, P, T0, dur, pad=None):
    """masked re-detrend (+ optional edge trim applied BEFORE detrending)"""
    if pad is not None:
        k = trim_edges(t, pad)
        t, f, fe = t[k], f[k], fe[k]
    ph = phase_fold(t, P, T0)
    inmask = np.abs(ph * P) < 0.75 * dur
    fd, fed, keep, mad = detrend(t, f, fe, 1.0, mask=inmask)
    keep = keep | inmask
    return t[keep], fd[keep], fed[keep]


# ------------------------------------------------------------------ A. raw flux around the three big "events"
big = [4, 36, 39]
normal = [59, 63]
fig, axes = plt.subplots(len(big) + len(normal), 2, figsize=(13, 2.6 * (len(big) + len(normal))))
for r, n in enumerate(big + normal):
    tm = T0 + n * P
    for c, (tt, ff, lab) in enumerate(((t, f, "PSF FLUX"), (ts, fs, "SAP_FLUX"))):
        ax = axes[r, c]
        s = np.abs(tt - tm) < 3.0
        ax.plot(tt[s] - tm, (ff[s] - 1) * 1e3, ".", ms=3)
        ax.axvspan(-PIPELINE_DUR_H / 48, PIPELINE_DUR_H / 48, color="C3", alpha=0.15)
        ax.set_title(f"epoch {n}  t_mid={tm:.2f} BKJD  {lab} (per-quarter normalised, raw)", fontsize=9)
        ax.set_xlabel("days from predicted mid-transit"); ax.set_ylabel("flux-1 [ppt]")
        ax.set_xlim(-3, 3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "events_raw_psf_vs_sap.png"), dpi=110)

# ------------------------------------------------------------------ B. gap-edge geometry of each epoch
segs = segments(t)
seg_bounds = [(t[a], t[b - 1]) for (a, b) in segs]
edge_info = []
for n in range(-1, 71):
    tm = T0 + n * P
    if tm < t.min() - 1 or tm > t.max() + 1:
        continue
    dists = [min(abs(tm - g0), abs(tm - g1)) for (g0, g1) in seg_bounds]
    k = int(np.argmin(dists))
    g0, g1 = seg_bounds[k]
    side = "start" if abs(tm - g0) < abs(tm - g1) else "end"
    edge_info.append({"epoch": n, "t_mid": float(tm), "d_edge": float(dists[k]), "edge_side": side})
out["epoch_edge_geometry"] = edge_info
near_edge_epochs = [e["epoch"] for e in edge_info if e["d_edge"] < 0.5]
out["epochs_within_0.5d_of_gap_edge"] = near_edge_epochs

# ------------------------------------------------------------------ C. depth with different treatments and box widths
treatments = {}
for dur_h in (3.0, 3.4, 8.0):
    dur = dur_h / 24.0
    row = {}
    # (1) all data, masked re-detrend, T0 re-optimised for this box
    t2, f2, fe2 = prep(t, f, fe, P, T0, dur)
    # re-centre T0 for this duration
    best = (-np.inf, T0)
    for T0b in T0 + np.linspace(-0.15, 0.15, 61):
        d, e, _ = box_depth(t2, f2, fe2, P, T0b, dur)
        if np.isfinite(d) and e > 0 and d / e > best[0]:
            best = (d / e, T0b)
    T0d = best[1]
    d, e, nin = box_depth(t2, f2, fe2, P, T0d, dur)
    s, _, _ = empirical_depth_error(t2, f2, fe2, P, T0d, dur, n=300)
    row["all"] = {"t0": float(T0d), "depth_ppm": float(d * 1e6), "emp_err_ppm": float(s * 1e6), "snr": float(d / s), "n_in": int(nin)}
    # (2) drop the epochs within 0.5 d of a gap edge
    drop = np.zeros_like(t2, bool)
    for n in near_edge_epochs:
        drop |= np.abs(t2 - (T0 + n * P)) < 3 * dur
    d2, _, nin2 = box_depth(t2[~drop], f2[~drop], fe2[~drop], P, T0d, dur)
    s2, _, _ = empirical_depth_error(t2[~drop], f2[~drop], fe2[~drop], P, T0d, dur, n=300)
    row["drop_gap_edge_epochs"] = {"depth_ppm": float(d2 * 1e6), "emp_err_ppm": float(s2 * 1e6), "snr": float(d2 / s2), "n_in": int(nin2),
                                   "epochs_dropped": near_edge_epochs}
    # (3) edge-trimmed (0.5 d) before detrending
    t3, f3, fe3 = prep(t, f, fe, P, T0, dur, pad=0.5)
    d3, _, nin3 = box_depth(t3, f3, fe3, P, T0d, dur)
    s3, _, _ = empirical_depth_error(t3, f3, fe3, P, T0d, dur, n=300)
    row["edge_trimmed_0.5d"] = {"depth_ppm": float(d3 * 1e6), "emp_err_ppm": float(s3 * 1e6), "snr": float(d3 / s3), "n_in": int(nin3)}
    # (4) drop only the two strongest epochs (4 and 39)
    drop2 = np.zeros_like(t2, bool)
    for n in (4, 39):
        drop2 |= np.abs(t2 - (T0 + n * P)) < 3 * dur
    d4, _, nin4 = box_depth(t2[~drop2], f2[~drop2], fe2[~drop2], P, T0d, dur)
    s4, _, _ = empirical_depth_error(t2[~drop2], f2[~drop2], fe2[~drop2], P, T0d, dur, n=300)
    row["drop_epochs_4_39"] = {"depth_ppm": float(d4 * 1e6), "emp_err_ppm": float(s4 * 1e6), "snr": float(d4 / s4), "n_in": int(nin4)}
    # (5) SAP at the same ephemeris / box
    ts2, fs2, fes2 = prep(ts, fs, fes, P, T0, dur)
    ds, _, nins = box_depth(ts2, fs2, fes2, P, T0d, dur)
    ss, _, _ = empirical_depth_error(ts2, fs2, fes2, P, T0d, dur, n=300)
    row["sap_all"] = {"depth_ppm": float(ds * 1e6), "emp_err_ppm": float(ss * 1e6), "snr": float(ds / ss), "n_in": int(nins)}
    # SAP with gap-edge epochs dropped
    drops = np.zeros_like(ts2, bool)
    for n in near_edge_epochs:
        drops |= np.abs(ts2 - (T0 + n * P)) < 3 * dur
    ds2, _, _ = box_depth(ts2[~drops], fs2[~drops], fes2[~drops], P, T0d, dur)
    ss2, _, _ = empirical_depth_error(ts2[~drops], fs2[~drops], fes2[~drops], P, T0d, dur, n=300)
    row["sap_drop_gap_edge_epochs"] = {"depth_ppm": float(ds2 * 1e6), "emp_err_ppm": float(ss2 * 1e6), "snr": float(ds2 / ss2)}
    # per-year PSF and SAP (all data)
    yr_of_q = {1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 2, 7: 2, 8: 2, 9: 3, 10: 3, 11: 3, 12: 3, 13: 4, 14: 4, 15: 4, 16: 4, 17: 4}
    yrs_psf = []
    yrs_sap = []
    for y in (1, 2, 3, 4):
        lo = {1: 131, 2: 460, 3: 820, 4: 1180}[y]; hi = {1: 460, 2: 820, 3: 1180, 4: 1600}[y]
        sy = (t2 >= lo) & (t2 < hi)
        dpy, _, _ = box_depth(t2[sy], f2[sy], fe2[sy], P, T0d, dur)
        epy, _, _ = empirical_depth_error(t2[sy], f2[sy], fe2[sy], P, T0d, dur, n=150)
        syd = sy & ~drop
        dpyd, _, _ = box_depth(t2[syd], f2[syd], fe2[syd], P, T0d, dur)
        epyd, _, _ = empirical_depth_error(t2[syd], f2[syd], fe2[syd], P, T0d, dur, n=150)
        yrs_psf.append({"year": y, "depth_ppm": float(dpy * 1e6), "err_ppm": float(epy * 1e6),
                        "depth_no_edge_ppm": float(dpyd * 1e6), "err_no_edge_ppm": float(epyd * 1e6)})
        ss_ = (ts2 >= lo) & (ts2 < hi)
        dsy, _, _ = box_depth(ts2[ss_], fs2[ss_], fes2[ss_], P, T0d, dur)
        esy, _, _ = empirical_depth_error(ts2[ss_], fs2[ss_], fes2[ss_], P, T0d, dur, n=150)
        yrs_sap.append({"year": y, "depth_ppm": float(dsy * 1e6), "err_ppm": float(esy * 1e6)})
    row["per_year_psf"] = yrs_psf
    row["per_year_sap"] = yrs_sap
    # per-event PSF vs SAP
    evs = []
    for (n, tm, nin_, nexp, nloc) in transit_epochs(t2, P, T0d, dur):
        dp, ep, _ = per_event_depth(t2, f2, fe2, tm, dur)
        dsp, esp, _ = per_event_depth(ts2, fs2, fes2, tm, dur)
        evs.append({"epoch": n, "psf_ppm": None if not np.isfinite(dp) else float(dp * 1e6),
                    "psf_err": None if not np.isfinite(ep) else float(ep * 1e6),
                    "sap_ppm": None if not np.isfinite(dsp) else float(dsp * 1e6),
                    "sap_err": None if not np.isfinite(esp) else float(esp * 1e6),
                    "cov": round(nin_ / nexp, 2)})
    row["per_event"] = evs
    treatments[f"box_{dur_h}h"] = row
    print(f"box {dur_h} h: all {d*1e6:.0f}+-{s*1e6:.0f} | no-edge-epochs {d2*1e6:.0f}+-{s2*1e6:.0f} | trimmed {d3*1e6:.0f}+-{s3*1e6:.0f} | drop4,39 {d4*1e6:.0f}+-{s4*1e6:.0f} | SAP {ds*1e6:.0f}+-{ss*1e6:.0f} | SAP no-edge {ds2*1e6:.0f}+-{ss2*1e6:.0f}", flush=True)
    print("   per-year PSF:", ", ".join(f"{y['depth_ppm']:.0f}+-{y['err_ppm']:.0f} (no-edge {y['depth_no_edge_ppm']:.0f})" for y in yrs_psf))
    print("   per-year SAP:", ", ".join(f"{y['depth_ppm']:.0f}+-{y['err_ppm']:.0f}" for y in yrs_sap))
out["treatments"] = treatments

# ------------------------------------------------------------------ D. blind BLS again on edge-trimmed data (0.5 d) and on SAP
k = trim_edges(t, 0.5)
fd, fed, keep, mad = detrend(t[k], f[k], fe[k], 1.0)
tt, ff, ffe = t[k][keep], fd[keep], fed[keep]
print(f"edge-trimmed PSF BLS: {len(tt)} points", flush=True)
_, periods, res, dt = run_bls(tt, ff, ffe)
pw = np.asarray(res.power)
pk = top_peaks(periods, pw, n=3)
out["bls_psf_edge_trimmed"] = {"n_points": int(len(tt)), "sde_top": float(sde_of(pw)),
                               "top3": [{"period_d": float(res.period[i]), "sde": float((p - np.median(pw)) / np.std(pw)),
                                         "depth_ppm": float(res.depth[i] * 1e6), "duration_h": float(res.duration[i] * 24)} for (P_, p, i) in pk]}
j = int(np.argmin(np.abs(periods - PIPELINE_P)))
jj = slice(max(j - 200, 0), j + 200)
jbest = jj.start + int(np.argmax(pw[jj]))
out["bls_psf_edge_trimmed"]["near_21.1d"] = {"period_d": float(periods[jbest]), "sde": float((pw[jbest] - np.median(pw)) / np.std(pw)),
                                             "depth_ppm": float(res.depth[jbest] * 1e6), "duration_h": float(res.duration[jbest] * 24),
                                             "rank_note": "local max within +-200 grid steps of 21.106 d"}
# rank of the 21.1 d peak in the trimmed periodogram
order = np.argsort(pw)[::-1]
rank = int(np.where(np.abs(periods[order] - periods[jbest]) < 1e-9)[0][0]) + 1 if np.any(np.abs(periods[order] - periods[jbest]) < 1e-9) else None
# count distinct peaks above it
pk_many = top_peaks(periods, pw, n=30)
n_above = sum(1 for (P_, p, i) in pk_many if p > pw[jbest] and abs(P_ - periods[jbest]) / periods[jbest] > 0.02)
out["bls_psf_edge_trimmed"]["near_21.1d"]["n_distinct_peaks_above"] = int(n_above)
print("edge-trimmed PSF top3:", [(round(x["period_d"], 4), round(x["sde"], 2)) for x in out["bls_psf_edge_trimmed"]["top3"]],
      " 21.1 d local:", out["bls_psf_edge_trimmed"]["near_21.1d"], flush=True)

fds, feds, keeps, mads = detrend(ts, fs, fes, 1.0)
_, periods_s, res_s, dt_s = run_bls(ts[keeps], fds[keeps], feds[keeps])
pws = np.asarray(res_s.power)
pks = top_peaks(periods_s, pws, n=3)
js = int(np.argmin(np.abs(periods_s - PIPELINE_P)))
jjs = slice(max(js - 200, 0), js + 200)
jbs = jjs.start + int(np.argmax(pws[jjs]))
out["bls_sap_blind"] = {"n_points": int(keeps.sum()), "sde_top": float(sde_of(pws)),
                        "top3": [{"period_d": float(res_s.period[i]), "sde": float((p - np.median(pws)) / np.std(pws)),
                                  "depth_ppm": float(res_s.depth[i] * 1e6), "duration_h": float(res_s.duration[i] * 24)} for (P_, p, i) in pks],
                        "near_21.1d": {"period_d": float(periods_s[jbs]), "sde": float((pws[jbs] - np.median(pws)) / np.std(pws)),
                                       "depth_ppm": float(res_s.depth[jbs] * 1e6)}}
print("SAP blind top3:", [(round(x["period_d"], 4), round(x["sde"], 2)) for x in out["bls_sap_blind"]["top3"]],
      " 21.1 d local:", out["bls_sap_blind"]["near_21.1d"], flush=True)

# ------------------------------------------------------------------ E. figure: periodograms + per-event PSF vs SAP
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
ax = axes[0]
ax.plot(periods, pw, lw=0.5, color="0.3"); ax.axvline(PIPELINE_P, color="C3", ls=":", label="21.106 d")
ax.set_xscale("log"); ax.set_title(f"PSF, edges trimmed 0.5 d: SDE_top={sde_of(pw):.1f}"); ax.set_xlabel("period [d]"); ax.legend(fontsize=7)
ax = axes[1]
ax.plot(periods_s, pws, lw=0.5, color="0.3"); ax.axvline(PIPELINE_P, color="C3", ls=":", label="21.106 d")
ax.set_xscale("log"); ax.set_title(f"SAP_FLUX blind BLS: SDE_top={sde_of(pws):.1f}"); ax.set_xlabel("period [d]"); ax.legend(fontsize=7)
ax = axes[2]
ev = treatments["box_3.4h"]["per_event"]
x = [e["psf_ppm"] for e in ev if e["psf_ppm"] is not None and e["sap_ppm"] is not None and e["cov"] >= 0.5]
y = [e["sap_ppm"] for e in ev if e["psf_ppm"] is not None and e["sap_ppm"] is not None and e["cov"] >= 0.5]
ax.plot(x, y, "o", ms=4); ax.plot([-1000, 8000], [-1000, 8000], "k:", lw=0.8, label="1:1")
ax.axhline(0, color="0.6", lw=0.5); ax.axvline(0, color="0.6", lw=0.5)
ax.set_xlabel("per-event PSF depth [ppm]"); ax.set_ylabel("per-event SAP depth [ppm]"); ax.set_title("3.4 h box, same epochs")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "followup_periodograms_sap_vs_psf.png"), dpi=110)


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


with open(os.path.join(OUT_DIR, "rederive_followup.json"), "w") as fh:
    json.dump(_clean(out), fh, indent=2)
print("done")
