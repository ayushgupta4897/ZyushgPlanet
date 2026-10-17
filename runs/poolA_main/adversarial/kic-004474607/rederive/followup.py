#!/usr/bin/env python
"""
Follow-up after comparing with the pipeline's vet.json (run AFTER results.json existed):
  * box depth at the pipeline's exact ephemeris (P=32.94272, t0=754.5747, dur=1.95 h) on MY
    detrending, for running-median windows 0.5 / 1.0 / 2.0 d, masked and unmasked, plus a
    detrend-free local-baseline estimate;
  * same at my ephemeris (P=32.98173, t0=754.113, dur=2.16 h);
  * gap / quarter-edge proximity of every covered event;
  * fine-peak period width (HWHM of delta log L) for a period uncertainty;
  * an event gallery (raw per-quarter-normalised flux around each covered event).
Writes comparison.json and events_gallery.png. Uses only my own module, never pipeline/.
"""
import sys, json
import numpy as np
sys.path.insert(0, "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-004474607/rederive")
import rederive_kic004474607 as R
from astropy.timeseries import BoxLeastSquares
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = R.OUT
res = json.load(open(f"{OUT}/results.json"))
PIPE = dict(P=32.94271710658947, t0=754.5746921888925, dur_h=1.95, depth_ppm=2140.678999968948,
            snr_emp=11.585667931163934, box_noise_ppm=184.7695801992392, n_in=64, n_events=17,
            odd_ppm=1960.0818450152824, even_ppm=2345.3557755830357, odd_even_sigma=1.0425794390838252,
            secondary_ppm=-12.088131066789742, secondary_sigma=-0.06542273383830266,
            event_weighted_depth_ppm=1025.2208677211386, event_chi2_const=45.30022853192067, event_dof=16,
            depth_by_quartile_ppm=[1042.9, 888.3, 1269.9, 490.5], tls_sde=5.676751608394032,
            tls_depth_ppm=(1 - 0.9983819211175828) * 1e6, tls_period=32.939176623992466,
            source="runs/poolA_main/vet/kic-004474607/vet.json (read after my results.json was written)")
MINE = dict(P=res["ephemeris"]["period_d"], t0=res["ephemeris"]["t0_bkjd"], dur_h=res["ephemeris"]["duration_box_h"])

t, f, e, q = R.load("FLUX")

def depth_at(P, t0, dur_d, window, masked, label):
    ph = ((t - t0) / P + 0.5) % 1.0 - 0.5
    intr = np.abs(ph * P) < 1.5 * dur_d
    flat, ef, tr, ok, mad = R.detrend(t, f, e, q, window=window, in_transit=intr if masked else None)
    d, n = R.box_depth(t[ok], flat[ok], P, t0, dur_d)
    err, bias, _ = R.random_phase_boxes(t[ok], flat[ok], P, t0, dur_d, n=1500, avoid=(0.0, 0.5), seed=7)
    return dict(label=label, window_d=window, masked=masked, depth_ppm=d * 1e6, err_ppm=err * 1e6, snr=d / err, n_in=n)

grid = []
for name, eph in [("pipeline", PIPE), ("mine", MINE)]:
    for w in (0.5, 1.0, 2.0):
        for masked in (False, True):
            grid.append(dict(ephemeris=name, **depth_at(eph["P"], eph["t0"], eph["dur_h"] / 24, w, masked, f"{name} w={w} masked={masked}")))
            print(grid[-1])

# detrend-free: per-quarter normalise only, local baseline +-1 d around each event
fn = f.copy()
for Q in np.unique(q):
    fn[q == Q] /= np.median(f[q == Q])
local = {}
for name, eph in [("pipeline", PIPE), ("mine", MINE)]:
    ev = R.per_event(t, fn, eph["P"], eph["t0"], eph["dur_h"] / 24)
    cov = [x for x in ev if x["covered"]]
    d = np.array([x["depth"] for x in cov]); er = np.array([x["err"] for x in cov]); w = 1 / er**2
    local[name] = dict(weighted_depth_ppm=float(np.sum(d * w) / np.sum(w) * 1e6), err_ppm=float(1 / np.sqrt(np.sum(w)) * 1e6),
                       n_events=len(cov), chi2_const=float(np.sum((d - np.sum(d * w) / np.sum(w))**2 * w)),
                       events=[dict(epoch=x["epoch"], tc=x["tc"], depth_ppm=x["depth"] * 1e6, err_ppm=x["err"] * 1e6, n_in=x["n_in"]) for x in cov])
    print(name, "local-baseline weighted depth", local[name]["weighted_depth_ppm"], "+-", local[name]["err_ppm"], "chi2", local[name]["chi2_const"])

# ---- gap proximity of each covered event (my ephemeris and pipeline's)
segs = R.segments(t, q)
edges = np.array(sorted(set([t[a] for a, b in segs] + [t[b - 1] for a, b in segs])))
qedges = np.array(sorted(set([t[q == Q].min() for Q in np.unique(q)] + [t[q == Q].max() for Q in np.unique(q)])))
gap = {}
for name, eph in [("pipeline", PIPE), ("mine", MINE)]:
    rows = []
    for x in local[name]["events"]:
        rows.append(dict(epoch=x["epoch"], tc=x["tc"], d_nearest_gap_d=float(np.min(np.abs(edges - x["tc"]))),
                         d_nearest_quarter_edge_d=float(np.min(np.abs(qedges - x["tc"]))), depth_ppm=x["depth_ppm"]))
    n_near = sum(r["d_nearest_gap_d"] < 0.5 for r in rows)
    gap[name] = dict(events=rows, n_within_0p5d_of_gap=n_near, n_events=len(rows))
    print(name, "events within 0.5 d of a gap/segment edge:", n_near, "/", len(rows))

# ---- fine peak width around my period (own detrend, single duration)
D = np.load(f"{OUT}/detrended_psf_masked.npz"); tt, yy, ye = D["t"], D["y"], D["ye"]
fine_P = np.linspace(MINE["P"] * (1 - 0.004), MINE["P"] * (1 + 0.004), 8001)
bls = BoxLeastSquares(tt, yy, dy=ye)
rf = bls.power(fine_P, [MINE["dur_h"] / 24], objective="likelihood", method="fast", oversample=20)
pwf = np.asarray(rf.power); j = int(np.argmax(pwf))
def width(drop):
    m = pwf >= pwf[j] - drop
    # contiguous region around j
    lo = j
    while lo > 0 and m[lo - 1]: lo -= 1
    hi = j
    while hi < len(m) - 1 and m[hi + 1]: hi += 1
    return float(0.5 * (fine_P[hi] - fine_P[lo]))
k_pipe = int(np.argmin(np.abs(fine_P - PIPE["P"])))
peak = dict(P_fine=float(fine_P[j]), power_max=float(pwf[j]), hwhm_dlogL_0p5_d=width(0.5), hwhm_dlogL_2_d=width(2.0),
            power_at_pipeline_P=float(pwf[k_pipe]), dlogL_mine_minus_pipeline=float(pwf[j] - pwf[k_pipe]),
            n_local_maxima_within_2_of_top=int(np.sum((pwf[1:-1] > pwf[:-2]) & (pwf[1:-1] > pwf[2:]) & (pwf[1:-1] > pwf[j] - 2))))
print("fine peak:", peak)

# ---- event gallery (my ephemeris), raw per-quarter-normalised flux
flat1, ef1, tr1, ok1, mad1 = R.detrend(t, f, e, q, window=1.0,
                                       in_transit=np.abs((((t - MINE["t0"]) / MINE["P"] + 0.5) % 1 - 0.5) * MINE["P"]) < 1.5 * MINE["dur_h"] / 24)
cov = local["mine"]["events"]
n = len(cov); ncol = 4; nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(16, 3.0 * nrow), sharey=True)
axes = axes.ravel()
for ax, x in zip(axes, cov):
    tc = x["tc"]; m = np.abs(t - tc) < 1.5
    ax.plot(t[m] - tc, (fn[m] - 1) * 1e6, ".", ms=3, color="#b8b7b1")
    ax.plot(t[m] - tc, (tr1[m] - 1) * 1e6, "-", lw=1.2, color="#2a78d6", label="running-median trend (1 d, masked)")
    ax.axvspan(-MINE["dur_h"] / 48, MINE["dur_h"] / 48, color="#e34948", alpha=0.15, label="my box")
    tp = PIPE["t0"] + PIPE["P"] * round((tc - PIPE["t0"]) / PIPE["P"])
    ax.axvline(tp - tc, color="#eb6834", ls="--", lw=1, label="pipeline tc")
    for qe in qedges:
        if abs(qe - tc) < 1.5: ax.axvline(qe - tc, color="k", ls=":", lw=1)
    ax.set_title(f"epoch {x['epoch']}  BKJD {tc:.2f}  local depth {x['depth_ppm']:.0f}+-{x['err_ppm']:.0f} ppm", fontsize=9)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-9000, 6000)
for ax in axes[n:]: ax.axis("off")
axes[0].legend(fontsize=7, loc="lower left")
fig.supxlabel("days from my mid-transit (P=%.5f d)" % MINE["P"]); fig.supylabel("ppm (per-quarter normalised, no detrend)")
fig.suptitle("kic-004474607 event gallery: raw flux around each covered epoch; dotted black = quarter edge", y=1.0)
fig.tight_layout(); fig.savefig(f"{OUT}/events_gallery.png", dpi=110, bbox_inches="tight")

comp = dict(pipeline_values=PIPE, my_values=dict(**MINE, depth_ppm=res["depth"]["depth_ppm"], depth_err_ppm=res["depth"]["depth_err_empirical_ppm"],
                                                 empirical_snr=res["depth"]["empirical_snr"], n_transits=res["transits"]["n_transits_covered"]),
            depth_vs_detrend_window=grid, local_baseline_no_detrend=local, gap_proximity=gap, fine_peak=peak,
            ephemeris_offset=dict(dP_d=MINE["P"] - PIPE["P"], dt0_d=MINE["t0"] - PIPE["t0"],
                                  epoch_where_ephemerides_cross=float(-(MINE["t0"] - PIPE["t0"]) / (MINE["P"] - PIPE["P"])),
                                  note="the two ephemerides coincide near epochs 11-13 (BKJD ~1150-1180) and drift apart by ~0.46 d at both ends"))
json.dump(comp, open(f"{OUT}/comparison.json", "w"), indent=2)
print("written comparison.json, events_gallery.png")
