#!/usr/bin/env python
"""Stage 2 of the independent re-derivation: evaluate at the pipeline period 29.433 d,
check the two deep single events, and test sensitivity to detrending window / event masking.
Reuses only my own functions from rederive_kic005866176.py (no pipeline/ code)."""
import sys, json, time
sys.path.insert(0, "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005866176/rederive")
import numpy as np
from astropy.timeseries import BoxLeastSquares
import rederive_kic005866176 as R

OUT = R.OUT; log = []
def L(s): print(s); log.append(s)
T0 = time.time()
t, f, fe, qu = R.load("FLUX")

def evaluate_at(P0, label, t, f, fe, qu, frac=0.003):
    fl, keep = R.detrend(t, f); t1, f1, fe1 = t[keep], fl[keep], fe[keep]
    fine = np.linspace(P0*(1-frac), P0*(1+frac), 3001); fdurs = np.linspace(1.0, 6.0, 21)/24
    rf = BoxLeastSquares(t1, f1, fe1).power(fine, fdurs, objective="snr"); k = np.argmax(rf.power)
    P, t0, dur = float(fine[k]), float(rf.transit_time[k]), float(rf.duration[k])
    inm = np.abs(R.phase(t, P, t0)) < dur/P
    fl2, keep2 = R.detrend(t, f, mask=inm); t2, f2, qu2 = t[keep2], fl2[keep2], qu[keep2]
    d, ef, nin = R.box_depth(t2, f2, P, t0, dur); ee, _ = R.empirical_depth_err(t2, f2, P, t0, dur)
    n = np.round((t2-t0)/P).astype(int); odd = n%2==1
    do, _, _ = R.box_depth(t2[odd], f2[odd], P, t0, dur); de, _, _ = R.box_depth(t2[~odd], f2[~odd], P, t0, dur)
    eo, _ = R.empirical_depth_err(t2[odd], f2[odd], P, t0, dur, 300); ev, _ = R.empirical_depth_err(t2[~odd], f2[~odd], P, t0, dur, 300)
    ds, _, _ = R.box_depth(t2, f2, P, t0, dur, 0.5)
    eps = R.transit_epochs(t2, P, t0, dur); pev = []
    for ep in eps:
        tc = t0+ep*P; s = np.abs(t2-tc) < 3*dur; dd, e1, nn = R.box_depth(t2[s], f2[s], P, t0, dur)
        # gap proximity: nearest cadence gap > 0.5 d to this event (in days), using all good cadences
        gaps = np.where(np.diff(t) > 0.5)[0]; edges = np.concatenate([[t.min()], t[gaps], t[gaps+1], [t.max()]])
        pev.append(dict(epoch=int(ep), tc_bkjd=float(tc), depth_ppm=float(dd*1e6), err_ppm=float(e1*1e6), n_in=nn, quarter=int(np.median(qu2[s])), days_to_nearest_gap_or_edge=float(np.min(np.abs(edges-tc)))))
    yr_of = lambda q: 1 if q<=4 else 2 if q<=8 else 3 if q<=12 else 4
    py = []
    for Y in [1,2,3,4]:
        s = np.array([yr_of(q)==Y for q in qu2])
        if s.sum() < 100: py.append(dict(year=Y, depth_ppm=None)); continue
        dd, _, _ = R.box_depth(t2[s], f2[s], P, t0, dur); e2, _ = R.empirical_depth_err(t2[s], f2[s], P, t0, dur, 250)
        py.append(dict(year=Y, depth_ppm=float(dd*1e6), err_ppm=float(e2*1e6), n_transits=len(R.transit_epochs(t2[s], P, t0, dur))))
    # leave-one-out: depth without the deepest single event
    deps = np.array([p["depth_ppm"] for p in pev]); imax = int(np.argmax(deps))
    s_lo = np.abs(t2 - pev[imax]["tc_bkjd"]) > 3*dur
    d_lo, _, _ = R.box_depth(t2[s_lo], f2[s_lo], P, t0, dur); e_lo, _ = R.empirical_depth_err(t2[s_lo], f2[s_lo], P, t0, dur, 300)
    # two deepest removed
    order = np.argsort(deps)[::-1]; s_lo2 = s_lo & (np.abs(t2 - pev[order[1]]["tc_bkjd"]) > 3*dur)
    d_lo2, _, _ = R.box_depth(t2[s_lo2], f2[s_lo2], P, t0, dur); e_lo2, _ = R.empirical_depth_err(t2[s_lo2], f2[s_lo2], P, t0, dur, 300)
    res = dict(label=label, period_d=P, t0_bkjd=t0, t0_bjd=t0+2454833.0, duration_h=dur*24, bls_depth_ppm=float(rf.depth[k]*1e6), bls_snr_power=float(rf.power[k]),
               depth_ppm=float(d*1e6), depth_err_empirical_ppm=float(ee*1e6), empirical_snr=float(d/ee), n_in=nin, n_transits=len(eps),
               odd_ppm=float(do*1e6), odd_err=float(eo*1e6), even_ppm=float(de*1e6), even_err=float(ev*1e6), odd_even_sigma=float(abs(do-de)/np.hypot(eo,ev)),
               secondary_ppm=float(ds*1e6), secondary_sigma=float(ds/ee), per_event=pev, per_year=py,
               deepest_event=pev[imax], depth_without_deepest_ppm=float(d_lo*1e6), err_without_deepest_ppm=float(e_lo*1e6), snr_without_deepest=float(d_lo/e_lo),
               depth_without_two_deepest_ppm=float(d_lo2*1e6), snr_without_two_deepest=float(d_lo2/e_lo2),
               single_event_fraction=float(np.max(deps)*pev[imax]["n_in"]/max(np.sum(deps*np.array([p["n_in"] for p in pev])),1e-9)))
    L(f"[{label}] P={P:.5f} t0={t0:.4f} dur={dur*24:.2f}h depth={d*1e6:.0f}+-{ee*1e6:.0f} ppm SNR={d/ee:.2f} n_tr={len(eps)} odd/even {do*1e6:.0f}/{de*1e6:.0f} ({res['odd_even_sigma']:.1f}s) sec {ds*1e6:.0f} ({ds/ee:.1f}s)")
    L(f"   deepest event epoch {pev[imax]['epoch']} tc={pev[imax]['tc_bkjd']:.3f} depth {deps[imax]:.0f} ppm, {pev[imax]['days_to_nearest_gap_or_edge']:.2f} d from gap/edge; without it depth {d_lo*1e6:.0f}+-{e_lo*1e6:.0f} (SNR {d_lo/e_lo:.2f}); without two deepest {d_lo2*1e6:.0f} (SNR {d_lo2/e_lo2:.2f})")
    L("   per-year: " + ", ".join(f"Y{p['year']} {p['depth_ppm']:.0f}+-{p['err_ppm']:.0f} ({p['n_transits']} tr)" if p['depth_ppm'] is not None else f"Y{p['year']} n/a" for p in py))
    for p in pev:
        if abs(p['depth_ppm']) > 1000: L(f"   deep event: epoch {p['epoch']} Q{p['quarter']} tc={p['tc_bkjd']:.3f} depth={p['depth_ppm']:.0f}+-{p['err_ppm']:.0f} n_in={p['n_in']} gap_dist={p['days_to_nearest_gap_or_edge']:.2f} d")
    return res

at_pipe = evaluate_at(29.433, "pipeline_P_29.433", t, f, fe, qu)
at_top = evaluate_at(44.148, "blind_top_P_44.148", t, f, fe, qu)

# ---- look at the raw (un-detrended, normalised) flux around the two deep events
L("raw normalised PSF flux around deep events (+-1.5 d, 6-h bins):")
raw_ctx = {}
for p in at_top["per_event"]:
    if p["depth_ppm"] > 1000:
        tc = p["tc_bkjd"]; s = np.abs(t-tc) < 1.5
        bins = np.arange(-1.5, 1.51, 0.25); bi = np.digitize(t[s]-tc, bins)
        prof = [float(np.mean(f[s][bi==i])) if np.sum(bi==i) else None for i in range(1, len(bins))]
        raw_ctx[str(p["epoch"])] = dict(tc=tc, n_pts_pm1p5d=int(s.sum()), profile_6h_bins=prof)
        L(f"   tc={tc:.3f}: n={s.sum()} pts in +-1.5 d; 6h-binned raw flux: " + " ".join("--" if v is None else f"{(v-1)*1e6:+.0f}" for v in prof))

# ---- sensitivity: detrend window and masking of the two deep events, blind BLS 1-50 d each time
sens = []
durs = np.array([1.5, 2.5, 4.0, 6.0, 9.0])/24
deep_tcs = [p["tc_bkjd"] for p in at_top["per_event"] if p["depth_ppm"] > 1000]
for window, mask_deep in [(0.5, False), (2.0, False), (1.0, True)]:
    fl, keep = R.detrend(t, f, window=window)
    if mask_deep:
        for tc in deep_tcs: keep &= np.abs(t - tc) > 0.5
    t1, f1, fe1 = t[keep], fl[keep], fe[keep]
    periods, res = R.bls_scan(t1, f1, fe1, 1.0, 50.0, durs); S = R.sde(res.power)
    pk = R.top_peaks(periods, res.power, 3)
    win = np.abs(periods/29.433-1) < 0.003
    row = dict(window_d=window, deep_events_masked=mask_deep, n_pts=int(len(t1)),
               top3=[dict(period_d=float(periods[i]), sde=float(S[i]), depth_ppm=float(res.depth[i]*1e6), dur_h=float(res.duration[i]*24), alias=R.is_alias(periods[i], 29.433)) for i in pk],
               sde_near_29p433=float(S[win].max()))
    sens.append(row)
    L(f"sensitivity window={window} d masked_deep={mask_deep}: top3 " + "; ".join(f"{r['period_d']:.4f} d SDE {r['sde']:.2f} ({r['alias']})" for r in row["top3"]) + f"; SDE@29.433={row['sde_near_29p433']:.2f}")

# ---- inverted-flux null: full 1-50 d blind BLS on 2-f with the default detrend
fl, keep = R.detrend(t, f); t1, f1, fe1 = t[keep], fl[keep], fe[keep]
periods, resi = R.bls_scan(t1, 2-f1, fe1, 1.0, 50.0, durs); Si = R.sde(resi.power); pk = R.top_peaks(periods, resi.power, 3)
inv = [dict(period_d=float(periods[i]), sde=float(Si[i]), depth_ppm=float(resi.depth[i]*1e6)) for i in pk]
L("inverted-flux null blind BLS top3: " + "; ".join(f"{r['period_d']:.4f} d SDE {r['sde']:.2f} depth {r['depth_ppm']:.0f}" for r in inv))

out = dict(at_pipeline_period=at_pipe, at_blind_top_period=at_top, raw_context_deep_events=raw_ctx, sensitivity=sens, inverted_null_top3=inv, log=log, runtime_s=time.time()-T0)
json.dump(out, open(f"{OUT}/results_stage2.json", "w"), indent=1)
L(f"wrote {OUT}/results_stage2.json in {time.time()-T0:.0f}s")
