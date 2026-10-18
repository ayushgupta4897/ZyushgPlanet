#!/usr/bin/env python
"""Stage 3: apples-to-apples at the pipeline's exact ephemeris, SAP depth at 29.43 d, final fold figure."""
import sys, json, time
sys.path.insert(0, "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005866176/rederive")
import numpy as np
import rederive_kic005866176 as R
OUT = R.OUT; log = []
def L(s): print(s); log.append(s)
Pp, t0p, durp = 29.433016127730333, 172.57296158611774, 1.95/24     # pipeline ephemeris (vet.json)
s2 = json.load(open(f"{OUT}/results_stage2.json")); mine = s2["at_pipeline_period"]
Pm, t0m, durm = mine["period_d"], mine["t0_bkjd"], mine["duration_h"]/24

res = {}
for col in ["FLUX", "SAP_FLUX"]:
    t, f, fe, qu = R.load(col)
    for lab, (P, t0, dur) in {"pipeline_ephemeris": (Pp, t0p, durp), "my_ephemeris": (Pm, t0m, durm)}.items():
        inm = np.abs(R.phase(t, P, t0)) < dur/P
        fl, keep = R.detrend(t, f, mask=inm); t2, f2, qu2 = t[keep], fl[keep], qu[keep]
        d, ef, nin = R.box_depth(t2, f2, P, t0, dur); ee, _ = R.empirical_depth_err(t2, f2, P, t0, dur)
        n = np.round((t2-t0)/P).astype(int); odd = n%2==1
        do, _, _ = R.box_depth(t2[odd], f2[odd], P, t0, dur); de, _, _ = R.box_depth(t2[~odd], f2[~odd], P, t0, dur)
        eo, _ = R.empirical_depth_err(t2[odd], f2[odd], P, t0, dur, 300); ev, _ = R.empirical_depth_err(t2[~odd], f2[~odd], P, t0, dur, 300)
        ds, _, _ = R.box_depth(t2, f2, P, t0, dur, 0.5)
        eps = R.transit_epochs(t2, P, t0, dur)
        # without the two gap-adjacent deep events
        far = (np.abs(t2-231.48) > 0.5) & (np.abs(t2-937.85) > 0.5)
        d_lo, _, _ = R.box_depth(t2[far], f2[far], P, t0, dur); e_lo, _ = R.empirical_depth_err(t2[far], f2[far], P, t0, dur, 300)
        yr_of = lambda q: 1 if q<=4 else 2 if q<=8 else 3 if q<=12 else 4
        py = []
        for Y in [1,2,3,4]:
            s = np.array([yr_of(q)==Y for q in qu2]); dd, _, _ = R.box_depth(t2[s], f2[s], P, t0, dur); e2, _ = R.empirical_depth_err(t2[s], f2[s], P, t0, dur, 250)
            py.append(dict(year=Y, depth_ppm=float(dd*1e6), err_ppm=float(e2*1e6)))
        res[f"{col}/{lab}"] = dict(P=P, t0_bkjd=t0, dur_h=dur*24, depth_ppm=float(d*1e6), err_emp_ppm=float(ee*1e6), snr=float(d/ee), n_in=nin, n_transits=len(eps),
                                   odd_ppm=float(do*1e6), even_ppm=float(de*1e6), odd_even_sigma=float(abs(do-de)/np.hypot(eo,ev)), secondary_ppm=float(ds*1e6), secondary_sigma=float(ds/ee),
                                   depth_without_gap_events_ppm=float(d_lo*1e6), err_without_gap_events_ppm=float(e_lo*1e6), snr_without_gap_events=float(d_lo/e_lo), per_year=py)
        L(f"{col:8s} {lab:18s} P={P:.5f} dur={dur*24:.2f}h depth={d*1e6:.0f}+-{ee*1e6:.0f} (SNR {d/ee:.2f}, n_in {nin}, {len(eps)} tr) odd/even {do*1e6:.0f}/{de*1e6:.0f} ({res[f'{col}/{lab}']['odd_even_sigma']:.1f}s) sec {ds*1e6:.0f} ({ds/ee:.1f}s); w/o gap events {d_lo*1e6:.0f}+-{e_lo*1e6:.0f} (SNR {d_lo/e_lo:.2f}); years " + " ".join(f"{p['depth_ppm']:.0f}+-{p['err_ppm']:.0f}" for p in py))
        if col == "FLUX" and lab == "pipeline_ephemeris":
            keepfig = (t2, f2, qu2, P, t0, dur, d, ee, do, de, eo, ev, ds, eps, d_lo, e_lo)
psf = res["FLUX/pipeline_ephemeris"]; sap = res["SAP_FLUX/pipeline_ephemeris"]
L(f"SAP/PSF depth ratio at pipeline ephemeris: {sap['depth_ppm']/psf['depth_ppm']:.2f}; at my ephemeris: {res['SAP_FLUX/my_ephemeris']['depth_ppm']/res['FLUX/my_ephemeris']['depth_ppm']:.2f}")
json.dump(dict(results=res, log=log), open(f"{OUT}/results_stage3.json", "w"), indent=1)

# ---------------- final figure at the pipeline period
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
t2, f2, qu2, P, t0, dur, d, ee, do, de, eo, ev, ds, eps, d_lo, e_lo = keepfig
pg = np.load(f"{OUT}/bls_periodogram.npy"); periods, S = pg[0], pg[2]
BLUE, ORANGE, GRAY, INK = "#2a78d6", "#eb6834", "#b8b7b2", "#52514e"
fig, ax = plt.subplots(2, 2, figsize=(12, 8.5), facecolor="#fcfcfb")
for a in ax.ravel():
    a.set_facecolor("#fcfcfb"); a.grid(True, color="#e8e7e3", lw=0.6); a.tick_params(colors=INK, labelsize=9)
    for s in a.spines.values(): s.set_color("#d8d7d2")
a = ax[0,0]; a.plot(periods, S, color=BLUE, lw=0.7, label="PSF flux, 1.0 d running-median detrend")
a.axvline(29.433, color=ORANGE, lw=1.0, ls="--", label="29.433 d (pipeline)"); a.axvline(44.144, color=INK, lw=0.8, ls=":", label="44.144 d = 3/2 x 29.433 (blind top peak)")
a.set_xscale("log"); a.set_xlabel("period [d]"); a.set_ylabel("BLS SDE (own grid, SNR objective)"); a.legend(fontsize=8, frameon=False, loc="upper left")
a.set_title(f"blind BLS 1-50 d: top peak SDE {S.max():.1f} at 44.14 d; SDE at 29.433 d = 6.8", fontsize=10, color=INK)
a = ax[0,1]; phh = R.phase(t2, P, t0)*P*24; s = np.abs(phh) < 12
a.plot(phh[s], (f2[s]-1)*1e6, ".", color=GRAY, ms=2, alpha=0.6, label="PSF flux, masked re-detrend")
bins = np.arange(-12, 12.01, 0.5); bi = np.digitize(phh[s], bins)
bm = np.array([np.mean(f2[s][bi==i]) if np.sum(bi==i)>3 else np.nan for i in range(1,len(bins))])
a.plot(0.5*(bins[1:]+bins[:-1]), (bm-1)*1e6, "-o", color=BLUE, ms=4, lw=1.5, label="30-min bins")
far = (np.abs(t2-231.48) > 0.5) & (np.abs(t2-937.85) > 0.5); bi2 = np.digitize(phh[s & far], bins)
bm2 = np.array([np.mean(f2[s & far][bi2==i]) if np.sum(bi2==i)>3 else np.nan for i in range(1,len(bins))])
a.plot(0.5*(bins[1:]+bins[:-1]), (bm2-1)*1e6, "-", color=ORANGE, lw=1.3, label="30-min bins, 2 gap-adjacent events removed")
a.axvspan(-dur*12, dur*12, color=ORANGE, alpha=0.10, lw=0)
a.set_ylim(-4000, 3000); a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm"); a.legend(fontsize=8, loc="lower left", frameon=False)
a.set_title(f"fold at P={P:.5f} d, {dur*24:.2f} h box: {d*1e6:.0f}±{ee*1e6:.0f} ppm (SNR {d/ee:.1f}); w/o 2 gap events {d_lo*1e6:.0f}±{e_lo*1e6:.0f} (SNR {d_lo/e_lo:.1f})", fontsize=9, color=INK)
a = ax[1,0]
n = np.round((t2-t0)/P).astype(int); odd = n%2==1
for sel, lab, col in [(odd, "odd", BLUE), (~odd, "even", ORANGE)]:
    pp = R.phase(t2[sel], P, t0)*P*24; ss = np.abs(pp) < 12; bi = np.digitize(pp[ss], bins)
    bm = np.array([np.mean(f2[sel][ss][bi==i]) if np.sum(bi==i)>3 else np.nan for i in range(1,len(bins))])
    a.plot(0.5*(bins[1:]+bins[:-1]), (bm-1)*1e6, "-o", color=col, ms=3.5, lw=1.3, label=lab)
ppsec = (((t2-t0)/P) % 1 - 0.5)*P*24; ss = np.abs(ppsec) < 12; bi = np.digitize(ppsec[ss], bins)
bm = np.array([np.mean(f2[ss][bi==i]) if np.sum(bi==i)>3 else np.nan for i in range(1,len(bins))])
a.plot(0.5*(bins[1:]+bins[:-1]), (bm-1)*1e6, "-", color=INK, lw=1.0, label="phase 0.5 (secondary)")
a.set_ylim(-1500, 800); a.set_xlabel("hours from mid-transit / mid-secondary"); a.set_ylabel("ppm (30-min bins)"); a.legend(fontsize=8, frameon=False)
a.set_title(f"odd {do*1e6:.0f}±{eo*1e6:.0f}, even {de*1e6:.0f}±{ev*1e6:.0f} ppm ({abs(do-de)/np.hypot(eo,ev):.1f}σ); secondary {ds*1e6:.0f} ppm ({ds/ee:.1f}σ)", fontsize=10, color=INK)
a = ax[1,1]
# raw context of the two deep events
tr, fr, _, _ = R.load("FLUX")
for tc, col, lab in [(231.48, BLUE, "event near BKJD 231.5 (Q2)"), (937.85, ORANGE, "event near BKJD 937.9 (Q10)")]:
    s = np.abs(tr-tc) < 2.0; a.plot((tr[s]-tc)*24, (fr[s]-1)*1e6, ".", color=col, ms=3, label=lab)
a.axvline(0, color=INK, lw=0.6); a.axvspan(-dur*12, dur*12, color=INK, alpha=0.08, lw=0)
a.set_xlabel("hours from predicted mid-transit"); a.set_ylabel("raw per-quarter-normalised PSF flux [ppm]"); a.legend(fontsize=8, frameon=False)
a.set_title("the two deepest 'transits' (~2-3 ppt) both start <0.1 d after a data gap (raw flux)", fontsize=10, color=INK)
fig.suptitle("kic-005866176 independent re-derivation (own code, no pipeline/ imports) — first-screen signal, not a planet", fontsize=12, color="#0b0b0b")
fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig(f"{OUT}/fold_rederive.png", dpi=130)
L(f"wrote {OUT}/fold_rederive.png and results_stage3.json")
