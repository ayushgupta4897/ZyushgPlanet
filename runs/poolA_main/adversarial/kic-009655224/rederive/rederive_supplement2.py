#!/usr/bin/env python
"""Supplement 2: (a) bootstrap-over-epochs period uncertainty (replaces the spiky dlogL half-width
in supplement 1, which is a cadence-quantisation artefact); (b) masked re-detrend depth at the
known ephemeris (post-detection refinement, labelled as such); (c) absolute flux deficit in e/s in
SAP vs PSF; (d) split-half period consistency."""
import sys, json, os, numpy as np
from astropy.timeseries import BoxLeastSquares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rederive_kic009655224 as R
OUT = R.OUT
res = json.load(open(os.path.join(OUT, "rederive_results.json")))
P = res["ephemeris"]["period_d"]; t0 = res["ephemeris"]["t0_bjd"] - 2454833.0; DUR = res["ephemeris"]["duration_h"] / 24
S = {}
arr = np.load(os.path.join(OUT, "flat_lc.npy")); t, y, dy, ysap, dysap, q = arr
rng = np.random.default_rng(7)

# (a) bootstrap over transit epochs
n_ep = np.round((t - t0) / P).astype(int)
eps = np.unique(n_ep)
pgrid = np.linspace(P * (1 - 0.0004), P * (1 + 0.0004), 1601)   # step 3.4e-6 d
Ps = []
for b in range(120):
    pick = rng.choice(eps, size=len(eps), replace=True)
    idx = np.concatenate([np.where(n_ep == e)[0] for e in pick])
    bls = BoxLeastSquares(t[idx], y[idx], dy[idx])
    r = bls.power(pgrid, [DUR], objective="likelihood", oversample=10)
    Ps.append(pgrid[int(np.argmax(r.power))])
Ps = np.array(Ps)
S["period_bootstrap"] = dict(n_boot=len(Ps), median_d=float(np.median(Ps)), std_d=float(np.std(Ps)),
                             p16=float(np.percentile(Ps, 16)), p84=float(np.percentile(Ps, 84)),
                             pipeline_P_d=R.PIPE_P, pipeline_offset_sigma=float(abs(np.median(Ps) - R.PIPE_P) / np.std(Ps)),
                             main_refined_P_d=P)
print("bootstrap P:", S["period_bootstrap"])

# (d) split-half period
half = t < np.median(t)
sh = {}
for name, s in (("first_half", half), ("second_half", ~half)):
    bls = BoxLeastSquares(t[s], y[s], dy[s])
    pg = np.linspace(P * (1 - 0.003), P * (1 + 0.003), 3001)
    r = bls.power(pg, [DUR], objective="likelihood", oversample=10)
    i = int(np.argmax(r.power))
    d, nin = R.box_depth(t[s], y[s], P, t0, DUR); rp = R.random_phase_depths(t[s], y[s], P, t0, DUR, n=500)
    sh[name] = dict(best_P_d=float(pg[i]), dlogL=float(r.power[i]), depth_at_global_ephemeris_ppm=float(d), err_ppm=float(np.std(rp)), n_in=int(nin))
    print(name, sh[name])
S["split_half"] = sh

# (b) masked re-detrend: exclude |dt|<1.0*DUR from the running-median trend estimate (1.0 d window), then box depth
D = R.load(R.FITS)
ta = D["t"]
ph = R.phase_of(ta, P, t0); intr = np.abs(ph * P) < 1.0 * DUR
mk = {}
for win in (0.5, 1.0, 2.0):
    mask = ~intr
    # iterate outlier mask like detrend() but keep transit excluded from trend
    for _ in range(3):
        tr = R.running_median_time(ta, D["fn"], win, mask)
        r = D["fn"] / tr - 1
        mad = 1.4826 * np.median(np.abs(r[mask] - np.median(r[mask])))
        mask = (~intr) & (r < 3 * mad) & (r > -8 * mad)
    tr = R.running_median_time(ta, D["fn"], win, mask)
    flat = D["fn"] / tr; r = flat - 1
    keep = (r < 4 * mad) & (r > -8 * mad)
    d, nin = R.box_depth(ta[keep], flat[keep], P, t0, DUR)
    rp = R.random_phase_depths(ta[keep], flat[keep], P, t0, DUR, n=800)
    mk[f"window_{win}d"] = dict(depth_ppm=float(d), err_ppm=float(np.std(rp)), n_in=int(nin), snr=float(d / np.std(rp)))
    print(f"masked re-detrend window {win} d: depth {d:.0f} +/- {np.std(rp):.0f} ppm (SNR {d/np.std(rp):.1f})")
    if win == 1.0:
        # odd/even and secondary on the masked-detrend flux too
        tt, yy = ta[keep], flat[keep]
        ne = np.round((tt - t0) / P).astype(int); odd = ne % 2 == 1
        do, _ = R.box_depth(tt[odd], yy[odd], P, t0, DUR); de, _ = R.box_depth(tt[~odd], yy[~odd], P, t0, DUR)
        eo = np.std(R.random_phase_depths(tt[odd], yy[odd], P, t0, DUR, n=400)); ee = np.std(R.random_phase_depths(tt[~odd], yy[~odd], P, t0, DUR, n=400))
        sec, _ = R.box_depth(tt, yy, P, t0, DUR, phase_center=0.5); es = np.std(R.random_phase_depths(tt, yy, P, t0, DUR, n=800, avoid=[0.0, 0.5]))
        mk["window_1.0d"].update(odd_ppm=float(do), odd_err=float(eo), even_ppm=float(de), even_err=float(ee), odd_even_sigma=float(abs(do - de) / np.hypot(eo, ee)),
                                 secondary_ppm=float(sec), secondary_err=float(es), secondary_sigma=float(sec / es))
        # SAP on masked detrend
        mask_s = ~intr
        for _ in range(3):
            trs = R.running_median_time(ta, D["sapn"], win, mask_s); rs = D["sapn"] / trs - 1
            mads = 1.4826 * np.median(np.abs(rs[mask_s] - np.median(rs[mask_s]))); mask_s = (~intr) & (rs < 3 * mads) & (rs > -8 * mads)
        trs = R.running_median_time(ta, D["sapn"], win, mask_s); flats = D["sapn"] / trs; rs = flats - 1
        keeps = (rs < 4 * mads) & (rs > -8 * mads)
        ds, nins = R.box_depth(ta[keeps], flats[keeps], P, t0, DUR); rps = R.random_phase_depths(ta[keeps], flats[keeps], P, t0, DUR, n=800)
        mk["window_1.0d"].update(sap_depth_ppm=float(ds), sap_err_ppm=float(np.std(rps)), sap_over_psf=float(ds / d))
        print(f"  masked 1.0 d: odd {do:.0f}+/-{eo:.0f} even {de:.0f}+/-{ee:.0f} ({abs(do-de)/np.hypot(eo,ee):.2f} sig); sec {sec:.0f}+/-{es:.0f}; SAP {ds:.0f}+/-{np.std(rps):.0f} ratio {ds/d:.2f}")
S["masked_redetrend"] = mk

# (c) absolute flux deficit e/s
med_psf = res["data"]["median_psf_flux_e_s"]; med_sap = res["data"]["median_sap_flux_e_s"]
dp, ep = res["depth"]["psf_depth_ppm"], res["depth"]["psf_depth_err_ppm"]; ds_, es_ = res["sap"]["depth_ppm"], res["sap"]["err_ppm"]
dF_psf, dF_sap = dp * 1e-6 * med_psf, ds_ * 1e-6 * med_sap
S["absolute_deficit"] = dict(psf_e_s=float(dF_psf), psf_err=float(ep * 1e-6 * med_psf), sap_e_s=float(dF_sap), sap_err=float(es_ * 1e-6 * med_sap),
                             ratio_sap_over_psf=float(dF_sap / dF_psf), ratio_err=float(dF_sap / dF_psf * np.hypot(ep / dp, es_ / ds_)),
                             note="stitched FLUX and SAP_FLUX are both rescaled per quarter to common levels (see supplement raw_per_quarter); e/s here are in the stitched scale")
print("absolute deficit", S["absolute_deficit"])
json.dump(S, open(os.path.join(OUT, "rederive_supplement2.json"), "w"), indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
print("wrote supplement2")
