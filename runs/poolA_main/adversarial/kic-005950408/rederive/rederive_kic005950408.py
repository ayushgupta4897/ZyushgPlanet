#!/usr/bin/env python
"""
Independent re-derivation for kic-005950408 (KBONUS-BKG PSF light curve).

Written from scratch with numpy / scipy / astropy only. Does NOT import or read
anything under pipeline/.  Nothing here is a planet; this is a re-derivation of
a first-screen signal's period, depth, duration and consistency statistics.

Steps
 1. Load LIGHTCURVE_STITCHED, keep SAP_QUALITY==0 and finite values,
    normalise per quarter, detrend with a time-windowed running median (~1 d)
    with an upward-only sigma clip (own implementation, not wotan).
 2. Blind BLS 1-50 d (astropy BoxLeastSquares, own log-uniform frequency grid),
    on the real flux and on the inverted flux (2 - f) as a null.
 3. At the best ephemeris: box depth with empirical random-phase error,
    odd/even, secondary at phase 0.5, per-year depths, N transits, duration.
 4. Same depth measurement on SAP_FLUX (dilution check).
Outputs: results JSON + fold figure in the same directory as this script.
"""
import json, os, sys, time
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FNAME = "kic-005950408"
FITS = "/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005950408_kepler_v1.0_lc.fits"
OUTDIR = os.path.dirname(os.path.abspath(__file__))
RNG = np.random.default_rng(20261017)

# pipeline numbers given in the task text (used ONLY for comparison at the end)
PIPE_P, PIPE_DEPTH_PPM, PIPE_DUR_H = 38.60670, 4104.0, 0.9

# ----------------------------------------------------------------------------- 1
def load(col="FLUX"):
    with fits.open(FITS) as f:
        d = f["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float)
        y = np.asarray(d[col], float)
        e = np.asarray(d[col + "_ERR"], float)
        q = np.asarray(d["SAP_QUALITY"], int)
        qu = np.asarray(d["QUARTER"], int)
    m = (q == 0) & np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (y > 0)
    t, y, e, qu = t[m], y[m], e[m], qu[m]
    o = np.argsort(t)
    return t[o], y[o], e[o], qu[o]

def running_median(t, y, window=1.0):
    """Median of y over |t_j - t_i| <= window/2, for each i (sorted t)."""
    half = window / 2.0
    lo = np.searchsorted(t, t - half, side="left")
    hi = np.searchsorted(t, t + half, side="right")
    out = np.empty_like(y)
    for i in range(len(t)):
        out[i] = np.median(y[lo[i]:hi[i]])
    return out

def detrend(t, y, e, qu, window=1.0, nclip=3.0):
    """Per-quarter normalisation, then running-median flattening.
    Upward-only clipping (positive outliers) so dips are never removed;
    additionally, points more than 5 sigma *below* are kept but excluded
    from the trend estimate on the second pass (so the trend is not dragged
    down into deep events).  Returns normalised, flattened flux and errors."""
    yn, en = np.empty_like(y), np.empty_like(e)
    for Q in np.unique(qu):
        s = qu == Q
        med = np.median(y[s])
        yn[s] = y[s] / med
        en[s] = e[s] / med
    # pass 1
    trend = np.ones_like(yn)
    for Q in np.unique(qu):
        s = qu == Q
        trend[s] = running_median(t[s], yn[s], window)
    r = yn / trend - 1.0
    mad = 1.4826 * np.median(np.abs(r - np.median(r)))
    good = r < nclip * mad                 # drop flares / cosmic rays (upward)
    # pass 2: trend from points that are not strong outliers in either direction,
    # evaluated at all points (interpolated across excluded ones)
    use = good & (r > -5 * mad)
    trend2 = np.ones_like(yn)
    for Q in np.unique(qu):
        s = qu == Q
        su = s & use
        tr = running_median(t[su], yn[su], window)
        trend2[s] = np.interp(t[s], t[su], tr)
    f = yn / trend2
    ef = en / trend2
    keep = good
    return t[keep], f[keep], ef[keep], qu[keep], mad

# ----------------------------------------------------------------------------- 2
def make_grid(t, pmin=1.0, pmax=50.0, dur_min=0.04, oversample=2.0):
    T = t.max() - t.min()
    dlnf = dur_min / (oversample * T)      # phase drift over baseline <= dur_min/oversample
    lnf = np.arange(np.log(1.0 / pmax), np.log(1.0 / pmin), dlnf)
    return 1.0 / np.exp(lnf)[::-1]         # periods ascending

DURATIONS = np.array([0.03, 0.045, 0.06, 0.08, 0.11, 0.15, 0.20, 0.27])

def run_bls(t, f, ef, periods, durations=DURATIONS):
    bls = BoxLeastSquares(t, f, dy=ef)
    res = bls.power(periods, durations, objective="likelihood")
    p = np.asarray(res.power)
    sde = (p.max() - np.mean(p)) / np.std(p)
    return bls, res, p, sde

def top_peaks(periods, power, n=3, exclude_frac=0.02):
    """Top n peaks, each at least exclude_frac (fractional period) from the others."""
    idx = np.argsort(power)[::-1]
    peaks = []
    for i in idx:
        P = periods[i]
        if all(abs(P - q) / q > exclude_frac for q, _ in peaks):
            peaks.append((P, i))
        if len(peaks) == n:
            break
    return peaks

def alias_of(P, P0, tol=0.01):
    for k in (1, 2, 3, 4, 1/2, 1/3, 1/4, 3/2, 2/3):
        if abs(P / P0 - k) < tol:
            return k
    return None

# ----------------------------------------------------------------------------- 3
def phase_fold(t, P, t0):
    ph = ((t - t0) / P + 0.5) % 1.0 - 0.5
    return ph

def box_depth(t, f, P, t0, dur, phase_center=0.0, n_random=2000, rng=RNG):
    """Depth = median(out) - mean(in) for a box of width dur centred at phase_center.
    Empirical error: std of the same statistic for n_random random phases
    (excluding windows overlapping the primary and the tested box)."""
    ph = phase_fold(t, P, t0)
    hw = 0.5 * dur / P
    inb = np.abs(((ph - phase_center + 0.5) % 1.0) - 0.5) < hw
    oot = np.abs(ph) > 2 * hw
    if inb.sum() == 0:
        return np.nan, np.nan, 0
    base = np.median(f[oot])
    depth = base - np.mean(f[inb])
    # random-phase boxes
    rs = []
    cands = rng.uniform(-0.5, 0.5, n_random * 2)
    for c in cands:
        if abs(c) < 3 * hw or abs(((c - phase_center + 0.5) % 1.0) - 0.5) < 3 * hw:
            continue
        m = np.abs(((ph - c + 0.5) % 1.0) - 0.5) < hw
        if m.sum() < max(3, 0.3 * inb.sum()):
            continue
        rs.append(base - np.mean(f[m]))
        if len(rs) >= n_random:
            break
    rs = np.array(rs)
    err = np.std(rs) if len(rs) > 10 else np.nan
    return depth, err, int(inb.sum())

def per_event(t, f, P, t0, dur):
    """Per-transit local depth with local out-of-transit baseline (+-3 dur)."""
    n = np.round((t - t0) / P)
    ev = []
    for k in np.unique(n):
        tc = t0 + k * P
        inb = np.abs(t - tc) < dur / 2
        loc = (np.abs(t - tc) > dur) & (np.abs(t - tc) < 3 * dur)
        if inb.sum() >= 1 and loc.sum() >= 5:
            d = np.median(f[loc]) - np.mean(f[inb])
            scat = np.std(f[loc]) / np.sqrt(inb.sum())
            ev.append(dict(epoch=int(k), tc_bkjd=float(tc), n_in=int(inb.sum()),
                           n_loc=int(loc.sum()), depth_ppm=float(d * 1e6),
                           err_ppm=float(scat * 1e6)))
    return ev

def refine_ephemeris(t, f, ef, P0, dur0):
    """Fine BLS around P0 to refine P, t0 and duration."""
    bls = BoxLeastSquares(t, f, dy=ef)
    T = t.max() - t.min()
    dP = 0.02 * dur0 * P0 / T
    periods = np.arange(P0 * 0.995, P0 * 1.005, dP)
    durs = np.array([0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.07, 0.08, 0.10, 0.12])
    res = bls.power(periods, durs, objective="likelihood", oversample=20)
    i = int(np.argmax(res.power))
    return (float(res.period[i]), float(res.transit_time[i]), float(res.duration[i]),
            float(res.depth[i]), float(res.depth_err[i]), float(res.depth_snr[i]), bls, res)

def duration_scan(t, f, P, t0):
    """Chi2-like scan of box duration at fixed P, t0 (max depth*sqrt(n_in) ~ max S/N)."""
    durs = np.arange(0.01, 0.20, 0.0025)
    ph = phase_fold(t, P, t0) * P
    oot = np.abs(ph) > 0.3
    base, sig = np.median(f[oot]), np.std(f[oot])
    snr = []
    for d in durs:
        m = np.abs(ph) < d / 2
        if m.sum() < 2:
            snr.append(0); continue
        snr.append((base - f[m].mean()) * np.sqrt(m.sum()) / sig)
    snr = np.array(snr)
    return durs, snr, float(durs[np.argmax(snr)])

def year_of_quarter(q):
    if q <= 4: return 1
    if q <= 8: return 2
    if q <= 12: return 3
    return 4

# ----------------------------------------------------------------------------- main
def main():
    t_start = time.time()
    out = dict(fname=FNAME, fits=FITS, script=os.path.abspath(__file__),
               date="2026-10-17", method={})
    t, y, e, qu = load("FLUX")
    out["n_cadences_quality0"] = int(len(t))
    t, f, ef, qu, mad = detrend(t, y, e, qu, window=1.0)
    out["n_cadences_after_clip"] = int(len(t))
    out["point_scatter_mad_ppm"] = float(mad * 1e6)
    out["baseline_d"] = float(t.max() - t.min())
    out["method"]["detrend"] = ("per-quarter median normalisation; 1.0 d time-windowed running median, "
                                "two-pass (upward 3-MAD clip removed; >5-MAD dips excluded from trend)")
    print(f"loaded {len(t)} pts, MAD={mad*1e6:.0f} ppm, baseline {out['baseline_d']:.1f} d")

    # --- blind BLS, real and inverted
    periods = make_grid(t)
    out["method"]["bls_grid"] = dict(pmin=1.0, pmax=50.0, n_periods=int(len(periods)),
                                    rule="log-uniform in f, dlnf = 0.04 d / (2 * baseline)",
                                    durations_d=DURATIONS.tolist(), objective="likelihood")
    print(f"BLS grid: {len(periods)} periods")
    t1 = time.time()
    bls, res, pw, sde = run_bls(t, f, ef, periods)
    print(f"real BLS done in {time.time()-t1:.0f}s, SDE={sde:.2f}")
    peaks = top_peaks(periods, pw, 3)
    out["bls_real"] = dict(sde=float(sde), top3=[])
    for P, i in peaks:
        out["bls_real"]["top3"].append(dict(period_d=float(P), power=float(pw[i]),
            sde=float((pw[i] - pw.mean()) / pw.std()),
            depth_ppm=float(res.depth[i] * 1e6), duration_d=float(res.duration[i]),
            t0_bkjd=float(res.transit_time[i]),
            alias_of_pipeline_P=alias_of(P, PIPE_P)))
        print(f"  peak P={P:.5f}  pow={pw[i]:.1f}  depth={res.depth[i]*1e6:.0f}ppm dur={res.duration[i]:.3f} alias={alias_of(P, PIPE_P)}")
    Pbest, ibest = peaks[0]
    out["recovered_blind"] = bool(alias_of(Pbest, PIPE_P, tol=0.002) == 1)
    out["pipeline_P_in_top3_or_alias"] = any(p["alias_of_pipeline_P"] is not None for p in out["bls_real"]["top3"])

    t1 = time.time()
    _, res_i, pw_i, sde_i = run_bls(t, 2.0 - f, ef, periods)
    print(f"inverted BLS done in {time.time()-t1:.0f}s, SDE={sde_i:.2f}")
    pk_i = top_peaks(periods, pw_i, 3)
    out["bls_inverted_null"] = dict(sde=float(sde_i),
        top3=[dict(period_d=float(P), power=float(pw_i[i]),
                   sde=float((pw_i[i]-pw_i.mean())/pw_i.std()),
                   depth_ppm=float(res_i.depth[i]*1e6)) for P, i in pk_i])
    # power of the real spectrum exactly at the pipeline period (nearest grid point)
    j = int(np.argmin(np.abs(periods - PIPE_P)))
    out["bls_real"]["power_at_pipeline_P"] = dict(period_d=float(periods[j]), power=float(pw[j]),
                                                 sde=float((pw[j]-pw.mean())/pw.std()))

    # --- refine at the best peak
    P1, t01, dur1, dep1, deperr1, dsnr1, bls_f, res_f = refine_ephemeris(t, f, ef, Pbest, float(res.duration[ibest]))
    print(f"refined: P={P1:.6f} t0={t01:.4f} dur={dur1*24:.2f}h depth={dep1*1e6:.0f}+-{deperr1*1e6:.0f}ppm (formal BLS snr {dsnr1:.1f})")
    durs, dsnr, dur_scan_best = duration_scan(t, f, P1, t01)
    dur_use = dur1
    out["ephemeris"] = dict(period_d=P1, t0_bkjd=t01, t0_bjd=t01 + 2454833.0,
                            duration_d_bls=dur1, duration_h_bls=dur1 * 24,
                            duration_h_snr_scan=dur_scan_best * 24,
                            bls_formal_depth_ppm=dep1 * 1e6, bls_formal_depth_err_ppm=deperr1 * 1e6,
                            bls_formal_snr=dsnr1)

    # --- depth with empirical error
    depth, derr, nin = box_depth(t, f, P1, t01, dur_use)
    out["depth"] = dict(depth_ppm=depth * 1e6, depth_err_ppm_random_phase=derr * 1e6,
                        empirical_snr=depth / derr, n_in_transit_points=nin,
                        statistic="median(out) - mean(in), box width = BLS duration; error = std over random-phase boxes")
    print(f"box depth {depth*1e6:.0f} +- {derr*1e6:.0f} ppm (empirical SNR {depth/derr:.1f}), n_in={nin}")

    # --- odd / even
    ph = phase_fold(t, P1, t01)
    n = np.round((t - t01) / P1).astype(int)
    hw = 0.5 * dur_use / P1
    inb = np.abs(ph) < hw
    oot = np.abs(ph) > 2 * hw
    base = np.median(f[oot])
    sig = np.std(f[oot])
    odd = inb & (n % 2 == 1); even = inb & (n % 2 == 0)
    d_odd = base - f[odd].mean(); d_even = base - f[even].mean()
    e_odd = sig / np.sqrt(odd.sum()); e_even = sig / np.sqrt(even.sum())
    oe_sigma = abs(d_odd - d_even) / np.hypot(e_odd, e_even)
    out["odd_even"] = dict(depth_odd_ppm=d_odd * 1e6, err_odd_ppm=e_odd * 1e6, n_odd=int(odd.sum()),
                           depth_even_ppm=d_even * 1e6, err_even_ppm=e_even * 1e6, n_even=int(even.sum()),
                           difference_sigma=float(oe_sigma))
    print(f"odd {d_odd*1e6:.0f}+-{e_odd*1e6:.0f} (n={odd.sum()})  even {d_even*1e6:.0f}+-{e_even*1e6:.0f} (n={even.sum()})  -> {oe_sigma:.2f} sigma")

    # --- secondary at phase 0.5 (same box) and a wider box
    sec, sec_err, nsec = box_depth(t, f, P1, t01, dur_use, phase_center=0.5)
    sec_w, sec_w_err, nsec_w = box_depth(t, f, P1, t01, 3 * dur_use, phase_center=0.5)
    out["secondary"] = dict(depth_ppm=sec * 1e6, err_ppm=sec_err * 1e6, sigma=sec / sec_err, n_in=nsec,
                            wide_box_depth_ppm=sec_w * 1e6, wide_box_err_ppm=sec_w_err * 1e6,
                            wide_box_sigma=sec_w / sec_w_err, wide_box_n_in=nsec_w)
    print(f"secondary @0.5: {sec*1e6:.0f}+-{sec_err*1e6:.0f} ppm ({sec/sec_err:.2f} sigma); 3x box {sec_w*1e6:.0f}+-{sec_w_err*1e6:.0f}")

    # --- scan phase 0.1-0.9 for the deepest secondary anywhere
    scan_ph = np.arange(0.08, 0.93, 0.005)
    scan = [box_depth(t, f, P1, t01, dur_use, phase_center=c, n_random=200)[0] for c in scan_ph]
    k = int(np.nanargmax(scan))
    out["secondary"]["max_anywhere_phase"] = float(scan_ph[k])
    out["secondary"]["max_anywhere_depth_ppm"] = float(scan[k] * 1e6)

    # --- per-event and per-year
    ev = per_event(t, f, P1, t01, dur_use)
    out["events"] = ev
    out["n_transits_with_data"] = len(ev)
    out["n_epochs_in_baseline"] = int(np.floor((t.max() - t01) / P1) - np.ceil((t.min() - t01) / P1) + 1)
    yrs = {}
    yr = np.array([year_of_quarter(q) for q in qu])
    for Y in (1, 2, 3, 4):
        s = yr == Y
        if inb[s].sum() == 0:
            yrs[f"year{Y}"] = dict(depth_ppm=None, err_ppm=None, n_in=0); continue
        dY = np.median(f[s & oot]) - f[s & inb].mean()
        eY = np.std(f[s & oot]) / np.sqrt(inb[s].sum())
        yrs[f"year{Y}"] = dict(depth_ppm=float(dY * 1e6), err_ppm=float(eY * 1e6), n_in=int(inb[s].sum()),
                               quarters=sorted(int(x) for x in np.unique(qu[s])))
    out["per_year"] = yrs
    for Y, v in yrs.items():
        print(f"  {Y}: depth {v['depth_ppm']} +- {v['err_ppm']} ppm (n_in={v['n_in']})")
    # single-event dominance: drop the deepest event, re-measure
    if ev:
        dd = np.array([x["depth_ppm"] for x in ev]); ee = np.array([x["err_ppm"] for x in ev])
        kmax = int(np.argmax(dd / np.maximum(ee, 1)))
        drop = np.abs(t - ev[kmax]["tc_bkjd"]) > dur_use
        d_wo, e_wo, _ = box_depth(t[drop], f[drop], P1, t01, dur_use, n_random=800)
        out["single_event_check"] = dict(strongest_event_epoch=ev[kmax]["epoch"],
                                         strongest_event_depth_ppm=float(dd[kmax]),
                                         depth_without_it_ppm=d_wo * 1e6, err_without_it_ppm=e_wo * 1e6,
                                         snr_without_it=d_wo / e_wo,
                                         per_event_snr=sorted((dd / np.maximum(ee, 1)).round(2).tolist(), reverse=True))
        print(f"strongest event epoch {ev[kmax]['epoch']} ({dd[kmax]:.0f} ppm); depth without it {d_wo*1e6:.0f}+-{e_wo*1e6:.0f}")

    # --- V-shape proxy: compare depth in the central third vs the full box
    core = np.abs(ph) < hw / 3
    d_core = base - f[core].mean() if core.sum() > 0 else np.nan
    out["shape"] = dict(depth_core_third_ppm=float(d_core * 1e6), depth_full_box_ppm=float(depth * 1e6),
                        core_over_full=float(d_core / depth) if depth > 0 else None, n_core=int(core.sum()))

    # ------------------------------------------------------------------------- 4 SAP
    ts, ys, es, qus = load("SAP_FLUX")
    ts, fs, efs, qus, mad_s = detrend(ts, ys, es, qus, window=1.0)
    d_sap, e_sap, n_sap = box_depth(ts, fs, P1, t01, dur_use)
    with fits.open(FITS) as F:
        d = F["LIGHTCURVE_STITCHED"].data
        m = d["SAP_QUALITY"] == 0
        med_psf = float(np.nanmedian(d["FLUX"][m])); med_sap = float(np.nanmedian(d["SAP_FLUX"][m]))
    flfrcsap, psffrac = 0.7948817533, 0.5660038863112153   # from targets CSV row
    expected_ratio_if_on_target = flfrcsap * med_psf / med_sap
    out["sap_check"] = dict(sap_depth_ppm=d_sap * 1e6, sap_err_ppm=e_sap * 1e6, sap_snr=d_sap / e_sap,
                            sap_point_scatter_mad_ppm=float(mad_s * 1e6),
                            psf_depth_ppm=depth * 1e6, ratio_sap_over_psf=float(d_sap / depth),
                            median_psf_flux_e_s=med_psf, median_sap_flux_e_s=med_sap,
                            flfrcsap=flfrcsap, psffrac=psffrac,
                            expected_ratio_if_signal_on_target=float(expected_ratio_if_on_target),
                            note=("If the dip is on this star, SAP depth ~ depth_psf * flfrcsap*F_psf/F_sap "
                                  "(target flux captured by the 2-pixel aperture divided by total aperture flux). "
                                  "A SAP depth well above that expectation means the aperture's contaminating flux "
                                  "also dips, i.e. the signal is at least partly from a neighbour."))
    print(f"SAP depth {d_sap*1e6:.0f}+-{e_sap*1e6:.0f} ppm; ratio SAP/PSF = {d_sap/depth:.2f}; expected if on target ~{expected_ratio_if_on_target:.2f}")

    # ------------------------------------------------------------------------- comparison
    comp = dict(pipeline=dict(period_d=PIPE_P, depth_ppm=PIPE_DEPTH_PPM, duration_h=PIPE_DUR_H),
                mine=dict(period_d=P1, depth_ppm=depth * 1e6, duration_h=dur_use * 24),
                period_diff_d=P1 - PIPE_P, depth_ratio_mine_over_pipeline=depth * 1e6 / PIPE_DEPTH_PPM,
                duration_ratio=dur_use * 24 / PIPE_DUR_H)
    out["comparison_with_pipeline_task_numbers"] = comp
    out["runtime_s"] = time.time() - t_start

    with open(os.path.join(OUTDIR, "rederive_results.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else float(o))

    # ------------------------------------------------------------------------- figure
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    a = ax[0, 0]
    a.plot(periods, pw, lw=0.5, color="0.3")
    a.axvline(PIPE_P, color="tab:red", ls="--", lw=1, label=f"pipeline P={PIPE_P}")
    a.set_xscale("log"); a.set_xlabel("period [d]"); a.set_ylabel("BLS power (dlogL)")
    a.set_title(f"blind BLS 1-50 d, real flux: SDE={sde:.1f} (inverted null SDE={sde_i:.1f})")
    a.legend(fontsize=8)
    a = ax[0, 1]
    phd = ph * P1 * 24
    s = np.abs(phd) < 12
    a.plot(phd[s], (f[s] - 1) * 1e6, ".", ms=2, color="0.6", alpha=0.6)
    # binned
    bins = np.arange(-12, 12.01, 0.25)
    bc = 0.5 * (bins[1:] + bins[:-1])
    bm = np.array([np.mean(f[s][(phd[s] >= bins[i]) & (phd[s] < bins[i + 1])]) if np.any((phd[s] >= bins[i]) & (phd[s] < bins[i + 1])) else np.nan for i in range(len(bc))])
    a.plot(bc, (bm - 1) * 1e6, "o-", ms=4, color="tab:blue", label="15-min bins")
    a.axvspan(-dur_use * 12, dur_use * 12, color="tab:orange", alpha=0.15, label=f"box {dur_use*24:.2f} h")
    a.axhline(-depth * 1e6, color="tab:red", lw=1, label=f"depth {depth*1e6:.0f}+-{derr*1e6:.0f} ppm")
    a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm"); a.set_ylim(-12000, 6000)
    a.set_title(f"PSF flux folded at P={P1:.5f} d, t0={t01:.3f} BKJD"); a.legend(fontsize=8)
    a = ax[1, 0]
    for lab, m_, c in (("odd", odd, "tab:purple"), ("even", even, "tab:green")):
        a.plot(phd[m_], (f[m_] - 1) * 1e6, "o", ms=4, color=c, alpha=0.7, label=f"{lab} (n={m_.sum()})")
    a.axhline(-d_odd * 1e6, color="tab:purple", lw=1); a.axhline(-d_even * 1e6, color="tab:green", lw=1)
    a.set_xlim(-dur_use * 24, dur_use * 24); a.set_ylim(-15000, 6000)
    a.set_title(f"odd/even in-box points: {oe_sigma:.2f} sigma difference"); a.set_xlabel("hours"); a.legend(fontsize=8)
    a = ax[1, 1]
    if ev:
        ep = np.array([x["epoch"] for x in ev]); dd = np.array([x["depth_ppm"] for x in ev]); ee = np.array([x["err_ppm"] for x in ev])
        a.errorbar(ep, dd, ee, fmt="o", color="k", ms=4)
        a.axhline(depth * 1e6, color="tab:red", lw=1, label="global box depth")
        a.axhline(0, color="0.5", lw=0.5)
    a.set_xlabel("epoch"); a.set_ylabel("local depth [ppm]"); a.set_title(f"per-event local depths (N={len(ev)})"); a.legend(fontsize=8)
    fig.suptitle(f"{FNAME}  independent re-derivation (first-screen signal; not a planet)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "rederive_fold.png"), dpi=130)
    print("done in %.0fs" % out["runtime_s"])

if __name__ == "__main__":
    main()
