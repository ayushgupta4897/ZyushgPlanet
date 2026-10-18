#!/usr/bin/env python
"""Independent re-derivation for kic-005866176 (KBONUS-BKG PSF light curve).

Written from scratch (numpy / scipy / astropy only). Does NOT import anything from
pipeline/ and does not use wotan. Steps:
  1. load FITS, SAP_QUALITY==0, finite, per-quarter median normalisation,
     time-windowed running-median detrend (own implementation, 1.0 d window,
     gap-aware, upper-outlier clipped, then a masked re-detrend at the final ephemeris)
  2. blind BLS 1-50 d (astropy BoxLeastSquares, own frequency grid), top-3 peaks, SDE
  3. depth with random-phase empirical error, odd/even, secondary (phase 0.5),
     per-Kepler-year depths, n_transits, duration
  4. same depth on SAP_FLUX (dilution check)
  5. results JSON + one fold figure
"""
import json, sys, os, time
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares

FITS = "/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005866176_kepler_v1.0_lc.fits"
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-005866176/rederive"
PIPE_P, PIPE_DEPTH, PIPE_DUR = 29.43300, 481.0, 2.2
RNG = np.random.default_rng(20261017)
WINDOW_D = 1.0

# ---------------------------------------------------------------- loading
def load(col):
    with fits.open(FITS) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float); f = np.asarray(d[col], float)
        fe = np.asarray(d[col + "_ERR"], float); q = np.asarray(d["SAP_QUALITY"]); qu = np.asarray(d["QUARTER"])
    m = (q == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (f > 0)
    t, f, fe, qu = t[m], f[m], fe[m], qu[m]
    o = np.argsort(t); t, f, fe, qu = t[o], f[o], fe[o], qu[o]
    # per-quarter normalisation by median
    fn = np.empty_like(f); fen = np.empty_like(f)
    for Q in np.unique(qu):
        s = qu == Q; med = np.median(f[s]); fn[s] = f[s] / med; fen[s] = fe[s] / med
    return t, fn, fen, qu

def running_median(t, f, window, mask=None, gap=0.5):
    """time-windowed running median, computed separately on gap-split segments.
    mask: boolean, True = exclude from the trend estimate (but trend still evaluated there)."""
    trend = np.full_like(f, np.nan)
    use = np.ones_like(f, bool) if mask is None else ~mask
    # split on gaps
    brk = np.where(np.diff(t) > gap)[0] + 1
    for seg in np.split(np.arange(len(t)), brk):
        ts, fs, us = t[seg], f[seg], use[seg]
        tu, fu = ts[us], fs[us]
        if len(tu) < 10:
            trend[seg] = np.median(fs) if len(fs) else np.nan; continue
        lo = np.searchsorted(tu, ts - window / 2); hi = np.searchsorted(tu, ts + window / 2)
        tr = np.empty(len(ts))
        for i in range(len(ts)):
            a, b = lo[i], hi[i]
            if b - a < 5:  # widen at edges / sparse regions
                a = max(0, a - 5); b = min(len(tu), b + 5)
            tr[i] = np.median(fu[a:b])
        trend[seg] = tr
    return trend

def detrend(t, f, mask=None, window=WINDOW_D):
    # pass 1: trend, clip upward (flare/outlier) and strong downward outliers from the trend estimate
    tr = running_median(t, f, window, mask)
    r = f / tr - 1; mad = 1.4826 * np.median(np.abs(r - np.median(r)))
    bad = (r > 3 * mad) | (r < -6 * mad)
    m2 = bad if mask is None else (bad | mask)
    tr = running_median(t, f, window, m2)
    fl = f / tr
    r = fl - 1; mad = 1.4826 * np.median(np.abs(r - np.median(r)))
    keep = (r < 4 * mad) & (r > -10 * mad)   # drop upward outliers only aggressively; keep dips
    return fl, keep

# ---------------------------------------------------------------- transit tools
def phase(t, P, t0):
    return ((t - t0) / P + 0.5) % 1.0 - 0.5  # in [-0.5,0.5), transit at 0

def box_depth(t, f, P, t0, dur, ph_center=0.0):
    ph = ((t - t0) / P - ph_center + 0.5) % 1.0 - 0.5
    intr = np.abs(ph) < 0.5 * dur / P
    oot = (np.abs(ph) > 1.0 * dur / P) & (np.abs(ph) < 3.0 * dur / P)   # local out-of-transit annulus
    if intr.sum() < 3: return np.nan, np.nan, int(intr.sum())
    d = (np.median(f[oot]) if oot.sum() > 10 else 1.0) - np.mean(f[intr])
    e = np.std(f[oot]) / np.sqrt(intr.sum()) if oot.sum() > 10 else np.nan
    return d, e, int(intr.sum())

def empirical_depth_err(t, f, P, t0, dur, n=600):
    """depth measured in boxes of the same duration at random phases, excluding +-2 dur of transit & secondary"""
    vals = []
    while len(vals) < n:
        pc = RNG.uniform(-0.5, 0.5)
        if abs(pc) < 2 * dur / P or abs(abs(pc) - 0.5) < 2 * dur / P: continue
        d, _, nin = box_depth(t, f, P, t0, dur, pc)
        if np.isfinite(d) and nin >= 3: vals.append(d)
    vals = np.array(vals)
    return float(np.std(vals)), float(np.median(vals))

def transit_epochs(t, P, t0, dur):
    n = np.round((t - t0) / P).astype(int)
    ph = t - (t0 + n * P)
    intr = np.abs(ph) < 0.5 * dur
    eps = np.unique(n[intr])
    out = []
    for e in eps:
        s = intr & (n == e)
        if s.sum() >= 3: out.append(int(e))
    return out

def sde(power):
    return (power - np.mean(power)) / np.std(power)

def bls_scan(t, f, fe, pmin, pmax, durations_d, drift_frac=1/3.0, oversample=1.0):
    base = t.max() - t.min()
    df = (min(durations_d) * drift_frac) / base / oversample
    freqs = np.arange(1.0 / pmax, 1.0 / pmin, df)
    periods = 1.0 / freqs
    bls = BoxLeastSquares(t, f, fe)
    res = bls.power(periods, durations_d, objective="snr")
    return periods, res

def top_peaks(periods, power, k=3, excl=0.02):
    """top k distinct peaks, excluding +-excl fractional neighbourhood of chosen ones"""
    idx = np.argsort(power)[::-1]; chosen = []
    for i in idx:
        if all(abs(periods[i] / periods[j] - 1) > excl for j in chosen):
            chosen.append(i)
        if len(chosen) == k: break
    return chosen

def is_alias(P, Pref, tol=0.01):
    for a, b in [(1,1),(1,2),(2,1),(1,3),(3,1),(2,3),(3,2),(1,4),(4,1)]:
        if abs(P * a / b / Pref - 1) < tol: return f"{a}/{b}"
    return None

# ================================================================ main
def main():
    T0 = time.time()
    log = []
    def L(s): print(s); log.append(s)

    t, f, fe, qu = load("FLUX")
    L(f"loaded PSF FLUX: {len(t)} good cadences, {t.min():.2f}-{t.max():.2f} BKJD, quarters {sorted(set(qu.tolist()))}")
    fl, keep = detrend(t, f)
    t1, f1, fe1, qu1 = t[keep], fl[keep], fe[keep] , qu[keep]
    fe1 = fe1 / 1.0
    rms = 1.4826 * np.median(np.abs(f1 - np.median(f1)))
    L(f"detrended (running median {WINDOW_D} d): {len(t1)} pts, robust point-to-point scatter {rms*1e6:.0f} ppm")

    # ---------------- 2. blind BLS 1-50 d
    durs_h = np.array([1.5, 2.5, 4.0, 6.0, 9.0]); durs = durs_h / 24
    periods, res = bls_scan(t1, f1, fe1, 1.0, 50.0, durs)
    pw = res.power; S = sde(pw)
    L(f"BLS grid: {len(periods)} periods, 1-50 d, durations {durs_h.tolist()} h; {time.time()-T0:.0f}s")
    pk = top_peaks(periods, pw, k=3)
    peaks = []
    for r, i in enumerate(pk):
        peaks.append(dict(rank=r+1, period_d=float(periods[i]), sde=float(S[i]), snr_power=float(pw[i]),
                          depth_ppm=float(res.depth[i]*1e6), duration_h=float(res.duration[i]*24),
                          t0_bkjd=float(res.transit_time[i]), alias_of_29p433=is_alias(periods[i], PIPE_P)))
        L(f"  peak {r+1}: P={periods[i]:.5f} d SDE={S[i]:.2f} depth={res.depth[i]*1e6:.0f} ppm dur={res.duration[i]*24:.2f} h alias={peaks[-1]['alias_of_29p433']}")
    # SDE at the pipeline period specifically
    j = np.argmin(np.abs(periods - PIPE_P)); win = np.abs(periods/PIPE_P - 1) < 0.003
    L(f"  SDE near 29.433 d (best within 0.3%): {S[win].max():.2f} at P={periods[win][np.argmax(S[win])]:.5f}")
    recovered = peaks[0]["alias_of_29p433"] is not None
    # also list the top peaks in SDE excluding 1/1 of 29.433 window for context
    # ---------------- refine top peak (or pipeline peak if top is an alias)
    Pbest0 = peaks[0]["period_d"]
    fine = np.linspace(Pbest0 * (1 - 0.002), Pbest0 * (1 + 0.002), 4001)
    fdurs = np.linspace(1.0, 5.0, 17) / 24
    bls = BoxLeastSquares(t1, f1, fe1); rf = bls.power(fine, fdurs, objective="snr")
    k = np.argmax(rf.power)
    P, t0, dur = float(fine[k]), float(rf.transit_time[k]), float(rf.duration[k])
    L(f"refined: P={P:.5f} d t0={t0:.4f} BKJD dur={dur*24:.2f} h BLS depth={rf.depth[k]*1e6:.0f}+-{rf.depth_err[k]*1e6:.0f} ppm (formal)")
    stats = bls.compute_stats(P, dur, t0)

    # ---------------- 3. masked re-detrend then depth statistics
    ph = phase(t, P, t0); inmask = np.abs(ph) < 1.0 * dur / P   # mask +-1 duration around transit
    fl2, keep2 = detrend(t, f, mask=inmask)
    t2, f2, qu2 = t[keep2], fl2[keep2], qu[keep2]
    d, e_formal, nin = box_depth(t2, f2, P, t0, dur)
    e_emp, med_rand = empirical_depth_err(t2, f2, P, t0, dur)
    L(f"masked re-detrend depth: {d*1e6:.0f} ppm, formal err {e_formal*1e6:.0f}, empirical (random-phase box) err {e_emp*1e6:.0f} ppm -> empirical SNR {d/e_emp:.2f}; n_in={nin}")
    # odd / even
    n = np.round((t2 - t0) / P).astype(int)
    odd = n % 2 == 1; even = ~odd
    d_odd, e_odd, _ = box_depth(t2[odd], f2[odd], P, t0, dur); d_even, e_even, _ = box_depth(t2[even], f2[even], P, t0, dur)
    e_odd_emp, _ = empirical_depth_err(t2[odd], f2[odd], P, t0, dur, n=300)
    e_even_emp, _ = empirical_depth_err(t2[even], f2[even], P, t0, dur, n=300)
    oe_sig = abs(d_odd - d_even) / np.hypot(e_odd_emp, e_even_emp)
    L(f"odd {d_odd*1e6:.0f}+-{e_odd_emp*1e6:.0f}  even {d_even*1e6:.0f}+-{e_even_emp*1e6:.0f} ppm  -> {oe_sig:.2f} sigma")
    # secondary at phase 0.5
    d_sec, _, nsec = box_depth(t2, f2, P, t0, dur, 0.5)
    sec_sig = d_sec / e_emp
    L(f"secondary (phase 0.5, same box): {d_sec*1e6:.0f} ppm = {sec_sig:.2f} sigma (n_in={nsec})")
    # scan of box depth over phase (for the figure + the strongest non-primary dip)
    pgrid = np.linspace(-0.5, 0.5, 400, endpoint=False)
    dscan = np.array([box_depth(t2, f2, P, t0, dur, pc)[0] for pc in pgrid])
    far = np.abs(pgrid) > 2 * dur / P
    imax = np.argmax(np.where(far, dscan, -np.inf))
    L(f"strongest other dip in phase: {dscan[imax]*1e6:.0f} ppm at phase {pgrid[imax]:.3f} ({dscan[imax]/e_emp:.2f} sigma)")
    # per-event
    eps = transit_epochs(t2, P, t0, dur); per_event = []
    for ep in eps:
        tc = t0 + ep * P; s = np.abs(t2 - tc) < 3 * dur
        dd, ee, nn = box_depth(t2[s], f2[s], P, t0, dur)
        per_event.append(dict(epoch=ep, tc_bkjd=float(tc), depth_ppm=float(dd*1e6), err_ppm=float(ee*1e6) if np.isfinite(ee) else None, n_in=nn, quarter=int(np.median(qu2[s]))))
    L(f"n_transits with >=3 in-transit points: {len(eps)}; epochs {eps}")
    for pe in per_event: L(f"   epoch {pe['epoch']:3d} Q{pe['quarter']:2d} tc={pe['tc_bkjd']:.3f} depth={pe['depth_ppm']:.0f}+-{pe['err_ppm'] if pe['err_ppm'] is None else round(pe['err_ppm'])} n={pe['n_in']}")
    # per Kepler year (by quarter block: Y1 Q0-4, Y2 Q5-8, Y3 Q9-12, Y4 Q13-17)
    yr_of = lambda q: 1 if q <= 4 else 2 if q <= 8 else 3 if q <= 12 else 4
    per_year = []
    for Y in [1,2,3,4]:
        s = np.array([yr_of(q) == Y for q in qu2])
        if s.sum() < 100: per_year.append(dict(year=Y, depth_ppm=None, err_ppm=None, n_transits=0)); continue
        dd, _, nn = box_depth(t2[s], f2[s], P, t0, dur)
        ee, _ = empirical_depth_err(t2[s], f2[s], P, t0, dur, n=250)
        ne = len(transit_epochs(t2[s], P, t0, dur))
        per_year.append(dict(year=Y, quarters=sorted(set(qu2[s].tolist())), depth_ppm=float(dd*1e6), err_ppm=float(ee*1e6), n_transits=ne, n_in=nn))
        L(f"year {Y} (Q{per_year[-1]['quarters']}): depth {dd*1e6:.0f}+-{ee*1e6:.0f} ppm, {ne} transits")
    # transit shape: trapezoid-ish check -- depth in the central half vs full box
    phc = phase(t2, P, t0); cen = np.abs(phc) < 0.25 * dur / P
    d_cen = np.median(f2[(np.abs(phc) > dur/P) & (np.abs(phc) < 3*dur/P)]) - np.mean(f2[cen])
    L(f"central-half depth {d_cen*1e6:.0f} ppm vs full-box {d*1e6:.0f} ppm")

    # ---------------- 4. SAP_FLUX
    ts, fs, fes, qus = load("SAP_FLUX")
    phs = phase(ts, P, t0); ms = np.abs(phs) < 1.0 * dur / P
    fls, keeps = detrend(ts, fs, mask=ms)
    ts2, fs2 = ts[keeps], fls[keeps]
    d_sap, _, nin_sap = box_depth(ts2, fs2, P, t0, dur)
    e_sap, _ = empirical_depth_err(ts2, fs2, P, t0, dur)
    rms_sap = 1.4826 * np.median(np.abs(fs2 - np.median(fs2)))
    L(f"SAP_FLUX: depth {d_sap*1e6:.0f}+-{e_sap*1e6:.0f} ppm (SNR {d_sap/e_sap:.2f}), scatter {rms_sap*1e6:.0f} ppm, n_in={nin_sap}; ratio SAP/PSF = {d_sap/d:.2f}")
    # BLS SDE in SAP at same P (narrow window) for completeness
    fine_s = np.linspace(P*(1-0.002), P*(1+0.002), 801)
    rs = BoxLeastSquares(ts2, fs2, fes[keeps]).power(fine_s, fdurs, objective="snr")
    L(f"SAP narrow-window BLS best: P={fine_s[np.argmax(rs.power)]:.5f} depth={rs.depth[np.argmax(rs.power)]*1e6:.0f} ppm")

    # ---------------- inverted-flux control on the same grid near P (cheap sanity)
    rinv = BoxLeastSquares(t1, 2 - f1, fe1).power(fine, fdurs, objective="snr")
    L(f"inverted flux narrow-window max SNR-power {rinv.power.max():.2f} vs real {rf.power.max():.2f}")

    # ---------------- physical plausibility
    Rstar = 0.9458  # Rsun, from targets csv
    rp_re = np.sqrt(max(d,0)) * Rstar * 109.1
    # expected duration for central transit, M~R for a ~0.95 Rsun dwarf
    a_au = (P/365.25)**(2/3) * (Rstar)**(1/3)
    T_central_h = (P/np.pi) * (Rstar*0.00465/a_au) * 24
    L(f"implied Rp ~ {rp_re:.2f} Re (undiluted, Rstar={Rstar}); central-transit duration expectation ~{T_central_h:.1f} h vs measured {dur*24:.2f} h")

    results = dict(
        fname="kic-005866176", kic=5866176, run_label="adversarial_rederive_2026-10-17",
        method=dict(quality="SAP_QUALITY==0", normalisation="per-quarter median", detrend=f"time-windowed running median, {WINDOW_D} d, gap-split, outlier-clipped; masked re-detrend (+-1 dur) for depth stats",
                    bls="astropy BoxLeastSquares objective=snr, 1-50 d, df=(1.5h/3)/baseline, durations 1.5-9 h; refinement +-0.2% with 1-5 h durations",
                    depth="box mean in-transit vs median of local annulus (1-3 dur); error = std of same-size boxes at 600 random phases"),
        n_points_good=int(len(t)), n_points_detrended=int(len(t1)), scatter_ppm=float(rms*1e6),
        blind_bls=dict(n_periods=int(len(periods)), top_peaks=peaks, sde_at_29p433=float(S[win].max()), recovered_blind=bool(recovered)),
        best=dict(period_d=P, t0_bkjd=t0, t0_bjd=t0 + 2454833.0, duration_h=dur*24, bls_depth_ppm=float(rf.depth[k]*1e6), bls_depth_err_formal_ppm=float(rf.depth_err[k]*1e6),
                  depth_ppm=float(d*1e6), depth_err_empirical_ppm=float(e_emp*1e6), depth_err_formal_ppm=float(e_formal*1e6), empirical_snr=float(d/e_emp),
                  random_phase_median_ppm=float(med_rand*1e6), central_half_depth_ppm=float(d_cen*1e6), n_in_transit=nin,
                  n_transits=len(eps), epochs=eps, per_event=per_event,
                  odd_depth_ppm=float(d_odd*1e6), odd_err_ppm=float(e_odd_emp*1e6), even_depth_ppm=float(d_even*1e6), even_err_ppm=float(e_even_emp*1e6), odd_even_sigma=float(oe_sig),
                  secondary_ppm=float(d_sec*1e6), secondary_sigma=float(sec_sig), strongest_other_dip_ppm=float(dscan[imax]*1e6), strongest_other_dip_phase=float(pgrid[imax]), strongest_other_dip_sigma=float(dscan[imax]/e_emp),
                  per_year=per_year, astropy_stats_depth_odd_even=[float(x) for x in stats["depth_odd"]] + [float(x) for x in stats["depth_even"]],
                  inverted_narrow_max_power=float(rinv.power.max()), real_narrow_max_power=float(rf.power.max())),
        sap=dict(depth_ppm=float(d_sap*1e6), depth_err_empirical_ppm=float(e_sap*1e6), snr=float(d_sap/e_sap), scatter_ppm=float(rms_sap*1e6), ratio_sap_over_psf=float(d_sap/d), n_in_transit=nin_sap),
        plausibility=dict(Rstar_rsun=Rstar, implied_rp_re=float(rp_re), central_duration_expect_h=float(T_central_h)),
        pipeline_reference=dict(period_d=PIPE_P, depth_ppm=PIPE_DEPTH, duration_h=PIPE_DUR),
        log=log, runtime_s=time.time()-T0)
    with open(f"{OUT}/results.json", "w") as fh: json.dump(results, fh, indent=1)
    np.save(f"{OUT}/bls_periodogram.npy", np.vstack([periods, pw, S]))

    # ---------------- figure
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    BLUE, ORANGE, GRAY, INK = "#2a78d6", "#eb6834", "#b8b7b2", "#52514e"
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5), facecolor="#fcfcfb")
    for a in ax.ravel():
        a.set_facecolor("#fcfcfb"); a.grid(True, color="#e8e7e3", lw=0.6); a.tick_params(colors=INK, labelsize=9)
        for s in a.spines.values(): s.set_color("#d8d7d2")
    # periodogram
    a = ax[0,0]; a.plot(periods, S, color=BLUE, lw=0.8)
    for pe_ in peaks: a.axvline(pe_["period_d"], color=ORANGE, lw=0.8, ls=":")
    a.axvline(PIPE_P, color=INK, lw=0.6, ls="--"); a.set_xscale("log"); a.set_xlabel("period [d]"); a.set_ylabel("BLS SDE (own grid, SNR objective)")
    a.set_title(f"blind BLS 1-50 d: top {peaks[0]['period_d']:.4f} d, SDE {peaks[0]['sde']:.1f}; dashed = 29.433 d", fontsize=10, color=INK)
    # fold
    a = ax[0,1]; phh = phase(t2, P, t0) * P * 24
    s = np.abs(phh) < 12
    a.plot(phh[s], (f2[s]-1)*1e6, ".", color=GRAY, ms=2, alpha=0.6, label="PSF flux, masked re-detrend")
    bins = np.arange(-12, 12.01, 0.5); bi = np.digitize(phh[s], bins)
    bm = np.array([np.mean(f2[s][bi==i]) if np.sum(bi==i)>3 else np.nan for i in range(1,len(bins))])
    a.plot(0.5*(bins[1:]+bins[:-1]), (bm-1)*1e6, "-o", color=BLUE, ms=4, lw=1.5, label="30-min bins")
    a.axvspan(-dur*12, dur*12, color=ORANGE, alpha=0.12, lw=0, label=f"box {dur*24:.2f} h, depth {d*1e6:.0f}±{e_emp*1e6:.0f} ppm")
    a.set_ylim(-4000, 3000); a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm"); a.legend(fontsize=8, loc="lower left", frameon=False)
    a.set_title(f"fold at P={P:.5f} d, t0={t0:.3f} BKJD, empirical SNR {d/e_emp:.1f}", fontsize=10, color=INK)
    # odd/even + secondary
    a = ax[1,0]
    for sel, lab, col, off in [(odd, "odd", BLUE, 0), (even, "even", ORANGE, 0)]:
        pp = phase(t2[sel], P, t0) * P * 24; ss = np.abs(pp) < 12; bi = np.digitize(pp[ss], bins)
        bm = np.array([np.mean(f2[sel][ss][bi==i]) if np.sum(bi==i)>3 else np.nan for i in range(1,len(bins))])
        a.plot(0.5*(bins[1:]+bins[:-1]), (bm-1)*1e6, "-o", color=col, ms=3.5, lw=1.3, label=lab)
    ppsec = (((t2 - t0)/P - 0.5 + 0.5) % 1 - 0.5) * P * 24; ss = np.abs(ppsec) < 12; bi = np.digitize(ppsec[ss], bins)
    bm = np.array([np.mean(f2[ss][bi==i]) if np.sum(bi==i)>3 else np.nan for i in range(1,len(bins))])
    a.plot(0.5*(bins[1:]+bins[:-1]), (bm-1)*1e6, "-", color=INK, lw=1.0, label="phase 0.5 (secondary)")
    a.set_ylim(-1500, 800); a.set_xlabel("hours from mid-transit / mid-secondary"); a.set_ylabel("ppm (30-min bins)"); a.legend(fontsize=8, frameon=False)
    a.set_title(f"odd {d_odd*1e6:.0f}±{e_odd_emp*1e6:.0f}, even {d_even*1e6:.0f}±{e_even_emp*1e6:.0f} ppm ({oe_sig:.1f}σ); secondary {d_sec*1e6:.0f} ppm ({sec_sig:.1f}σ)", fontsize=10, color=INK)
    # per-event depths
    a = ax[1,1]
    xs = [pe_["tc_bkjd"] for pe_ in per_event]; ys = [pe_["depth_ppm"] for pe_ in per_event]; es = [pe_["err_ppm"] or 0 for pe_ in per_event]
    a.errorbar(xs, ys, yerr=es, fmt="o", color=BLUE, ms=5, lw=1, capsize=2, label="per-event box depth")
    a.axhline(d*1e6, color=ORANGE, lw=1, label=f"global {d*1e6:.0f} ppm"); a.axhline(0, color=INK, lw=0.6)
    a.axhline(d_sap*1e6, color=INK, lw=1, ls="--", label=f"SAP_FLUX global {d_sap*1e6:.0f} ppm")
    a.set_xlabel("mid-transit time [BKJD]"); a.set_ylabel("depth [ppm]"); a.legend(fontsize=8, frameon=False)
    a.set_title(f"{len(eps)} transits with data; per-year: " + ", ".join(f"Y{py['year']} {py['depth_ppm']:.0f}" if py['depth_ppm'] is not None else f"Y{py['year']} n/a" for py in per_year) + " ppm", fontsize=10, color=INK)
    fig.suptitle("kic-005866176 independent re-derivation (no pipeline code) — signal, not a planet", fontsize=12, color="#0b0b0b")
    fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig(f"{OUT}/fold_rederive.png", dpi=130)
    L(f"done in {time.time()-T0:.0f}s; wrote {OUT}/results.json, fold_rederive.png")

if __name__ == "__main__":
    main()
