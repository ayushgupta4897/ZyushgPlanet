#!/usr/bin/env python
"""
Independent re-derivation for kic-009655224 (KBONUS-BKG PSF light curve).

Written from scratch with numpy / scipy / astropy only. Does NOT import or read
anything under pipeline/.  No wotan, no lightkurve, no TLS.

Steps
 1. load stitched LC, SAP_QUALITY==0, per-quarter normalisation, own detrend
    (time-windowed running median, 1.0 d window, iterated with outlier masking)
 2. blind BLS 1-50 d (astropy BoxLeastSquares, own period grid), top-3 peaks,
    identical run on the inverted flux (2 - f) as the null
 3. at the best ephemeris: box depth with random-phase empirical error,
    odd/even, secondary at phase 0.5, per-year depths, n_transits, duration scan,
    per-event depths (single-event dominance)
 4. same depth measurement on SAP_FLUX (dilution check)
 5. results JSON + one multi-panel fold figure

Run id: adversarial/kic-009655224/rederive  (date 2026-10-17)
"""
import json, sys, time, os
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FNAME = "kic-009655224"
FITS = "/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-009655224_kepler_v1.0_lc.fits"
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-009655224/rederive"
PIPE_P, PIPE_DEPTH, PIPE_DUR_H = 6.70350, 638.0, 2.4   # pipeline values to compare against (given in task)
RNG = np.random.default_rng(20261017)
CADENCE = 29.4244 / 60.0 / 24.0   # d

# ----------------------------------------------------------------------------- 1. load
def load(path):
    with fits.open(path) as f:
        d = f["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float)
        fl = np.asarray(d["FLUX"], float)
        fe = np.asarray(d["FLUX_ERR"], float)
        sap = np.asarray(d["SAP_FLUX"], float)
        sape = np.asarray(d["SAP_FLUX_ERR"], float)
        q = np.asarray(d["SAP_QUALITY"], int)
        qu = np.asarray(d["QUARTER"], int)
    m = (q == 0) & np.isfinite(t) & np.isfinite(fl) & np.isfinite(sap) & (fl > 0) & (sap > 0)
    o = np.argsort(t[m])
    out = dict(t=t[m][o], f=fl[m][o], fe=fe[m][o], sap=sap[m][o], sape=sape[m][o], q=qu[m][o])
    # per-quarter normalisation
    for key, ekey in (("f", "fe"), ("sap", "sape")):
        x = out[key].copy(); e = out[ekey].copy()
        for qq in np.unique(out["q"]):
            s = out["q"] == qq
            med = np.nanmedian(x[s]); x[s] /= med; e[s] /= med
        out[key + "n"] = x; out[ekey + "n"] = e
    return out


def running_median_time(t, y, window=1.0, mask=None):
    """Median of y within +/- window/2 in time around each point (own implementation)."""
    if mask is None:
        mask = np.ones_like(y, bool)
    tm, ym = t[mask], y[mask]
    lo = np.searchsorted(tm, t - window / 2.0)
    hi = np.searchsorted(tm, t + window / 2.0)
    trend = np.empty_like(y)
    for i in range(len(t)):
        a, b = lo[i], hi[i]
        if b - a < 5:   # too few neighbours: widen
            a = max(0, a - 5); b = min(len(tm), b + 5)
        trend[i] = np.median(ym[a:b])
    return trend


def detrend(t, y, window=1.0, nsig_hi=3.0, nsig_lo=8.0, niter=3):
    """Divide by iterated running median; upward outliers (>nsig_hi) and extreme
    downward outliers (>nsig_lo) excluded from the trend estimate, then clipped."""
    mask = np.ones_like(y, bool)
    for _ in range(niter):
        tr = running_median_time(t, y, window, mask)
        r = y / tr - 1.0
        mad = 1.4826 * np.median(np.abs(r[mask] - np.median(r[mask])))
        mask = (r < nsig_hi * mad) & (r > -nsig_lo * mad)
    tr = running_median_time(t, y, window, mask)
    flat = y / tr
    r = flat - 1.0
    mad = 1.4826 * np.median(np.abs(r - np.median(r)))
    keep = (r < 4.0 * mad) & (r > -8.0 * mad)
    return flat, keep, mad


# ----------------------------------------------------------------------------- 2. BLS
def period_grid(pmin, pmax, baseline, dur_min, oversample=3.0):
    """dP = dur_min * P^2 / (baseline * oversample): phase drift over baseline < dur/oversample."""
    fmin, fmax = 1.0 / pmax, 1.0 / pmin
    df = dur_min / (baseline * oversample)   # in frequency: d(1/P) = dP/P^2
    n = int(np.ceil((fmax - fmin) / df))
    freqs = fmin + df * np.arange(n + 1)
    return np.sort(1.0 / freqs)


def sde_of(power):
    p = np.asarray(power, float)
    return (p.max() - p.mean()) / p.std()


def detrended_sde(periods, power, win=501):
    """SDE after removing a running median (in grid index) from the power spectrum."""
    from scipy.ndimage import median_filter
    base = median_filter(power, size=win, mode="nearest")
    r = power - base
    return (r.max() - r.mean()) / r.std(), r


def top_peaks(periods, power, n=3, excl_frac=0.02):
    """Top-n distinct peaks: exclude +/- excl_frac in period around each selected peak."""
    p = power.copy()
    out = []
    for _ in range(n):
        i = int(np.argmax(p))
        if not np.isfinite(p[i]):
            break
        out.append((float(periods[i]), float(power[i]), i))
        sel = np.abs(periods / periods[i] - 1.0) < excl_frac
        p[sel] = -np.inf
    return out


def run_bls(t, y, dy, periods, durations, objective="likelihood"):
    bls = BoxLeastSquares(t, y, dy)
    res = bls.power(periods, durations, objective=objective)
    return bls, res


# ----------------------------------------------------------------------------- 3. statistics at the ephemeris
def phase_of(t, P, t0):
    return ((t - t0) / P + 0.5) % 1.0 - 0.5   # [-0.5,0.5), transit at 0


def box_depth(t, y, P, t0, dur, phase_center=0.0, oot_excl=2.0):
    """Depth (ppm) = median(out-of-transit) - mean(in-box); in-box: |phase*P - c| < dur/2."""
    ph = phase_of(t, P, t0)
    dt = (ph - phase_center) * P
    dt = (dt + P / 2) % P - P / 2
    inb = np.abs(dt) < dur / 2
    oot = np.abs(dt) > oot_excl * dur / 2
    if inb.sum() < 3:
        return np.nan, inb.sum()
    return (np.median(y[oot]) - np.mean(y[inb])) * 1e6, int(inb.sum())


def random_phase_depths(t, y, P, t0, dur, n=1000, avoid=None):
    """Depths of the same box placed at random phases, avoiding +/-1.5 dur around
    phases in `avoid` (fractions). Returns array of depths (ppm)."""
    if avoid is None:
        avoid = [0.0]
    w = 1.5 * dur / P
    deps = []
    tries = 0
    while len(deps) < n and tries < 20 * n:
        tries += 1
        c = RNG.uniform(-0.5, 0.5)
        if any(abs(((c - a) + 0.5) % 1.0 - 0.5) < w for a in avoid):
            continue
        d, nin = box_depth(t, y, P, t0, dur, phase_center=c)
        if np.isfinite(d):
            deps.append(d)
    return np.array(deps)


def refine_ephemeris(t, y, dy, P0, dur0, frac=0.004, ndur=25):
    """Fine BLS around P0; returns P, t0, dur, depth, Delta-loglike."""
    periods = np.linspace(P0 * (1 - frac), P0 * (1 + frac), 4001)
    durations = np.linspace(0.5 / 24, 6.0 / 24, ndur)
    bls = BoxLeastSquares(t, y, dy)
    res = bls.power(periods, durations, objective="likelihood", oversample=20)
    i = int(np.argmax(res.power))
    return (float(res.period[i]), float(res.transit_time[i]), float(res.duration[i]),
            float(res.depth[i]) * 1e6, float(res.power[i]), res)


def duration_scan(t, y, P, t0, durs_h):
    """Delta chi2 (vs flat) of a box of each duration at fixed P,t0; depth free."""
    out = []
    for dh in durs_h:
        dur = dh / 24.0
        ph = phase_of(t, P, t0); dt = ph * P
        inb = np.abs(dt) < dur / 2
        oot = ~inb
        if inb.sum() < 3:
            out.append((dh, np.nan, np.nan)); continue
        sig = 1.4826 * np.median(np.abs(y[oot] - np.median(y[oot])))
        d = np.median(y[oot]) - np.mean(y[inb])
        dchi2 = d ** 2 * inb.sum() / sig ** 2
        out.append((dh, d * 1e6, dchi2))
    return out


def per_event(t, y, P, t0, dur, local_win=1.0):
    """Per-epoch depth: local OOT median (+/- local_win d, excluding +/- dur) minus in-box mean."""
    n = np.round((t - t0) / P).astype(int)
    rows = []
    for ep in np.unique(n):
        tc = t0 + ep * P
        inb = np.abs(t - tc) < dur / 2
        loc = (np.abs(t - tc) < local_win) & (np.abs(t - tc) > dur)
        if inb.sum() == 0:
            continue
        if loc.sum() < 10:
            rows.append((int(ep), float(tc), int(inb.sum()), int(loc.sum()), np.nan, np.nan)); continue
        sig = 1.4826 * np.median(np.abs(y[loc] - np.median(y[loc])))
        d = (np.median(y[loc]) - np.mean(y[inb])) * 1e6
        e = sig / np.sqrt(inb.sum()) * 1e6
        # gap proximity: distance to nearest data gap > 0.5 d edge
        rows.append((int(ep), float(tc), int(inb.sum()), int(loc.sum()), float(d), float(e)))
    return rows


# ============================================================================= main
def main():
    T0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    D = load(FITS)
    t = D["t"]
    print(f"loaded {len(t)} good cadences; span {t.max()-t.min():.1f} d; quarters {np.unique(D['q']).tolist()}")

    # ---- detrend PSF flux and SAP flux with the same method
    WIN = 1.0
    flat, keep, mad = detrend(t, D["fn"], window=WIN)
    flat_sap, keep_sap, mad_sap = detrend(t, D["sapn"], window=WIN)
    k = keep & keep_sap
    t, y, ysap, q = t[k], flat[k], flat_sap[k], D["q"][k]
    # per-quarter robust scatter as dy
    dy = np.empty_like(y); dysap = np.empty_like(y)
    for qq in np.unique(q):
        s = q == qq
        dy[s] = 1.4826 * np.median(np.abs(y[s] - np.median(y[s])))
        dysap[s] = 1.4826 * np.median(np.abs(ysap[s] - np.median(ysap[s])))
    print(f"after clipping {len(t)} cadences; PSF MAD-scatter {mad*1e6:.0f} ppm; SAP {mad_sap*1e6:.0f} ppm")
    np.save(os.path.join(OUT, "flat_lc.npy"), np.vstack([t, y, dy, ysap, dysap, q]))

    # ---- 2. blind BLS 1-50 d
    baseline = t.max() - t.min()
    DUR_MIN = 0.04
    durations = np.array([0.04, 0.06, 0.08, 0.10, 0.13, 0.17, 0.22, 0.30])
    periods = period_grid(1.0, 50.0, baseline, DUR_MIN, oversample=3.0)
    print(f"BLS grid: {len(periods)} periods, {len(durations)} durations; baseline {baseline:.1f} d")
    tb = time.time()
    bls, res = run_bls(t, y, dy, periods, durations)
    print(f"BLS real done in {time.time()-tb:.0f} s")
    sde_raw = sde_of(res.power)
    sde_det, pdet = detrended_sde(res.period, res.power)
    peaks = top_peaks(res.period, pdet, n=3)
    peaks_raw = top_peaks(res.period, res.power, n=3)
    ibest = peaks[0][2]
    P_best = float(res.period[ibest]); t0_best = float(res.transit_time[ibest])
    dur_best = float(res.duration[ibest]); depth_best = float(res.depth[ibest]) * 1e6
    print("top-3 (detrended power):", [(round(p, 5), round(v, 2)) for p, v, _ in peaks])
    print("top-3 (raw power):     ", [(round(p, 5), round(v, 2)) for p, v, _ in peaks_raw])
    print(f"best P={P_best:.5f} t0={t0_best:.4f} dur={dur_best*24:.2f} h depth={depth_best:.0f} ppm SDE raw={sde_raw:.2f} det={sde_det:.2f}")

    # alias check vs pipeline period
    def alias_of(p, pref=PIPE_P):
        for num, den in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (4, 1)]:
            if abs(p / (pref * num / den) - 1.0) < 0.002:
                return f"{num}/{den}xP_pipe"
        return None
    peak_info = []
    for p, v, i in peaks:
        peak_info.append(dict(period_d=p, power_detrended=v, power_raw=float(res.power[i]),
                              sde_local=float((v - pdet.mean()) / pdet.std()),
                              duration_h=float(res.duration[i] * 24), depth_ppm=float(res.depth[i] * 1e6),
                              t0_bkjd=float(res.transit_time[i]), alias_of_pipeline=alias_of(p)))

    # ---- inverted-flux null (same grid, same code)
    tb = time.time()
    _, res_inv = run_bls(t, 2.0 - y, dy, periods, durations)
    sde_inv_raw = sde_of(res_inv.power)
    sde_inv_det, pdet_inv = detrended_sde(res_inv.period, res_inv.power)
    pk_inv = top_peaks(res_inv.period, pdet_inv, n=3)
    print(f"BLS inverted done in {time.time()-tb:.0f} s: SDE raw={sde_inv_raw:.2f} det={sde_inv_det:.2f}; top P={pk_inv[0][0]:.4f}")

    # ---- 3. refine at the best peak (and, separately, at the pipeline period if different)
    recovered_blind = abs(P_best / PIPE_P - 1.0) < 0.002 or alias_of(P_best) is not None
    P, t0, dur, depth_bls, dll, res_fine = refine_ephemeris(t, y, dy, P_best, dur_best)
    print(f"refined: P={P:.6f} t0={t0:.5f} dur={dur*24:.2f} h depth(BLS)={depth_bls:.0f} ppm dlogL={dll:.1f}")
    # put t0 near the data middle-ish first transit for reporting
    n0 = np.floor((t.min() - t0) / P)
    t0 = t0 + n0 * P
    if t0 < t.min(): t0 += P

    # duration scan (box, depth free)
    dscan = duration_scan(t, y, P, t0, np.arange(0.5, 6.01, 0.25))
    dsc = np.array([(a, b, c) for a, b, c in dscan], float)
    j = int(np.nanargmax(dsc[:, 2]))
    dur_scan_h = float(dsc[j, 0])
    # use the finer of the two duration estimates (BLS refine) for the box statistics
    DUR = dur

    # box depth + empirical error (PSF flux)
    depth, nin = box_depth(t, y, P, t0, DUR)
    rp = random_phase_depths(t, y, P, t0, DUR, n=1500, avoid=[0.0])
    depth_err = float(np.std(rp)); rp_mean = float(np.mean(rp))
    snr_emp = depth / depth_err
    ph = phase_of(t, P, t0)
    oot = np.abs(ph * P) > DUR
    formal_err = float(np.std(y[oot]) / np.sqrt(nin) * 1e6)
    print(f"PSF box depth {depth:.0f} +/- {depth_err:.0f} ppm (random-phase; mean {rp_mean:.0f}); formal {formal_err:.0f}; N_in={nin}; SNR_emp={snr_emp:.2f}")

    # odd / even
    n_ep = np.round((t - t0) / P).astype(int)
    odd = n_ep % 2 == 1; even = ~odd
    d_odd, nin_odd = box_depth(t[odd], y[odd], P, t0, DUR)
    d_even, nin_even = box_depth(t[even], y[even], P, t0, DUR)
    rp_odd = random_phase_depths(t[odd], y[odd], P, t0, DUR, n=600)
    rp_even = random_phase_depths(t[even], y[even], P, t0, DUR, n=600)
    e_odd, e_even = float(np.std(rp_odd)), float(np.std(rp_even))
    oe_sigma = abs(d_odd - d_even) / np.hypot(e_odd, e_even)
    print(f"odd {d_odd:.0f}+/-{e_odd:.0f} (N={nin_odd})  even {d_even:.0f}+/-{e_even:.0f} (N={nin_even})  diff {oe_sigma:.2f} sigma")

    # secondary at phase 0.5
    sec, nin_sec = box_depth(t, y, P, t0, DUR, phase_center=0.5)
    rp_sec = random_phase_depths(t, y, P, t0, DUR, n=1500, avoid=[0.0, 0.5])
    sec_err = float(np.std(rp_sec))
    sec_sigma = sec / sec_err
    # also scan phase 0.3-0.7 for the deepest secondary-like box
    ph_scan = np.linspace(0.2, 0.8, 241)
    sec_scan = np.array([box_depth(t, y, P, t0, DUR, phase_center=c)[0] for c in ph_scan])
    print(f"secondary at 0.5: {sec:.0f} +/- {sec_err:.0f} ppm ({sec_sigma:.2f} sigma); max box 0.2-0.8: {np.nanmax(sec_scan):.0f} ppm at phase {ph_scan[np.nanargmax(sec_scan)]:.3f}")

    # per-year depths (Kepler years by quarter group)
    years = {"Y1_Q3-4": [3, 4], "Y2_Q5-8": [5, 6, 7, 8], "Y3_Q9-12": [9, 10, 11, 12], "Y4_Q13-17": [13, 14, 15, 16, 17]}
    per_year = {}
    for name, qs in years.items():
        s = np.isin(q, qs)
        if s.sum() < 100:
            continue
        d_y, n_y = box_depth(t[s], y[s], P, t0, DUR)
        rp_y = random_phase_depths(t[s], y[s], P, t0, DUR, n=500)
        per_year[name] = dict(depth_ppm=float(d_y), err_ppm=float(np.std(rp_y)), n_in=int(n_y),
                              quarters=sorted(np.unique(q[s]).tolist()))
        print(f"  {name}: {d_y:.0f} +/- {np.std(rp_y):.0f} ppm (N_in={n_y})")
    # per-quarter depths too
    per_quarter = {}
    for qq in np.unique(q):
        s = q == qq
        d_q, n_q = box_depth(t[s], y[s], P, t0, DUR)
        rp_q = random_phase_depths(t[s], y[s], P, t0, DUR, n=300)
        per_quarter[int(qq)] = dict(depth_ppm=float(d_q), err_ppm=float(np.std(rp_q)), n_in=int(n_q))

    # number of transits and per-event depths
    ev = per_event(t, y, P, t0, DUR)
    exp_n = DUR / CADENCE
    n_any = len(ev)
    n_half = sum(1 for r in ev if r[2] >= 0.5 * exp_n)
    n_full = sum(1 for r in ev if r[2] >= int(np.floor(exp_n)))
    ev_d = np.array([r[4] for r in ev]); ev_e = np.array([r[5] for r in ev]); ev_ok = np.isfinite(ev_d) & np.isfinite(ev_e)
    # weighted mean of per-event depths; single-event dominance
    w = 1 / ev_e[ev_ok] ** 2
    ev_mean = float(np.sum(w * ev_d[ev_ok]) / np.sum(w)); ev_mean_err = float(1 / np.sqrt(np.sum(w)))
    snr_contrib = (ev_d[ev_ok] / ev_e[ev_ok])
    order = np.argsort(-snr_contrib)
    tot = np.sum(snr_contrib ** 2)
    top1_frac = float(snr_contrib[order[0]] ** 2 / tot); top3_frac = float(np.sum(snr_contrib[order[:3]] ** 2) / tot)
    # depth with the 3 most significant events removed
    drop = set(int(ev[i][0]) for i in np.where(ev_ok)[0][order[:3]])
    keep3 = ~np.isin(n_ep, list(drop))
    d_no3, _ = box_depth(t[keep3], y[keep3], P, t0, DUR)
    rp_no3 = random_phase_depths(t[keep3], y[keep3], P, t0, DUR, n=600)
    print(f"events: any={n_any} >=half={n_half} full={n_full}; weighted per-event mean {ev_mean:.0f}+/-{ev_mean_err:.0f}; top1 {top1_frac:.2f} top3 {top3_frac:.2f} of sum(snr^2); depth w/o top3 events {d_no3:.0f}+/-{np.std(rp_no3):.0f}")

    # ---- 4. SAP flux depth at the same ephemeris
    depth_sap, nin_sap = box_depth(t, ysap, P, t0, DUR)
    rp_sap = random_phase_depths(t, ysap, P, t0, DUR, n=1500)
    sap_err = float(np.std(rp_sap))
    ratio = depth_sap / depth
    ratio_err = ratio * np.hypot(sap_err / depth_sap if depth_sap != 0 else np.inf, depth_err / depth)
    med_psf = float(np.median(D["f"])); med_sap = float(np.median(D["sap"]))
    exp_ratio_fluxes = med_psf / med_sap   # if PSF flux == target flux and all extra SAP flux is from neighbours
    # header-based expectation: CROWDSAP = target fraction of SAP flux
    with fits.open(FITS) as f:
        crowd = [f[i].header.get("CROWDSAP") for i in range(2, len(f), 2)]
        flfrc = [f[i].header.get("FLFRCSAP") for i in range(2, len(f), 2)]
        psffrac = [f[i].header.get("PSFFRAC") for i in range(2, len(f), 2)]
    exp_ratio_crowd = float(np.mean(crowd))
    print(f"SAP box depth {depth_sap:.0f} +/- {sap_err:.0f} ppm; ratio SAP/PSF = {ratio:.2f} +/- {ratio_err:.2f}; expected if on-target: {exp_ratio_fluxes:.2f} (flux medians) or {exp_ratio_crowd:.2f} (mean CROWDSAP)")
    # SAP odd/even + secondary for completeness
    sec_sap, _ = box_depth(t, ysap, P, t0, DUR, phase_center=0.5)

    # ---- also: statistics at the pipeline ephemeris exactly, if different from best
    at_pipe = None
    if abs(P / PIPE_P - 1) > 1e-4:
        Pp, t0p, durp, dbp, dllp, _ = refine_ephemeris(t, y, dy, PIPE_P, PIPE_DUR_H / 24, frac=0.002)
        dp, ninp = box_depth(t, y, Pp, t0p, durp)
        rpp = random_phase_depths(t, y, Pp, t0p, durp, n=800)
        at_pipe = dict(period_d=Pp, t0_bkjd=t0p, duration_h=durp * 24, depth_ppm=float(dp), depth_err_ppm=float(np.std(rpp)), n_in=int(ninp))
        print("at pipeline period:", at_pipe)

    # ---- 5. results JSON
    R = dict(
        fname=FNAME, run_id="poolA_main/adversarial/kic-009655224/rederive", date="2026-10-17",
        method=dict(quality="SAP_QUALITY==0", normalisation="per-quarter median", detrend=f"iterated time-windowed running median, window {WIN} d, 3-sigma upward / 8-sigma downward mask, then 4/8-sigma clip",
                    bls="astropy BoxLeastSquares, objective=likelihood, own grid dP = dur_min P^2/(T*3), dur_min=0.04 d",
                    durations_d=durations.tolist(), period_range_d=[1.0, 50.0], n_periods=int(len(periods)),
                    depth="median(OOT, |dt|>dur) - mean(in-box, |dt|<dur/2); error = std of the same box at random phases (avoiding transit)"),
        data=dict(n_good=int(len(t)), span_d=float(baseline), quarters=sorted(np.unique(q).tolist()), psf_mad_ppm=float(mad * 1e6), sap_mad_ppm=float(mad_sap * 1e6),
                  median_psf_flux_e_s=med_psf, median_sap_flux_e_s=med_sap, crowdsap_per_quarter=crowd, flfrcsap_per_quarter=flfrc, psffrac_per_quarter=psffrac),
        blind_bls=dict(top3_peaks=peak_info, sde_raw=float(sde_raw), sde_detrended=float(sde_det), best_period_d=P_best, recovered_blind=bool(recovered_blind),
                       pipeline_period_d=PIPE_P, best_vs_pipeline_frac=float(P_best / PIPE_P - 1)),
        inverted_null=dict(sde_raw=float(sde_inv_raw), sde_detrended=float(sde_inv_det), top3=[dict(period_d=p, power_detrended=v) for p, v, _ in pk_inv]),
        ephemeris=dict(period_d=float(P), t0_bkjd=float(t0), t0_bjd=float(t0 + 2454833.0), duration_h=float(DUR * 24), duration_h_boxscan=dur_scan_h,
                       depth_bls_ppm=float(depth_bls), delta_loglike=float(dll), duration_scan=[dict(dur_h=a, depth_ppm=b, dchi2=c) for a, b, c in dscan]),
        depth=dict(psf_depth_ppm=float(depth), psf_depth_err_ppm=depth_err, random_phase_mean_ppm=rp_mean, formal_err_ppm=formal_err, n_in=int(nin), empirical_snr=float(snr_emp)),
        odd_even=dict(odd_ppm=float(d_odd), odd_err_ppm=e_odd, even_ppm=float(d_even), even_err_ppm=e_even, sigma=float(oe_sigma), n_in_odd=int(nin_odd), n_in_even=int(nin_even)),
        secondary=dict(phase=0.5, depth_ppm=float(sec), err_ppm=sec_err, sigma=float(sec_sigma), max_box_0p2_0p8_ppm=float(np.nanmax(sec_scan)), max_box_phase=float(ph_scan[np.nanargmax(sec_scan)]),
                       sap_secondary_ppm=float(sec_sap)),
        per_year=per_year, per_quarter=per_quarter,
        transits=dict(n_epochs_with_data=n_any, n_epochs_ge_half_coverage=n_half, n_epochs_full_coverage=n_full, expected_in_box_cadences=float(exp_n),
                      weighted_event_mean_ppm=ev_mean, weighted_event_mean_err_ppm=ev_mean_err, top1_event_frac_of_snr2=top1_frac, top3_event_frac_of_snr2=top3_frac,
                      depth_without_top3_events_ppm=float(d_no3), depth_without_top3_err_ppm=float(np.std(rp_no3)),
                      events=[dict(epoch=r[0], tc_bkjd=r[1], n_in=r[2], n_local=r[3], depth_ppm=r[4], err_ppm=r[5]) for r in ev]),
        sap=dict(depth_ppm=float(depth_sap), err_ppm=sap_err, n_in=int(nin_sap), ratio_sap_over_psf=float(ratio), ratio_err=float(ratio_err),
                 expected_ratio_if_on_target_from_flux_medians=float(exp_ratio_fluxes), expected_ratio_if_on_target_from_mean_crowdsap=exp_ratio_crowd),
        at_pipeline_period=at_pipe,
        runtime_s=float(time.time() - T0),
    )
    with open(os.path.join(OUT, "rederive_results.json"), "w") as fh:
        json.dump(R, fh, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
    # events csv
    with open(os.path.join(OUT, "events_rederive.csv"), "w") as fh:
        fh.write("epoch,tc_bkjd,n_in,n_local,depth_ppm,err_ppm\n")
        for r in ev:
            fh.write(",".join(str(x) for x in r) + "\n")
    np.save(os.path.join(OUT, "bls_power.npy"), np.vstack([res.period, res.power, pdet, res_inv.power]))

    # ---- figure
    fig, ax = plt.subplots(3, 2, figsize=(13, 12))
    a = ax[0, 0]
    a.plot(res.period, pdet, lw=0.4, color="k")
    a.plot(res_inv.period, pdet_inv - pdet.max() * 1.1, lw=0.4, color="0.6", label="inverted (offset)")
    for p, v, i in peaks:
        a.axvline(p, color="r", alpha=0.3)
    a.axvline(PIPE_P, color="b", ls="--", alpha=0.6, label=f"pipeline P={PIPE_P}")
    a.set_xscale("log"); a.set_xlabel("period [d]"); a.set_ylabel("BLS dlogL (median-detrended)")
    a.set_title(f"blind BLS 1-50 d: SDE={sde_det:.1f} (inv {sde_inv_det:.1f}); best P={P_best:.4f} d"); a.legend(fontsize=8)

    def binned(ph, yy, nb=120, lo=-0.5, hi=0.5):
        edges = np.linspace(lo, hi, nb + 1); idx = np.digitize(ph, edges) - 1
        bm = np.array([np.mean(yy[idx == i]) if np.sum(idx == i) > 0 else np.nan for i in range(nb)])
        be = np.array([np.std(yy[idx == i]) / np.sqrt(np.sum(idx == i)) if np.sum(idx == i) > 1 else np.nan for i in range(nb)])
        return 0.5 * (edges[1:] + edges[:-1]), bm, be

    a = ax[0, 1]
    bc, bm, be = binned(ph, y, nb=200)
    a.plot(ph, (y - 1) * 1e6, ",", color="0.8")
    a.errorbar(bc, (bm - 1) * 1e6, be * 1e6, fmt="o", ms=3, color="k")
    a.set_ylim(-1500, 1500); a.set_xlabel("phase"); a.set_ylabel("ppm")
    a.set_title(f"PSF flux folded at P={P:.5f} d (200 bins); secondary at 0.5: {sec:.0f}+/-{sec_err:.0f} ppm")
    a.axhline(0, color="0.5", lw=0.5); a.axvline(0.5, color="r", ls=":", alpha=0.5); a.axvline(-0.5, color="r", ls=":", alpha=0.5)

    a = ax[1, 0]
    zoom = np.abs(ph * P) < 0.4
    bc, bm, be = binned(ph[zoom] * P * 24, y[zoom], nb=48, lo=-0.4 * 24, hi=0.4 * 24)
    a.plot(ph[zoom] * P * 24, (y[zoom] - 1) * 1e6, ".", ms=2, color="0.8", alpha=0.5)
    a.errorbar(bc, (bm - 1) * 1e6, be * 1e6, fmt="o", ms=4, color="k", label="PSF")
    bc2, bm2, be2 = binned(ph[zoom] * P * 24, ysap[zoom], nb=48, lo=-0.4 * 24, hi=0.4 * 24)
    a.errorbar(bc2 + 0.1, (bm2 - 1) * 1e6, be2 * 1e6, fmt="s", ms=3, color="tab:orange", alpha=0.8, label="SAP")
    a.axvspan(-DUR * 12, DUR * 12, color="b", alpha=0.08)
    a.axhline(-depth, color="k", ls="--", lw=0.8); a.axhline(-depth_sap, color="tab:orange", ls="--", lw=0.8)
    a.set_ylim(-1500, 1000); a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm"); a.legend(fontsize=8)
    a.set_title(f"zoom: PSF {depth:.0f}+/-{depth_err:.0f} ppm (SNR {snr_emp:.1f}); SAP {depth_sap:.0f}+/-{sap_err:.0f}; dur {DUR*24:.2f} h")

    a = ax[1, 1]
    for sel, lab, col, off in ((odd, f"odd {d_odd:.0f}+/-{e_odd:.0f}", "tab:blue", -0.1), (even, f"even {d_even:.0f}+/-{e_even:.0f}", "tab:red", 0.1)):
        z = sel & zoom
        bc, bm, be = binned(ph[z] * P * 24, y[z], nb=32, lo=-0.4 * 24, hi=0.4 * 24)
        a.errorbar(bc + off, (bm - 1) * 1e6, be * 1e6, fmt="o", ms=3, color=col, label=lab)
    a.axvspan(-DUR * 12, DUR * 12, color="b", alpha=0.08); a.set_ylim(-1500, 1000)
    a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm"); a.legend(fontsize=8); a.set_title(f"odd/even: {oe_sigma:.2f} sigma")

    a = ax[2, 0]
    tcs = np.array([r[1] for r in ev]); dd = np.array([r[4] for r in ev]); ee = np.array([r[5] for r in ev])
    a.errorbar(tcs, dd, ee, fmt="o", ms=3, color="k", alpha=0.6)
    a.axhline(depth, color="b", ls="--"); a.axhline(0, color="0.5", lw=0.5)
    a.set_xlabel("transit mid-time [BKJD]"); a.set_ylabel("per-event depth [ppm]"); a.set_ylim(-6000, 6000)
    a.set_title(f"per-event depths ({n_any} epochs with data, {n_half} >=half coverage); top3 share of sum(snr^2)={top3_frac:.2f}")

    a = ax[2, 1]
    names = list(per_year.keys()); vals = [per_year[n]["depth_ppm"] for n in names]; errs = [per_year[n]["err_ppm"] for n in names]
    a.errorbar(range(len(names)), vals, errs, fmt="o", color="k", capsize=3)
    a.set_xticks(range(len(names))); a.set_xticklabels(names)
    a.axhline(depth, color="b", ls="--"); a.axhline(0, color="0.5", lw=0.5); a.set_ylabel("depth [ppm]"); a.set_title("per-year box depth (PSF)")
    fig.suptitle(f"{FNAME} independent re-derivation (adversarial/rederive, 2026-10-17)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fold_rederive.png"), dpi=130)
    print(f"done in {time.time()-T0:.0f} s; wrote {OUT}")


if __name__ == "__main__":
    main()
