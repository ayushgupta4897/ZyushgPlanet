#!/usr/bin/env python
"""
Independent re-derivation for kic-010213563 (Gaia DR3 2128591018876020864).

Written from scratch with numpy / scipy / astropy only. Deliberately shares NO
code with /pipeline/ (not imported, not read before these numbers were produced).

Steps
 1. load KBONUS-BKG stitched light curve, SAP_QUALITY==0, finite values
 2. per-quarter median normalisation
 3. detrend with my own time-windowed running median (1.0 d), segments split at
    gaps > 0.5 d, asymmetric sigma clipping (upper 3.5 sigma, lower 8 sigma so
    dips survive)
 4. blind BLS 1-50 d (astropy BoxLeastSquares, my own log-spaced frequency grid)
    on real flux AND on inverted flux (2 - f) as a null
 5. refine at the top peak, box depth with random-phase empirical error,
    odd/even, secondary at phase 0.5, per-year depths, per-event local depths,
    number of transits, duration
 6. same box depth on SAP_FLUX (dilution check)
Outputs: results.json, fold.png, bls_periodogram.csv, events.csv
"""
import json
import os
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- settings
FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        "hlsp_kbonus-bkg_kepler_kepler_kic-010213563_kepler_v1.0_lc.fits")
OUT = "/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-010213563/rederive"
FNAME = "kic-010213563"
PIPE_P = 36.76710          # pipeline value we are testing against (period only)
RNG = np.random.default_rng(20261017)

WIN_D = 1.0                # running-median window (days)
GAP_D = 0.5                # split segments at gaps longer than this
PMIN, PMAX = 1.0, 50.0     # blind search range (days)
DURS = np.array([0.05, 0.08, 0.12, 0.17, 0.24, 0.33])   # days
OVERSAMPLE = 3.0
N_RANDOM = 3000            # random-phase boxes for the empirical error
CADENCE = 0.0204340        # d, Kepler long cadence

# colours (dataviz reference palette, light mode)
C_BLUE, C_ORANGE, C_AQUA, C_GREY = "#2a78d6", "#eb6834", "#1baf7a", "#9a9a94"

if os.environ.get("REDERIVE_QUICK"):      # smoke test only: narrow grid, scratch output
    PMIN, PMAX, OVERSAMPLE, N_RANDOM = 30.0, 40.0, 1.0, 200
    OUT = os.environ["REDERIVE_QUICK"]
os.makedirs(OUT, exist_ok=True)
log_lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log_lines.append(s)


# ---------------------------------------------------------------- 1. load
def load(column):
    with fits.open(FITS) as f:
        d = f["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float)
        y = np.asarray(d[column], float)
        e = np.asarray(d[column + "_ERR"], float)
        q = np.asarray(d["QUARTER"], int)
        qual = np.asarray(d["SAP_QUALITY"], int)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (y > 0)
    t, y, e, q = t[m], y[m], e[m], q[m]
    o = np.argsort(t)
    return t[o], y[o], e[o], q[o]


# ---------------------------------------------------------------- 2. normalise
def normalise_per_quarter(t, y, e, q):
    yn = np.empty_like(y)
    en = np.empty_like(e)
    for qq in np.unique(q):
        s = q == qq
        med = np.median(y[s])
        yn[s] = y[s] / med
        en[s] = e[s] / med
    return yn, en


# ---------------------------------------------------------------- 3. detrend
def running_median(t, y, win, mask=None, min_pts=8):
    """Time-windowed running median. Points with mask==False are excluded from
    the median (used to mask transits) but still receive a trend value."""
    if mask is None:
        mask = np.ones_like(y, bool)
    tm, ym = t[mask], y[mask]
    lo = np.searchsorted(tm, t - win / 2.0)
    hi = np.searchsorted(tm, t + win / 2.0)
    trend = np.full_like(y, np.nan)
    for i in range(len(t)):
        if hi[i] - lo[i] >= min_pts:
            trend[i] = np.median(ym[lo[i]:hi[i]])
    return trend


def segments(t, gap):
    br = np.where(np.diff(t) > gap)[0] + 1
    return np.split(np.arange(len(t)), br)


def detrend(t, y, win=WIN_D, transit_mask=None, n_iter=3):
    """Returns flattened flux, trend, and a boolean 'keep' array (outliers
    removed). transit_mask (True = in transit) is excluded from the trend."""
    keep = np.ones_like(y, bool)
    use = np.ones_like(y, bool) if transit_mask is None else ~transit_mask
    trend = np.full_like(y, np.nan)
    for seg in segments(t, GAP_D):
        ts, ys = t[seg], y[seg]
        good = use[seg].copy()
        for _ in range(n_iter):
            tr = running_median(ts, ys, win, mask=good)
            res = ys - tr
            ok = np.isfinite(res)
            mad = 1.4826 * np.median(np.abs(res[ok & good] - np.median(res[ok & good])))
            if not np.isfinite(mad) or mad <= 0:
                break
            new_good = use[seg] & ok & (res < 3.5 * mad) & (res > -8.0 * mad)
            if np.array_equal(new_good, good):
                break
            good = new_good
        tr = running_median(ts, ys, win, mask=good)
        trend[seg] = tr
        res = ys - tr
        ok = np.isfinite(res)
        mad = 1.4826 * np.median(np.abs(res[ok] - np.median(res[ok])))
        keep[seg] = ok & (res < 3.5 * mad) & (res > -8.0 * mad)
    flat = y / trend
    return flat, trend, keep


# ---------------------------------------------------------------- 4. BLS grid
def frequency_grid(t, pmin, pmax, dmin, oversample):
    """Log-spaced: df_i = f_i * dmin / (oversample * baseline), i.e. the phase
    drift over the baseline is < dmin/oversample at every trial period."""
    T = t.max() - t.min()
    fmin, fmax = 1.0 / pmax, 1.0 / pmin
    step = dmin / (oversample * T)                  # d ln f
    n = int(np.ceil(np.log(fmax / fmin) / step)) + 1
    f = fmin * np.exp(np.arange(n) * step)
    return f[f <= fmax]


def sde_of(power):
    p = np.asarray(power, float)
    ok = np.isfinite(p)
    return (np.nanmax(p) - np.mean(p[ok])) / np.std(p[ok])


def spectrum_detrend(periods, power, nwin=2001):
    """Median-filter the periodogram in index space (grid is log-uniform in f)
    and return (power - median)/std -> SDE-like statistic per period."""
    from scipy.ndimage import median_filter
    med = median_filter(power, size=nwin, mode="nearest")
    resid = power - med
    ok = np.isfinite(resid)
    return resid / np.std(resid[ok])


def top_peaks(periods, stat, n=3, excl_frac=0.03):
    """Greedy peak picking; exclude +-excl_frac in period around each pick."""
    s = np.array(stat, float)
    out = []
    for _ in range(n):
        if not np.any(np.isfinite(s)):
            break
        i = int(np.nanargmax(s))
        out.append((float(periods[i]), float(s[i]), i))
        s[np.abs(periods - periods[i]) < excl_frac * periods[i]] = np.nan
    return out


def is_alias(p, pref, tol=0.01):
    for nmr in range(1, 5):
        for den in range(1, 5):
            pa = pref * nmr / den
            if abs(p - pa) / pa < tol:
                return f"{nmr}/{den}"
    return None


# ---------------------------------------------------------------- 5. box stats
def phase_days(t, P, t0):
    return ((t - t0 + 0.5 * P) % P) - 0.5 * P


def box_depth(t, f, P, t0, dur, out_lo=1.0, out_hi=3.0):
    """depth = mean(out-of-transit local) - mean(in-transit), fractional.
    in: |phi| < dur/2 ; out: out_lo*dur < |phi| < out_hi*dur."""
    ph = np.abs(phase_days(t, P, t0))
    inm = ph < dur / 2.0
    outm = (ph > out_lo * dur) & (ph < out_hi * dur)
    if inm.sum() < 3 or outm.sum() < 10:
        return np.nan, int(inm.sum())
    return float(np.mean(f[outm]) - np.mean(f[inm])), int(inm.sum())


def random_phase_error(t, f, P, t0, dur, n=N_RANDOM, avoid=2.0):
    """Same statistic at n random phases that avoid +-avoid*dur around the
    transit and around phase 0.5. Returns (robust std, mean, all values)."""
    vals = []
    tries = 0
    while len(vals) < n and tries < 20 * n:
        tries += 1
        shift = RNG.uniform(0, P)
        dphi = min(shift % P, P - shift % P)
        dsec = min(abs(shift - P / 2) % P, P - abs(shift - P / 2) % P)
        if dphi < avoid * dur or dsec < avoid * dur:
            continue
        v, nin = box_depth(t, f, P, t0 + shift, dur)
        if np.isfinite(v):
            vals.append(v)
    vals = np.array(vals)
    mad = 1.4826 * np.median(np.abs(vals - np.median(vals)))
    return float(mad), float(np.std(vals)), float(np.mean(vals)), vals


def per_event(t, f, P, t0, dur, gap_d=GAP_D):
    """Local depth per epoch: baseline = median of |dt| in (dur/2+0.05, dur/2+1.0)."""
    # gap edges for proximity: segment boundaries + data ends
    segs = segments(t, gap_d)
    edges = [t[0], t[-1]]
    for s in segs:
        edges += [t[s[0]], t[s[-1]]]
    edges = np.array(sorted(set(edges)))
    n_lo = int(np.floor((t.min() - t0) / P)) - 1
    n_hi = int(np.ceil((t.max() - t0) / P)) + 1
    rows = []
    for n in range(n_lo, n_hi + 1):
        tm = t0 + n * P
        if tm < t.min() - dur or tm > t.max() + dur:
            continue
        dt = np.abs(t - tm)
        inm = dt < dur / 2.0
        outm = (dt > dur / 2.0 + 0.05) & (dt < dur / 2.0 + 1.0)
        n_in, n_out = int(inm.sum()), int(outm.sum())
        exp_in = dur / CADENCE
        if n_in == 0:
            rows.append(dict(epoch=n, tmid=tm, n_in=0, n_out=n_out, frac_in=0.0,
                             depth_ppm=np.nan, err_ppm=np.nan,
                             gap_dist_d=float(np.min(np.abs(edges - tm)))))
            continue
        if n_out < 10:
            dep, err = np.nan, np.nan
        else:
            base = np.median(f[outm])
            sig = 1.4826 * np.median(np.abs(f[outm] - base))
            dep = base - np.mean(f[inm])
            err = sig * np.sqrt(1.0 / n_in + 1.0 / n_out)
        rows.append(dict(epoch=n, tmid=float(tm), n_in=n_in, n_out=n_out,
                         frac_in=float(n_in / exp_in),
                         depth_ppm=float(dep * 1e6) if np.isfinite(dep) else np.nan,
                         err_ppm=float(err * 1e6) if np.isfinite(err) else np.nan,
                         gap_dist_d=float(np.min(np.abs(edges - tm)))))
    return rows


def wmean(d, e):
    d, e = np.asarray(d, float), np.asarray(e, float)
    ok = np.isfinite(d) & np.isfinite(e) & (e > 0)
    if ok.sum() == 0:
        return np.nan, np.nan, 0
    w = 1.0 / e[ok] ** 2
    return float(np.sum(w * d[ok]) / np.sum(w)), float(1.0 / np.sqrt(np.sum(w))), int(ok.sum())


def kepler_year(q):
    # Q1-4 -> 1, Q5-8 -> 2, Q9-12 -> 3, Q13-17 -> 4
    return np.clip((np.asarray(q) - 1) // 4 + 1, 1, 4)


# ================================================================= main
T0 = time.time()
res = {"fname": FNAME, "date": "2026-10-17",
       "method": {"quality": "SAP_QUALITY==0 & finite", "normalise": "per-quarter median",
                  "detrend": f"time-windowed running median {WIN_D} d, segments split at gaps>{GAP_D} d, "
                             "clip upper 3.5 MAD / lower 8 MAD; second pass with transit masked",
                  "bls": f"astropy BoxLeastSquares, likelihood objective, {PMIN}-{PMAX} d, "
                         f"log-uniform frequency grid dlnf = dmin/({OVERSAMPLE}*baseline), "
                         f"durations {DURS.tolist()} d",
                  "depth": "box: mean(1*dur<|phi|<3*dur) - mean(|phi|<dur/2); error = robust std of "
                           f"{N_RANDOM} random-phase boxes avoiding transit and phase 0.5",
                  "no_pipeline_code": True}}

# ---- PSF flux
t, y_raw, e_raw, q = load("FLUX")
yn, en = normalise_per_quarter(t, y_raw, e_raw, q)
flat0, trend0, keep0 = detrend(t, yn)
m0 = keep0 & np.isfinite(flat0)
log(f"[load] N quality==0: {len(t)}  kept after detrend/clip: {m0.sum()}  "
    f"baseline {t.min():.2f}-{t.max():.2f} BKJD ({t.max()-t.min():.1f} d)")
scat = 1.4826 * np.median(np.abs(flat0[m0] - np.median(flat0[m0])))
log(f"[noise] per-cadence robust scatter of flattened PSF flux: {scat*1e6:.0f} ppm; "
    f"median quoted error {np.median(en[m0])*1e6:.0f} ppm")
res["n_points_quality0"] = int(len(t))
res["n_points_used"] = int(m0.sum())
res["baseline_d"] = float(t.max() - t.min())
res["scatter_ppm_per_cadence"] = float(scat * 1e6)

# ---- blind BLS, real
tt, ff, ee = t[m0], flat0[m0], en[m0]
freq = frequency_grid(tt, PMIN, PMAX, DURS.min(), OVERSAMPLE)
periods = np.sort(1.0 / freq)
log(f"[bls] {len(periods)} trial periods, {len(DURS)} durations")
bls = BoxLeastSquares(tt, ff, dy=ee)
t1 = time.time()
r = bls.power(periods, DURS, objective="likelihood")
log(f"[bls] real done in {time.time()-t1:.0f} s")
power = np.asarray(r.power, float)
sde_stat = spectrum_detrend(periods, power)
peaks_raw = top_peaks(periods, power, 3)
peaks_det = top_peaks(periods, sde_stat, 3)
raw_sde = sde_of(power)
log(f"[bls] top peaks (raw power): " + "; ".join(f"P={p:.5f} pow={s:.1f}" for p, s, _ in peaks_raw))
log(f"[bls] top peaks (median-detrended SDE): " + "; ".join(f"P={p:.5f} SDE={s:.2f}" for p, s, _ in peaks_det))
log(f"[bls] raw SDE of top peak = {raw_sde:.2f}")
ibest = peaks_raw[0][2]
P_bls = float(r.period[ibest]); t0_bls = float(r.transit_time[ibest])
dur_bls = float(r.duration[ibest]); dep_bls = float(r.depth[ibest])
log(f"[bls] best: P={P_bls:.5f} t0={t0_bls:.4f} dur={dur_bls*24:.2f} h depth={dep_bls*1e6:.0f} ppm")
alias = is_alias(P_bls, PIPE_P)
log(f"[bls] top peak vs pipeline {PIPE_P}: alias ratio = {alias}")

res["bls_blind"] = {
    "n_periods": int(len(periods)),
    "top3_raw_power": [{"period_d": p, "power": s, "alias_of_pipeline": is_alias(p, PIPE_P)} for p, s, _ in peaks_raw],
    "top3_detrended_sde": [{"period_d": p, "sde": s, "alias_of_pipeline": is_alias(p, PIPE_P)} for p, s, _ in peaks_det],
    "sde_raw": float(raw_sde),
    "sde_detrended": float(peaks_det[0][1]),
    "best_period_d": P_bls, "best_t0_bkjd": t0_bls, "best_duration_h": dur_bls * 24,
    "best_depth_ppm": dep_bls * 1e6,
    "top_peak_is_pipeline_period_or_alias": alias is not None,
    "alias_ratio": alias,
}

# ---- blind BLS, inverted null (same grid)
t1 = time.time()
rinv = bls.__class__(tt, 2.0 - ff, dy=ee).power(periods, DURS, objective="likelihood")
log(f"[bls] inverted done in {time.time()-t1:.0f} s")
pinv = np.asarray(rinv.power, float)
sde_inv = spectrum_detrend(periods, pinv)
peaks_inv = top_peaks(periods, sde_inv, 3)
log(f"[bls] inverted top peaks (detrended SDE): " + "; ".join(f"P={p:.5f} SDE={s:.2f}" for p, s, _ in peaks_inv))
res["bls_inverted_null"] = {"sde_raw": float(sde_of(pinv)),
                            "top3_detrended_sde": [{"period_d": p, "sde": s} for p, s, _ in peaks_inv]}

np.savetxt(os.path.join(OUT, "bls_periodogram.csv"),
           np.column_stack([periods, power, sde_stat, pinv, sde_inv])[::5],
           delimiter=",", header="period_d,power_real,sde_real,power_inverted,sde_inverted", comments="")

# ---- refine at the top peak: fine period + duration grid
Pc = P_bls
dP = 0.004 * Pc
fine_P = np.linspace(Pc - dP, Pc + dP, 1601)
fine_D = np.arange(0.06, 0.45, 0.01)
rf = bls.power(fine_P, fine_D, objective="likelihood")
j = int(np.argmax(rf.power))
P_ref, t0_ref, dur_ref, dep_ref = float(rf.period[j]), float(rf.transit_time[j]), float(rf.duration[j]), float(rf.depth[j])
log(f"[refine] P={P_ref:.5f} t0={t0_ref:.4f} dur={dur_ref*24:.2f} h depth={dep_ref*1e6:.0f} ppm")

# ---- re-detrend with transit masked (my own second pass)
in_tr = np.abs(phase_days(t, P_ref, t0_ref)) < 0.75 * dur_ref
flat, trend, keep = detrend(t, yn, transit_mask=in_tr)
m = keep & np.isfinite(flat)
# never clip in-transit points on the lower side (they were excluded from the trend)
m |= (in_tr & np.isfinite(flat))
tt, ff, ee, qq = t[m], flat[m], en[m], q[m]
log(f"[refine] masked re-detrend kept {m.sum()} points")

# fine duration scan on the masked-detrended flux: maximise box SNR
best = None
for dd in np.arange(0.08, 0.40, 0.005):
    for sh in np.linspace(-0.05, 0.05, 11):
        dep, nin = box_depth(tt, ff, P_ref, t0_ref + sh, dd)
        if not np.isfinite(dep):
            continue
        snr = dep * np.sqrt(nin)
        if best is None or snr > best[0]:
            best = (snr, dd, sh, dep)
dur_fin = best[1]; t0_fin = t0_ref + best[2]
log(f"[refine] box-SNR duration scan: dur={dur_fin*24:.2f} h, t0 shift {best[2]*24*60:.1f} min")
# keep BLS duration as the primary number, report the scan as a check
duration_h = dur_ref * 24.0

# ---- primary box depth + empirical error
depth, n_in = box_depth(tt, ff, P_ref, t0_ref, dur_ref)
err_mad, err_std, mean_rand, rand_vals = random_phase_error(tt, ff, P_ref, t0_ref, dur_ref)
snr_emp = depth / err_mad
log(f"[depth] PSF box depth = {depth*1e6:.0f} +- {err_mad*1e6:.0f} ppm (random-phase MAD; std {err_std*1e6:.0f}); "
    f"n_in={n_in}; empirical SNR={snr_emp:.1f}; random-box mean {mean_rand*1e6:.0f} ppm")

# ---- secondary at phase 0.5
sec, n_sec = box_depth(tt, ff, P_ref, t0_ref + P_ref / 2.0, dur_ref)
sec_sigma = sec / err_mad
log(f"[secondary] phase-0.5 box depth = {sec*1e6:.0f} ppm ({sec_sigma:.2f} sigma), n_in={n_sec}")

# ---- per-event
events = per_event(tt, ff, P_ref, t0_ref, dur_ref)
obs = [ev for ev in events if ev["frac_in"] >= 0.5 and np.isfinite(ev["depth_ppm"])]
n_transits = len(obs)
log(f"[events] {len(events)} epochs in span, {n_transits} with >=50% in-transit coverage")
for ev in events:
    log("   epoch %3d tmid %9.3f n_in %2d frac %.2f depth %7.0f +- %5.0f ppm  gap_dist %.2f d"
        % (ev["epoch"], ev["tmid"], ev["n_in"], ev["frac_in"],
           ev["depth_ppm"] if np.isfinite(ev["depth_ppm"]) else -1,
           ev["err_ppm"] if np.isfinite(ev["err_ppm"]) else -1, ev["gap_dist_d"]))
ev_mean, ev_err, ev_n = wmean([e["depth_ppm"] for e in obs], [e["err_ppm"] for e in obs])
log(f"[events] weighted mean of per-event depths = {ev_mean:.0f} +- {ev_err:.0f} ppm (n={ev_n})")
with open(os.path.join(OUT, "events.csv"), "w") as fh:
    fh.write("epoch,tmid_bkjd,n_in,n_out,frac_in,depth_ppm,err_ppm,gap_dist_d\n")
    for ev in events:
        fh.write("%d,%.5f,%d,%d,%.3f,%.1f,%.1f,%.3f\n" % (ev["epoch"], ev["tmid"], ev["n_in"], ev["n_out"],
                                                         ev["frac_in"], ev["depth_ppm"], ev["err_ppm"], ev["gap_dist_d"]))

# ---- odd / even (per-event weighted means, and box statistic on subsets)
odd = [e for e in obs if e["epoch"] % 2 == 1]
even = [e for e in obs if e["epoch"] % 2 == 0]
d_odd, e_odd, n_odd = wmean([e["depth_ppm"] for e in odd], [e["err_ppm"] for e in odd])
d_even, e_even, n_even = wmean([e["depth_ppm"] for e in even], [e["err_ppm"] for e in even])
oe_sigma = abs(d_odd - d_even) / np.sqrt(e_odd ** 2 + e_even ** 2)
log(f"[odd/even] odd {d_odd:.0f}+-{e_odd:.0f} (n={n_odd})  even {d_even:.0f}+-{e_even:.0f} (n={n_even})  -> {oe_sigma:.2f} sigma")
# box version with 2P folding
ep = np.round((tt - t0_ref) / P_ref).astype(int)
sel_odd = (ep % 2 == 1); sel_even = (ep % 2 == 0)
d_odd_box, _ = box_depth(tt[sel_odd], ff[sel_odd], P_ref, t0_ref, dur_ref)
d_even_box, _ = box_depth(tt[sel_even], ff[sel_even], P_ref, t0_ref, dur_ref)
e_odd_box = random_phase_error(tt[sel_odd], ff[sel_odd], P_ref, t0_ref, dur_ref, n=1000)[0]
e_even_box = random_phase_error(tt[sel_even], ff[sel_even], P_ref, t0_ref, dur_ref, n=1000)[0]
oe_sigma_box = abs(d_odd_box - d_even_box) / np.sqrt(e_odd_box ** 2 + e_even_box ** 2)
log(f"[odd/even] box: odd {d_odd_box*1e6:.0f}+-{e_odd_box*1e6:.0f}  even {d_even_box*1e6:.0f}+-{e_even_box*1e6:.0f} -> {oe_sigma_box:.2f} sigma")

# ---- per-year
yr = kepler_year(qq)
per_year = []
for Y in (1, 2, 3, 4):
    s = yr == Y
    dY, nY = box_depth(tt[s], ff[s], P_ref, t0_ref, dur_ref)
    eY = random_phase_error(tt[s], ff[s], P_ref, t0_ref, dur_ref, n=800)[0] if nY >= 3 else np.nan
    nev = sum(1 for e in obs if Y == kepler_year(np.array([qq[np.argmin(np.abs(tt - e['tmid']))]]))[0])
    per_year.append({"year": Y, "depth_ppm": float(dY * 1e6) if np.isfinite(dY) else None,
                     "err_ppm": float(eY * 1e6) if np.isfinite(eY) else None, "n_in": int(nY), "n_events": int(nev)})
    log(f"[year {Y}] depth {dY*1e6 if np.isfinite(dY) else float('nan'):.0f} +- {eY*1e6 if np.isfinite(eY) else float('nan'):.0f} ppm  n_in={nY}  events={nev}")

# ---- depth with transit-unmasked detrending for reference
depth_unmasked, _ = box_depth(t[m0], flat0[m0], P_ref, t0_ref, dur_ref)
log(f"[depth] with unmasked 1-d running median: {depth_unmasked*1e6:.0f} ppm (masked: {depth*1e6:.0f})")

# ---- SAP flux
ts, ys_raw, es_raw, qs = load("SAP_FLUX")
ysn, esn = normalise_per_quarter(ts, ys_raw, es_raw, qs)
in_tr_s = np.abs(phase_days(ts, P_ref, t0_ref)) < 0.75 * dur_ref
flat_s, trend_s, keep_s = detrend(ts, ysn, transit_mask=in_tr_s)
ms = (keep_s & np.isfinite(flat_s)) | (in_tr_s & np.isfinite(flat_s))
tts, ffs = ts[ms], flat_s[ms]
scat_s = 1.4826 * np.median(np.abs(ffs - np.median(ffs)))
dep_s, nin_s = box_depth(tts, ffs, P_ref, t0_ref, dur_ref)
err_s = random_phase_error(tts, ffs, P_ref, t0_ref, dur_ref, n=2000)[0]
ratio = dep_s / depth
ratio_err = ratio * np.sqrt((err_s / dep_s) ** 2 + (err_mad / depth) ** 2) if dep_s != 0 else np.nan
log(f"[SAP] scatter {scat_s*1e6:.0f} ppm; SAP box depth = {dep_s*1e6:.0f} +- {err_s*1e6:.0f} ppm (SNR {dep_s/err_s:.1f}); "
    f"SAP/PSF = {ratio:.2f} +- {ratio_err:.2f}")
# aperture flux budget: how much flux does the SAP aperture hold relative to the PSF flux of the target?
fr_q = []
for qq_ in np.unique(qs):
    a = qs == qq_
    b = q == qq_
    fr_q.append(float(np.median(ys_raw[a]) / np.median(y_raw[b])))
fr_q = np.array(fr_q)
FLFRCSAP_CSV = 0.8577   # from targets/poolA_kbonus_bg_dwarfs_ranked.csv (KBONUS flux fraction of target in aperture)
exp_ratio_on_target = FLFRCSAP_CSV / np.median(fr_q)
log(f"[SAP] median(SAP_FLUX)/median(PSF FLUX) per quarter: min {fr_q.min():.2f} med {np.median(fr_q):.2f} max {fr_q.max():.2f}")
log(f"[SAP] if the dip is on the target and the aperture holds flfrcsap={FLFRCSAP_CSV} of its PSF flux, "
    f"expected SAP/PSF depth ratio ~ {exp_ratio_on_target:.2f} (KBONUS crowdsap in CSV: 0.8307)")
res["sap_flux_budget"] = {"sap_over_psf_median_flux_per_quarter": fr_q.tolist(),
                          "expected_depth_ratio_if_on_target_flfrcsap_over_fluxratio": float(exp_ratio_on_target),
                          "csv_crowdsap": 0.8307, "csv_flfrcsap": FLFRCSAP_CSV}
# SAP per-event depths too
ev_s = per_event(tts, ffs, P_ref, t0_ref, dur_ref)
obs_s = [e for e in ev_s if e["frac_in"] >= 0.5 and np.isfinite(e["depth_ppm"])]
evs_mean, evs_err, _ = wmean([e["depth_ppm"] for e in obs_s], [e["err_ppm"] for e in obs_s])
log(f"[SAP] weighted per-event mean depth {evs_mean:.0f} +- {evs_err:.0f} ppm")

# ---- summary numbers
res["ephemeris"] = {"period_d": P_ref, "t0_bkjd": t0_ref, "t0_bjd": t0_ref + 2454833.0,
                    "duration_h_bls": duration_h, "duration_h_boxsnr_scan": dur_fin * 24.0,
                    "bls_depth_ppm": dep_ref * 1e6}
res["psf_depth"] = {"depth_ppm": depth * 1e6, "err_ppm_random_phase_mad": err_mad * 1e6,
                    "err_ppm_random_phase_std": err_std * 1e6, "empirical_snr": snr_emp,
                    "n_in_transit_cadences": n_in, "random_box_mean_ppm": mean_rand * 1e6,
                    "depth_ppm_unmasked_detrend": depth_unmasked * 1e6,
                    "per_event_weighted_mean_ppm": ev_mean, "per_event_weighted_err_ppm": ev_err}
res["secondary"] = {"depth_ppm": sec * 1e6, "sigma": sec_sigma, "n_in": n_sec}
res["odd_even"] = {"odd_ppm": d_odd, "odd_err_ppm": e_odd, "n_odd": n_odd,
                   "even_ppm": d_even, "even_err_ppm": e_even, "n_even": n_even,
                   "sigma_per_event": oe_sigma,
                   "odd_box_ppm": d_odd_box * 1e6, "even_box_ppm": d_even_box * 1e6,
                   "odd_box_err_ppm": e_odd_box * 1e6, "even_box_err_ppm": e_even_box * 1e6,
                   "sigma_box": oe_sigma_box}
res["per_year"] = per_year
res["n_transits_observed"] = n_transits
res["n_epochs_in_span"] = len(events)
res["events"] = events
res["sap"] = {"depth_ppm": dep_s * 1e6, "expected_sap_over_psf_if_on_target": float(exp_ratio_on_target), "err_ppm": err_s * 1e6, "snr": dep_s / err_s,
              "scatter_ppm": scat_s * 1e6, "per_event_weighted_mean_ppm": evs_mean,
              "per_event_weighted_err_ppm": evs_err, "sap_over_psf": ratio, "sap_over_psf_err": ratio_err}

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(3, 2, figsize=(12, 12))
fig.patch.set_facecolor("#fcfcfb")
# (a) periodogram
ax = axes[0, 0]
ax.plot(periods, sde_stat, color=C_BLUE, lw=0.6)
for p, s, _ in peaks_det:
    ax.annotate(f"{p:.3f} d", (p, s), textcoords="offset points", xytext=(4, 2), fontsize=8)
ax.axvline(PIPE_P, color=C_ORANGE, lw=0.8, ls="--", label=f"pipeline {PIPE_P} d")
ax.set_xscale("log"); ax.set_xlabel("period (d)"); ax.set_ylabel("BLS SDE (median-detrended)")
ax.set_title("Blind BLS 1-50 d, PSF flux"); ax.legend(fontsize=8, frameon=False)
# (b) inverted
ax = axes[0, 1]
ax.plot(periods, sde_inv, color=C_GREY, lw=0.6)
ax.set_xscale("log"); ax.set_xlabel("period (d)"); ax.set_ylabel("BLS SDE (median-detrended)")
ax.set_title("Inverted-flux null (2 - f), same grid")
ax.set_ylim(axes[0, 0].get_ylim())
# (c) fold
ax = axes[1, 0]
ph = phase_days(tt, P_ref, t0_ref) * 24.0
s = np.abs(ph) < 24.0
ax.plot(ph[s], (ff[s] - 1) * 1e6, ".", color=C_GREY, ms=2, alpha=0.5)
bins = np.arange(-24, 24.01, 0.5)
idx = np.digitize(ph[s], bins)
bx = 0.5 * (bins[1:] + bins[:-1])
by = np.array([np.mean(ff[s][idx == i + 1]) if np.any(idx == i + 1) else np.nan for i in range(len(bx))])
ax.plot(bx, (by - 1) * 1e6, "o-", color=C_BLUE, ms=4, lw=1.2, label="30-min bins")
ax.axvspan(-dur_ref * 12, dur_ref * 12, color=C_ORANGE, alpha=0.12, lw=0)
ax.axhline(-depth * 1e6, color=C_ORANGE, lw=1, ls="--", label=f"box depth {depth*1e6:.0f} ppm")
ax.set_ylim(-8000, 6000); ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("flux - 1 (ppm)")
ax.set_title(f"PSF flux folded at P={P_ref:.5f} d"); ax.legend(fontsize=8, frameon=False)
# (d) odd / even
ax = axes[1, 1]
for sel, col, lab in ((sel_odd, C_BLUE, "odd"), (sel_even, C_ORANGE, "even")):
    phs = phase_days(tt[sel], P_ref, t0_ref) * 24.0
    ss = np.abs(phs) < 24.0
    idx = np.digitize(phs[ss], bins)
    byy = np.array([np.mean(ff[sel][ss][idx == i + 1]) if np.any(idx == i + 1) else np.nan for i in range(len(bx))])
    ax.plot(bx, (byy - 1) * 1e6, "o-", color=col, ms=3.5, lw=1.1, label=lab)
ax.set_ylim(-6000, 4000); ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("flux - 1 (ppm)")
ax.set_title(f"odd vs even, {oe_sigma:.1f} sigma (per-event)"); ax.legend(fontsize=8, frameon=False)
# (e) secondary + SAP fold
ax = axes[2, 0]
phsec = phase_days(tt, P_ref, t0_ref + P_ref / 2) * 24.0
s2 = np.abs(phsec) < 24.0
idx = np.digitize(phsec[s2], bins)
by2 = np.array([np.mean(ff[s2][idx == i + 1]) if np.any(idx == i + 1) else np.nan for i in range(len(bx))])
ax.plot(bx, (by2 - 1) * 1e6, "o-", color=C_BLUE, ms=3.5, lw=1.1, label=f"phase 0.5: {sec*1e6:.0f} ppm")
phS = phase_days(tts, P_ref, t0_ref) * 24.0
s3 = np.abs(phS) < 24.0
idx = np.digitize(phS[s3], bins)
by3 = np.array([np.mean(ffs[s3][idx == i + 1]) if np.any(idx == i + 1) else np.nan for i in range(len(bx))])
ax.plot(bx, (by3 - 1) * 1e6, "s-", color=C_AQUA, ms=3.5, lw=1.1, label=f"SAP flux at transit: {dep_s*1e6:.0f} ppm")
ax.set_ylim(-6000, 4000); ax.set_xlabel("hours from phase 0.5 / mid-transit"); ax.set_ylabel("flux - 1 (ppm)")
ax.set_title("secondary (PSF) and SAP-flux primary"); ax.legend(fontsize=8, frameon=False)
# (f) per-event depths
ax = axes[2, 1]
tm = np.array([e["tmid"] for e in events]); dpp = np.array([e["depth_ppm"] for e in events])
epp = np.array([e["err_ppm"] for e in events]); fr = np.array([e["frac_in"] for e in events])
okk = np.isfinite(dpp)
ax.errorbar(tm[okk & (fr >= 0.5)], dpp[okk & (fr >= 0.5)], yerr=epp[okk & (fr >= 0.5)], fmt="o", color=C_BLUE, ms=5, lw=1, label="covered events")
ax.errorbar(tm[okk & (fr < 0.5)], dpp[okk & (fr < 0.5)], yerr=epp[okk & (fr < 0.5)], fmt="o", mfc="none", color=C_GREY, ms=5, lw=1, label="partial (<50%)")
ax.axhline(depth * 1e6, color=C_ORANGE, lw=1, ls="--", label="global box depth")
ax.axhline(0, color=C_GREY, lw=0.6)
ax.set_xlabel("BKJD (d)"); ax.set_ylabel("local depth (ppm)"); ax.set_title(f"per-event depths ({n_transits} covered)")
ax.legend(fontsize=8, frameon=False)
fig.suptitle(f"{FNAME}: independent re-derivation (no pipeline code)  P={P_ref:.5f} d  depth={depth*1e6:.0f}+-{err_mad*1e6:.0f} ppm  dur={duration_h:.2f} h", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(os.path.join(OUT, "fold.png"), dpi=130)
log(f"[fig] wrote fold.png")

res["runtime_s"] = time.time() - T0
res["log"] = log_lines


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


with open(os.path.join(OUT, "results.json"), "w") as fh:
    json.dump(_clean(res), fh, indent=1)
log(f"[done] {res['runtime_s']:.0f} s; results.json written")
