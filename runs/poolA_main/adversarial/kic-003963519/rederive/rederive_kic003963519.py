#!/usr/bin/env python
"""
Independent re-derivation for kic-003963519 (KBONUS-BKG PSF light curve).

Written from scratch with numpy / scipy / astropy only. Deliberately shares NO code
with /home/user/ZyushgPlanet/pipeline/ (not imported, not read before this script
produced its numbers). No wotan, no TLS: detrending is a time-based running median.

Steps
 1. load LIGHTCURVE_STITCHED, SAP_QUALITY == 0, per-quarter median normalisation,
    running-median detrend (window W d) with upward-outlier clipping only
 2. blind astropy BoxLeastSquares 1-50 d on my own period grid
 3. at the best ephemeris: box depth + empirical (random-phase) error, odd/even,
    secondary at phase 0.5, per-year depths, number of transits, duration
 4. same depth measurement on SAP_FLUX (dilution check)
 5. inverted-flux (2 - f) null on the same grid (cheap sanity check on the SDE)

Usage: python rederive_kic003963519.py [--quick]   (--quick = coarse grid, for timing)
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares

FITS_PATH = ('/home/user/ZyushgPlanet/data/kbonus_batch/'
             'hlsp_kbonus-bkg_kepler_kepler_kic-003963519_kepler_v1.0_lc.fits')
OUT_DIR = '/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-003963519/rederive'
FNAME = 'kic-003963519'
BKJD0 = 2454833.0
PIPE_P = 6.81560          # pipeline value, used ONLY to ask "is it the top peak?"
RNG = np.random.default_rng(20261017)


# ----------------------------------------------------------------------------- loading
def load(path, flux_col='FLUX'):
    with fits.open(path) as h:
        d = h['LIGHTCURVE_STITCHED'].data
        t = np.asarray(d['TIME'], float)
        f = np.asarray(d[flux_col], float)
        fe = np.asarray(d[flux_col + '_ERR'], float)
        q = np.asarray(d['SAP_QUALITY'], int)
        quarter = np.asarray(d['QUARTER'], int)
    good = (q == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (f > 0)
    return t[good], f[good], fe[good], quarter[good]


def running_median(t, y, window):
    """Time-based running median: for each point, median of y within +-window/2 in time."""
    order = np.argsort(t)
    ts, ys = t[order], y[order]
    lo = np.searchsorted(ts, ts - window / 2.0, side='left')
    hi = np.searchsorted(ts, ts + window / 2.0, side='right')
    out = np.empty_like(ys)
    for i in range(len(ts)):
        out[i] = np.median(ys[lo[i]:hi[i]])
    res = np.empty_like(out)
    res[order] = out
    return res


def detrend(t, f, fe, quarter, window=1.0, clip_sigma=4.0):
    """Per-quarter median normalisation, running-median flattening, upward-outlier clip.
    Returns t, flat flux (around 1), flat error, quarter, and the trend (for reference)."""
    tt, ff, ee, qq = [], [], [], []
    for Q in np.unique(quarter):
        m = quarter == Q
        tq, fq, eq = t[m], f[m], fe[m]
        med = np.median(fq)
        fq = fq / med
        eq = eq / med
        # two-pass: first trend, drop upward outliers, second trend
        tr = running_median(tq, fq, window)
        r = fq / tr - 1.0
        mad = 1.4826 * np.median(np.abs(r - np.median(r)))
        keep = r < clip_sigma * mad          # only clip UPWARD outliers (never dips)
        tr2 = running_median(tq[keep], fq[keep], window)
        tr_full = np.interp(tq, tq[keep], tr2)
        flat = fq / tr_full
        tt.append(tq[keep]); ff.append(flat[keep]); ee.append(eq[keep] / tr_full[keep]); qq.append(np.full(keep.sum(), Q))
    t = np.concatenate(tt); f = np.concatenate(ff); e = np.concatenate(ee); q = np.concatenate(qq)
    o = np.argsort(t)
    return t[o], f[o], e[o], q[o]


# ----------------------------------------------------------------------------- search
def period_grid(pmin, pmax, baseline, dur_ref=0.08, k=2.0):
    """Geometric grid: successive periods shift the last transit by dur_ref/k over the baseline."""
    ps = [pmin]
    while ps[-1] < pmax:
        p = ps[-1]
        ps.append(p * (1.0 + dur_ref / (baseline * k)))
    return np.array(ps)


def sde_from_power(power, periods, smooth_frac=0.05):
    """SDE two ways: raw (max-mean)/std and after removing a running-median baseline in log-period."""
    pw = np.asarray(power, float)
    raw = (pw.max() - pw.mean()) / pw.std()
    # baseline: running median over a window in index space (grid is ~geometric)
    n = len(pw)
    w = max(51, int(smooth_frac * n) | 1)
    idx = np.arange(n)
    lo = np.clip(idx - w // 2, 0, n); hi = np.clip(idx + w // 2 + 1, 0, n)
    base = np.array([np.median(pw[lo[i]:hi[i]]) for i in range(0, n)])
    resid = pw - base
    det = (resid.max() - resid.mean()) / resid.std()
    return float(raw), float(det), resid


def top_peaks(periods, power, n=3, sep=0.02):
    """Greedy distinct peaks: the next peak must differ in period by > sep (fractional)."""
    order = np.argsort(power)[::-1]
    picked = []
    for i in order:
        if all(abs(periods[i] / periods[j] - 1) > sep for j in picked):
            picked.append(i)
        if len(picked) >= n:
            break
    return picked


def alias_of(p, pref, tol=0.002):
    for num, den in [(1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (3, 2), (2, 3), (4, 1), (1, 4), (5, 1), (1, 5)]:
        if abs(p / (pref * num / den) - 1) < tol:
            return f'{num}/{den} x P_pipe' if (num, den) != (1, 1) else 'same'
    return None


def run_bls(t, f, e, periods, durations, objective='likelihood'):
    bls = BoxLeastSquares(t, f, dy=e)
    res = bls.power(periods, durations, objective=objective)
    return bls, res


# ----------------------------------------------------------------------------- statistics at an ephemeris
def phase_fold(t, P, t0):
    ph = ((t - t0) / P + 0.5) % 1.0 - 0.5
    return ph


def box_depth(t, f, P, t0, dur, oot_halfwidth=None):
    """Depth = mean(out-of-transit) - mean(in-transit) in relative flux (positive = dip).
    Out-of-transit = everything with |phase*P| > dur (i.e. outside 2x box), or a local window."""
    ph = phase_fold(t, P, t0) * P
    inm = np.abs(ph) < dur / 2.0
    if oot_halfwidth is None:
        oot = np.abs(ph) > dur
    else:
        oot = (np.abs(ph) > dur) & (np.abs(ph) < oot_halfwidth)
    if inm.sum() < 2 or oot.sum() < 10:
        return np.nan, np.nan, int(inm.sum())
    d = np.mean(f[oot]) - np.mean(f[inm])
    # formal error from point scatter
    s = np.std(f[oot])
    err = s * np.sqrt(1.0 / inm.sum() + 1.0 / oot.sum())
    return d, err, int(inm.sum())


def empirical_depth_error(t, f, P, t0, dur, n=3000, avoid=(0.0, 0.5), avoid_width=None):
    """Std of box depths at random phases, avoiding the transit and the secondary windows."""
    if avoid_width is None:
        avoid_width = 2.0 * dur / P
    ph = phase_fold(t, P, t0)      # in [-0.5, 0.5)
    vals = []
    tries = 0
    while len(vals) < n and tries < 20 * n:
        tries += 1
        c = RNG.uniform(-0.5, 0.5)
        if any(min(abs(c - a), 1 - abs(c - a)) < avoid_width for a in avoid):
            continue
        dph = (ph - c + 0.5) % 1.0 - 0.5
        inm = np.abs(dph) * P < dur / 2.0
        oot = np.abs(dph) * P > dur
        if inm.sum() < 5:
            continue
        vals.append(np.mean(f[oot]) - np.mean(f[inm]))
    vals = np.array(vals)
    return float(np.std(vals)), float(np.mean(vals)), len(vals)


def per_event(t, f, q, P, t0, dur, local_halfwidth_factor=3.0):
    """Per-transit local depth: baseline from +-3 dur around mid-transit excluding +-dur."""
    n_min = int(np.floor((t.min() - t0) / P)) - 1
    n_max = int(np.ceil((t.max() - t0) / P)) + 1
    events = []
    for n in range(n_min, n_max + 1):
        tc = t0 + n * P
        if tc < t.min() - dur or tc > t.max() + dur:
            continue
        dt = t - tc
        inm = np.abs(dt) < dur / 2.0
        loc = (np.abs(dt) > dur) & (np.abs(dt) < local_halfwidth_factor * dur)
        n_in = int(inm.sum())
        if n_in == 0:
            events.append(dict(epoch=n, tc_bkjd=tc, n_in=0, depth=np.nan, err=np.nan, quarter=None))
            continue
        if loc.sum() < 5:
            events.append(dict(epoch=n, tc_bkjd=tc, n_in=n_in, depth=np.nan, err=np.nan,
                               quarter=int(np.median(q[inm]))))
            continue
        base = np.median(f[loc])
        d = base - np.mean(f[inm])
        s = np.std(f[loc])
        err = s * np.sqrt(1.0 / n_in + 1.0 / loc.sum())
        events.append(dict(epoch=n, tc_bkjd=tc, n_in=n_in, depth=d, err=err, quarter=int(np.median(q[inm]))))
    return events


def trapezoid(ph_days, depth, t14, t23_frac):
    """Symmetric trapezoid in time-from-midtransit (days). t23 = t23_frac * t14."""
    t23 = t23_frac * t14
    x = np.abs(ph_days)
    m = np.ones_like(x)
    full = x <= t23 / 2
    ingress = (x > t23 / 2) & (x < t14 / 2)
    m[full] = 1 - depth
    m[ingress] = 1 - depth * (t14 / 2 - x[ingress]) / ((t14 - t23) / 2 + 1e-12)
    return m


def fit_trapezoid(t, f, P, t0, dur_guess, depth_guess):
    """Least-squares trapezoid fit near transit; returns depth, T14, T23/T14, dt0."""
    from scipy.optimize import least_squares
    ph = phase_fold(t, P, t0) * P
    win = np.abs(ph) < 4 * dur_guess
    x, y = ph[win], f[win]

    def resid(p):
        depth, t14, t23f, dt0 = p
        return trapezoid(x - dt0, depth, t14, t23f) - y

    p0 = [depth_guess, dur_guess, 0.6, 0.0]
    lb = [0.0, 0.3 * dur_guess, 0.0, -0.5 * dur_guess]
    ub = [10 * depth_guess + 1e-3, 3.0 * dur_guess, 1.0, 0.5 * dur_guess]
    r = least_squares(resid, p0, bounds=(lb, ub))
    return dict(depth=float(r.x[0]), t14_d=float(r.x[1]), t23_over_t14=float(r.x[2]), dt0_d=float(r.x[3]),
                chi2=float(np.sum(r.fun ** 2)), n=int(win.sum()))


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--window', type=float, default=1.0, help='running-median window (d)')
    ap.add_argument('--no-null', action='store_true')
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    T0 = time.time()
    log = lambda *a: print(f'[{time.time()-T0:7.1f}s]', *a, flush=True)

    # ---------- 1. load + detrend (PSF flux)
    t_raw, f_raw, e_raw, q_raw = load(FITS_PATH, 'FLUX')
    log(f'loaded {len(t_raw)} quality==0 cadences, quarters {sorted(int(x) for x in set(q_raw))}')
    t, f, e, q = detrend(t_raw, f_raw, e_raw, q_raw, window=args.window)
    baseline = t.max() - t.min()
    rms = float(np.std(f))
    mad_rms = float(1.4826 * np.median(np.abs(f - np.median(f))))
    log(f'detrended: {len(t)} pts, baseline {baseline:.1f} d, rms {rms*1e6:.0f} ppm, robust rms {mad_rms*1e6:.0f} ppm')

    # ---------- 2. blind BLS 1-50 d
    pmin, pmax = 1.0, 50.0
    if args.quick:
        periods = period_grid(pmin, pmax, baseline, dur_ref=0.08, k=0.25)
        durations = np.array([0.05, 0.1, 0.2])
    else:
        periods = period_grid(pmin, pmax, baseline, dur_ref=0.08, k=2.0)
        durations = np.array([0.04, 0.06, 0.08, 0.11, 0.15, 0.20, 0.28])
    log(f'BLS grid: {len(periods)} periods, durations {durations}')
    bls, res = run_bls(t, f, e, periods, durations)
    log('BLS done')
    power = np.asarray(res.power)
    sde_raw, sde_det, resid = sde_from_power(power, periods)
    picks = top_peaks(periods, power, n=3)
    picks_det = top_peaks(periods, resid, n=3)
    peaks = []
    for rank, i in enumerate(picks, 1):
        peaks.append(dict(rank=rank, period_d=float(periods[i]), power=float(power[i]),
                          depth_ppm=float(res.depth[i] * 1e6), duration_d=float(res.duration[i]),
                          t0_bkjd=float(res.transit_time[i]), depth_snr=float(res.depth_snr[i]),
                          relation_to_pipeline_P=alias_of(periods[i], PIPE_P)))
        log(f'peak {rank}: P={periods[i]:.5f} d power={power[i]:.2f} depth={res.depth[i]*1e6:.0f} ppm '
            f'dur={res.duration[i]*24:.2f} h  rel={alias_of(periods[i], PIPE_P)}')
    peaks_det = [dict(rank=r + 1, period_d=float(periods[i]), resid_power=float(resid[i]),
                      relation_to_pipeline_P=alias_of(periods[i], PIPE_P)) for r, i in enumerate(picks_det)]
    log(f'SDE raw={sde_raw:.2f}  SDE baseline-subtracted={sde_det:.2f}')
    best = picks[0]
    P_coarse = float(periods[best])
    recovered = alias_of(P_coarse, PIPE_P) == 'same'
    recovered_alias = alias_of(P_coarse, PIPE_P)

    # ---------- 2b. refine around the top peak AND its 1/2 and 2x aliases (fine period, fine duration)
    # The coarse grid cannot separate P from 2P when depths are equal; decide with a fine grid.
    fine_dur = np.linspace(0.03, 0.30, 28)
    alias_tests = []
    for mult in (1.0, 0.5, 2.0):
        Pc = P_coarse * mult
        if not (pmin <= Pc <= pmax):
            continue
        fP = np.linspace(Pc * (1 - 0.003), Pc * (1 + 0.003), 1500)
        _, rr = run_bls(t, f, e, fP, fine_dur)
        jj = int(np.argmax(rr.power))
        alias_tests.append(dict(mult=mult, period_d=float(fP[jj]), power=float(rr.power[jj]),
                                depth_ppm=float(rr.depth[jj] * 1e6), duration_h=float(rr.duration[jj] * 24),
                                t0_bkjd=float(rr.transit_time[jj]), _fP=fP, _rr=rr, _j=jj))
        log(f'fine grid at {mult:g} x top peak: P={fP[jj]:.6f} power={rr.power[jj]:.2f} depth={rr.depth[jj]*1e6:.0f} ppm dur={rr.duration[jj]*24:.2f} h')
    bestA = max(alias_tests, key=lambda a: a['power'])
    fine_P, rfine, j = bestA['_fP'], bestA['_rr'], bestA['_j']
    adopted_mult = bestA['mult']
    for a in alias_tests:
        for kk in ('_fP', '_rr', '_j'):
            a.pop(kk)
    P = float(fine_P[j]); dur = float(rfine.duration[j]); t0 = float(rfine.transit_time[j])
    # move t0 to the epoch nearest the data midpoint (minimises P-t0 covariance)
    n_mid = np.round((0.5 * (t.min() + t.max()) - t0) / P)
    t0_mid = t0 + n_mid * P
    log(f'refined: P={P:.6f} d, dur={dur*24:.2f} h, t0={t0_mid:.4f} BKJD, BLS depth={rfine.depth[j]*1e6:.0f} ppm')
    # period uncertainty: half-width of the fine peak at power - 1 (delta log-likelihood = 1)
    above = fine_P[rfine.power >= rfine.power[j] - 1.0]
    P_err = float(0.5 * (above.max() - above.min())) if len(above) > 1 else float(fine_P[1] - fine_P[0])

    # ---------- 3. statistics at the ephemeris (PSF flux)
    d, d_err_formal, n_in = box_depth(t, f, P, t0_mid, dur)
    d_err_emp, d_mean_emp, n_emp = empirical_depth_error(t, f, P, t0_mid, dur)
    snr_emp = d / d_err_emp
    log(f'box depth {d*1e6:.0f} +- {d_err_emp*1e6:.0f} ppm (empirical, {n_emp} random phases; formal +- {d_err_formal*1e6:.0f}); '
        f'n_in={n_in}; SNR_emp={snr_emp:.2f}')
    # scaled-empirical errors for subsets: err_sub = err_emp * sqrt(n_in_all / n_in_sub)
    ph = phase_fold(t, P, t0_mid) * P
    epoch = np.round((t - t0_mid) / P).astype(int)
    inm = np.abs(ph) < dur / 2
    oot = np.abs(ph) > dur
    odd = inm & (epoch % 2 != 0); even = inm & (epoch % 2 == 0)
    base = np.mean(f[oot])
    d_odd = base - np.mean(f[odd]); d_even = base - np.mean(f[even])
    e_odd = d_err_emp * np.sqrt(n_in / max(odd.sum(), 1)); e_even = d_err_emp * np.sqrt(n_in / max(even.sum(), 1))
    oe_sigma = abs(d_odd - d_even) / np.sqrt(e_odd ** 2 + e_even ** 2)
    log(f'odd {d_odd*1e6:.0f}+-{e_odd*1e6:.0f} ({odd.sum()} pts)  even {d_even*1e6:.0f}+-{e_even*1e6:.0f} ({even.sum()} pts)  -> {oe_sigma:.2f} sigma')
    # secondary at phase 0.5
    d_sec, _, n_sec = box_depth(t, f, P, t0_mid + P / 2, dur)
    e_sec = d_err_emp * np.sqrt(n_in / max(n_sec, 1))
    sec_sigma = d_sec / e_sec
    log(f'secondary at phase 0.5: {d_sec*1e6:.0f} +- {e_sec*1e6:.0f} ppm ({n_sec} pts) -> {sec_sigma:.2f} sigma')
    # if the phase-0.5 'secondary' is significant AND consistent with the primary depth, P/2 is the likely true period
    half_period_flag = bool(sec_sigma > 3 and abs(d_sec - d) / np.sqrt(e_sec ** 2 + d_err_emp ** 2) < 2)
    log(f'half-period flag (secondary significant and equal to primary): {half_period_flag}')
    # strongest secondary anywhere in phase 0.2-0.8 (post-hoc, descriptive)
    sec_scan = []
    for c in np.linspace(0.15, 0.85, 141):
        dd, _, nn = box_depth(t, f, P, t0_mid + c * P, dur)
        sec_scan.append((c, dd, nn))
    sec_scan = np.array(sec_scan)
    k = int(np.nanargmax(sec_scan[:, 1]))
    # per-year depths (Kepler years: Q1-4, Q5-8, Q9-12, Q13-17)
    year_of = lambda Q: 1 if Q <= 4 else 2 if Q <= 8 else 3 if Q <= 12 else 4
    yrs = np.array([year_of(Q) for Q in q])
    per_year = []
    for Y in (1, 2, 3, 4):
        m = yrs == Y
        if m.sum() == 0:
            continue
        dy_, _, ny = box_depth(t[m], f[m], P, t0_mid, dur)
        ey = d_err_emp * np.sqrt(n_in / max(ny, 1))
        per_year.append(dict(year=Y, quarters=sorted(set(int(x) for x in q[m])), depth_ppm=float(dy_ * 1e6),
                             err_ppm=float(ey * 1e6), n_in=int(ny)))
        log(f'year {Y}: depth {dy_*1e6:.0f} +- {ey*1e6:.0f} ppm ({ny} in-transit pts)')
    # per-event
    events = per_event(t, f, q, P, t0_mid, dur)
    expected_in = dur / np.median(np.diff(t[:200]))      # ~ dur / 29.4 min
    observed = [ev for ev in events if ev['n_in'] >= max(2, 0.5 * expected_in)]
    n_transits = len(observed)
    n_epochs_in_baseline = int(np.floor((t.max() - t0_mid) / P) - np.ceil((t.min() - t0_mid) / P)) + 1
    ev_d = np.array([ev['depth'] for ev in observed if np.isfinite(ev['depth'])])
    ev_e = np.array([ev['err'] for ev in observed if np.isfinite(ev['depth'])])
    # single-event dominance: fraction of the summed (depth/err^2) carried by the strongest event
    w = ev_d / ev_e ** 2
    dom = float(np.max(w) / np.sum(w)) if np.sum(w) > 0 else np.nan
    chi2_const = float(np.sum(((ev_d - np.average(ev_d, weights=1 / ev_e ** 2)) / ev_e) ** 2))
    log(f'{n_transits} transits observed of {n_epochs_in_baseline} epochs in baseline; '
        f'strongest-event share of weight {dom:.2f}; chi2 of constant depth {chi2_const:.1f} / {len(ev_d)-1} dof')
    # trapezoid duration
    trap = fit_trapezoid(t, f, P, t0_mid, dur, d)
    log(f'trapezoid: depth {trap["depth"]*1e6:.0f} ppm, T14 {trap["t14_d"]*24:.2f} h, T23/T14 {trap["t23_over_t14"]:.2f}, dt0 {trap["dt0_d"]*24*60:.1f} min')

    # ---------- 4. SAP flux at the same ephemeris
    ts, fs_raw, es_raw, qs_raw = load(FITS_PATH, 'SAP_FLUX')
    ts, fs, es, qs = detrend(ts, fs_raw, es_raw, qs_raw, window=args.window)
    d_sap, d_sap_formal, n_in_sap = box_depth(ts, fs, P, t0_mid, dur)
    d_sap_err, _, _ = empirical_depth_error(ts, fs, P, t0_mid, dur, n=1500)
    sap_ratio = d_sap / d
    log(f'SAP box depth {d_sap*1e6:.0f} +- {d_sap_err*1e6:.0f} ppm; SAP/PSF depth ratio {sap_ratio:.2f}; '
        f'SAP rms {np.std(fs)*1e6:.0f} ppm')
    # median raw flux levels (e/s) -- crude dilution context
    with fits.open(FITS_PATH) as h:
        hdr0 = h[0].header
        dd = h['LIGHTCURVE_STITCHED'].data
        med_psf = float(np.nanmedian(dd['FLUX'])); med_sap = float(np.nanmedian(dd['SAP_FLUX']))
        tpforg = hdr0.get('TPFORG')
    # also a short BLS on the SAP flux around the same period (does SAP independently see it?)
    _, rsap = run_bls(ts, fs, es, fine_P, fine_dur)
    jsap = int(np.argmax(rsap.power))

    # ---------- 5. inverted-flux null (same grid, same code path)
    null = None
    if not args.no_null:
        _, rnull = run_bls(t, 2.0 - f, e, periods, durations)
        pn = np.asarray(rnull.power)
        n_raw, n_det, _ = sde_from_power(pn, periods)
        kn = int(np.argmax(pn))
        null = dict(sde_raw=n_raw, sde_detrended=n_det, best_period_d=float(periods[kn]),
                    best_depth_ppm=float(rnull.depth[kn] * 1e6), best_duration_d=float(rnull.duration[kn]))
        log(f'inverted null: SDE raw {n_raw:.2f} det {n_det:.2f}, best P {periods[kn]:.4f} d depth {rnull.depth[kn]*1e6:.0f} ppm')

    # ---------- results JSON
    results = dict(
        fname=FNAME, date='2026-10-17', role='independent re-derivation (no pipeline code)',
        data=dict(fits=FITS_PATH, n_quality0=int(len(t_raw)), n_used=int(len(t)), baseline_d=float(baseline),
                  quarters=sorted(int(x) for x in set(q)), rms_ppm=rms * 1e6, robust_rms_ppm=mad_rms * 1e6,
                  detrend='per-quarter median normalisation + time-based running median, window %.2f d, '
                          'upward outliers > 4 MAD-sigma clipped (two-pass)' % args.window),
        bls=dict(pmin_d=pmin, pmax_d=pmax, n_periods=int(len(periods)), durations_d=durations.tolist(),
                 grid='geometric, dP = P * 0.08 d / (baseline * %.2f)' % (0.25 if args.quick else 2.0),
                 objective='likelihood', sde_raw=sde_raw, sde_baseline_subtracted=sde_det,
                 top3_peaks=peaks, top3_peaks_baseline_subtracted=peaks_det,
                 pipeline_period_d=PIPE_P, top_peak_is_pipeline_period=bool(recovered),
                 top_peak_relation_to_pipeline=recovered_alias,
                 fine_grid_alias_tests=alias_tests, adopted_multiple_of_top_peak=adopted_mult,
                 adopted_period_relation_to_pipeline=alias_of(P, PIPE_P)),
        ephemeris=dict(period_d=P, period_err_d=P_err, t0_bkjd=t0_mid, t0_bjd=t0_mid + BKJD0,
                       duration_bls_d=dur, duration_bls_h=dur * 24,
                       duration_trapezoid_t14_h=trap['t14_d'] * 24, trapezoid=trap),
        depth=dict(box_depth_ppm=d * 1e6, depth_err_empirical_ppm=d_err_emp * 1e6,
                   depth_err_formal_ppm=d_err_formal * 1e6, empirical_snr=snr_emp,
                   n_random_phases=n_emp, random_phase_mean_ppm=d_mean_emp * 1e6, n_in_transit=n_in,
                   bls_depth_ppm=float(rfine.depth[j] * 1e6)),
        odd_even=dict(odd_ppm=d_odd * 1e6, odd_err_ppm=e_odd * 1e6, n_odd=int(odd.sum()),
                      even_ppm=d_even * 1e6, even_err_ppm=e_even * 1e6, n_even=int(even.sum()), sigma=float(oe_sigma)),
        secondary=dict(phase=0.5, depth_ppm=d_sec * 1e6, err_ppm=e_sec * 1e6, sigma=float(sec_sigma), n_in=n_sec,
                       half_period_flag=half_period_flag,
                       strongest_in_phase_0p15_0p85=dict(phase=float(sec_scan[k, 0]), depth_ppm=float(sec_scan[k, 1] * 1e6),
                                                         sigma=float(sec_scan[k, 1] / (d_err_emp * np.sqrt(n_in / max(sec_scan[k, 2], 1)))))),
        per_year=per_year,
        transits=dict(n_observed=n_transits, n_epochs_in_baseline=n_epochs_in_baseline,
                      expected_in_transit_points_per_event=float(expected_in),
                      strongest_event_weight_share=dom, chi2_constant_depth=chi2_const, dof=int(len(ev_d) - 1),
                      events=[{k_: (float(v) if isinstance(v, (float, np.floating)) else v) for k_, v in ev.items()} for ev in events]),
        sap=dict(box_depth_ppm=d_sap * 1e6, depth_err_empirical_ppm=d_sap_err * 1e6, n_in=n_in_sap,
                 sap_over_psf_depth_ratio=float(sap_ratio), sap_rms_ppm=float(np.std(fs) * 1e6),
                 median_psf_flux_e_s=med_psf, median_sap_flux_e_s=med_sap, psf_over_sap_flux=med_psf / med_sap,
                 sap_bls_best_period_d=float(fine_P[jsap]), sap_bls_depth_ppm=float(rsap.depth[jsap] * 1e6),
                 sap_bls_duration_h=float(rsap.duration[jsap] * 24), tpforg=tpforg),
        inverted_null=null,
        runtime_s=time.time() - T0, quick=args.quick,
    )
    with open(os.path.join(OUT_DIR, 'rederive_results.json'), 'w') as fh:
        json.dump(results, fh, indent=2, default=float)
    np.savez_compressed(os.path.join(OUT_DIR, 'bls_spectrum.npz'), periods=periods, power=power)
    with open(os.path.join(OUT_DIR, 'events.csv'), 'w') as fh:
        fh.write('epoch,tc_bkjd,quarter,n_in,depth_ppm,err_ppm\n')
        for ev in events:
            fh.write(f"{ev['epoch']},{ev['tc_bkjd']:.5f},{ev['quarter']},{ev['n_in']},"
                     f"{ev['depth']*1e6 if np.isfinite(ev['depth']) else ''},{ev['err']*1e6 if np.isfinite(ev['err']) else ''}\n")

    # ---------- figure
    make_figure(t, f, q, P, t0_mid, dur, d, d_err_emp, periods, power, resid, picks, events, observed,
                ts, fs, d_sap, trap, results)
    log('done')
    return results


def make_figure(t, f, q, P, t0, dur, d, derr, periods, power, resid, picks, events, observed, ts, fs, d_sap, trap, R):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ink = '#1f2328'; muted = '#7a828c'; grid = '#e6e8eb'
    c_pts = '#b8c0c8'; c_bin = '#1d4ed8'; c_odd = '#1d4ed8'; c_even = '#b45309'; c_model = '#111827'; c_sap = '#6b7280'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': muted, 'axes.labelcolor': ink, 'xtick.color': ink,
                         'ytick.color': ink, 'axes.grid': True, 'grid.color': grid, 'grid.linewidth': 0.6,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axs = plt.subplots(3, 2, figsize=(11, 11))
    fig.suptitle(f'{FNAME}: independent re-derivation (PSF flux, own running-median detrend, own BLS grid)\n'
                 f'P = {P:.5f} d, T0 = {t0 + BKJD0:.4f} BJD, box depth = {d*1e6:.0f} +- {derr*1e6:.0f} ppm (empirical), '
                 f'duration (BLS) = {dur*24:.2f} h', fontsize=10, color=ink)

    def binned(x, y, nb, lo, hi):
        edges = np.linspace(lo, hi, nb + 1)
        idx = np.digitize(x, edges) - 1
        bx, by, be = [], [], []
        for i in range(nb):
            m = idx == i
            if m.sum() >= 3:
                bx.append(0.5 * (edges[i] + edges[i + 1])); by.append(np.mean(y[m])); be.append(np.std(y[m]) / np.sqrt(m.sum()))
        return np.array(bx), np.array(by), np.array(be)

    # (0,0) periodogram
    ax = axs[0, 0]
    ax.plot(periods, power, color=c_pts, lw=0.6)
    for r, i in enumerate(picks):
        ax.plot(periods[i], power[i], 'o', ms=5, color=c_bin if r == 0 else c_even, mec='white')
        ax.annotate(f'{periods[i]:.4f} d', (periods[i], power[i]), textcoords='offset points', xytext=(4, 2), fontsize=8, color=ink)
    ax.axvline(PIPE_P, color=c_model, ls=':', lw=0.8)
    ax.set_xscale('log'); ax.set_xlabel('period (d)'); ax.set_ylabel('BLS log-likelihood power')
    ax.set_title(f'blind BLS 1-50 d: SDE raw {R["bls"]["sde_raw"]:.1f}, baseline-subtracted {R["bls"]["sde_baseline_subtracted"]:.1f}', fontsize=9)

    # (0,1) full fold
    ax = axs[0, 1]
    ph = phase_fold(t, P, t0)
    ax.plot(ph, (f - 1) * 1e6, '.', ms=1.5, color=c_pts, rasterized=True)
    bx, by, be = binned(ph, (f - 1) * 1e6, 200, -0.5, 0.5)
    ax.errorbar(bx, by, be, fmt='o', ms=3, color=c_bin, ecolor=c_bin, elinewidth=0.6, capsize=0)
    ax.set_ylim(-6000, 6000); ax.set_xlabel('phase'); ax.set_ylabel('relative flux (ppm)')
    ax.set_title('full phase fold (200 bins); secondary at phase 0.5 = '
                 f'{R["secondary"]["depth_ppm"]:.0f} +- {R["secondary"]["err_ppm"]:.0f} ppm', fontsize=9)
    ax.axvline(0.5 - 1, color=muted, lw=0.5); ax.axvline(0.5, color=muted, lw=0.5)

    # (1,0) zoom with trapezoid
    ax = axs[1, 0]
    phd = ph * P * 24
    win = np.abs(phd) < 12
    ax.plot(phd[win], (f[win] - 1) * 1e6, '.', ms=2, color=c_pts, rasterized=True)
    bx, by, be = binned(phd[win], (f[win] - 1) * 1e6, 48, -12, 12)
    ax.errorbar(bx, by, be, fmt='o', ms=3.5, color=c_bin, ecolor=c_bin, elinewidth=0.6, capsize=0, label='PSF flux, 30-min bins')
    xx = np.linspace(-12, 12, 600)
    ax.plot(xx, (trapezoid(xx / 24 - trap['dt0_d'], trap['depth'], trap['t14_d'], trap['t23_over_t14']) - 1) * 1e6,
            color=c_model, lw=1.2, label=f'trapezoid: {trap["depth"]*1e6:.0f} ppm, T14 {trap["t14_d"]*24:.2f} h')
    ax.axvspan(-dur * 12, dur * 12, color=c_bin, alpha=0.06, lw=0)
    ax.set_ylim(-4000, 3000); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('relative flux (ppm)')
    ax.legend(loc='lower right', fontsize=8, frameon=False); ax.set_title('zoom on transit (shaded: BLS box)', fontsize=9)

    # (1,1) odd / even and SAP
    ax = axs[1, 1]
    epoch = np.round((t - t0) / P).astype(int)
    for lab, m, col, off in (('odd', epoch % 2 != 0, c_odd, 0.0), ('even', epoch % 2 == 0, c_even, 0.0)):
        bx, by, be = binned(phd[win & m], (f[win & m] - 1) * 1e6, 36, -12, 12)
        ax.errorbar(bx + off, by, be, fmt='o', ms=3.5, color=col, ecolor=col, elinewidth=0.6, capsize=0, label=f'{lab} epochs')
    phs = phase_fold(ts, P, t0) * P * 24
    ws = np.abs(phs) < 12
    bx, by, be = binned(phs[ws], (fs[ws] - 1) * 1e6, 36, -12, 12)
    ax.plot(bx, by, '-', color=c_sap, lw=1.0, label=f'SAP flux (binned): depth {d_sap*1e6:.0f} ppm')
    ax.set_ylim(-4000, 3000); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('relative flux (ppm)')
    ax.legend(loc='lower right', fontsize=8, frameon=False)
    ax.set_title(f'odd {R["odd_even"]["odd_ppm"]:.0f} vs even {R["odd_even"]["even_ppm"]:.0f} ppm '
                 f'({R["odd_even"]["sigma"]:.1f} sigma); SAP/PSF depth ratio {R["sap"]["sap_over_psf_depth_ratio"]:.2f}', fontsize=9)

    # (2,0) per-event depths
    ax = axs[2, 0]
    ev_t = np.array([ev['tc_bkjd'] for ev in observed]); ev_d = np.array([ev['depth'] for ev in observed]) * 1e6
    ev_e = np.array([ev['err'] for ev in observed]) * 1e6
    ok = np.isfinite(ev_d)
    ax.errorbar(ev_t[ok], ev_d[ok], ev_e[ok], fmt='o', ms=3.5, color=c_bin, ecolor=c_bin, elinewidth=0.6, capsize=0)
    ax.axhline(d * 1e6, color=c_model, lw=1.0); ax.axhline(0, color=muted, lw=0.6)
    for Y, x0 in ((1, 131), (2, 539), (3, 906), (4, 1274)):
        ax.axvline(x0, color=grid, lw=0.8)
    ax.set_xlabel('mid-transit time (BKJD)'); ax.set_ylabel('local depth (ppm)')
    ax.set_title(f'{len(observed)} observed events of {R["transits"]["n_epochs_in_baseline"]} epochs; '
                 f'max single-event weight share {R["transits"]["strongest_event_weight_share"]:.2f}', fontsize=9)

    # (2,1) per-year depths
    ax = axs[2, 1]
    yrs = [py['year'] for py in R['per_year']]; yd = [py['depth_ppm'] for py in R['per_year']]; ye = [py['err_ppm'] for py in R['per_year']]
    ax.errorbar(yrs, yd, ye, fmt='o', ms=6, color=c_bin, ecolor=c_bin, elinewidth=1, capsize=3)
    ax.axhline(d * 1e6, color=c_model, lw=1.0, label='all-data box depth')
    ax.axhline(0, color=muted, lw=0.6)
    for x_, y_, py in zip(yrs, yd, R['per_year']):
        ax.annotate(f'{y_:.0f} ppm\n{py["n_in"]} pts', (x_, y_), textcoords='offset points', xytext=(8, -4), fontsize=8, color=ink)
    ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(['Y1 (Q1-4)', 'Y2 (Q5-8)', 'Y3 (Q9-12)', 'Y4 (Q13-17)'])
    ax.set_ylabel('box depth (ppm)'); ax.set_title('per-year depths (empirical error scaled by sqrt(N))', fontsize=9)
    ax.legend(loc='lower right', fontsize=8, frameon=False)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(OUT_DIR, 'rederive_fold.png'), dpi=140)
    plt.close(fig)


if __name__ == '__main__':
    main()
