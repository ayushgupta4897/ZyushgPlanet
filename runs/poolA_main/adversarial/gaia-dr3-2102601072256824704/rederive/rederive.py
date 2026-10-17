#!/usr/bin/env python
"""
Independent re-derivation for gaia-dr3-2102601072256824704 (KBONUS-BKG).
Written from scratch with numpy/astropy only. No code from pipeline/ is imported
or read. Detrending is a gap-aware, iteratively sigma-clipped running median
(window WIN_D days) applied per quarter after per-quarter median normalisation.

Outputs (same directory as this script):
  results.json   all numbers quoted in the write-up
  fold.png       phase-folded light curve, odd/even, secondary, per-year
  bls_periodogram.png
  rederive.log   stdout copy
"""
import json, sys, time, os
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares

OUT = os.path.dirname(os.path.abspath(__file__))
FITS = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2102601072256824704_kepler_v1.0_lc.fits'
FNAME = 'gaia-dr3-2102601072256824704'
RNG = np.random.default_rng(20261017)
WIN_D = 1.0          # running-median window (days)
PIPE_P = 32.55710    # pipeline's claimed period (for alias comparison only; search is blind)
PMIN, PMAX = 1.0, 50.0

def log(*a):
    print(*a, flush=True)

# ----------------------------------------------------------------- loading
def load(col):
    with fits.open(FITS) as h:
        d = h['LIGHTCURVE_STITCHED'].data
        t = np.asarray(d['TIME'], float)
        f = np.asarray(d[col], float)
        fe = np.asarray(d[col + '_ERR'], float)
        q = np.asarray(d['SAP_QUALITY'], int)
        quarter = np.asarray(d['QUARTER'], int)
    m = (q == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (f > 0)
    return t[m], f[m], fe[m], quarter[m]

# ----------------------------------------------------------------- detrending
def running_median(t, f, win, nsig=3.0, niter=3, mask=None):
    """Gap-aware sliding-window median. Points flagged by `mask` (True = exclude,
    e.g. in-transit) and iteratively clipped low/high outliers do not enter the
    median but still receive a trend value."""
    n = len(t)
    good = np.ones(n, bool) if mask is None else ~mask
    trend = np.full(n, np.nan)
    for it in range(niter):
        tg, fg = t[good], f[good]
        lo = np.searchsorted(tg, t - win / 2)
        hi = np.searchsorted(tg, t + win / 2)
        for i in range(n):
            if hi[i] - lo[i] >= 5:
                trend[i] = np.median(fg[lo[i]:hi[i]])
        # fill any NaN (tiny isolated chunks) with nearest finite value
        bad = ~np.isfinite(trend)
        if bad.any():
            idx = np.arange(n)
            trend[bad] = np.interp(idx[bad], idx[~bad], trend[~bad])
        resid = f - trend
        sig = 1.4826 * np.median(np.abs(resid[good] - np.median(resid[good])))
        # clip only positive outliers hard (flares/cosmic rays); clip negative
        # outliers only beyond nsig*2 so that shallow transits are not eaten.
        newgood = (resid < nsig * sig) & (resid > -2 * nsig * sig)
        if mask is not None:
            newgood &= ~mask
        if np.array_equal(newgood, good):
            break
        good = newgood
    return trend

def detrend(t, f, fe, quarter, win=WIN_D, mask=None):
    fn = np.empty_like(f); en = np.empty_like(f)
    for q in np.unique(quarter):
        s = quarter == q
        med = np.median(f[s])
        fq = f[s] / med
        eq = fe[s] / med
        # split on gaps > 0.5 d so medians do not straddle gaps
        tq = t[s]
        gaps = np.where(np.diff(tq) > 0.5)[0]
        edges = np.concatenate([[0], gaps + 1, [len(tq)]])
        tr = np.empty_like(fq)
        mq = None if mask is None else mask[s]
        for a, b in zip(edges[:-1], edges[1:]):
            if b - a < 10:
                tr[a:b] = np.median(fq[a:b])
                continue
            tr[a:b] = running_median(tq[a:b], fq[a:b], win,
                                     mask=None if mq is None else mq[a:b])
        fn[s] = fq / tr
        en[s] = eq / tr
    # final global 5-sigma clip of positive outliers only
    sig = 1.4826 * np.median(np.abs(fn - np.median(fn)))
    keep = fn < 1 + 5 * sig
    return fn, en, keep

# ----------------------------------------------------------------- BLS grid
def make_freq_grid(t, pmin, pmax, dur_min, oversample=2.0):
    """Frequency grid where the accumulated phase drift over the baseline is
    < dur_min/oversample: df = dur_min * f / (oversample * T)  -> log-uniform."""
    T = t.max() - t.min()
    fmin, fmax = 1 / pmax, 1 / pmin
    dlnf = dur_min / (oversample * T)
    n = int(np.ceil(np.log(fmax / fmin) / dlnf))
    return np.exp(np.linspace(np.log(fmin), np.log(fmax), n))

def bls_search(t, f, e, periods, durations):
    bls = BoxLeastSquares(t, f, dy=e)
    res = bls.power(periods, durations, objective='likelihood')
    return bls, res

def sde_from_power(periods, power, smooth_frac=0.02):
    """SDE after removing a slowly varying trend in power (running median in
    log-period, window = smooth_frac of the log range)."""
    lp = np.log(periods)
    order = np.argsort(lp)
    lp, pw = lp[order], power[order]
    win = smooth_frac * (lp[-1] - lp[0])
    lo = np.searchsorted(lp, lp - win / 2); hi = np.searchsorted(lp, lp + win / 2)
    # coarse running median evaluated on a subgrid then interpolated (speed)
    sub = np.arange(0, len(lp), 200)
    med = np.array([np.median(pw[lo[i]:hi[i]]) for i in sub])
    trend = np.interp(lp, lp[sub], med)
    r = pw - trend
    sde = (r - np.median(r)) / np.std(r)
    out = np.empty_like(sde); out[order] = sde
    return out

def top_peaks(periods, sde, power, n=3, excl=0.03):
    idx = np.argsort(sde)[::-1]
    picks = []
    for i in idx:
        p = periods[i]
        if all(abs(p / q[0] - 1) > excl for q in picks):
            picks.append((p, sde[i], power[i], i))
        if len(picks) == n:
            break
    return picks

def alias_of(p, pref, tol=0.01):
    for num, den in [(1,1),(2,1),(1,2),(3,1),(1,3),(3,2),(2,3),(4,1),(1,4)]:
        if abs(p / (pref * num / den) - 1) < tol:
            return f'{num}/{den}'
    return None

# ----------------------------------------------------------------- depth tools
def phase(t, P, t0):
    return ((t - t0) / P + 0.5) % 1.0 - 0.5   # in [-0.5,0.5), transit at 0

def box_depth(t, f, P, t0, dur, oot_excl=2.0):
    ph = phase(t, P, t0) * P   # days from mid-transit
    inn = np.abs(ph) < dur / 2
    oot = (np.abs(ph) > oot_excl * dur / 2)
    if inn.sum() < 3:
        return np.nan, np.nan, int(inn.sum())
    d = np.mean(f[oot]) - np.mean(f[inn])
    err = np.std(f[oot]) / np.sqrt(inn.sum())  # formal white-noise error
    return d, err, int(inn.sum())

def empirical_depth_err(t, f, P, t0, dur, nrand=2000):
    """Place boxes of the same duration at random phases (avoiding +-2 dur of the
    transit and of phase 0.5), measure the 'depth' each time; the std of that
    distribution is the empirical 1-sigma error on a box depth at this P/dur,
    including red noise."""
    ph = phase(t, P, t0) * P
    oot_all = np.abs(ph) > 2 * dur
    base = np.mean(f[oot_all])
    vals = []
    tries = 0
    while len(vals) < nrand and tries < nrand * 5:
        tries += 1
        c = RNG.uniform(-0.5, 0.5) * P
        if abs(c) < 1.5 * dur or abs(abs(c) - P / 2) < 1.5 * dur:
            continue
        sel = np.abs(ph - c) < dur / 2
        if sel.sum() < 3:
            continue
        vals.append(base - np.mean(f[sel]))
    vals = np.array(vals)
    return np.std(vals), np.mean(vals), len(vals)

def segment_edges(t, gap=0.5):
    g = np.where(np.diff(t) > gap)[0]
    return np.concatenate([[t[0]], t[g + 1]]), np.concatenate([t[g], [t[-1]]])

def event_table(t, f, P, t0, dur):
    """Per-transit local depth (in - local out within +-1 d), count of in-transit points,
    nearest intra-window gap and distance to the nearest data-segment edge (gaps>0.5 d)."""
    n_lo = int(np.floor((t.min() - t0) / P)); n_hi = int(np.ceil((t.max() - t0) / P))
    seg_s, seg_e = segment_edges(t)
    rows = []
    for n in range(n_lo, n_hi + 1):
        tc = t0 + n * P
        loc = np.abs(t - tc) < 1.0
        inn = np.abs(t - tc) < dur / 2
        oot = loc & (np.abs(t - tc) > dur)
        if inn.sum() == 0:
            continue
        cover = inn.sum() / (dur / (29.4 / 1440))
        d = np.mean(f[oot]) - np.mean(f[inn]) if oot.sum() > 5 else np.nan
        de = np.std(f[oot]) / np.sqrt(inn.sum()) if oot.sum() > 5 else np.nan
        # gap proximity: nearest cadence gap > 0.1 d within 1 d of tc
        tt = t[loc]
        dg = np.diff(tt)
        gap_near = np.nan
        if len(tt) > 1 and (dg > 0.1).any():
            gi = np.where(dg > 0.1)[0]
            gap_near = float(np.min(np.minimum(np.abs(tt[gi] - tc), np.abs(tt[gi + 1] - tc))))
        edge = float(min(np.min(np.abs(seg_s - tc)), np.min(np.abs(seg_e - tc))))
        rows.append(dict(epoch=n, tc=float(tc), n_in=int(inn.sum()), coverage=float(cover), segment_edge_dist_d=edge,
                         depth_ppm=float(d * 1e6) if np.isfinite(d) else None,
                         depth_err_ppm=float(de * 1e6) if np.isfinite(de) else None,
                         nearest_gap_d=None if not np.isfinite(gap_near) else gap_near))
    return rows

# ===================================================================== MAIN
def main():
    T0 = time.time()
    out = dict(fname=FNAME, date='2026-10-17', method=dict(
        quality='SAP_QUALITY==0, finite TIME/FLUX/FLUX_ERR, FLUX>0',
        normalise='per-quarter median',
        detrend=f'gap-aware sliding running median, window {WIN_D} d, 3 iterations sigma-clip (+3/-6 sigma), segments split at gaps>0.5 d',
        bls='astropy BoxLeastSquares, likelihood objective, own log-uniform frequency grid with phase drift < dur_min/2 over baseline',
        depth_err='std of box means at random phases of the same duration (red-noise aware), out-of-transit excludes +-2 durations'))

    # ---------- PSF flux
    t, f, fe, quarter = load('FLUX')
    log(f'loaded {len(t)} good cadences, {t.min():.2f}-{t.max():.2f} BKJD, quarters {np.unique(quarter).tolist()}')
    fn, en, keep = detrend(t, f, fe, quarter)
    t_, f_, e_, q_ = t[keep], fn[keep], en[keep], quarter[keep]
    log(f'after detrend/clip: {len(t_)} cadences; point-to-point rms {1e6*np.std(np.diff(f_))/np.sqrt(2):.0f} ppm')
    out['n_cadences_used'] = int(len(t_))
    out['rms_ppm_per_cadence'] = float(1e6 * np.std(np.diff(f_)) / np.sqrt(2))

    # ---------- blind BLS 1-50 d
    durations = np.array([0.06, 0.09, 0.13, 0.18, 0.25, 0.35])
    freqs = make_freq_grid(t_, PMIN, PMAX, durations.min(), oversample=2.0)
    periods = np.sort(1 / freqs)
    log(f'BLS grid: {len(periods)} periods, durations {durations.tolist()}')
    t1 = time.time()
    bls, res = bls_search(t_, f_, e_, periods, durations)
    log(f'BLS done in {time.time()-t1:.0f}s')
    sde = sde_from_power(periods, res.power)
    peaks = top_peaks(periods, sde, res.power, n=3)
    out['bls'] = dict(pmin=PMIN, pmax=PMAX, n_periods=int(len(periods)), durations_d=durations.tolist(),
                      top_peaks=[dict(period_d=float(p), sde=float(s), power=float(pw),
                                      duration_d=float(res.duration[i]), t0=float(res.transit_time[i]),
                                      depth_ppm=float(res.depth[i] * 1e6),
                                      alias_of_pipeline_P=alias_of(p, PIPE_P)) for p, s, pw, i in peaks])
    for pk in out['bls']['top_peaks']:
        log('  peak', pk)
    best = peaks[0]
    out['recovered_blind'] = alias_of(best[0], PIPE_P) == '1/1'
    out['top_peak_is_alias_of_pipeline'] = alias_of(best[0], PIPE_P)

    # also: SDE at the pipeline period in this blind search (nearest grid point)
    j = np.argmin(np.abs(periods - PIPE_P))
    out['bls']['sde_at_pipeline_period'] = float(sde[j])
    out['bls']['rank_of_pipeline_period'] = int(np.sum(sde > sde[j]) + 1)

    # ---------- inverted-flux null with identical grid
    t1 = time.time()
    _, resi = bls_search(t_, 2 - f_, e_, periods, durations)
    sdei = sde_from_power(periods, resi.power)
    pki = top_peaks(periods, sdei, resi.power, n=3)
    out['inverted_null'] = dict(top_peaks=[dict(period_d=float(p), sde=float(s)) for p, s, pw, i in pki],
                                max_sde=float(np.max(sdei)))
    log(f'inverted null done in {time.time()-t1:.0f}s; max SDE {np.max(sdei):.2f} at P={pki[0][0]:.4f}')

    # ---------- refine around top peak
    P0 = best[0]
    fine_P = np.linspace(P0 * (1 - 0.002), P0 * (1 + 0.002), 4001)
    fine_dur = np.linspace(0.04, 0.4, 37)
    rf = bls.power(fine_P, fine_dur, objective='likelihood', oversample=20)
    k = np.argmax(rf.power)
    P, t0, dur = float(rf.period[k]), float(rf.transit_time[k]), float(rf.duration[k])
    log(f'refined: P={P:.5f} t0={t0:.4f} dur={dur*24:.2f} h bls_depth={rf.depth[k]*1e6:.0f}+-{rf.depth_err[k]*1e6:.0f} ppm')
    # duration scan at fixed P/t0 for a duration estimate and its ~1-sigma range
    out['refined'] = dict(period_d=P, t0_bkjd=t0, t0_bjd=t0 + 2454833.0, duration_h=dur * 24,
                          bls_depth_ppm=float(rf.depth[k] * 1e6), bls_depth_err_ppm=float(rf.depth_err[k] * 1e6),
                          bls_snr=float(rf.depth_snr[k]), bls_loglike=float(rf.log_likelihood[k]))

    # ---------- masked re-detrend (so the 1-d median does not eat the transit)
    ph_d = phase(t, P, t0) * P
    tmask = np.abs(ph_d) < 1.5 * dur
    fn2, en2, keep2 = detrend(t, f, fe, quarter, mask=tmask)
    keep2 &= np.isfinite(fn2)
    t2, f2, e2, q2 = t[keep2], fn2[keep2], en2[keep2], quarter[keep2]

    def measure(tt, ff, label):
        d, derr_formal, nin = box_depth(tt, ff, P, t0, dur)
        emp, bias, nr = empirical_depth_err(tt, ff, P, t0, dur)
        # odd/even
        ep = np.round((tt - t0) / P).astype(int)
        odd = ep % 2 == 1
        d_odd, e_odd, n_odd = box_depth(tt[odd], ff[odd], P, t0, dur)
        d_even, e_even, n_even = box_depth(tt[~odd], ff[~odd], P, t0, dur)
        # empirical per-half errors: scale emp by sqrt(ratio of points)
        e_odd_emp = emp * np.sqrt(nin / max(n_odd, 1)); e_even_emp = emp * np.sqrt(nin / max(n_even, 1))
        oe_sig = abs(d_odd - d_even) / np.hypot(e_odd_emp, e_even_emp)
        # secondary at phase 0.5
        d_sec, e_sec_f, n_sec = box_depth(tt, ff, P, t0 + P / 2, dur)
        sec_sig = d_sec / emp
        # per-year depths (Kepler years by BKJD)
        yrs = [(131, 500), (500, 865), (865, 1230), (1230, 1600)]
        py = []
        for a, b in yrs:
            s = (tt >= a) & (tt < b)
            dy, ey, ny = box_depth(tt[s], ff[s], P, t0, dur)
            py.append(dict(bkjd_range=[a, b], depth_ppm=None if not np.isfinite(dy) else float(dy * 1e6),
                           err_ppm_emp=None if ny == 0 else float(emp * np.sqrt(nin / ny) * 1e6), n_in=ny))
        ev = event_table(tt, ff, P, t0, dur)
        n_tr = sum(1 for r in ev if r['coverage'] > 0.5)
        r = dict(label=label, depth_ppm=float(d * 1e6), depth_err_formal_ppm=float(derr_formal * 1e6),
                 depth_err_empirical_ppm=float(emp * 1e6), random_box_mean_ppm=float(bias * 1e6),
                 empirical_snr=float(d / emp), n_in_transit=nin,
                 odd=dict(depth_ppm=float(d_odd * 1e6), err_ppm=float(e_odd_emp * 1e6), n_in=n_odd),
                 even=dict(depth_ppm=float(d_even * 1e6), err_ppm=float(e_even_emp * 1e6), n_in=n_even),
                 odd_even_sigma=float(oe_sig),
                 secondary=dict(depth_ppm=float(d_sec * 1e6), err_ppm=float(emp * 1e6), sigma=float(sec_sig), n_in=n_sec),
                 per_year=py, n_transits_covered=int(n_tr), n_epochs_with_data=len(ev), events=ev)
        log(f'[{label}] depth {d*1e6:.0f} +- {emp*1e6:.0f} ppm (formal {derr_formal*1e6:.0f}), SNR {d/emp:.1f}, '
            f'odd {d_odd*1e6:.0f} even {d_even*1e6:.0f} ({oe_sig:.1f} sig), sec {d_sec*1e6:.0f} ppm ({sec_sig:.1f} sig), '
            f'n_tr {n_tr}, per-year {[round(x["depth_ppm"]) if x["depth_ppm"] is not None else None for x in py]}')
        return r

    out['psf_unmasked_detrend'] = measure(t_, f_, 'PSF, unmasked 1-d median detrend')
    out['psf_masked_detrend'] = measure(t2, f2, 'PSF, transit-masked 1-d median detrend')

    # duration scan on the masked light curve: chi2 of box vs duration at fixed P, t0
    durs = np.linspace(0.03, 0.45, 43)
    bls2 = BoxLeastSquares(t2, f2, dy=e2)
    ll = []
    for dd in durs:
        st = bls2.compute_stats(P, dd, t0)
        ll.append(st['depth'][0] / st['depth'][1])
    ll = np.array(ll)
    kbest = int(np.argmax(ll))
    out['duration_scan'] = dict(durations_h=(durs * 24).tolist(), depth_snr=ll.tolist(),
                                best_duration_h=float(durs[kbest] * 24))
    # also a proper fine re-fit of P,t0,dur on masked data
    rf2 = bls2.power(fine_P, fine_dur, objective='likelihood', oversample=20)
    k2 = np.argmax(rf2.power)
    out['refined_masked'] = dict(period_d=float(rf2.period[k2]), t0_bkjd=float(rf2.transit_time[k2]),
                                 duration_h=float(rf2.duration[k2] * 24), bls_depth_ppm=float(rf2.depth[k2] * 1e6),
                                 bls_depth_err_ppm=float(rf2.depth_err[k2] * 1e6))
    log('refined on masked detrend:', out['refined_masked'])

    # ---------- SAP flux
    ts, fs, fes, qs = load('SAP_FLUX')
    fns, ens, keeps = detrend(ts, fs, fes, qs)
    keeps &= np.isfinite(fns)
    out['sap_unmasked_detrend'] = measure(ts[keeps], fns[keeps], 'SAP, unmasked 1-d median detrend')
    dp = out['psf_unmasked_detrend']['depth_ppm']; ds = out['sap_unmasked_detrend']['depth_ppm']
    ep_ = out['psf_unmasked_detrend']['depth_err_empirical_ppm']; es_ = out['sap_unmasked_detrend']['depth_err_empirical_ppm']
    ratio = ds / dp
    ratio_err = abs(ratio) * np.hypot(ep_ / dp, es_ / ds)
    out['sap_vs_psf'] = dict(psf_depth_ppm=dp, sap_depth_ppm=ds, ratio_sap_over_psf=float(ratio), ratio_err=float(ratio_err),
                             note='both on unmasked 1-d running-median detrend, full data, same P/t0/dur')
    log(f'SAP/PSF depth ratio {ratio:.2f} +- {ratio_err:.2f}')

    # ---------- edge-dominated events (objective rule, fixed before looking at which epochs qualify):
    # |local depth| > 5 x median |local depth| of all covered events AND within 0.15 d of a segment edge
    ev = out['psf_unmasked_detrend']['events']
    dd = np.array([e['depth_ppm'] for e in ev if e['depth_ppm'] is not None])
    thr = 5 * np.median(np.abs(dd))
    dom = [e for e in ev if e['depth_ppm'] is not None and abs(e['depth_ppm']) > thr]
    edge_dom = [e for e in dom if e['segment_edge_dist_d'] < 0.15]
    log(f'dominant events (|depth|>{thr:.0f} ppm): {[(e["epoch"], round(e["depth_ppm"]), round(e["segment_edge_dist_d"],3)) for e in dom]}')
    nall = out['psf_unmasked_detrend']['n_in_transit']
    contrib = sum(e['depth_ppm'] * e['n_in'] for e in edge_dom) / nall if edge_dom else 0.0
    out['edge_dominated_events'] = dict(rule='|local depth| > 5*median|local depth| and segment_edge_dist < 0.15 d',
                                        threshold_ppm=float(thr), events=edge_dom,
                                        contribution_to_mean_depth_ppm=float(contrib),
                                        fraction_of_mean_depth=float(contrib / dp))
    # excise them (+-1 d around each) and re-measure / re-search
    exc_ep = [e['epoch'] for e in edge_dom]
    ep_all = np.round((t_ - t0) / P).astype(int)
    excl = np.isin(ep_all, exc_ep) & (np.abs(phase(t_, P, t0) * P) < 1.0)
    tx, fx, ex = t_[~excl], f_[~excl], e_[~excl]
    out['psf_excised'] = measure(tx, fx, f'PSF, unmasked detrend, epochs {exc_ep} excised')
    for k in ('events',):
        out['psf_excised'].pop(k)
    t1 = time.time()
    rx = BoxLeastSquares(tx, fx, dy=ex).power(periods, durations, objective='likelihood')
    sdex = sde_from_power(periods, rx.power); pkx = top_peaks(periods, sdex, rx.power, n=3)
    jx = np.argmin(np.abs(periods - P))
    out['psf_excised']['bls'] = dict(top_peaks=[dict(period_d=float(p), sde=float(s), alias_of_pipeline_P=alias_of(p, PIPE_P)) for p, s, pw, i in pkx],
                                     sde_at_pipeline_period=float(sdex[jx]), rank_of_pipeline_period=int(np.sum(sdex > sdex[jx]) + 1))
    log(f'excised BLS ({time.time()-t1:.0f}s): top {[(round(p,4), round(s,2)) for p,s,pw,i in pkx]}; SDE at P_pipe {sdex[jx]:.2f} rank {out["psf_excised"]["bls"]["rank_of_pipeline_period"]}')

    # ---------- principled alternative: trim TRIM d from every segment edge BEFORE detrending
    TRIM = 0.25
    def trim_edges(tt, trim=TRIM):
        s_, e_s = segment_edges(tt)
        ds_ = np.min(np.abs(tt[:, None] - s_[None, :]), axis=1); de_ = np.min(np.abs(tt[:, None] - e_s[None, :]), axis=1)
        return (ds_ > trim) & (de_ > trim)
    okp = trim_edges(t)
    fnT, enT, keepT = detrend(t[okp], f[okp], fe[okp], quarter[okp])
    tT, fT, eT = t[okp][keepT], fnT[keepT], enT[keepT]
    out['psf_edge_trimmed'] = measure(tT, fT, f'PSF, {TRIM} d trimmed at every segment edge, then 1-d median detrend')
    out['psf_edge_trimmed']['n_cadences'] = int(len(tT)); out['psf_edge_trimmed']['trim_d'] = TRIM
    t1 = time.time()
    bT = BoxLeastSquares(tT, fT, dy=eT)
    rT = bT.power(periods, durations, objective='likelihood'); sdeT = sde_from_power(periods, rT.power); pkT = top_peaks(periods, sdeT, rT.power, n=3)
    rTi = BoxLeastSquares(tT, 2 - fT, dy=eT).power(periods, durations, objective='likelihood'); sdeTi = sde_from_power(periods, rTi.power); pkTi = top_peaks(periods, sdeTi, rTi.power, n=3)
    jT = np.argmin(np.abs(periods - P))
    out['psf_edge_trimmed']['bls'] = dict(top_peaks=[dict(period_d=float(p), sde=float(s), alias_of_pipeline_P=alias_of(p, PIPE_P)) for p, s, pw, i in pkT],
                                          sde_at_pipeline_period=float(sdeT[jT]), rank_of_pipeline_period=int(np.sum(sdeT > sdeT[jT]) + 1),
                                          inverted_top_peaks=[dict(period_d=float(p), sde=float(s)) for p, s, pw, i in pkTi], inverted_max_sde=float(np.max(sdeTi)))
    log(f'edge-trimmed BLS+null ({time.time()-t1:.0f}s): top {[(round(p,4), round(s,2)) for p,s,pw,i in pkT]}; SDE at P_pipe {sdeT[jT]:.2f} rank {out["psf_edge_trimmed"]["bls"]["rank_of_pipeline_period"]}; inverted max {np.max(sdeTi):.2f}')
    oks = trim_edges(ts)
    fnsT, ensT, keepsT = detrend(ts[oks], fs[oks], fes[oks], qs[oks])
    out['sap_edge_trimmed'] = measure(ts[oks][keepsT], fnsT[keepsT], f'SAP, {TRIM} d edge-trimmed, 1-d median detrend')
    out['sap_edge_trimmed'].pop('events')

    # ---------- comparison with the pipeline (read only now, after all numbers above were produced)
    try:
        vet = json.load(open('/home/user/ZyushgPlanet/runs/poolA_main/vet/gaia-dr3-2102601072256824704/vet.json'))
        pipe_ev = np.genfromtxt('/home/user/ZyushgPlanet/runs/poolA_main/vet/gaia-dr3-2102601072256824704/events.csv', delimiter=',', names=True)
        big = pipe_ev[np.abs(pipe_ev['depth']) > 5e-4]
        out['pipeline_vet'] = dict(P=vet['P'], dur_h=vet['dur_h'], t0_bkjd=vet['t0_bkjd'], depth_pass2_ppm=vet['depth_pass2_ppm'],
                                   snr_emp_pass2=vet['snr_emp_pass2'], odd_even_sigma=vet['odd_even_sigma'], secondary_depth_ppm=vet['secondary_depth_ppm'],
                                   n_events=vet['n_events'], event_chi2_const=vet['event_chi2_const'], event_dof=vet['event_dof'],
                                   event_weighted_depth_ppm=vet['event_weighted_depth_ppm'], tls_sde=vet['tls']['SDE'], tls_period=vet['tls']['period'],
                                   stage2_pass=vet['stage2_pass'], stage2_rules=vet['stage2_rules'], stage2_pass_v14=vet.get('stage2_pass_v14'),
                                   events_csv_deep_events=[dict(n=int(r['n']), tc=float(r['tc_bkjd']), depth_ppm=float(r['depth'] * 1e6), gap_dist_d=float(r['gap_dist_d'])) for r in big])
        out['pipeline_vet']['pixel_json_exists'] = os.path.exists('/home/user/ZyushgPlanet/runs/poolA_main/pixel/gaia-dr3-2102601072256824704/pixel.json')
    except Exception as exc:
        out['pipeline_vet'] = dict(error=str(exc))

    # ---------- figure
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 2, figsize=(13, 16))
    def binned(x, y, w):
        edges = np.arange(x.min(), x.max() + w, w)
        idx = np.digitize(x, edges)
        sel = [i for i in range(1, len(edges)) if (idx == i).sum() > 2]
        xb = np.array([x[idx == i].mean() for i in sel]); yb = np.array([y[idx == i].mean() for i in sel])
        eb = np.array([y[idx == i].std() / np.sqrt((idx == i).sum()) for i in sel])
        return xb, yb, eb
    r = out['psf_unmasked_detrend']; rT_ = out['psf_edge_trimmed']
    ph = phase(t_, P, t0) * P * 24
    a = ax[0, 0]; s = np.abs(ph) < 24
    a.plot(ph[s], (f_[s] - 1) * 1e6, '.', color='0.75', ms=2)
    xb, yb, eb = binned(ph[s], (f_[s] - 1) * 1e6, 0.5); a.errorbar(xb, yb, eb, fmt='o', color='C3', ms=4)
    a.axvspan(-dur * 12, dur * 12, color='C0', alpha=0.15); a.set_ylim(-1500, 1500); a.set_xlabel('hours from mid-transit'); a.set_ylabel('ppm')
    a.set_title(f'{FNAME} PSF flux, all data (own detrend)\nP={P:.5f} d t0={t0:.4f} BKJD dur={dur*24:.2f} h depth={dp:.0f}+-{ep_:.0f} ppm (30-min bins)', fontsize=10)
    a = ax[0, 1]; phT = phase(tT, P, t0) * P * 24; s = np.abs(phT) < 24
    a.plot(phT[s], (fT[s] - 1) * 1e6, '.', color='0.75', ms=2)
    xb, yb, eb = binned(phT[s], (fT[s] - 1) * 1e6, 0.5); a.errorbar(xb, yb, eb, fmt='o', color='C3', ms=4)
    a.axvspan(-dur * 12, dur * 12, color='C0', alpha=0.15); a.set_ylim(-1500, 1500); a.set_xlabel('hours'); a.set_ylabel('ppm')
    a.set_title(f'same ephemeris, {TRIM} d trimmed at every segment edge before detrending\ndepth={rT_["depth_ppm"]:.0f}+-{rT_["depth_err_empirical_ppm"]:.0f} ppm (SNR {rT_["empirical_snr"]:.1f})', fontsize=10)
    a = ax[1, 0]
    phs = phase(ts[keeps], P, t0) * P * 24; s = np.abs(phs) < 24
    xb, yb, eb = binned(phs[s], (fns[keeps][s] - 1) * 1e6, 1.0); a.errorbar(xb, yb, eb, fmt='o', color='C2', ms=4, label=f'SAP {ds:.0f}+-{es_:.0f} ppm')
    xb, yb, eb = binned(ph[np.abs(ph) < 24], (f_[np.abs(ph) < 24] - 1) * 1e6, 1.0); a.errorbar(xb + 0.15, yb, eb, fmt='s', color='C3', ms=3, alpha=0.6, label=f'PSF {dp:.0f}+-{ep_:.0f} ppm')
    a.axvspan(-dur * 12, dur * 12, color='C0', alpha=0.15); a.legend(); a.set_ylim(-800, 800); a.set_title('SAP vs PSF, all data (1-h bins)'); a.set_xlabel('hours'); a.set_ylabel('ppm')
    a = ax[1, 1]
    epk = np.round((t_ - t0) / P).astype(int)
    for par, col, lab in [(1, 'C1', 'odd'), (0, 'C4', 'even')]:
        sel = (epk % 2 == par) & (np.abs(ph) < 24)
        xb, yb, eb = binned(ph[sel], (f_[sel] - 1) * 1e6, 1.0); a.errorbar(xb + (0.15 if par else -0.15), yb, eb, fmt='o', ms=4, color=col, label=lab)
    a.set_title(f"odd {r['odd']['depth_ppm']:.0f} / even {r['even']['depth_ppm']:.0f} ppm ({r['odd_even_sigma']:.1f} sigma), all data; edge-trimmed: {rT_['odd']['depth_ppm']:.0f}/{rT_['even']['depth_ppm']:.0f} ({rT_['odd_even_sigma']:.1f} sigma)", fontsize=9)
    a.axvspan(-dur * 12, dur * 12, color='C0', alpha=0.15); a.legend(); a.set_ylim(-800, 800); a.set_xlabel('hours'); a.set_ylabel('ppm')
    a = ax[2, 0]
    ph5 = phase(t_, P, t0 + P / 2) * P * 24; s = np.abs(ph5) < 24
    xb, yb, eb = binned(ph5[s], (f_[s] - 1) * 1e6, 1.0); a.errorbar(xb, yb, eb, fmt='o', color='C5', ms=4)
    a.axvspan(-dur * 12, dur * 12, color='C0', alpha=0.15); a.set_ylim(-800, 800)
    a.set_title(f"secondary (phase 0.5): {r['secondary']['depth_ppm']:.0f} ppm ({r['secondary']['sigma']:.1f} sigma)"); a.set_xlabel('hours'); a.set_ylabel('ppm')
    a = ax[2, 1]
    evs = [e for e in r['events'] if e['depth_ppm'] is not None]
    cols = ['r' if e['segment_edge_dist_d'] < 0.15 else 'C0' for e in evs]
    a.errorbar([e['tc'] for e in evs], [e['depth_ppm'] for e in evs], [e['depth_err_ppm'] for e in evs], fmt='none', ecolor='0.6')
    a.scatter([e['tc'] for e in evs], [e['depth_ppm'] for e in evs], c=cols, s=25, zorder=3)
    for e in edge_dom:
        a.annotate(f"n={e['epoch']}\nedge {e['segment_edge_dist_d']:.2f} d", (e['tc'], e['depth_ppm']), fontsize=8, textcoords='offset points', xytext=(5, -25))
    a.axhline(dp, color='C3', label=f'mean depth {dp:.0f}'); a.axhline(out['psf_excised']['depth_ppm'], color='C2', ls='--', label=f"excised {out['psf_excised']['depth_ppm']:.0f}"); a.axhline(0, color='k', lw=0.5)
    a.legend(fontsize=8); a.set_title('per-event local depth (red = within 0.15 d of a segment edge)'); a.set_xlabel('BKJD'); a.set_ylabel('ppm')
    a = ax[3, 0]
    a.plot(periods, sde, lw=0.5, color='C0', label='real, all data')
    a.plot(periods, sdei, lw=0.5, color='0.6', alpha=0.7, label='inverted (2-f), all data')
    a.axvline(PIPE_P, color='g', ls='--', alpha=0.6, label='pipeline P')
    a.set_xscale('log'); a.set_xlabel('period (d)'); a.set_ylabel('SDE'); a.legend(fontsize=8)
    a.set_title(f'blind BLS 1-50 d: top real SDE {sde.max():.1f} @ {periods[np.argmax(sde)]:.3f} d; inverted max {np.max(sdei):.1f}', fontsize=10)
    a = ax[3, 1]
    a.plot(periods, sdeT, lw=0.5, color='C0', label='real, edge-trimmed'); a.plot(periods, sdeTi, lw=0.5, color='0.6', alpha=0.7, label='inverted, edge-trimmed')
    a.axvline(PIPE_P, color='g', ls='--', alpha=0.6, label='pipeline P')
    a.set_xscale('log'); a.set_xlabel('period (d)'); a.set_ylabel('SDE'); a.legend(fontsize=8)
    a.set_title(f'edge-trimmed BLS: SDE at pipeline P = {sdeT[jT]:.1f} (rank {out["psf_edge_trimmed"]["bls"]["rank_of_pipeline_period"]}); top {pkT[0][0]:.3f} d SDE {pkT[0][1]:.1f}', fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fold.png'), dpi=110)

    out['runtime_s'] = time.time() - T0
    with open(os.path.join(OUT, 'results.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    log('wrote results.json; runtime', round(out['runtime_s']), 's')

if __name__ == '__main__':
    main()
