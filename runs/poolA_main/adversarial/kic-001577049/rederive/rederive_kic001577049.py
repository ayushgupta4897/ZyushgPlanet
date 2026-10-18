#!/usr/bin/env python
"""
Independent re-derivation of the kic-001577049 first-screen signal (ZyushgPlanet poolA_main).

Written from scratch with numpy / astropy / matplotlib only. No code under pipeline/ was imported,
read or consulted before these numbers were produced. Detrending is a time-windowed running median
(not wotan). BLS is astropy.timeseries.BoxLeastSquares on a period grid built here (Ofir 2014 spacing).

Outputs (all under runs/poolA_main/adversarial/kic-001577049/rederive/):
  rederive_results.json   every number quoted in the write-up
  fold_rederive.png       fold figure
  rederive.log            console log

Run:  . .venv/bin/activate && python -I rederive_kic001577049.py
"""
import json, os, sys, time
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares, LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FITS = '/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-001577049_kepler_v1.0_lc.fits'
OUT = '/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-001577049/rederive'
PIPE_P, PIPE_DEPTH_PPM, PIPE_DUR_H = 0.53790, 1324.0, 1.1      # pipeline values quoted in the task text (for comparison only)
RNG = np.random.default_rng(20261017)
WINDOW_D = 1.0          # running-median window (days)
N_RANDOM = 2000         # random-phase boxes for the empirical depth error
T_START = time.time()
LOG = []

def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True); LOG.append(s)

def J(x):
    """json-safe"""
    if isinstance(x, dict): return {str(k): J(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [J(v) for v in x]
    if isinstance(x, np.ndarray): return [J(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    return x

# ----------------------------------------------------------------------------- loading / detrending
def load(col):
    with fits.open(FITS) as f:
        h0 = f[0].header
        d = f['LIGHTCURVE_STITCHED'].data
        t = np.asarray(d['TIME'], float); fl = np.asarray(d[col], float); fe = np.asarray(d[col + '_ERR'], float)
        q = np.asarray(d['SAP_QUALITY'], int); qu = np.asarray(d['QUARTER'], int)
    m = (q == 0) & np.isfinite(t) & np.isfinite(fl) & np.isfinite(fe) & (fl > 0)
    o = np.argsort(t[m])
    meta = {k: h0.get(k) for k in ('OBJECT', 'KEPLERID', 'GAIAID', 'RA_OBJ', 'DEC_OBJ', 'KEPMAG', 'GMAG', 'TPFORG')}
    return t[m][o], fl[m][o], fe[m][o], qu[m][o], meta, int(len(t))

def fold(t, P, T0):
    """phase in days, in [-P/2, P/2)"""
    return (t - T0 + 0.5 * P) % P - 0.5 * P

def running_median(t, y, window=1.0, mask=None, min_pts=10):
    """Time-windowed running median per contiguous segment (gaps > 0.5 d split segments).
    Points with mask==False are excluded from the median but still get a trend value."""
    n = len(t); tr = np.full(n, np.nan)
    mask = np.ones(n, bool) if mask is None else mask
    gaps = np.where(np.diff(t) > 0.5)[0]
    starts = np.r_[0, gaps + 1]; ends = np.r_[gaps + 1, n]; h = window / 2
    for s, e in zip(starts, ends):
        ts = t[s:e]; ys = y[s:e]; ms = mask[s:e]; tg = ts[ms]; yg = ys[ms]
        lo = np.searchsorted(tg, ts - h, 'left'); hi = np.searchsorted(tg, ts + h, 'right')
        for i in range(e - s):
            if hi[i] - lo[i] >= min_pts:
                tr[s + i] = np.median(yg[lo[i]:hi[i]])
    return tr

def prepare(col, mask_ephem=None, window=WINDOW_D):
    """quality==0 -> per-quarter median normalisation -> running median -> clip (+4 / -10 robust sigma)."""
    t, fl, fe, qu, meta, n_raw = load(col)
    fn = np.empty_like(fl); medq = {}
    for Q in np.unique(qu):
        m = qu == Q; mq = np.median(fl[m]); medq[int(Q)] = float(mq); fn[m] = fl[m] / mq
    mask = None
    if mask_ephem is not None:
        P, T0, dur = mask_ephem
        mask = np.abs(fold(t, P, T0)) > 1.0 * dur
    tr = running_median(t, fn, window, mask=mask)
    y = fn / tr
    ok = np.isfinite(y)
    r = y - 1
    mad = 1.4826 * np.median(np.abs(r[ok] - np.median(r[ok])))
    ok &= (r < 4 * mad) & (r > -10 * mad)
    sig = np.empty_like(y)
    for Q in np.unique(qu):
        m = (qu == Q) & ok
        sig[qu == Q] = 1.4826 * np.median(np.abs(y[m] - np.median(y[m])))
    return dict(col=col, t=t[ok], y=y[ok], sig=sig[ok], qu=qu[ok], meta=meta, n_raw=n_raw, n_quality0=int(len(t)),
                n_used=int(ok.sum()), n_clipped=int((~ok).sum()), per_quarter_median_eps=medq,
                median_flux_eps=float(np.median(fl)), global_robust_sigma_ppm=float(1e6 * mad),
                per_quarter_sigma_ppm={int(Q): float(1e6 * np.median(sig[(qu == Q) & ok])) for Q in np.unique(qu)})

# ----------------------------------------------------------------------------- BLS helpers
def ofir_grid(tspan_d, Pmin, Pmax, M_sun=0.75, R_sun=0.77, OS=3):
    """Ofir (2014) optimal frequency sampling, uniform in f^(1/3); density set by the transit duration a
    star of (M,R) implies at each period, with oversampling OS over the baseline tspan_d."""
    G, Msun, Rsun = 6.674e-11, 1.989e30, 6.957e8
    T = tspan_d * 86400.
    A = (2 * np.pi) ** (2 / 3) / np.pi * R_sun * Rsun / (G * M_sun * Msun) ** (1 / 3) / (T * OS)
    fmin, fmax = 1 / (Pmax * 86400), 1 / (Pmin * 86400)
    N = int(np.ceil((fmax ** (1 / 3) - fmin ** (1 / 3) + A / 3) * 3 / A))
    f = (fmin ** (1 / 3) + A * np.arange(N + 1) / 3) ** 3
    f = f[(f <= fmax) & (f >= fmin)]
    return np.sort(1 / f / 86400.), float(A), N

def top_peaks(res, k=3, excl=0.02):
    pw = np.asarray(res.power); order = np.argsort(pw)[::-1]; picked = []
    for i in order:
        if all(abs(res.period[i] - res.period[j]) / res.period[j] > excl for j in picked):
            picked.append(i)
        if len(picked) == k: break
    sde = (pw - pw.mean()) / pw.std()
    return [dict(period_d=float(res.period[i]), power=float(pw[i]), sde=float(sde[i]), depth_ppm=float(1e6 * res.depth[i]),
                 depth_snr_formal=float(res.depth_snr[i]), duration_h=float(24 * res.duration[i]),
                 t0_bkjd=float(res.transit_time[i])) for i in picked]

def alias_of(P, Pref, tol=0.003):
    r = P / Pref
    for n in range(1, 500):
        if abs(r - n) < tol * n: return f'{n} x Pref'
        if abs(r - 1 / n) < tol / n: return f'Pref / {n}'
    for n in range(1, 12):
        for m in range(2, 12):
            if abs(r - n / m) < tol: return f'{n}/{m} x Pref'
    return None

def run_bls(d, periods, durations, oversample=5, label=''):
    t1 = time.time()
    bls = BoxLeastSquares(d['t'], d['y'], dy=d['sig'])
    res = bls.power(periods, durations, objective='likelihood', oversample=oversample)
    pw = np.asarray(res.power)
    log(f'  BLS[{label}] {len(periods)} periods x {len(durations)} durations in {time.time()-t1:.0f} s; '
        f'top P={res.period[np.argmax(pw)]:.5f} SDE={(pw.max()-pw.mean())/pw.std():.1f}')
    return res

# ----------------------------------------------------------------------------- box measurements
def box_depth(t, y, P, T0, dur, center_phase=0.0):
    ph = fold(t, P, T0 + center_phase * P)
    inm = np.abs(ph) < dur / 2; out = np.abs(ph) > dur
    if inm.sum() < 3 or out.sum() < 10: return np.nan, int(inm.sum())
    return float(np.median(y[out]) - np.mean(y[inm])), int(inm.sum())

def empirical_err(t, y, P, T0, dur, n=N_RANDOM, avoid=((0.0, 1.5), (0.5, 1.5))):
    """std of box 'depths' at random phases (same P, dur), avoiding windows around the primary and phase 0.5."""
    deps = []; tries = 0
    while len(deps) < n and tries < 20 * n:
        tries += 1
        c = RNG.uniform(0, 1)
        if any(abs(((c - pc + 0.5) % 1) - 0.5) * P < k * dur for pc, k in avoid): continue
        dd, nin = box_depth(t, y, P, T0, dur, center_phase=c)
        if np.isfinite(dd): deps.append(dd)
    deps = np.array(deps)
    return float(np.std(deps)), float(np.mean(deps)), int(len(deps))

def measure(d, P, T0, dur, P_fund=None, tag=''):
    t, y, qu = d['t'], d['y'], d['qu']
    out = dict(period_d=P, t0_bkjd=T0, duration_h=24 * dur)
    dep, nin = box_depth(t, y, P, T0, dur)
    err, mu, nr = empirical_err(t, y, P, T0, dur)
    out.update(depth_ppm=1e6 * dep, depth_err_ppm_empirical=1e6 * err, n_in_transit_points=nin,
               random_phase_mean_ppm=1e6 * mu, n_random_boxes=nr, empirical_snr=dep / err if err > 0 else np.nan)
    # secondary at phase 0.5
    sdep, snin = box_depth(t, y, P, T0, dur, 0.5)
    out.update(secondary_ppm=1e6 * sdep, secondary_sigma=sdep / err if err > 0 else np.nan, secondary_n_points=snin)
    # odd / even
    ep = np.round((t - T0) / P).astype(int)
    oe = {}
    for name, sel in (('odd', ep % 2 == 1), ('even', ep % 2 == 0)):
        dd, nn = box_depth(t[sel], y[sel], P, T0, dur)
        ee, _, _ = empirical_err(t[sel], y[sel], P, T0, dur, n=N_RANDOM // 2)
        oe[name] = dict(depth_ppm=1e6 * dd, err_ppm=1e6 * ee, n_points=nn)
    out['odd_even'] = oe
    out['odd_even_sigma'] = abs(oe['odd']['depth_ppm'] - oe['even']['depth_ppm']) / np.hypot(oe['odd']['err_ppm'], oe['even']['err_ppm'])
    # per year (Kepler year by quarter: Q1-4 -> 1, Q5-8 -> 2, Q9-12 -> 3, Q13-17 -> 4) and per quarter
    year = (qu - 1) // 4 + 1
    py = {}
    for Y in np.unique(year):
        sel = year == Y
        dd, nn = box_depth(t[sel], y[sel], P, T0, dur); ee, _, _ = empirical_err(t[sel], y[sel], P, T0, dur, n=N_RANDOM // 4)
        py[int(Y)] = dict(depth_ppm=1e6 * dd, err_ppm=1e6 * ee, n_points=nn, quarters=[int(q) for q in np.unique(qu[sel])])
    out['per_year'] = py
    pq = {}
    for Q in np.unique(qu):
        sel = qu == Q; dd, nn = box_depth(t[sel], y[sel], P, T0, dur)
        pq[int(Q)] = dict(depth_ppm=1e6 * dd, n_points=nn)
    out['per_quarter'] = pq
    # number of transits (epochs with >=1 / >=2 in-box cadences) and expected count over the baseline
    inm = np.abs(fold(t, P, T0)) < dur / 2
    ep_in = ep[inm]; u, c = np.unique(ep_in, return_counts=True)
    out.update(n_transits_ge1=int(len(u)), n_transits_ge2=int((c >= 2).sum()),
               n_epochs_spanned=int(np.ptp(ep) + 1), median_points_per_transit=float(np.median(c)) if len(c) else 0)
    # depths at the sub-harmonic phases (k/m) if a fundamental is supplied: a genuine single transit shows
    # nothing there; a signal whose true period is P/m shows the same depth at every k/m.
    if P_fund:
        m_ = int(round(P / P_fund)); sub = {}
        for k in range(m_):
            dd, nn = box_depth(t, y, P, T0, dur, k / m_)
            sub[f'{k}/{m_}'] = dict(depth_ppm=1e6 * dd, sigma=dd / err if err > 0 else np.nan, n_points=nn)
        out['depth_at_subharmonic_phases'] = sub; out['harmonic_ratio_m'] = m_
    # out-of-box harmonic content at P (ellipsoidal / reflection proxy)
    ph = fold(t, P, T0) / P; o = np.abs(ph) > dur / P
    X = np.c_[np.ones(o.sum()), np.cos(2 * np.pi * ph[o]), np.sin(2 * np.pi * ph[o]), np.cos(4 * np.pi * ph[o]), np.sin(4 * np.pi * ph[o])]
    cf = np.linalg.lstsq(X, y[o], rcond=None)[0]
    out['out_of_box_harmonics_ppm'] = dict(A1=1e6 * np.hypot(cf[1], cf[2]), A2=1e6 * np.hypot(cf[3], cf[4]))
    log(f'  [{tag}] P={P:.6f} T0={T0:.4f} dur={24*dur:.2f} h: depth={1e6*dep:.0f}+/-{1e6*err:.0f} ppm (emp SNR {dep/err:.1f}), '
        f'nin={nin}, secondary={1e6*sdep:.0f} ppm ({sdep/err:.1f} sig), odd/even {oe["odd"]["depth_ppm"]:.0f}/{oe["even"]["depth_ppm"]:.0f} '
        f'({out["odd_even_sigma"]:.1f} sig), n_transits>=1 {len(u)}')
    return out

def refine(d, P0, dur_grid, half_width_frac=0.004, label=''):
    """fine BLS around P0; step chosen so the accumulated phase drift over the baseline is < dur/20."""
    T = np.ptp(d['t']); dmin = dur_grid.min()
    step = P0 * dmin / (20 * T)
    per = np.arange(P0 * (1 - half_width_frac), P0 * (1 + half_width_frac), step)
    res = run_bls(d, per, dur_grid, oversample=10, label=f'refine {label}')
    i = int(np.argmax(res.power))
    stats = BoxLeastSquares(d['t'], d['y'], dy=d['sig']).compute_stats(res.period[i], res.duration[i], res.transit_time[i])
    return dict(period_d=float(res.period[i]), duration_d=float(res.duration[i]), t0_bkjd=float(res.transit_time[i]),
                depth_ppm=float(1e6 * res.depth[i]), depth_err_ppm_formal=float(1e6 * res.depth_err[i]),
                depth_snr_formal=float(res.depth_snr[i]), power=float(res.power[i]), n_grid=int(len(per)), step_d=float(step),
                astropy_stats=dict(depth_odd_ppm=float(1e6 * stats['depth_odd'][0]), depth_even_ppm=float(1e6 * stats['depth_even'][0]),
                                   depth_half_ppm=float(1e6 * stats['depth_half'][0]), depth_phased_ppm=float(1e6 * stats['depth_phased'][0]),
                                   harmonic_delta_log_likelihood=float(stats['harmonic_delta_log_likelihood']))), res

def binned(ph, y, nb):
    edges = np.linspace(-0.5, 0.5, nb + 1); idx = np.clip(np.digitize(ph, edges) - 1, 0, nb - 1)
    c = 0.5 * (edges[1:] + edges[:-1]); mu = np.full(nb, np.nan); se = np.full(nb, np.nan); n = np.zeros(nb, int)
    for b in range(nb):
        m = idx == b; n[b] = m.sum()
        if n[b] > 2: mu[b] = np.mean(y[m]); se[b] = np.std(y[m]) / np.sqrt(n[b])
    return c, mu, se, n

def shape_fit(d, P, T0, nb=50):
    """On the fold at P: chi2 of (a) box, (b) trapezoid, (c) sinusoid (fundamental + 1st harmonic). V-shape proxy = T23/T14."""
    ph = fold(d['t'], P, T0) / P
    c, mu, se, n = binned(ph, d['y'], nb)
    ok = np.isfinite(mu); c, mu, se = c[ok], mu[ok], se[ok]
    # sinusoid
    X = np.c_[np.ones_like(c), np.cos(2 * np.pi * c), np.sin(2 * np.pi * c), np.cos(4 * np.pi * c), np.sin(4 * np.pi * c)]
    cf = np.linalg.lstsq(X / se[:, None], mu / se, rcond=None)[0]
    chi_sin = float(np.sum(((mu - X @ cf) / se) ** 2))
    # trapezoid grid (depth analytic for given T14, T23, centre=0)
    best = (np.inf, None)
    for T14 in np.linspace(0.05, 0.98, 60):          # total duration in phase units
        for r in np.linspace(0.0, 1.0, 21):          # T23/T14  (0 = V, 1 = box)
            T23 = r * T14
            x = np.abs(c)
            shape = np.where(x <= T23 / 2, 1.0, np.where(x >= T14 / 2, 0.0, (T14 / 2 - x) / ((T14 - T23) / 2 + 1e-12)))
            A = np.c_[np.ones_like(c), -shape]
            cc = np.linalg.lstsq(A / se[:, None], mu / se, rcond=None)[0]
            chi = float(np.sum(((mu - A @ cc) / se) ** 2))
            if chi < best[0]: best = (chi, dict(T14_phase=float(T14), T23_over_T14=float(r), depth_ppm=float(1e6 * cc[1]), chi2=chi))
    trap = best[1]
    # pure box
    bb = (np.inf, None)
    for T14 in np.linspace(0.05, 0.98, 120):
        shape = (np.abs(c) <= T14 / 2).astype(float); A = np.c_[np.ones_like(c), -shape]
        cc = np.linalg.lstsq(A / se[:, None], mu / se, rcond=None)[0]; chi = float(np.sum(((mu - A @ cc) / se) ** 2))
        if chi < bb[0]: bb = (chi, dict(T14_phase=float(T14), depth_ppm=float(1e6 * cc[1]), chi2=chi))
    nbins = int(ok.sum())
    return dict(n_bins=nbins, sinusoid=dict(chi2=chi_sin, k=5, A1_ppm=float(1e6 * np.hypot(cf[1], cf[2])), A2_ppm=float(1e6 * np.hypot(cf[3], cf[4])),
                                                 bic=chi_sin + 5 * np.log(nbins)),
                trapezoid=dict(**trap, k=4, bic=trap['chi2'] + 4 * np.log(nbins)), box=dict(**bb[1], k=3, bic=bb[1]['chi2'] + 3 * np.log(nbins)),
                binned=dict(phase=c, mean_ppm=1e6 * (mu - 1), se_ppm=1e6 * se))

# ============================================================================= main
def main():
    os.makedirs(OUT, exist_ok=True)
    R = dict(fname='kic-001577049', date='2026-10-17', fits=FITS, method=dict(
        quality='SAP_QUALITY == 0, finite, flux > 0', normalisation='per-quarter median', detrend=f'time-windowed running median, {WINDOW_D} d window, '
        'segments split at gaps > 0.5 d, >= 10 points per window', clipping='+4 / -10 robust (MAD) sigma after detrending',
        bls='astropy BoxLeastSquares, objective=likelihood, per-quarter robust sigma as dy, Ofir-2014 period grid (M=0.75 Msun, R=0.77 Rsun, OS=3)',
        empirical_error=f'std of box depths at {N_RANDOM} random phases (same P, duration), avoiding +/-1.5 dur around phase 0 and 0.5',
        sde='(peak - mean) / std of the BLS likelihood power over the grid'))

    # ---------------- 1. load + detrend
    log('1. loading / detrending PSF flux')
    d = prepare('FLUX')
    R['meta'] = d['meta']
    R['data'] = {k: d[k] for k in ('n_raw', 'n_quality0', 'n_used', 'n_clipped', 'per_quarter_median_eps', 'median_flux_eps',
                                   'global_robust_sigma_ppm', 'per_quarter_sigma_ppm')}
    R['data'].update(quarters=[int(q) for q in np.unique(d['qu'])], t_min=float(d['t'].min()), t_max=float(d['t'].max()),
                     baseline_d=float(np.ptp(d['t'])), days_with_data=float(len(d['t']) * 0.0204))
    log(f"  N raw {d['n_raw']}, quality0 {d['n_quality0']}, used {d['n_used']}, robust sigma {d['global_robust_sigma_ppm']:.0f} ppm, baseline {np.ptp(d['t']):.1f} d")
    tspan = float(np.ptp(d['t']))

    # ---------------- 2. blind BLS 1-50 d (+ inverted null)
    log('2. blind BLS 1-50 d')
    per, A, N = ofir_grid(tspan, 1.0, 50.0, OS=3)
    durs = np.array([0.025, 0.035, 0.05, 0.07, 0.1, 0.14, 0.2, 0.28])
    res = run_bls(d, per, durs, oversample=5, label='1-50 d real')
    peaks = top_peaks(res, 3)
    for p in peaks: p['alias_of_pipeline_P'] = alias_of(p['period_d'], PIPE_P)
    dinv = dict(d, y=2 - d['y'])
    resi = run_bls(dinv, per, durs, oversample=5, label='1-50 d inverted')
    peaks_inv = top_peaks(resi, 3)
    pw = np.asarray(res.power); pwi = np.asarray(resi.power)
    R['blind_bls_1_50d'] = dict(n_periods=int(len(per)), ofir_A=A, durations_d=durs, oversample=5, top_peaks=peaks,
                                 inverted_null_top_peaks=peaks_inv, sde_real=peaks[0]['sde'], sde_inverted=peaks_inv[0]['sde'],
                                 pipeline_period_in_range=bool(1.0 <= PIPE_P <= 50.0),
                                 note='the pipeline period 0.5379 d lies below the 1 d floor of this search; only its integer multiples can appear here')
    # power at the integer multiples of the pipeline period, for the alias comb
    comb = []
    for n in range(2, 20):
        Pn = PIPE_P * n
        if 1.0 <= Pn <= 50:
            w = np.abs(per - Pn) / Pn < 0.003
            if w.any():
                i = np.argmax(np.where(w, pw, -np.inf))
                comb.append(dict(n=n, period_d=float(per[i]), sde=float((pw[i] - pw.mean()) / pw.std()), depth_ppm=float(1e6 * res.depth[i]),
                                 snr_formal=float(res.depth_snr[i]), rank=int((pw > pw[i]).sum()) + 1))
    R['blind_bls_1_50d']['pipeline_period_multiples'] = comb
    log('  top-3 real:', json.dumps(J(peaks)))
    log('  top-3 inverted:', json.dumps(J(peaks_inv)))

    # ---------------- 2b. supplementary short-period BLS 0.08-1 d (to locate the fundamental)
    log('2b. supplementary BLS 0.08-1.0 d')
    per_s, A_s, _ = ofir_grid(tspan, 0.08, 1.0, OS=3)
    durs_s = np.array([0.01, 0.015, 0.02, 0.03, 0.045, 0.065])
    res_s = run_bls(d, per_s, durs_s, oversample=5, label='0.08-1 d real')
    peaks_s = top_peaks(res_s, 8, excl=0.02)
    pws = np.asarray(res_s.power)
    # fundamental = shortest period among the strong comb (peaks within 3x of the pipeline period ratio pattern)
    P_fund0 = min(p['period_d'] for p in peaks_s[:5])
    # check the comb: are the strong peaks integer multiples of P_fund0?
    for p in peaks_s:
        p['alias_of_pipeline_P'] = alias_of(p['period_d'], PIPE_P); p['ratio_to_shortest_strong_peak'] = p['period_d'] / P_fund0
    resi_s = run_bls(dinv, per_s, durs_s, oversample=5, label='0.08-1 d inverted')
    R['supplementary_bls_0p08_1d'] = dict(n_periods=int(len(per_s)), durations_d=durs_s, top_peaks=peaks_s, sde_real=peaks_s[0]['sde'],
                                           inverted_null_top_peaks=top_peaks(resi_s, 3), label='post-hoc supplementary search, outside the 1-50 d blind window')
    log('  top peaks 0.08-1 d:', json.dumps(J(peaks_s)))

    # ---------------- 2c. Lomb-Scargle (is the signal sinusoidal? where is the fundamental?)
    log('2c. Lomb-Scargle 0.2-48 c/d')
    fgrid = np.arange(0.2, 48.0, 1 / (tspan * 5))
    ls = LombScargle(d['t'], d['y'], d['sig'])
    lsp = ls.power(fgrid, method='fast')
    order = np.argsort(lsp)[::-1]; lspk = []
    for i in order:
        if all(abs(fgrid[i] - q['frequency_cpd']) > 0.05 for q in lspk):
            amp = ls.model_parameters(fgrid[i]); lspk.append(dict(frequency_cpd=float(fgrid[i]), period_d=float(1 / fgrid[i]), power=float(lsp[i]),
                                                                  semi_amplitude_ppm=float(1e6 * np.hypot(amp[1], amp[2])),
                                                                  fap=float(ls.false_alarm_probability(lsp[i], method='baluev', minimum_frequency=0.2, maximum_frequency=48.0))))
        if len(lspk) == 6: break
    f1 = lspk[0]['frequency_cpd']
    R['lomb_scargle'] = dict(n_frequencies=int(len(fgrid)), top_peaks=lspk, kepler_lc_sampling_cpd=1 / 0.0204340,
                             super_nyquist_alias_of_top_peak_cpd=float(1 / 0.0204340 - f1),
                             ratio_pipeline_period_to_top_LS_period=float(PIPE_P * f1))
    log('  LS top peaks:', json.dumps(J(lspk)))

    # ---------------- 3. refine ephemerides: (A) the pipeline-comparable 0.538 d alias, (B) the fundamental
    log('3. refining ephemerides')
    # (A) around the pipeline period
    refA, resA = refine(d, PIPE_P, np.linspace(0.015, 0.09, 16), label='0.538 d alias')
    log('  refined alias:', json.dumps(J({k: v for k, v in refA.items() if k != 'astropy_stats'})))
    # masked re-detrend at the refined alias ephemeris, then measure
    dA = prepare('FLUX', mask_ephem=(refA['period_d'], refA['t0_bkjd'], refA['duration_d']))
    # fundamental from LS (sanity: ratio to the alias should be an integer)
    P_fund = 1 / f1
    m_ratio = refA['period_d'] / P_fund
    log(f'  LS fundamental P={P_fund:.6f} d ({24*P_fund:.3f} h); alias/fundamental = {m_ratio:.4f}')
    m_int = int(round(m_ratio))
    P_fund_from_alias = refA['period_d'] / m_int if abs(m_ratio - m_int) < 0.02 else P_fund
    refB, resB = refine(d, P_fund_from_alias, np.linspace(0.008, 0.06, 14), half_width_frac=0.002, label='fundamental')
    log('  refined fundamental:', json.dumps(J({k: v for k, v in refB.items() if k != 'astropy_stats'})))
    R['refined'] = dict(alias_0p538=refA, fundamental=refB, alias_over_fundamental=float(refA['period_d'] / refB['period_d']),
                        fundamental_from_LS_d=P_fund, m=m_int)

    # ---------------- 3b. measurements at both ephemerides (PSF flux)
    log('3b. measurements on PSF flux')
    MA = measure(dA, refA['period_d'], refA['t0_bkjd'], refA['duration_d'], P_fund=refB['period_d'], tag='PSF @ 0.538 d alias, masked re-detrend')
    MA_nomask = measure(d, refA['period_d'], refA['t0_bkjd'], refA['duration_d'], tag='PSF @ 0.538 d alias, unmasked detrend')
    MB = measure(d, refB['period_d'], refB['t0_bkjd'], refB['duration_d'], tag='PSF @ fundamental')
    # shape on the folds
    SA = shape_fit(dA, refA['period_d'], refA['t0_bkjd'], nb=50)
    SB = shape_fit(d, refB['period_d'], refB['t0_bkjd'], nb=40)
    R['psf'] = dict(alias_0p538=MA, alias_0p538_unmasked_detrend=MA_nomask, fundamental=MB,
                    shape_alias_fold={k: v for k, v in SA.items() if k != 'binned'}, shape_fundamental_fold={k: v for k, v in SB.items() if k != 'binned'})
    log('  shape (fundamental fold): sinusoid chi2 %.1f (BIC %.1f, A1 %.0f ppm), trapezoid chi2 %.1f (BIC %.1f, T23/T14 %.2f), box chi2 %.1f (BIC %.1f)' % (
        SB['sinusoid']['chi2'], SB['sinusoid']['bic'], SB['sinusoid']['A1_ppm'], SB['trapezoid']['chi2'], SB['trapezoid']['bic'], SB['trapezoid']['T23_over_T14'],
        SB['box']['chi2'], SB['box']['bic']))

    # ---------------- 4. SAP flux
    log('4. SAP flux')
    s = prepare('SAP_FLUX', mask_ephem=(refA['period_d'], refA['t0_bkjd'], refA['duration_d']))
    MSA = measure(s, refA['period_d'], refA['t0_bkjd'], refA['duration_d'], P_fund=refB['period_d'], tag='SAP @ 0.538 d alias')
    MSB = measure(s, refB['period_d'], refB['t0_bkjd'], refB['duration_d'], tag='SAP @ fundamental')
    SSB = shape_fit(s, refB['period_d'], refB['t0_bkjd'], nb=40)
    Fp, Fs = d['median_flux_eps'], s['median_flux_eps']
    ratio = MSA['depth_ppm'] / MA['depth_ppm']
    R['sap'] = dict(data={k: s[k] for k in ('n_used', 'n_clipped', 'median_flux_eps', 'global_robust_sigma_ppm', 'per_quarter_sigma_ppm')},
                    alias_0p538=MSA, fundamental=MSB, shape_fundamental_fold={k: v for k, v in SSB.items() if k != 'binned'},
                    dilution=dict(median_psf_flux_eps=Fp, median_sap_flux_eps=Fs, flux_ratio_psf_over_sap=Fp / Fs,
                                  depth_ratio_sap_over_psf_alias=ratio, depth_ratio_err=ratio * np.hypot(MSA['depth_err_ppm_empirical'] / MSA['depth_ppm'], MA['depth_err_ppm_empirical'] / MA['depth_ppm']),
                                  depth_ratio_sap_over_psf_fundamental=MSB['depth_ppm'] / MB['depth_ppm'],
                                  absolute_deficit_psf_eps=1e-6 * MA['depth_ppm'] * Fp, absolute_deficit_sap_eps=1e-6 * MSA['depth_ppm'] * Fs,
                                  expected_depth_ratio_if_source_fully_in_both=Fp / Fs,
                                  interpretation=('If the eclipsing/variable source were fully captured in both products, the SAP depth should be diluted to '
                                                  'depth_psf x (F_psf/F_sap). A SAP/PSF depth ratio well above F_psf/F_sap means the SAP aperture contains '
                                                  'more of the varying source (in e/s) than the PSF flux attributes to the target, i.e. the variability is '
                                                  'dominated by light that the PSF fit only partly assigns to this star.')))
    log(f'  dilution: F_psf/F_sap = {Fp/Fs:.3f}; depth SAP/PSF (alias) = {ratio:.2f}; abs deficit PSF {1e-6*MA["depth_ppm"]*Fp:.1f} e/s vs SAP {1e-6*MSA["depth_ppm"]*Fs:.1f} e/s')

    # zero-point check (separate script in this directory)
    zp = os.path.join(OUT, 'zero_point_check.json')
    if os.path.exists(zp):
        R['zero_point_check'] = json.load(open(zp))['summary']

    # ---------------- 5. physical plausibility of the fundamental
    G, Msun, Rsun = 6.674e-11, 1.989e30, 6.957e8
    for key, P in (('fundamental', refB['period_d']), ('alias_0p538', refA['period_d'])):
        a = (G * 0.75 * Msun * (P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
        R.setdefault('orbit_if_planet_around_0p77Rsun_0p75Msun', {})[key] = dict(period_d=P, a_over_Rstar=a / (0.77 * Rsun),
                                                                                  central_duration_h=24 * P / np.pi / (a / (0.77 * Rsun)))

    # ---------------- 6. comparison with the pipeline numbers quoted in the task
    comp = dict(pipeline=dict(period_d=PIPE_P, depth_ppm=PIPE_DEPTH_PPM, duration_h=PIPE_DUR_H),
                mine_at_alias=dict(period_d=refA['period_d'], depth_ppm_box=MA['depth_ppm'], depth_ppm_bls=refA['depth_ppm'], duration_h=24 * refA['duration_d']),
                period_agrees=abs(refA['period_d'] - PIPE_P) < 1e-4, depth_ratio_mine_over_pipeline=MA['depth_ppm'] / PIPE_DEPTH_PPM,
                duration_agrees=abs(24 * refA['duration_d'] - PIPE_DUR_H) < 0.25)
    R['comparison_with_pipeline'] = comp
    R['runtime_s'] = time.time() - T_START

    # ---------------- figure
    make_figure(d, dA, s, refA, refB, MA, MSA, MB, res, resi, res_s, fgrid, lsp, SA, SB, SSB, R)
    R['log'] = LOG
    json.dump(J(R), open(os.path.join(OUT, 'rederive_results.json'), 'w'), indent=1)
    open(os.path.join(OUT, 'rederive.log'), 'w').write('\n'.join(LOG) + '\n')
    log(f'done in {time.time()-T_START:.0f} s')

# ----------------------------------------------------------------------------- figure
C_PSF, C_SAP, C_INV, C_TXT, C_TXT2, C_SURF, C_GRID = '#2a78d6', '#eb6834', '#52514e', '#0b0b0b', '#52514e', '#fcfcfb', '#e6e5e1'

def style(ax, title, xl, yl):
    ax.set_facecolor(C_SURF); ax.set_title(title, color=C_TXT, fontsize=9.5, loc='left')
    ax.set_xlabel(xl, color=C_TXT2, fontsize=8.5); ax.set_ylabel(yl, color=C_TXT2, fontsize=8.5)
    ax.tick_params(colors=C_TXT2, labelsize=8); ax.grid(True, color=C_GRID, lw=0.6)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'): ax.spines[sp].set_color(C_GRID)

def make_figure(d, dA, s, refA, refB, MA, MSA, MB, res, resi, res_s, fgrid, lsp, SA, SB, SSB, R):
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.6), facecolor=C_SURF)
    fig.suptitle('kic-001577049 independent re-derivation (PSF flux unless stated; 1 d running-median detrend; quality==0)', color=C_TXT, fontsize=11, x=0.01, ha='left')
    PA, T0A, durA = refA['period_d'], refA['t0_bkjd'], refA['duration_d']; PB, T0B, durB = refB['period_d'], refB['t0_bkjd'], refB['duration_d']
    # (a) fold at the 0.538 d alias
    ax = axs[0, 0]; ph = fold(dA['t'], PA, T0A) / PA; c, mu, se, n = binned(ph, dA['y'], 100)
    ax.plot(ph, 1e6 * (dA['y'] - 1), '.', ms=2, color=C_GRID, alpha=0.6, label='cadences')
    ax.errorbar(c, 1e6 * (mu - 1), 1e6 * se, fmt='o', ms=4, color=C_PSF, lw=1, label='100 bins')
    m_ = R['refined']['m']
    for k in range(m_): ax.axvline(((k / m_ + 0.5) % 1) - 0.5, color=C_SAP, lw=1, ls=':', alpha=0.9)
    ax.axvspan(-durA / PA / 2, durA / PA / 2, color=C_PSF, alpha=0.08, lw=0)
    ax.set_ylim(-2500, 2500); ax.legend(fontsize=7.5, loc='upper right', frameon=False)
    style(ax, f'(a) fold at P = {PA:.5f} d (pipeline alias); dotted = k/{m_} phases\nbox depth {MA["depth_ppm"]:.0f} +/- {MA["depth_err_ppm_empirical"]:.0f} ppm (random-phase err)', 'phase', 'relative flux [ppm]')
    # (b) fold at the fundamental, PSF and SAP
    ax = axs[0, 1]
    for dd, col, lab in ((d, C_PSF, 'PSF FLUX'), (s, C_SAP, 'SAP_FLUX')):
        ph = fold(dd['t'], PB, T0B) / PB; c, mu, se, n = binned(ph, dd['y'], 40)
        ax.errorbar(c, 1e6 * (mu - 1), 1e6 * se, fmt='o', ms=4, color=col, lw=1, label=lab)
    xx = np.linspace(-0.5, 0.5, 300)
    ax.plot(xx, SB['sinusoid']['A1_ppm'] * 0 + 1e6 * 0 + np.interp(xx, SB['binned']['phase'], SB['binned']['mean_ppm']), color=C_PSF, lw=0.8, alpha=0.3)
    ax.axvspan(-durB / PB / 2, durB / PB / 2, color=C_PSF, alpha=0.08, lw=0)
    ax.legend(fontsize=7.5, loc='upper right', frameon=False)
    style(ax, f'(b) fold at the fundamental P = {PB:.6f} d = {24*PB:.3f} h\nsinusoid A1 = {SB["sinusoid"]["A1_ppm"]:.0f} ppm (PSF), {SSB["sinusoid"]["A1_ppm"]:.0f} ppm (SAP); T23/T14 = {SB["trapezoid"]["T23_over_T14"]:.2f}', 'phase', 'relative flux [ppm]')
    # (c) Lomb-Scargle
    ax = axs[0, 2]; ax.plot(fgrid, lsp, color=C_PSF, lw=0.7)
    f1 = R['lomb_scargle']['top_peaks'][0]['frequency_cpd']
    for k in range(1, 6): ax.axvline(k * f1, color=C_SAP, lw=0.8, ls=':', alpha=0.8)
    ax.axvline(1 / PIPE_P, color=C_INV, lw=0.8, ls='--'); ax.text(1 / PIPE_P, ax.get_ylim()[1] * 0.9, ' 1/0.5379 d', fontsize=7.5, color=C_TXT2)
    ax.set_xlim(0.2, 48)
    style(ax, f'(c) Lomb-Scargle of detrended PSF flux; top peak f = {f1:.4f} c/d (P = {24/f1:.3f} h), dotted = harmonics', 'frequency [c/d]', 'LS power')
    # (d) blind BLS 1-50 d, real vs inverted
    ax = axs[1, 0]; pw = np.asarray(res.power); pwi = np.asarray(resi.power)
    ax.plot(res.period, (pw - pw.mean()) / pw.std(), color=C_PSF, lw=0.6, label='real flux')
    ax.plot(resi.period, (pwi - pwi.mean()) / pwi.std(), color=C_INV, lw=0.6, alpha=0.8, label='inverted flux (2 - f)')
    for n in range(2, 12): ax.axvline(PIPE_P * n, color=C_SAP, lw=0.7, ls=':', alpha=0.7)
    ax.set_xscale('log'); ax.legend(fontsize=7.5, loc='upper right', frameon=False)
    style(ax, f'(d) blind BLS 1-50 d (Ofir grid, {len(res.period)} periods); dotted = n x 0.5379 d\ntop real P = {R["blind_bls_1_50d"]["top_peaks"][0]["period_d"]:.4f} d SDE {R["blind_bls_1_50d"]["sde_real"]:.1f}; inverted SDE {R["blind_bls_1_50d"]["sde_inverted"]:.1f}', 'period [d]', 'SDE')
    # (e) short-period BLS 0.08-1 d
    ax = axs[1, 1]; pws = np.asarray(res_s.power)
    ax.plot(res_s.period, (pws - pws.mean()) / pws.std(), color=C_PSF, lw=0.6)
    for n in range(1, 10): ax.axvline(PB * n, color=C_SAP, lw=0.7, ls=':', alpha=0.7)
    ax.axvline(PIPE_P, color=C_INV, lw=0.8, ls='--')
    style(ax, f'(e) supplementary BLS 0.08-1 d; dotted = n x {PB:.5f} d, dashed = 0.5379 d', 'period [d]', 'SDE')
    # (f) per-year and odd/even depths at the alias, PSF vs SAP
    ax = axs[1, 2]
    labels = []; xs = 0
    for key, M, col, off in (('PSF', MA, C_PSF, -0.12), ('SAP', MSA, C_SAP, 0.12)):
        ys = []; es = []; xx = []
        for i, (Y, v) in enumerate(sorted(M['per_year'].items())):
            xx.append(i + off); ys.append(v['depth_ppm']); es.append(v['err_ppm'])
            if key == 'PSF': labels.append(f'Y{Y}')
        base = len(M['per_year'])
        for j, nm in enumerate(('odd', 'even')):
            xx.append(base + j + off); ys.append(M['odd_even'][nm]['depth_ppm']); es.append(M['odd_even'][nm]['err_ppm'])
            if key == 'PSF': labels.append(nm)
        xx.append(base + 2 + off); ys.append(M['secondary_ppm']); es.append(M['depth_err_ppm_empirical'])
        if key == 'PSF': labels.append('phase 0.5')
        ax.errorbar(xx, ys, es, fmt='o', ms=5, color=col, lw=1.2, label=f'{key} (all: {M["depth_ppm"]:.0f} +/- {M["depth_err_ppm_empirical"]:.0f} ppm)')
    ax.axhline(0, color=C_GRID, lw=1); ax.axhline(PIPE_DEPTH_PPM, color=C_INV, lw=0.8, ls='--'); ax.text(-0.3, PIPE_DEPTH_PPM * 1.03, 'pipeline 1324 ppm', fontsize=7.5, color=C_TXT2)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels); ax.legend(fontsize=7.5, loc='lower left', frameon=False)
    style(ax, f'(f) box depth at the 0.538 d alias ephemeris, by Kepler year / odd-even / phase 0.5\nodd-even {MA["odd_even_sigma"]:.1f} sigma (PSF); secondary {MA["secondary_sigma"]:.1f} sigma', '', 'depth [ppm]')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(OUT, 'fold_rederive.png'), dpi=130, facecolor=C_SURF)
    plt.close(fig)

if __name__ == '__main__':
    main()
