#!/usr/bin/env python
"""
Independent re-derivation for kic-004474607 (KBONUS-BKG, Kepler Bonus background dwarf).

Written from scratch with numpy / scipy / astropy only. Deliberately shares NO code with
ZyushgPlanet/pipeline (not imported, not read). This is step 4 of the "before anything goes
public" list in CLAUDE.md: a second, independent re-derivation of depth and period.

Stages (cached to .npz so the long search is resumable):
  search       load -> quality==0 -> per-quarter normalise -> running-median detrend (1 d window)
               -> blind BLS 1-50 d on my own geometric period grid, real and inverted (2-f).
  characterise fine BLS around the best peak, box depth with random-phase empirical error,
               odd/even, secondary at phase 0.5, per-year / per-quarter / per-event depths,
               SAP_FLUX depth at the same ephemeris, figure + results JSON.

Usage:
  python -I rederive_kic004474607.py search
  python -I rederive_kic004474607.py characterise
"""
import sys, os, json, time as _time
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from astropy.time import Time

FNAME = "kic-004474607"
FITS = ("/home/user/ZyushgPlanet/data/kbonus_batch/"
        f"hlsp_kbonus-bkg_kepler_kepler_{FNAME}_kepler_v1.0_lc.fits")
OUT = f"/home/user/ZyushgPlanet/runs/poolA_main/adversarial/{FNAME}/rederive"
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------- search set-up
PMIN, PMAX = 1.0, 50.0
DURATIONS = np.array([0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.20, 0.25, 0.30])   # days
DETREND_WINDOW = 1.0        # days, running median
CLIP_SIGMA_POS = 3.0        # positive outliers clipped (flares / cosmic rays)
CLIP_SIGMA_NEG = 6.0        # only extreme negative outliers clipped (protect transits)
PIPELINE_P = 32.94270       # value quoted in the task for comparison only (not used in search)

# ----------------------------------------------------------------------------- loader

def load(col="FLUX"):
    """Return time (BKJD), flux (e/s), flux_err, quarter for SAP_QUALITY==0 finite cadences."""
    with fits.open(FITS) as h:
        d = h[1].data   # LIGHTCURVE_STITCHED
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d[col], float)
        ecol = col + "_ERR" if col + "_ERR" in d.columns.names else "FLUX_ERR"
        e = np.asarray(d[ecol], float)
        q = np.asarray(d["QUARTER"], int)
        qual = np.asarray(d["SAP_QUALITY"], int)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (f > 0)
    return t[m], f[m], e[m], q[m]


def segments(t, q, gap=0.5):
    """Split into contiguous segments: new segment at quarter change or gap > `gap` days."""
    brk = np.where((np.diff(t) > gap) | (np.diff(q) != 0))[0] + 1
    edges = np.concatenate([[0], brk, [len(t)]])
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def running_median(t, y, window, mask=None):
    """Running median of y in a centred time window of `window` days.
    Points with mask==False are excluded from the median but still get a trend value."""
    if mask is None:
        mask = np.ones_like(y, bool)
    tt, yy = t[mask], y[mask]
    half = window / 2.0
    lo = np.searchsorted(tt, t - half, side="left")
    hi = np.searchsorted(tt, t + half, side="right")
    out = np.empty_like(y)
    for i in range(len(t)):
        a, b = lo[i], hi[i]
        if b - a < 5:   # too few points in window: widen
            a, b = max(0, a - 5), min(len(tt), b + 5)
        out[i] = np.median(yy[a:b]) if b > a else np.nan
    return out


def detrend(t, f, e, q, window=DETREND_WINDOW, in_transit=None, n_iter=2):
    """Per-quarter normalise, running-median flatten per segment, sigma-clip.
    in_transit: optional boolean mask of cadences to exclude from the trend estimate
    (used for the masked re-detrend at the found ephemeris)."""
    fn = f.copy(); en = e.copy()
    for Q in np.unique(q):
        m = q == Q
        med = np.median(f[m])
        fn[m] = f[m] / med
        en[m] = e[m] / med
    trend = np.ones_like(fn)
    keep = np.ones_like(fn, bool)
    if in_transit is not None:
        keep &= ~in_transit
    for it in range(n_iter):
        for a, b in segments(t, q):
            trend[a:b] = running_median(t[a:b], fn[a:b], window, mask=keep[a:b])
        resid = fn / trend - 1.0
        mad = 1.4826 * np.median(np.abs(resid[keep] - np.median(resid[keep])))
        good = (resid < CLIP_SIGMA_POS * mad) & (resid > -CLIP_SIGMA_NEG * mad)
        keep = good if in_transit is None else (good & ~in_transit)
    flat = fn / trend
    eflat = en / trend
    ok = np.isfinite(flat) & (np.abs(flat - 1) < CLIP_SIGMA_POS * mad) | (flat - 1 < 0)
    ok &= np.isfinite(flat) & (flat - 1 > -CLIP_SIGMA_NEG * mad)
    return flat, eflat, trend, ok, mad


# ----------------------------------------------------------------------------- period grid

def period_grid(t, pmin=PMIN, pmax=PMAX, dmin=DURATIONS.min(), frac=0.5):
    """Geometric grid: dP/P = frac * dmin / baseline, so that the accumulated phase drift over
    the full baseline for a neighbouring trial period is < frac * (shortest box)."""
    T = t.max() - t.min()
    step = frac * dmin / T
    n = int(np.ceil(np.log(pmax / pmin) / step)) + 1
    return pmin * np.exp(np.arange(n) * step)


def sde(power):
    p = np.asarray(power, float)
    return (p.max() - p.mean()) / p.std()


def top_peaks(periods, power, n=3, exclude_frac=0.03):
    """Top n distinct peaks: after taking a peak, exclude +-exclude_frac in period around it."""
    p = power.copy()
    out = []
    for _ in range(n):
        i = int(np.nanargmax(p))
        out.append((float(periods[i]), float(power[i]), i))
        m = np.abs(periods / periods[i] - 1) < exclude_frac
        p[m] = -np.inf
    return out


def alias_label(P, Pref):
    r = P / Pref
    for num, den in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (4, 1)]:
        if abs(r - num / den) < 0.01:
            return f"{num}/{den} x P_ref" if (num, den) != (1, 1) else "= P_ref"
    return "unrelated"


# ----------------------------------------------------------------------------- stage: search

def stage_search():
    t0_wall = _time.time()
    t, f, e, q = load("FLUX")
    flat, eflat, trend, ok, mad = detrend(t, f, e, q)
    print(f"loaded {len(t)} good cadences; after clip {ok.sum()}; robust scatter {mad*1e6:.0f} ppm")
    t, y, ye, q = t[ok], flat[ok], eflat[ok], q[ok]
    np.savez(f"{OUT}/detrended_psf.npz", t=t, y=y, ye=ye, q=q, mad=mad)

    periods = period_grid(t)
    print(f"period grid: {len(periods)} trial periods {periods[0]:.3f}-{periods[-1]:.3f} d, "
          f"{len(DURATIONS)} durations")
    res = {}
    for label, yy in [("real", y), ("inverted", 2.0 - y)]:
        tic = _time.time()
        bls = BoxLeastSquares(t, yy, dy=ye)
        r = bls.power(periods, DURATIONS, objective="likelihood", method="fast")
        pw = np.asarray(r.power)
        res[label] = dict(power=pw, depth=np.asarray(r.depth), duration=np.asarray(r.duration),
                          t0=np.asarray(r.transit_time), depth_snr=np.asarray(r.depth_snr))
        i = int(np.argmax(pw))
        print(f"[{label}] {_time.time()-tic:.0f}s  best P={periods[i]:.5f} d  dur={r.duration[i]:.3f} d "
              f"depth={r.depth[i]*1e6:.0f} ppm  SDE={sde(pw):.2f}")
        for P, pwv, idx in top_peaks(periods, pw):
            print(f"    peak P={P:.5f}  power={pwv:.1f}  dur={r.duration[idx]:.3f}  "
                  f"depth={r.depth[idx]*1e6:.0f} ppm  {alias_label(P, PIPELINE_P)}")
    np.savez(f"{OUT}/bls_search.npz", periods=periods,
             **{f"{k}_{kk}": vv for k, v in res.items() for kk, vv in v.items()})
    print(f"search stage done in {_time.time()-t0_wall:.0f}s")


# ----------------------------------------------------------------------------- stage: characterise

def box_depth(t, y, P, t0, dur, phase_centre=0.0):
    """Mean(out) - mean(in) for a box of width dur centred at phase_centre (in days: phase*P).
    Out-of-transit = everything farther than 1.5*dur from the box centre. Returns depth, n_in."""
    ph = ((t - t0) / P + 0.5) % 1.0 - 0.5          # phase in [-0.5, 0.5)
    dph = ((ph - phase_centre) + 0.5) % 1.0 - 0.5  # phase distance from box centre
    dt = dph * P
    inb = np.abs(dt) < dur / 2
    oot = np.abs(dt) > 1.5 * dur
    if inb.sum() == 0:
        return np.nan, 0
    return float(np.mean(y[oot]) - np.mean(y[inb])), int(inb.sum())


def random_phase_boxes(t, y, P, t0, dur, n=3000, avoid=(0.0,), avoid_width=None, seed=1):
    """Empirical depth-error: measure box 'depths' at n random phases that avoid the given
    phase centres by avoid_width (days)."""
    rng = np.random.default_rng(seed)
    if avoid_width is None:
        avoid_width = 2.0 * dur
    vals = []
    tries = 0
    while len(vals) < n and tries < 20 * n:
        tries += 1
        pc = rng.uniform(-0.5, 0.5)
        if any(abs(((pc - a) + 0.5) % 1.0 - 0.5) * P < avoid_width for a in avoid):
            continue
        d, nin = box_depth(t, y, P, t0, dur, phase_centre=pc)
        if nin > 0 and np.isfinite(d):
            vals.append(d)
    vals = np.asarray(vals)
    return float(np.std(vals)), float(np.mean(vals)), vals


def per_event(t, y, P, t0, dur, cadence=0.0204):
    """Per-transit-epoch local depth: in-transit mean vs local out-of-transit (within +-1 d
    but > 1.5 dur away). Returns list of dicts."""
    n_exp = max(1, int(round(dur / cadence)))
    epochs = np.arange(np.floor((t.min() - t0) / P) - 1, np.ceil((t.max() - t0) / P) + 2).astype(int)
    out = []
    for n in epochs:
        tc = t0 + n * P
        if tc < t.min() - dur or tc > t.max() + dur:
            continue
        dt = t - tc
        inb = np.abs(dt) < dur / 2
        loc = (np.abs(dt) > 1.5 * dur) & (np.abs(dt) < 1.0)
        if inb.sum() == 0 or loc.sum() < 10:
            out.append(dict(epoch=int(n), tc=float(tc), n_in=int(inb.sum()), n_loc=int(loc.sum()),
                            depth=None, err=None, covered=False))
            continue
        d = float(np.mean(y[loc]) - np.mean(y[inb]))
        s = float(np.std(y[loc]))
        err = s * np.sqrt(1.0 / inb.sum() + 1.0 / loc.sum())
        # distance to nearest data gap > 0.5 d (or segment end)
        out.append(dict(epoch=int(n), tc=float(tc), n_in=int(inb.sum()), n_loc=int(loc.sum()),
                        depth=d, err=float(err), covered=bool(inb.sum() >= 0.5 * n_exp)))
    return out


def stage_characterise(seed=None, tag=""):
    """seed: characterise the local BLS maximum within +-0.5% of this period instead of the global
    top peak (used to characterise the pipeline's ~32.94 d signal when it is not the top peak).
    tag: suffix for output files."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    S = np.load(f"{OUT}/bls_search.npz")
    D = np.load(f"{OUT}/detrended_psf.npz")
    periods = S["periods"]; pw = S["real_power"]; pw_inv = S["inverted_power"]
    t, y, ye, q = D["t"], D["y"], D["ye"], D["q"]

    peaks = top_peaks(periods, pw)
    if seed is None:
        P_best, _, ibest = peaks[0]
    else:
        msk = np.abs(periods / seed - 1) < 0.005
        ibest = int(np.where(msk)[0][np.argmax(pw[msk])]); P_best = float(periods[ibest])
    rank = int(np.sum(pw > pw[ibest]) )  # number of grid points with higher power
    dur0 = float(S["real_duration"][ibest]); t00 = float(S["real_t0"][ibest])
    print(f"coarse best: P={P_best:.5f} dur={dur0:.3f} t0={t00:.4f} (grid points with higher power: {rank})")

    # ---- fine BLS around the best peak (own grid, finer durations)
    fine_P = np.linspace(P_best * (1 - 0.004), P_best * (1 + 0.004), 4001)
    fine_dur = np.arange(0.03, 0.31, 0.01)
    bls = BoxLeastSquares(t, y, dy=ye)
    rf = bls.power(fine_P, fine_dur, objective="likelihood", method="fast", oversample=20)
    j = int(np.argmax(rf.power))
    P = float(fine_P[j]); dur = float(rf.duration[j]); t0 = float(rf.transit_time[j])
    # 1-sigma period uncertainty: width of the fine peak where power drops by 0.5 (delta-logL)
    above = fine_P[rf.power >= rf.power[j] - 0.5]
    P_err = float(0.5 * (above.max() - above.min())) if len(above) > 1 else float(np.nan)
    print(f"fine: P={P:.5f}+-{P_err:.5f} dur={dur:.3f} d ({dur*24:.2f} h) t0={t0:.4f}")

    # ---- masked re-detrend at this ephemeris (exclude +-1.5 dur around each transit from trend)
    tr, fr, er, qr = load("FLUX")
    ph = ((tr - t0) / P + 0.5) % 1.0 - 0.5
    intr = np.abs(ph * P) < 1.5 * dur
    flat2, eflat2, trend2, ok2, mad2 = detrend(tr, fr, er, qr, in_transit=intr)
    t2, y2, ye2, q2 = tr[ok2], flat2[ok2], eflat2[ok2], qr[ok2]
    np.savez(f"{OUT}/detrended_psf_masked{tag}.npz", t=t2, y=y2, ye=ye2, q=q2, mad=mad2)

    # ---- depth with empirical error (random-phase boxes), on masked re-detrend
    depth, n_in = box_depth(t2, y2, P, t0, dur)
    err_emp, bias, rvals = random_phase_boxes(t2, y2, P, t0, dur, avoid=(0.0, 0.5))
    oot = np.abs(((t2 - t0) / P + 0.5) % 1.0 - 0.5) * P > 1.5 * dur
    err_analytic = float(np.std(y2[oot]) / np.sqrt(n_in))
    snr_emp = depth / err_emp
    # same on the un-masked detrend (what the blind search saw)
    depth_unmasked, n_in_u = box_depth(t, y, P, t0, dur)

    # ---- odd / even
    ph2 = (t2 - t0) / P
    epoch = np.round(ph2).astype(int)
    dt2 = (ph2 - epoch) * P
    inb = np.abs(dt2) < dur / 2
    oot_all = np.abs(dt2) > 1.5 * dur
    base = np.mean(y2[oot_all])
    odd = inb & (epoch % 2 == 1); even = inb & (epoch % 2 == 0)
    s_oot = np.std(y2[oot_all])
    d_odd = base - np.mean(y2[odd]); e_odd = s_oot / np.sqrt(odd.sum())
    d_even = base - np.mean(y2[even]); e_even = s_oot / np.sqrt(even.sum())
    # empirical scaling: use ratio of empirical/analytic error from the full box
    scale = err_emp / err_analytic if err_analytic > 0 else 1.0
    e_odd *= scale; e_even *= scale
    odd_even_sigma = abs(d_odd - d_even) / np.hypot(e_odd, e_even)

    # ---- secondary at phase 0.5 (same box), also scan phases 0.3-0.7 for the deepest box
    sec, n_sec = box_depth(t2, y2, P, t0, dur, phase_centre=0.5)
    sec_err, _, _ = random_phase_boxes(t2, y2, P, t0, dur, avoid=(0.0, 0.5), seed=2)
    sec_sigma = sec / sec_err
    scan_ph = np.linspace(0.2, 0.8, 241)
    scan = np.array([box_depth(t2, y2, P, t0, dur, phase_centre=pc)[0] for pc in scan_ph])
    k = int(np.nanargmax(scan))
    sec_scan = dict(phase=float(scan_ph[k]), depth_ppm=float(scan[k] * 1e6),
                    sigma=float(scan[k] / sec_err))

    # ---- per-event
    ev = per_event(t2, y2, P, t0, dur)
    cov = [x for x in ev if x["covered"]]
    n_transits = len(cov)
    ev_d = np.array([x["depth"] for x in cov]); ev_e = np.array([x["err"] for x in cov])
    ev_snr = ev_d / ev_e
    # single-event dominance: SNR of the stacked signal with the strongest event removed
    w = 1 / ev_e**2
    stack_snr = float(np.sum(ev_d * w) / np.sqrt(np.sum(w)))
    k_top = int(np.argmax(ev_snr)); m_rest = np.ones(len(cov), bool); m_rest[k_top] = False
    stack_snr_minus_top = float(np.sum(ev_d[m_rest] * w[m_rest]) / np.sqrt(np.sum(w[m_rest])))
    chi2_const = float(np.sum((ev_d - np.average(ev_d, weights=w))**2 * w))

    # ---- per calendar year and per quarter depths
    years = Time(t2 + 2454833.0, format="jd").decimalyear.astype(int)
    per_year = {}
    for Y in np.unique(years):
        m = years == Y
        d_y, n_y = box_depth(t2[m], y2[m], P, t0, dur)
        e_y = np.std(y2[m & oot_all]) / np.sqrt(max(n_y, 1)) * scale
        n_ev = len([x for x in cov if int(Time(x["tc"] + 2454833.0, format="jd").decimalyear) == Y])
        per_year[int(Y)] = dict(depth_ppm=float(d_y * 1e6), err_ppm=float(e_y * 1e6),
                                n_in=int(n_y), n_transits=n_ev)
    per_quarter = {}
    for Q in np.unique(q2):
        m = q2 == Q
        d_q, n_q = box_depth(t2[m], y2[m], P, t0, dur)
        e_q = np.std(y2[m & oot_all]) / np.sqrt(max(n_q, 1)) * scale
        per_quarter[int(Q)] = dict(depth_ppm=float(d_q * 1e6), err_ppm=float(e_q * 1e6), n_in=int(n_q))

    # ---- duration from the folded, binned profile (FWHM) as a cross-check on the box
    fold_dt = dt2.copy()
    binw = 0.0204
    edges = np.arange(-0.5, 0.5 + binw, binw)
    idx = np.digitize(fold_dt, edges) - 1
    prof = np.array([np.mean(y2[idx == i]) if np.sum(idx == i) > 3 else np.nan for i in range(len(edges) - 1)])
    cen = 0.5 * (edges[1:] + edges[:-1])
    dip = base - prof
    half = 0.5 * np.nanmax(dip[np.abs(cen) < 0.2])
    below = cen[(dip > half) & (np.abs(cen) < 0.3)]
    fwhm = float(below.max() - below.min() + binw) if len(below) else float(np.nan)

    # ---- SAP_FLUX at the same ephemeris (dilution check)
    ts, fs, es, qs = load("SAP_FLUX")
    phs = ((ts - t0) / P + 0.5) % 1.0 - 0.5
    intr_s = np.abs(phs * P) < 1.5 * dur
    flat_s, eflat_s, _, ok_s, mad_s = detrend(ts, fs, es, qs, in_transit=intr_s)
    ts, ys = ts[ok_s], flat_s[ok_s]
    d_sap, n_sap = box_depth(ts, ys, P, t0, dur)
    e_sap, _, _ = random_phase_boxes(ts, ys, P, t0, dur, avoid=(0.0, 0.5), seed=3)
    # and PSF_FLAT_FLUX (KBONUS's own flattened PSF flux) as a third opinion
    tp, fp, ep, qp = load("PSF_FLAT_FLUX")
    php = ((tp - t0) / P + 0.5) % 1.0 - 0.5
    flat_p, _, _, ok_p, mad_p = detrend(tp, fp, ep, qp, in_transit=np.abs(php * P) < 1.5 * dur)
    d_pff, n_pff = box_depth(tp[ok_p], flat_p[ok_p], P, t0, dur)
    e_pff, _, _ = random_phase_boxes(tp[ok_p], flat_p[ok_p], P, t0, dur, avoid=(0.0, 0.5), seed=4)
    ratio = d_sap / depth
    ratio_err = ratio * np.hypot(e_sap / d_sap, err_emp / depth) if d_sap != 0 else np.nan

    # header crowding numbers for interpretation
    with fits.open(FITS) as h:
        crowd = {}
        for hdu in h[2:]:
            if hdu.name.startswith("LIGHTCURVE_Q"):
                crowd[hdu.name[-3:].strip("_Q")] = dict(
                    CROWDSAP=float(hdu.header.get("CROWDSAP", np.nan)),
                    FLFRCSAP=float(hdu.header.get("FLFRCSAP", np.nan)),
                    PSFFRAC=float(hdu.header.get("PSFFRAC", np.nan)),
                    NPIXSAP=int(hdu.header.get("NPIXSAP", -1)))

    # ---- physical plausibility (Rstar from target table, 0.87 Rsun)
    Rstar = 0.872151
    rp_rearth = np.sqrt(depth) * Rstar * 109.1
    # expected max duration for central transit on a 0.87 Rsun ~0.9 Msun star
    Mstar = 0.9
    a_au = (Mstar * (P / 365.25)**2)**(1 / 3)
    T14_central_h = 24 * P / np.pi * np.arcsin(Rstar * 0.00465 / a_au)

    # ---- SDE etc
    SDE = float(sde(pw)); SDE_inv = float(sde(pw_inv))
    inv_peaks = top_peaks(periods, pw_inv)

    res = dict(
        fname=FNAME, date="2026-10-17", role="independent re-derivation (no pipeline code read or imported)",
        data=dict(fits=os.path.basename(FITS), column_primary="FLUX (PSF)", quality_cut="SAP_QUALITY==0",
                  n_good=int(len(tr)), n_after_clip=int(len(t2)), quarters=[int(x) for x in np.unique(q2)],
                  baseline_d=float(t2.max() - t2.min()), robust_scatter_ppm=float(mad2 * 1e6)),
        method=dict(normalise="per-quarter median", detrend=f"running median, {DETREND_WINDOW} d window, per segment "
                    f"(quarter / >0.5 d gap), 2 iterations, clip +{CLIP_SIGMA_POS}/-{CLIP_SIGMA_NEG} MAD-sigma",
                    masked_redetrend="same, excluding +-1.5 dur around each transit at the found ephemeris",
                    search=f"astropy BoxLeastSquares, objective=likelihood, {PMIN}-{PMAX} d, geometric grid dP/P = 0.5*dmin/T "
                           f"({len(periods)} periods), durations {DURATIONS.tolist()} d",
                    sde="(max - mean)/std of the BLS log-likelihood power over the whole grid",
                    depth_error="std of box depths at 3000 random phases avoiding +-2 dur of phase 0 and 0.5"),
        blind_search=dict(
            top_peaks=[dict(period_d=P_, power=pw_, duration_d=float(S["real_duration"][i_]),
                            depth_ppm=float(S["real_depth"][i_] * 1e6), t0_bkjd=float(S["real_t0"][i_]),
                            relation_to_pipeline_P=alias_label(P_, PIPELINE_P)) for P_, pw_, i_ in peaks],
            sde=SDE, recovered_pipeline_period=bool(abs(P_best / PIPELINE_P - 1) < 0.002),
            characterised_peak=dict(seed=seed, coarse_period_d=P_best, coarse_power=float(pw[ibest]),
                                    n_grid_points_with_higher_power=rank,
                                    is_global_top_peak=bool(rank == 0)),
            sde_1_40d=float(sde(pw[periods <= 40])), sde_inverted_1_40d=float(sde(pw_inv[periods <= 40])),
            top_peak_1_40d=float(periods[periods <= 40][np.argmax(pw[periods <= 40])]),
            inverted_null=dict(sde=SDE_inv, top_peaks=[dict(period_d=P_, power=pw_,
                               depth_ppm=float(S["inverted_depth"][i_] * 1e6)) for P_, pw_, i_ in inv_peaks]),
            n_periods=int(len(periods))),
        ephemeris=dict(period_d=P, period_err_d=P_err, t0_bkjd=t0, t0_bjd=t0 + 2454833.0,
                       duration_box_d=dur, duration_box_h=dur * 24, duration_fwhm_h=fwhm * 24,
                       fine_grid="+-0.4% in P, 4001 steps; durations 0.03-0.30 d step 0.01"),
        depth=dict(depth_ppm=float(depth * 1e6), depth_err_empirical_ppm=float(err_emp * 1e6),
                   depth_err_analytic_ppm=float(err_analytic * 1e6), empirical_snr=float(snr_emp),
                   n_in_transit=int(n_in), random_box_mean_ppm=float(bias * 1e6),
                   depth_unmasked_detrend_ppm=float(depth_unmasked * 1e6),
                   bls_depth_at_peak_ppm=float(rf.depth[j] * 1e6), bls_depth_snr_at_peak=float(rf.depth_snr[j])),
        odd_even=dict(depth_odd_ppm=float(d_odd * 1e6), err_odd_ppm=float(e_odd * 1e6), n_odd=int(odd.sum()),
                      depth_even_ppm=float(d_even * 1e6), err_even_ppm=float(e_even * 1e6), n_even=int(even.sum()),
                      sigma=float(odd_even_sigma)),
        secondary=dict(phase=0.5, depth_ppm=float(sec * 1e6), err_ppm=float(sec_err * 1e6), sigma=float(sec_sigma),
                       n_in=int(n_sec), deepest_box_phase_0p2_0p8=sec_scan),
        transits=dict(n_transits_covered=n_transits, n_epochs_in_baseline=len(ev),
                      stacked_snr_from_events=stack_snr, stacked_snr_without_strongest=stack_snr_minus_top,
                      strongest_event=dict(epoch=cov[k_top]["epoch"], tc_bkjd=cov[k_top]["tc"],
                                           depth_ppm=float(ev_d[k_top] * 1e6), snr=float(ev_snr[k_top])),
                      chi2_constant_depth=chi2_const, dof=int(len(cov) - 1),
                      events=[dict(x, depth_ppm=None if x["depth"] is None else x["depth"] * 1e6,
                                   err_ppm=None if x["err"] is None else x["err"] * 1e6) for x in ev]),
        per_year=per_year, per_quarter=per_quarter,
        sap_vs_psf=dict(depth_sap_ppm=float(d_sap * 1e6), err_sap_ppm=float(e_sap * 1e6), n_in_sap=int(n_sap),
                        depth_psf_ppm=float(depth * 1e6), depth_psf_flat_flux_ppm=float(d_pff * 1e6),
                        err_psf_flat_flux_ppm=float(e_pff * 1e6),
                        ratio_sap_over_psf=float(ratio), ratio_err=float(ratio_err),
                        header_crowding_per_quarter=crowd,
                        target_table=dict(psffrac=0.609585, crowdsap=0.703728, flfrcsap=0.779184, npixsap=2)),
        physical=dict(Rstar_rsun=Rstar, rp_rearth_if_on_target=float(rp_rearth),
                      T14_central_h_for_0p9Msun=float(T14_central_h), a_au=float(a_au)),
    )
    with open(f"{OUT}/results{tag}.json", "w") as fh:
        json.dump(res, fh, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

    # ----------------------------------------------------------------------- figure
    import matplotlib.gridspec as gridspec
    C_DATA = "#b8b7b1"; C_BIN = "#2a78d6"; C_ODD = "#e34948"; C_EVEN = "#1baf7a"; C_SAP = "#eb6834"; C_INV = "#9a9995"  # dataviz reference palette slots 1,8,3,2 + neutrals
    fig = plt.figure(figsize=(14, 15))
    gs = gridspec.GridSpec(5, 2, figure=fig, hspace=0.42, wspace=0.22)

    ax = fig.add_subplot(gs[0, :])
    ax.plot(periods, pw, color=C_BIN, lw=0.6, label="PSF flux (real)")
    ax.plot(periods, pw_inv, color=C_INV, lw=0.5, alpha=0.8, label="inverted 2-f (null)")
    for k_, (P_, pw_, i_) in enumerate(peaks):
        ax.annotate(f"{P_:.3f} d", (P_, pw_), textcoords="offset points", xytext=(4, 2), fontsize=8, color=C_ODD)
    ax.axvline(PIPELINE_P, color=C_ODD, ls=":", lw=0.8, label=f"pipeline P={PIPELINE_P} d")
    ax.set_xscale("log"); ax.set_xlabel("period [d]"); ax.set_ylabel("BLS delta log L")
    ax.set_title(f"{FNAME}: blind BLS 1-50 d, own grid ({len(periods)} periods). SDE real={SDE:.1f}, inverted={SDE_inv:.1f}")
    ax.legend(fontsize=8, loc="upper left")

    ax = fig.add_subplot(gs[1, :])
    ax.plot(t2, (y2 - 1) * 1e6, ".", ms=1.5, color=C_DATA)
    for x in ev:
        ax.axvline(x["tc"], color=C_ODD if x["covered"] else "#dddddd", lw=0.6, alpha=0.7)
    ax.set_ylim(-9000, 6000); ax.set_xlabel("BKJD [d]"); ax.set_ylabel("ppm")
    ax.set_title(f"detrended PSF flux (masked re-detrend); red = {n_transits} covered transit epochs, grey = epochs in gaps")

    def fold_panel(ax, tt, yy, lab, color, win=0.6):
        phd = (((tt - t0) / P + 0.5) % 1.0 - 0.5) * P
        m = np.abs(phd) < win
        ax.plot(phd[m] * 24, (yy[m] - 1) * 1e6, ".", ms=2, color=C_DATA, alpha=0.6)
        bw = 0.0204 * 24
        be = np.arange(-win * 24, win * 24 + bw, bw)
        bi = np.digitize(phd[m] * 24, be) - 1
        bm = np.array([np.mean(yy[m][bi == i]) if np.sum(bi == i) > 2 else np.nan for i in range(len(be) - 1)])
        bs = np.array([np.std(yy[m][bi == i]) / np.sqrt(np.sum(bi == i)) if np.sum(bi == i) > 2 else np.nan for i in range(len(be) - 1)])
        bc = 0.5 * (be[1:] + be[:-1])
        ax.errorbar(bc, (bm - 1) * 1e6, yerr=bs * 1e6, fmt="o", ms=3.5, color=color, lw=1, label=lab)
        ax.axvspan(-dur * 12, dur * 12, color=color, alpha=0.08)
        ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("ppm"); ax.set_ylim(-5500, 3500)
        ax.legend(fontsize=8, loc="lower left")

    ax = fig.add_subplot(gs[2, 0])
    fold_panel(ax, t2, y2, f"PSF: depth {depth*1e6:.0f}+-{err_emp*1e6:.0f} ppm, SNR {snr_emp:.1f}", C_BIN)
    ax.set_title(f"fold at P={P:.5f} d, t0={t0:.3f} BKJD, box {dur*24:.2f} h")

    ax = fig.add_subplot(gs[2, 1])
    mo = (epoch % 2 == 1); me_ = ~mo
    phd = dt2
    for mm, col, lab in [(mo, C_ODD, f"odd {d_odd*1e6:.0f}+-{e_odd*1e6:.0f}"), (me_, C_EVEN, f"even {d_even*1e6:.0f}+-{e_even*1e6:.0f}")]:
        sel = mm & (np.abs(phd) < 0.6)
        bw = 0.0204 * 24 * 1.5
        be = np.arange(-0.6 * 24, 0.6 * 24 + bw, bw)
        bi = np.digitize(phd[sel] * 24, be) - 1
        bm = np.array([np.mean(y2[sel][bi == i]) if np.sum(bi == i) > 2 else np.nan for i in range(len(be) - 1)])
        bs = np.array([np.std(y2[sel][bi == i]) / np.sqrt(np.sum(bi == i)) if np.sum(bi == i) > 2 else np.nan for i in range(len(be) - 1)])
        ax.errorbar(0.5 * (be[1:] + be[:-1]), (bm - 1) * 1e6, yerr=bs * 1e6, fmt="o", ms=3, color=col, lw=1, label=lab + " ppm")
    ax.axvspan(-dur * 12, dur * 12, color="k", alpha=0.06)
    ax.set_ylim(-5500, 3500); ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("ppm")
    ax.set_title(f"odd vs even: |diff| = {odd_even_sigma:.1f} sigma"); ax.legend(fontsize=8, loc="lower left")

    ax = fig.add_subplot(gs[3, 0])
    phs2 = (((t2 - t0) / P) % 1.0)   # 0..1, secondary at 0.5
    m = np.abs(phs2 - 0.5) * P < 0.6
    ax.plot((phs2[m] - 0.5) * P * 24, (y2[m] - 1) * 1e6, ".", ms=2, color=C_DATA, alpha=0.6)
    bw = 0.0204 * 24
    be = np.arange(-0.6 * 24, 0.6 * 24 + bw, bw)
    bi = np.digitize((phs2[m] - 0.5) * P * 24, be) - 1
    bm = np.array([np.mean(y2[m][bi == i]) if np.sum(bi == i) > 2 else np.nan for i in range(len(be) - 1)])
    bs = np.array([np.std(y2[m][bi == i]) / np.sqrt(np.sum(bi == i)) if np.sum(bi == i) > 2 else np.nan for i in range(len(be) - 1)])
    ax.errorbar(0.5 * (be[1:] + be[:-1]), (bm - 1) * 1e6, yerr=bs * 1e6, fmt="o", ms=3.5, color=C_BIN, lw=1)
    ax.axvspan(-dur * 12, dur * 12, color="k", alpha=0.06)
    ax.set_ylim(-5500, 3500); ax.set_xlabel("hours from phase 0.5"); ax.set_ylabel("ppm")
    ax.set_title(f"secondary at phase 0.5: {sec*1e6:.0f}+-{sec_err*1e6:.0f} ppm ({sec_sigma:.1f} sigma)")

    ax = fig.add_subplot(gs[3, 1])
    fold_panel(ax, ts, ys, f"SAP: depth {d_sap*1e6:.0f}+-{e_sap*1e6:.0f} ppm (SAP/PSF = {ratio:.2f})", C_SAP)
    ax.set_title("same ephemeris on SAP_FLUX (dilution check)")

    ax = fig.add_subplot(gs[4, 0])
    tcs = np.array([x["tc"] for x in cov])
    ax.errorbar(tcs, ev_d * 1e6, yerr=ev_e * 1e6, fmt="o", ms=4, color=C_BIN, lw=1)
    ax.axhline(depth * 1e6, color=C_ODD, lw=1, label=f"stacked {depth*1e6:.0f} ppm")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("BKJD of event"); ax.set_ylabel("local depth [ppm]")
    ax.set_title(f"per-event depths ({n_transits} covered); stacked SNR {stack_snr:.1f}, without strongest {stack_snr_minus_top:.1f}")
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[4, 1])
    ys_ = sorted(per_year)
    ax.errorbar([str(Y) for Y in ys_], [per_year[Y]["depth_ppm"] for Y in ys_], yerr=[per_year[Y]["err_ppm"] for Y in ys_],
                fmt="s", ms=6, color=C_BIN, lw=1.2, label="per calendar year")
    qs_ = sorted(per_quarter)
    ax.errorbar([f"Q{Q}" for Q in qs_], [per_quarter[Q]["depth_ppm"] for Q in qs_], yerr=[per_quarter[Q]["err_ppm"] for Q in qs_],
                fmt="o", ms=4, color=C_SAP, lw=1, label="per quarter")
    ax.axhline(depth * 1e6, color=C_ODD, lw=1); ax.axhline(0, color="k", lw=0.5)
    ax.set_ylabel("depth [ppm]"); ax.set_title("depth consistency by year / quarter"); ax.legend(fontsize=8)
    ax.tick_params(axis="x", labelsize=8)

    fig.suptitle(f"{FNAME} independent re-derivation (no pipeline code): P={P:.5f} d, depth={depth*1e6:.0f} ppm, "
                 f"dur={dur*24:.2f} h, {n_transits} transits", y=0.995, fontsize=12)
    fig.savefig(f"{OUT}/fold_rederive{tag}.png", dpi=130, bbox_inches="tight")
    print(json.dumps({k: res[k] for k in ("ephemeris", "depth", "odd_even", "secondary", "sap_vs_psf", "per_year")},
                     indent=1, default=str)[:6000])
    print("transits:", n_transits, "stack SNR", stack_snr, "minus top", stack_snr_minus_top)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "search"
    if stage == "search":
        stage_search()
    elif stage == "characterise":
        # characterise [seed_period] [tag]
        seed = float(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != "top" else None
        tag = sys.argv[3] if len(sys.argv) > 3 else ""
        stage_characterise(seed=seed, tag=tag)
    else:
        raise SystemExit(f"unknown stage {stage}")
