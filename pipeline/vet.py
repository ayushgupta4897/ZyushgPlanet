"""Vetting diagnostics at a fixed ephemeris, plus null tests.

All depths are in fractional units (multiply by 1e6 for ppm) measured on the detrended flux (~1.0).
Box depth = 1 - mean(in-transit flux); its error combines the scatter of the in-transit points and the
scatter of the out-of-transit points (both robust), so it is NOT the optimistic formal white-noise error.
"""
from __future__ import annotations

import numpy as np

from .fetch import BKJD_OFFSET


def phase(t, period, t0):
    """Phase in [-0.5, 0.5) with transit at 0."""
    return ((t - t0 + 0.5 * period) % period) / period - 0.5


def epoch(t, period, t0):
    return np.round((t - t0) / period).astype(int)


def _robust_sigma(x):
    x = np.asarray(x, float)
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) > 2 else np.std(x)


def box_depth(t, f, period, t0, duration, phase_center=0.0, out_min_sep=1.5):
    """Depth of a box of width `duration` centred at `phase_center` (fraction of period) in the fold."""
    ph = phase(t, period, t0) - phase_center
    ph = (ph + 0.5) % 1.0 - 0.5
    half = 0.5 * duration / period
    inn = np.abs(ph) < half
    out = np.abs(phase(t, period, t0)) > out_min_sep * half  # away from the primary transit
    if inn.sum() < 2:
        return {"depth": np.nan, "depth_err": np.nan, "n_in": int(inn.sum()), "n_out": int(out.sum())}
    base = np.median(f[out])
    depth = base - np.mean(f[inn])
    err = np.sqrt(_robust_sigma(f[inn]) ** 2 / inn.sum() + _robust_sigma(f[out]) ** 2 / max(out.sum(), 1))
    return {"depth": float(depth), "depth_err": float(err), "n_in": int(inn.sum()), "n_out": int(out.sum()),
            "sigma": float(depth / err)}


def odd_even(t, f, period, t0, duration):
    ep = epoch(t, period, t0)
    res = {}
    for name, sel in (("odd", ep % 2 == 1), ("even", ep % 2 == 0)):
        res[name] = box_depth(t[sel], f[sel], period, t0, duration)
    d = res["odd"]["depth"] - res["even"]["depth"]
    e = np.hypot(res["odd"]["depth_err"], res["even"]["depth_err"])
    res["difference_sigma"] = float(d / e) if e > 0 else np.nan
    return res


def secondary_eclipse(t, f, period, t0, duration, n_phase=200):
    """Box depth at phase 0.5 and a scan of all phases away from the primary (most significant other dip)."""
    sec = box_depth(t, f, period, t0, duration, phase_center=0.5)
    half = 0.5 * duration / period
    grid = np.linspace(-0.5, 0.5, n_phase, endpoint=False)
    grid = grid[np.abs(grid) > 3 * half]
    scan = []
    for pc in grid:
        b = box_depth(t, f, period, t0, duration, phase_center=pc)
        scan.append((pc, b["depth"], b["depth_err"]))
    scan = np.array(scan)
    sig = scan[:, 1] / scan[:, 2]
    k = int(np.nanargmax(sig))
    return {"phase_0p5": sec, "scan_max_dip": {"phase": float(scan[k, 0]), "depth": float(scan[k, 1]),
                                                "depth_err": float(scan[k, 2]), "sigma": float(sig[k])},
            "scan_n_phases": int(len(scan)), "scan_sigma_std": float(np.nanstd(sig))}


def per_transit_depths(t, f, q, period, t0, duration, local_window_factor=4.0, min_in=2, min_out=8):
    """Local depth of every individual event: mean(local out-of-transit) - mean(in-transit)."""
    ep = epoch(t, period, t0)
    out = []
    for n in np.unique(ep):
        tc = t0 + n * period
        dt = t - tc
        inn = np.abs(dt) < 0.5 * duration
        loc = (np.abs(dt) > 0.75 * duration) & (np.abs(dt) < local_window_factor * duration)
        if inn.sum() < min_in or loc.sum() < min_out:
            continue
        base = np.median(f[loc]); sig = _robust_sigma(f[loc])
        depth = base - np.mean(f[inn])
        err = sig * np.sqrt(1.0 / inn.sum() + 1.0 / loc.sum())
        out.append({"epoch": int(n), "tc_bkjd": float(tc), "tc_bjd_tdb": float(tc + BKJD_OFFSET),
                    "quarter": int(np.median(q[inn])), "n_in": int(inn.sum()), "n_out": int(loc.sum()),
                    "depth_ppm": float(depth * 1e6), "depth_err_ppm": float(err * 1e6)})
    if out:
        d = np.array([o["depth_ppm"] for o in out]); e = np.array([o["depth_err_ppm"] for o in out])
        w = 1 / e ** 2
        summary = {"n_events": len(out), "weighted_mean_depth_ppm": float(np.sum(w * d) / np.sum(w)),
                   "weighted_mean_err_ppm": float(1 / np.sqrt(np.sum(w))), "median_depth_ppm": float(np.median(d)),
                   "scatter_of_event_depths_ppm": float(np.std(d)), "n_positive": int((d > 0).sum()),
                   "max_single_event_fraction": float(np.max(w * d) / np.sum(w * d)) if np.sum(w * d) > 0 else np.nan,
                   "chi2_constant_depth": float(np.sum(((d - np.sum(w * d) / np.sum(w)) / e) ** 2)), "dof": len(out) - 1}
    else:
        summary = {"n_events": 0}
    return {"events": out, "summary": summary}


def empirical_noise(t, f, period, t0, duration, n_random=400, seed=1):
    """Fold at random T0 offsets (>= 2 durations from the real transit) and measure the same box depth.
    The std of those depths is the empirical noise of a folded box depth at this period/duration, including
    red noise; empirical SNR = measured depth / that std."""
    rng = np.random.default_rng(seed)
    half = duration / period
    offs = rng.uniform(2 * half, 1 - 2 * half, n_random) * period
    depths = np.array([box_depth(t, f, period, t0 + o, duration, out_min_sep=0.0)["depth"] for o in offs])
    depths = depths[np.isfinite(depths)]
    real = box_depth(t, f, period, t0, duration)
    sd = float(np.std(depths))
    return {"n_random": int(len(depths)), "random_depth_mean_ppm": float(np.mean(depths) * 1e6),
            "random_depth_std_ppm": sd * 1e6, "random_depth_max_ppm": float(np.max(depths) * 1e6),
            "real_depth_ppm": real["depth"] * 1e6, "real_depth_formal_err_ppm": real["depth_err"] * 1e6,
            "empirical_snr": float(real["depth"] / sd) if sd > 0 else np.nan,
            "formal_snr": float(real["depth"] / real["depth_err"])}


def ephemeris_depth(t, f, period, t0, duration):
    """Depth of a box at a FIXED (external) ephemeris; used for the neighbour and for Kuber's ephemeris."""
    b = box_depth(t, f, period, t0, duration)
    return {k: (v * 1e6 if k in ("depth", "depth_err") else v) for k, v in b.items()} | {"units": "ppm"}


def compare_ephemeris(period, t0_bkjd, ref_period, ref_t0_bjd, baseline_mid_bkjd):
    """Compare a fitted (P, T0) with a reference ephemeris: period difference and T0 offset modulo P,
    evaluated at the epoch nearest the middle of the baseline (where the two are least covariant)."""
    ref_t0 = ref_t0_bjd - BKJD_OFFSET
    n_ref = np.round((baseline_mid_bkjd - ref_t0) / ref_period)
    tc_ref = ref_t0 + n_ref * ref_period
    n = np.round((tc_ref - t0_bkjd) / period)
    tc = t0_bkjd + n * period
    return {"ref_period_d": ref_period, "ref_t0_bjd_tdb": ref_t0_bjd, "ref_t0_bkjd": ref_t0,
            "period_diff_d": float(period - ref_period), "period_diff_ppm": float((period / ref_period - 1) * 1e6),
            "t0_offset_min_at_mid_baseline": float((tc - tc_ref) * 1440), "epoch_mid_baseline": int(n_ref),
            "fitted_tc_near_mid_bkjd": float(tc), "ref_tc_near_mid_bkjd": float(tc_ref)}
