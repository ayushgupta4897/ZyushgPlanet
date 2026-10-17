"""
Batch transit screen for Kepler Bonus (KBONUS-BKG) light curves.

One star in, one JSON-able dict out.  Designed to be called thousands of times
(locally with multiprocessing, or on Modal with .map).  No plotting here.

Pipeline per star
-----------------
1. download the KBONUS-BKG FITS (direct MAST URL derived from the catalogue `fname`)
2. read LIGHTCURVE_STITCHED: TIME (BKJD = BJD - 2454833), QUARTER, FLUX (PSF flux),
   FLUX_ERR, SAP_QUALITY
3. per-quarter median normalisation, wotan biweight flattening (window_length days),
   upward-outlier clipping only (a dip is never clipped)
4. astropy BoxLeastSquares on a log-spaced period grid; SDE from the periodogram,
   depth SNR from the fit, and an *empirical* depth SNR from random off-transit boxes
5. diagnostics: n_transits, odd/even, secondary at phase 0.5, max single-event share
6. iterate once: mask the first signal and search again (second signal)
7. the identical procedure on the INVERTED flux (2 - f): the null control

Everything that decides whether a star "passes" lives in prereg_*.json, not here.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
import warnings
from dataclasses import dataclass, asdict

import numpy as np

warnings.filterwarnings("ignore")

BKJD_OFFSET = 2454833.0
MAST_BASE = "https://archive.stsci.edu/hlsps/kbonus-bkg/lcs"


# --------------------------------------------------------------------------- #
# data access
# --------------------------------------------------------------------------- #
def kbonus_url(fname: str) -> str:
    """Direct MAST URL for a KBONUS-BKG light-curve file.

    fname is the catalogue column, e.g. 'kic-005689121' or 'gaia-dr3-2050192850525687424'.
    Verified 2026-10-17:
      kic-005689121            -> lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits
      gaia-dr3-2050192850525687424 -> lcs/2050/2050192850525687424/hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2050192850525687424_kepler_v1.0_lc.fits
    """
    if fname.startswith("kic-"):
        key = fname[4:]
    elif fname.startswith("gaia-dr3-"):
        key = fname[9:]
    else:
        raise ValueError(f"unrecognised fname {fname!r}")
    return f"{MAST_BASE}/{key[:4]}/{key}/hlsp_kbonus-bkg_kepler_kepler_{fname}_kepler_v1.0_lc.fits"


def fetch_kbonus(fname: str, cache_dir: str | None = None, keep: bool = True, timeout: int = 120) -> str:
    """Download (or reuse) the FITS file; returns the local path."""
    import requests

    url = kbonus_url(fname)
    cache_dir = cache_dir or os.path.join(os.getcwd(), "data", "kbonus_batch")
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, os.path.basename(url))
    if os.path.exists(path) and os.path.getsize(path) > 100_000:
        return path
    tmp = path + ".part"
    for attempt in range(4):
        try:
            with requests.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            os.replace(tmp, path)
            return path
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"download failed for {fname}: {last}")


def read_kbonus(path: str, flux_col: str = "FLUX"):
    """Return time (BKJD), flux, flux_err, quality, quarter, header-derived meta."""
    from astropy.io import fits

    with fits.open(path, memmap=False) as h:
        d = h["LIGHTCURVE_STITCHED"].data
        hdr = h[0].header
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d[flux_col], float)
        fe = np.asarray(d[flux_col + "_ERR"] if flux_col + "_ERR" in d.names else d["FLUX_ERR"], float)
        q = np.asarray(d["SAP_QUALITY"], int)
        quarter = np.asarray(d["QUARTER"], int)
        meta = {
            "object": hdr.get("OBJECT"),
            "gaiaid": hdr.get("GAIAID"),
            "tpforg": hdr.get("TPFORG"),
            "kepmag": hdr.get("KEPMAG"),
            "ra": hdr.get("RA_OBJ"),
            "dec": hdr.get("DEC_OBJ"),
        }
    return t, f, fe, q, quarter, meta


# --------------------------------------------------------------------------- #
# preparation
# --------------------------------------------------------------------------- #
def prepare(t, f, fe, q, quarter, window_length_d: float = 1.0, clip_sigma_up: float = 5.0,
            min_points_per_quarter: int = 200):
    """Quality mask, per-quarter normalisation, biweight flattening, upward clipping.

    Returns t, y (flattened, ~1), e (normalised errors), quarter, and info dict.
    window_length_d: ~3x the longest transit duration you care about.  For P <= 40 d
    around a <= 1 Rsun dwarf the duration is <= ~8 h, so 1.0 d is the conservative choice
    (Chopra's injection-recovery preferred 0.75 d for TESS short-period targets).
    """
    from wotan import flatten

    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (q == 0) & (f > 0)
    t, f, fe, quarter = t[ok], f[ok], fe[ok], quarter[ok]
    order = np.argsort(t)
    t, f, fe, quarter = t[order], f[order], fe[order], quarter[order]

    y = np.empty_like(f)
    e = np.empty_like(f)
    keep = np.ones(len(t), bool)
    for qq in np.unique(quarter):
        m = quarter == qq
        if m.sum() < min_points_per_quarter:
            keep[m] = False
            continue
        med = np.nanmedian(f[m])
        y[m] = f[m] / med
        e[m] = fe[m] / med
    t, y, e, quarter = t[keep], y[keep], e[keep], quarter[keep]
    if len(t) < 1000:
        return t, y, e, quarter, {"n_after_quality": int(len(t)), "too_few": True}

    # biweight per quarter via break_tolerance (gaps > 0.5 d start a new segment)
    flat, trend = flatten(t, y, method="biweight", window_length=window_length_d,
                          break_tolerance=0.5, return_trend=True, edge_cutoff=0.0)
    good = np.isfinite(flat)
    t, flat, e, quarter = t[good], flat[good], e[good], quarter[good]

    # robust scatter and *upward-only* clipping (never clip a dip)
    resid = flat - 1.0
    mad = np.nanmedian(np.abs(resid - np.nanmedian(resid))) * 1.4826
    up = resid > clip_sigma_up * mad
    t, flat, e, quarter = t[~up], flat[~up], e[~up], quarter[~up]
    p2p = np.nanmedian(np.abs(np.diff(flat))) * 1.4826 / np.sqrt(2)
    info = {
        "n_after_quality": int(len(t)),
        "n_quarters": int(len(np.unique(quarter))),
        "baseline_d": float(t.max() - t.min()),
        "mad_ppm": float(mad * 1e6),
        "p2p_ppm": float(p2p * 1e6),
        "cdpp6h_ppm_whitenoise": float(p2p * 1e6 / np.sqrt(12)),
        "n_clipped_up": int(up.sum()),
        "too_few": False,
    }
    return t, flat, e, quarter, info


# --------------------------------------------------------------------------- #
# search
# --------------------------------------------------------------------------- #
def period_grid(pmin: float, pmax: float, baseline: float, min_duration_d: float, oversample: float = 3.0):
    """Log-spaced period grid such that the phase drift across the baseline is
    <= min_duration/oversample at the LONGEST period (conservative everywhere else)."""
    # dP/P * (baseline/P) * P = dP * baseline / P  <= dur/os  ->  dlnP <= dur/(os*baseline)
    dlnp = min_duration_d / (oversample * baseline)
    n = int(np.ceil(np.log(pmax / pmin) / dlnp))
    n = int(min(max(n, 20_000), 400_000))  # hard cap to keep the per-star cost bounded
    return np.exp(np.linspace(np.log(pmin), np.log(pmax), n))


def _sde(power):
    p = np.asarray(power, float)
    p = p[np.isfinite(p)]
    s = np.std(p)
    return float((np.max(p) - np.median(p)) / s) if s > 0 else 0.0


def _empirical_box_noise(t, y, period, duration, t0, n_boxes: int = 60, rng=None):
    """Std of the mean of boxes of width `duration` placed at random phases away
    from the transit (>= 2 durations from phase 0).  Captures correlated noise that
    the formal depth error ignores (Chopra's fix)."""
    rng = rng or np.random.default_rng(42)
    phase = ((t - t0) / period + 0.5) % 1.0 - 0.5  # in [-0.5, 0.5)
    w = duration / period
    vals = []
    tries = 0
    while len(vals) < n_boxes and tries < 20 * n_boxes:
        tries += 1
        c = rng.uniform(-0.5, 0.5)
        if abs(c) < 2.5 * w or abs(abs(c) - 0.5) < 1.5 * w:
            continue
        m = np.abs(phase - c) < w / 2
        if m.sum() >= 5:
            vals.append(np.mean(y[m]))
    if len(vals) < 10:
        return np.nan
    return float(np.std(vals))


def bls_one(t, y, e, pmin, pmax, durations_h, oversample=3.0, objective="likelihood"):
    from astropy.timeseries import BoxLeastSquares

    durations = np.asarray(durations_h, float) / 24.0
    periods = period_grid(pmin, pmax, t.max() - t.min(), durations.min(), oversample)
    bls = BoxLeastSquares(t, y, e)
    res = bls.power(periods, durations, objective=objective)
    i = int(np.nanargmax(res.power))
    P, D, T0 = float(res.period[i]), float(res.duration[i]), float(res.transit_time[i])
    stats = bls.compute_stats(P, D, T0)
    # per-transit share of the log-likelihood (single-event dominance)
    ptl = np.asarray(stats["per_transit_log_likelihood"], float)
    ptc = np.asarray(stats["per_transit_count"], int)
    n_tr = int((ptc >= 2).sum())
    pos = np.clip(ptl, 0, None)
    max_single = float(pos.max() / pos.sum()) if pos.sum() > 0 else np.nan
    depth, depth_err = float(stats["depth"][0]), float(stats["depth"][1])
    dodd, dodd_e = float(stats["depth_odd"][0]), float(stats["depth_odd"][1])
    deven, deven_e = float(stats["depth_even"][0]), float(stats["depth_even"][1])
    # astropy naming: depth_phased = model shifted by half a period (the SECONDARY-eclipse test);
    # depth_half = model at half the period (NOT a secondary test; kept for information only)
    dsec, dsec_e = float(stats["depth_phased"][0]), float(stats["depth_phased"][1])
    dhalf, dhalf_e = float(stats["depth_half"][0]), float(stats["depth_half"][1])
    box_noise = _empirical_box_noise(t, y, P, D, T0)
    have_box = bool(box_noise and np.isfinite(box_noise) and box_noise > 0)
    snr_emp = depth / box_noise if have_box else np.nan
    # formal sigmas (white-noise errors) and empirical sigmas (folded-box scatter, which
    # includes correlated noise).  The screen uses the empirical ones when available.
    oe_sigma_formal = abs(dodd - deven) / np.sqrt(dodd_e**2 + deven_e**2) if (dodd_e > 0 and deven_e > 0) else np.nan
    oe_sigma_emp = abs(dodd - deven) / (2.0 * box_noise) if have_box else np.nan  # each half has ~sqrt(2) x noise
    sec_sigma_formal = dsec / dsec_e if dsec_e > 0 else np.nan
    sec_sigma_emp = dsec / box_noise if have_box else np.nan
    oe_sigma = oe_sigma_emp if np.isfinite(oe_sigma_emp) else oe_sigma_formal
    sec_sigma = sec_sigma_emp if np.isfinite(sec_sigma_emp) else sec_sigma_formal
    # top-5 distinct peaks for aliasing diagnostics
    order = np.argsort(res.power)[::-1]
    peaks = []
    for j in order:
        pj = float(res.period[j])
        if all(abs(pj - pk) / pk > 0.02 for pk, _ in peaks):
            peaks.append((pj, float(res.power[j])))
        if len(peaks) >= 5:
            break
    return {
        "period_d": P, "duration_h": D * 24.0, "t0_bkjd": T0, "t0_bjd": T0 + BKJD_OFFSET,
        "depth_ppm": depth * 1e6, "depth_err_ppm": depth_err * 1e6,
        "snr_fit": float(stats["depth_snr"]) if "depth_snr" in stats else (depth / depth_err if depth_err > 0 else np.nan),
        "snr_emp": snr_emp, "box_noise_ppm": (box_noise * 1e6) if np.isfinite(box_noise) else np.nan,
        "sde": _sde(res.power), "power": float(res.power[i]),
        "n_transits": n_tr, "max_single_event_frac": max_single,
        "depth_odd_ppm": dodd * 1e6, "depth_even_ppm": deven * 1e6,
        "odd_even_sigma": float(oe_sigma) if np.isfinite(oe_sigma) else None,
        "odd_even_sigma_formal": float(oe_sigma_formal) if np.isfinite(oe_sigma_formal) else None,
        "secondary_sigma": float(sec_sigma) if np.isfinite(sec_sigma) else None,
        "secondary_sigma_formal": float(sec_sigma_formal) if np.isfinite(sec_sigma_formal) else None,
        "secondary_depth_ppm": dsec * 1e6,
        "half_period_depth_ppm": dhalf * 1e6,
        "n_periods": int(len(periods)),
        "top_peaks": [[p, w] for p, w in peaks],
        "_mask": _phase_mask(t, P, min(2.0 * D, 0.5 * P), T0),
    }


def _phase_mask(t, period, width_d, t0):
    """In-transit mask of total width `width_d` centred on the ephemeris (no astropy
    validation, so it cannot fail for short periods)."""
    phase = ((t - t0) / period + 0.5) % 1.0 - 0.5
    return np.abs(phase) < 0.5 * width_d / period


def search_star(fname: str, config: dict, cache_dir: str | None = None, keep_fits: bool = False) -> dict:
    """Full per-star screen (real + inverted).  Never raises; errors are returned."""
    t_start = time.time()
    out = {"fname": fname, "config_sha": config.get("_sha"), "ok": False}
    try:
        path = fetch_kbonus(fname, cache_dir=cache_dir, keep=keep_fits)
        t, f, fe, q, quarter, meta = read_kbonus(path, flux_col=config.get("flux_col", "FLUX"))
        out.update(meta)
        out["n_raw"] = int(len(t))
        t, y, e, quarter, info = prepare(t, f, fe, q, quarter,
                                         window_length_d=config.get("window_length_d", 1.0),
                                         clip_sigma_up=config.get("clip_sigma_up", 5.0))
        out.update(info)
        if info.get("too_few"):
            out["ok"] = True
            out["skipped"] = "too few points"
            return out
        pmin, pmax = config.get("pmin_d", 0.5), config.get("pmax_d", 40.0)
        durs = config.get("durations_h", [1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
        osf = config.get("oversample", 3.0)

        def run(yy, label):
            r1 = bls_one(t, yy, e, pmin, pmax, durs, osf)
            m1 = r1.pop("_mask")
            r2 = None
            mode = config.get("second_signal", "on_pass")
            want_s2 = (mode is True) or (mode == "on_pass" and label == "real"
                                         and passes_first_screen(r1, config)[0])
            if want_s2 and (~m1).sum() > 1000:
                r2 = bls_one(t[~m1], yy[~m1], e[~m1], pmin, pmax, durs, osf)
                r2.pop("_mask")
            return {"s1": r1, "s2": r2}

        out["real"] = run(y, "real")
        if config.get("run_inverted", True):
            out["inverted"] = run(2.0 - y, "inverted")
        out["ok"] = True
    except Exception as ex:  # noqa: BLE001
        out["error"] = f"{type(ex).__name__}: {ex}"
    finally:
        out["runtime_s"] = round(time.time() - t_start, 1)
        if not keep_fits:
            try:
                p = os.path.join(cache_dir or os.path.join(os.getcwd(), "data", "kbonus_batch"),
                                 os.path.basename(kbonus_url(fname)))
                if os.path.exists(p):
                    os.remove(p)
            except Exception:  # noqa: BLE001
                pass
    return out


def load_config(path: str) -> dict:
    with open(path) as fh:
        cfg = json.load(fh)
    cfg["_sha"] = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]
    return cfg


def passes_first_screen(sig: dict, cfg: dict) -> tuple[bool, list[str]]:
    """Apply the pre-registered first-screen thresholds to one signal dict."""
    th = cfg["first_screen"]
    fails = []
    if sig is None:
        return False, ["no signal"]
    if not (sig.get("sde") or 0) >= th["sde_min"]:
        fails.append("sde")
    snr = sig.get("snr_emp")
    if snr is None or not np.isfinite(snr) or snr < th["snr_emp_min"]:
        fails.append("snr_emp")
    if (sig.get("n_transits") or 0) < th["n_transits_min"]:
        fails.append("n_transits")
    if sig.get("max_single_event_frac") is not None and sig["max_single_event_frac"] > th["max_single_event_frac"]:
        fails.append("single_event")
    oe = sig.get("odd_even_sigma")
    if oe is not None and oe > th["odd_even_sigma_max"]:
        fails.append("odd_even")
    sec = sig.get("secondary_sigma")
    if sec is not None and sec > th["secondary_sigma_max"]:
        fails.append("secondary")
    if sig.get("depth_ppm", 0) > th["depth_ppm_max"] or sig.get("depth_ppm", 0) <= 0:
        fails.append("depth")
    if sig.get("period_d", 0) < th.get("period_min_d", 0):
        fails.append("period")
    return len(fails) == 0, fails
