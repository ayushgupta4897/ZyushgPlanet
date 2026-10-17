"""
Stage-2 vetting for first-screen survivors (Kepler Bonus light curves).

For each survivor this module re-derives everything from the light curve with the
signal MASKED during detrending (the first-pass biweight biases depths and creates
shoulders), then runs the checks that separate planets from the two dominant
false-positive classes seen in the screen: periodic variables / eclipsing binaries
(both-polarity, Lomb-Scargle amplitude, V-shape, secondary) and Kepler's monthly
downlink / quarter-edge systematics (events piling up next to data gaps).

    python -m pipeline.vet_batch --run-id poolA_main              # all survivors in survivors.csv
    python -m pipeline.vet_batch --run-id poolA_main --fnames kic-012345678 --top 10

Outputs: runs/<run-id>/vet/<fname>/{vet.json, fold.png, events.png, periodograms.png}
and runs/<run-id>/vet_summary.csv.  The pass/fail criteria live in the prereg file
under "stage2_vetting" and are reported, never silently applied.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.batch_search import (BKJD_OFFSET, _empirical_box_noise, _phase_mask, fetch_kbonus,  # noqa: E402
                                   load_config, read_kbonus)
from pipeline.ephemeris_match import match_signal  # noqa: E402


def _gap_edges(t, min_gap_d=0.5):
    """Times of the start/end of every gap longer than min_gap_d, plus the global ends."""
    dt = np.diff(t)
    idx = np.where(dt > min_gap_d)[0]
    edges = [t[0], t[-1]]
    for i in idx:
        edges += [t[i], t[i + 1]]
    return np.sort(np.array(edges))


def _normalise(t, f, fe, q, quarter):
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (q == 0) & (f > 0)
    t, f, fe, quarter = t[ok], f[ok], fe[ok], quarter[ok]
    o = np.argsort(t)
    t, f, fe, quarter = t[o], f[o], fe[o], quarter[o]
    y = np.empty_like(f)
    e = np.empty_like(f)
    for qq in np.unique(quarter):
        m = quarter == qq
        med = np.nanmedian(f[m])
        y[m] = f[m] / med
        e[m] = fe[m] / med
    return t, y, e, quarter


def _flatten(t, y, window, mask=None):
    from wotan import flatten
    kw = dict(method="biweight", window_length=window, break_tolerance=0.5, return_trend=True, edge_cutoff=0.0)
    if mask is not None:
        flat, trend = flatten(t, y, mask=mask, **kw)
    else:
        flat, trend = flatten(t, y, **kw)
    return flat, trend


def _box_depth(t, y, P, D, T0, width_factor=1.0):
    m = _phase_mask(t, P, D * width_factor, T0)
    if m.sum() < 5:
        return np.nan, np.nan, int(m.sum())
    return float(1.0 - np.mean(y[m])), float(np.std(y[m]) / np.sqrt(m.sum())), int(m.sum())


def _refine_bls(t, y, e, P0, D0_h, frac=0.005, n=4000):
    from astropy.timeseries import BoxLeastSquares
    periods = np.linspace(P0 * (1 - frac), P0 * (1 + frac), n)
    durs = np.array([0.5, 0.75, 1.0, 1.25, 1.5]) * D0_h / 24.0
    durs = durs[durs < 0.4 * P0]
    bls = BoxLeastSquares(t, y, e)
    r = bls.power(periods, durs, objective="snr")
    i = int(np.nanargmax(r.power))
    return float(r.period[i]), float(r.duration[i]), float(r.transit_time[i]), bls


def _events(t, y, P, D, T0, edges):
    """Per-event local depth, in-transit count, and distance to the nearest gap edge."""
    n0 = int(np.floor((t.min() - T0) / P)) - 1
    n1 = int(np.ceil((t.max() - T0) / P)) + 1
    rows = []
    for n in range(n0, n1 + 1):
        tc = T0 + n * P
        if tc < t.min() - D or tc > t.max() + D:
            continue
        inm = np.abs(t - tc) < 0.5 * D
        outm = (np.abs(t - tc) < 3.0 * D) & (np.abs(t - tc) > 0.75 * D)
        if inm.sum() < 2 or outm.sum() < 6:
            continue
        depth = float(np.mean(y[outm]) - np.mean(y[inm]))
        err = float(np.std(y[outm]) * np.sqrt(1.0 / inm.sum() + 1.0 / outm.sum()))
        dgap = float(np.min(np.abs(edges - tc)))
        rows.append({"n": n, "tc_bkjd": tc, "depth": depth, "err": err, "n_in": int(inm.sum()), "n_out": int(outm.sum()), "gap_dist_d": dgap})
    return pd.DataFrame(rows)


def _lomb_scargle(t, y):
    from astropy.timeseries import LombScargle
    f = np.linspace(1 / 60.0, 1 / 0.2, 60000)
    ls = LombScargle(t, y)
    p = ls.power(f)
    i = int(np.argmax(p))
    model = ls.model(t, f[i])
    amp = float((np.max(model) - np.min(model)) / 2.0)
    return float(1.0 / f[i]), amp, float(p[i])


def _tls(t, y, e, P, rstar, mstar, frac=0.02):
    from transitleastsquares import transitleastsquares
    model = transitleastsquares(t, y, e)
    r = model.power(period_min=P * (1 - frac), period_max=P * (1 + frac), R_star=rstar, M_star=mstar,
                    R_star_min=max(0.1, rstar * 0.7), R_star_max=rstar * 1.3, M_star_min=max(0.1, mstar * 0.7),
                    M_star_max=mstar * 1.3, oversampling_factor=5, duration_grid_step=1.05, use_threads=2,
                    show_progress_bar=False)
    out = {k: (float(getattr(r, k)) if np.isscalar(getattr(r, k)) else None) for k in
           ("period", "SDE", "snr", "depth", "duration", "rp_rs", "odd_even_mismatch", "transit_count", "distinct_transit_count")}
    out["rp_earth"] = float(r.rp_rs * rstar * 109.1) if np.isfinite(r.rp_rs) else None
    out["transit_depths_mean"] = float(np.nanmean(r.transit_depths)) if hasattr(r, "transit_depths") else None
    return out


def vet_one(fname: str, P0: float, D0_h: float, T0_bkjd: float, star: dict, cfg: dict, out_dir: str, make_plots=True) -> dict:
    t_start = time.time()
    os.makedirs(out_dir, exist_ok=True)
    res = {"fname": fname, "P_screen": P0, "dur_screen_h": D0_h, "t0_screen_bkjd": T0_bkjd}
    try:
        path = fetch_kbonus(fname, cache_dir="data/kbonus_batch", keep=True)
        t, f, fe, q, quarter, meta = read_kbonus(path, flux_col=cfg.get("flux_col", "FLUX"))
        res.update(meta)
        t, y, e, quarter = _normalise(t, f, fe, q, quarter)
        window = cfg.get("window_length_d", 1.0)
        # pass 1: blind detrend (as in the screen) and refine the ephemeris
        y1, _ = _flatten(t, y, window)
        good = np.isfinite(y1)
        t, y, e, quarter, y1 = t[good], y[good], e[good], quarter[good], y1[good]
        P, D, T0, bls = _refine_bls(t, y1, e, P0, D0_h)
        res.update({"P": P, "dur_h": D * 24, "t0_bkjd": T0, "t0_bjd": T0 + BKJD_OFFSET})
        # pass 2: re-detrend with the transit masked, then measure everything on y2
        mask = _phase_mask(t, P, 2.0 * D, T0)
        y2, trend2 = _flatten(t, y, window, mask=mask)
        g2 = np.isfinite(y2)
        t2, y2, e2, q2 = t[g2], y2[g2], e[g2], quarter[g2]
        d1, _, _ = _box_depth(t, y1, P, D, T0)
        d2, d2e, nin = _box_depth(t2, y2, P, D, T0)
        box = _empirical_box_noise(t2, y2, P, D, T0, n_boxes=80)
        res.update({"depth_pass1_ppm": d1 * 1e6, "depth_pass2_ppm": d2 * 1e6, "depth_pass2_formal_err_ppm": d2e * 1e6,
                    "box_noise_ppm": box * 1e6 if np.isfinite(box) else None, "snr_emp_pass2": (d2 / box) if np.isfinite(box) and box > 0 else None,
                    "n_in_transit": nin, "n_points": int(len(t2))})
        # odd / even / secondary on y2
        ph = ((t2 - T0) / P)
        cyc = np.floor(ph + 0.5).astype(int)
        inm = _phase_mask(t2, P, D, T0)
        odd = inm & (cyc % 2 == 1)
        even = inm & (cyc % 2 == 0)
        dodd = 1 - np.mean(y2[odd]) if odd.sum() > 3 else np.nan
        deven = 1 - np.mean(y2[even]) if even.sum() > 3 else np.nan
        res["depth_odd_ppm"], res["depth_even_ppm"] = dodd * 1e6, deven * 1e6
        res["odd_even_sigma"] = float(abs(dodd - deven) / (2 * box)) if np.isfinite(box) and box > 0 else None
        sec_m = _phase_mask(t2, P, D, T0 + 0.5 * P)
        dsec = 1 - np.mean(y2[sec_m]) if sec_m.sum() > 3 else np.nan
        res["secondary_depth_ppm"] = dsec * 1e6
        res["secondary_sigma"] = float(dsec / box) if np.isfinite(box) and box > 0 else None
        # V-shape proxy: depth in the central third vs full box
        core = _phase_mask(t2, P, D / 3.0, T0)
        dcore = 1 - np.mean(y2[core]) if core.sum() > 3 else np.nan
        res["core_to_box_depth_ratio"] = float(dcore / d2) if d2 and np.isfinite(dcore) else None
        # events and gap proximity
        edges = _gap_edges(t2)
        ev = _events(t2, y2, P, D, T0, edges)
        res["n_events"] = int(len(ev))
        if len(ev):
            res["gap_frac_1d"] = float((ev.gap_dist_d < 1.0).mean())
            res["gap_frac_2d"] = float((ev.gap_dist_d < 2.0).mean())
            res["frac_events_positive"] = float((ev.depth > 0).mean())
            w = 1 / ev.err**2
            res["event_weighted_depth_ppm"] = float((ev.depth * w).sum() / w.sum() * 1e6)
            res["event_chi2_const"] = float((((ev.depth - (ev.depth * w).sum() / w.sum()) ** 2) * w).sum())
            res["event_dof"] = int(len(ev) - 1)
            ev.to_csv(os.path.join(out_dir, "events.csv"), index=False)
            # per-year depth consistency (4 bins of time)
            bins = np.quantile(ev.tc_bkjd, [0, 0.25, 0.5, 0.75, 1.0])
            yd = []
            for a, b in zip(bins[:-1], bins[1:]):
                s = ev[(ev.tc_bkjd >= a) & (ev.tc_bkjd <= b)]
                ww = 1 / s.err**2
                yd.append([float((s.depth * ww).sum() / ww.sum() * 1e6), float(np.sqrt(1 / ww.sum()) * 1e6), int(len(s))])
            res["depth_by_quartile_ppm"] = yd
        # variability
        lsP, lsA, lsPow = _lomb_scargle(t, y1)
        res.update({"ls_period_d": lsP, "ls_half_amplitude_ppm": lsA * 1e6, "ls_power": lsPow,
                    "ls_amp_over_depth": (lsA / d2) if d2 else None,
                    "ls_period_matches": bool(any(abs(lsP - k * P) / (k * P) < 0.02 for k in (1, 2, 0.5)))})
        # duration plausibility: T14max for circular, central transit  ~ 13 h (P/yr)^(1/3) (rho/rho_sun)^(-1/3)
        rstar = float(star.get("Rstar", 0.8) or 0.8)
        mstar = rstar if rstar < 1.0 else 1.0
        rho = mstar / rstar**3
        t14max_h = 13.0 * (P / 365.25) ** (1 / 3) * rho ** (-1 / 3) * (1 + 0.05)
        res.update({"rstar_used": rstar, "t14_max_h": t14max_h, "duration_ratio": (D * 24) / t14max_h})
        res["rp_earth_box"] = float(np.sqrt(max(d2, 0)) * rstar * 109.1)
        # TLS in a narrow window (shape-aware SDE and radius)
        try:
            res["tls"] = _tls(t2, y2, e2, P, rstar, mstar)
        except Exception as ex:  # noqa: BLE001
            res["tls"] = {"error": str(ex)[:200]}
        # known objects
        hits = match_signal(meta["ra"], meta["dec"], P, tpforg_kic=meta.get("tpforg"), kic=star.get("kic"),
                            radius_arcmin=3.0, t0_bkjd=T0)
        res["known_matches"] = hits.to_dict("records") if len(hits) else []
        # plots
        if make_plots:
            _plots(out_dir, t2, y2, P, D, T0, ev, edges, fname, res)
        res["ok"] = True
    except Exception as ex:  # noqa: BLE001
        res["ok"] = False
        res["error"] = f"{type(ex).__name__}: {ex}"
    res["runtime_s"] = round(time.time() - t_start, 1)
    ok2, fails2 = stage2_verdict(res, cfg, both_polarity=bool(star.get("_both_polarity", False)))
    res["stage2_pass"], res["stage2_fails"] = ok2, fails2
    with open(os.path.join(out_dir, "vet.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=_jd)
    return res


def stage2_verdict(r: dict, cfg: dict, both_polarity: bool = False) -> tuple[bool, list[str]]:
    """Apply the registered stage-2 criteria; returns (passes, list_of_failed_criteria)."""
    c = cfg.get("stage2_vetting", {})
    fails = []
    if not r.get("ok"):
        return False, ["error"]
    def g(k):
        v = r.get(k)
        return None if v is None or (isinstance(v, float) and not np.isfinite(v)) else v
    if g("gap_frac_1d") is not None and g("gap_frac_1d") > c.get("gap_frac_1d_max", 1):
        fails.append("gap_frac")
    if g("frac_events_positive") is not None and g("frac_events_positive") < c.get("frac_events_positive_min", 0):
        fails.append("events_positive")
    d1, d2 = g("depth_pass1_ppm"), g("depth_pass2_ppm")
    if d1 and d2 and d1 > 0:
        lo, hi = c.get("depth_pass2_over_pass1_range", [0, 99])
        if not (lo <= d2 / d1 <= hi):
            fails.append("depth_pass_ratio")
    if g("event_chi2_const") is not None and g("event_dof"):
        if g("event_chi2_const") / max(1, g("event_dof")) > c.get("event_reduced_chi2_max", 99):
            fails.append("event_chi2")
    if g("odd_even_sigma") is not None and g("odd_even_sigma") > c.get("odd_even_sigma_max", 99):
        fails.append("odd_even")
    if g("secondary_sigma") is not None and g("secondary_sigma") > c.get("secondary_sigma_max", 99):
        fails.append("secondary")
    if g("core_to_box_depth_ratio") is not None and g("core_to_box_depth_ratio") < c.get("core_to_box_depth_ratio_min", 0):
        fails.append("v_shape")
    if g("ls_amp_over_depth") is not None and g("ls_amp_over_depth") > c.get("ls_amp_over_depth_max", 99):
        fails.append("variability_amp")
    if c.get("ls_period_matches_forbidden") and r.get("ls_period_matches"):
        fails.append("variability_period")
    if c.get("both_polarity_forbidden") and both_polarity:
        fails.append("both_polarity")
    if g("duration_ratio") is not None and g("duration_ratio") > c.get("duration_ratio_max", 99):
        fails.append("duration")
    tls = r.get("tls") or {}
    if isinstance(tls, dict) and "SDE" in tls:
        if (tls.get("SDE") or 0) < c.get("tls_sde_min", 0):
            fails.append("tls_sde")
        if (tls.get("snr") or 0) < c.get("tls_snr_min", 0):
            fails.append("tls_snr")
    if c.get("known_match_same_star_or_tpf_forbidden"):
        if any(h.get("same_kic") or h.get("same_tpf") for h in r.get("known_matches", [])):
            fails.append("known_object")
    if g("rp_earth_box") is not None and g("rp_earth_box") > c.get("rp_earth_max_for_candidate", 99):
        fails.append("radius")
    return len(fails) == 0, fails


def _jd(o):
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _plots(out_dir, t, y, P, D, T0, ev, edges, fname, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ph = ((t - T0) / P + 0.5) % 1.0 - 0.5
    fig, ax = plt.subplots(1, 3, figsize=(16, 4))
    w = 4 * D / P
    m = np.abs(ph) < w
    ax[0].plot(ph[m] * P * 24, (y[m] - 1) * 1e6, ".", ms=2, alpha=0.3, color="0.6")
    bins = np.linspace(-w, w, 60)
    idx = np.digitize(ph[m], bins)
    bx = [(ph[m][idx == i]).mean() * P * 24 for i in range(1, len(bins)) if (idx == i).sum() > 3]
    by = [((y[m][idx == i]).mean() - 1) * 1e6 for i in range(1, len(bins)) if (idx == i).sum() > 3]
    ax[0].plot(bx, by, "o-", color="C1", ms=4)
    ax[0].axvspan(-D * 12, D * 12, color="C0", alpha=0.1)
    ax[0].set_xlabel("hours from mid-transit")
    ax[0].set_ylabel("ppm")
    ax[0].set_title(f"{fname}  P={P:.4f} d  depth={res.get('depth_pass2_ppm', 0):.0f} ppm (masked detrend)")
    ms = np.abs(((t - T0 - 0.5 * P) / P + 0.5) % 1.0 - 0.5) < w
    ph2 = ((t - T0 - 0.5 * P) / P + 0.5) % 1.0 - 0.5
    ax[1].plot(ph2[ms] * P * 24, (y[ms] - 1) * 1e6, ".", ms=2, alpha=0.3, color="0.6")
    ax[1].set_title(f"phase 0.5 (secondary): {res.get('secondary_depth_ppm', 0):.0f} ppm")
    ax[1].set_xlabel("hours from phase 0.5")
    if len(ev):
        ax[2].errorbar(ev.tc_bkjd, ev.depth * 1e6, ev.err * 1e6, fmt="o", ms=3, color="C0")
        for g in edges:
            ax[2].axvline(g, color="r", lw=0.4, alpha=0.5)
        ax[2].axhline(res.get("depth_pass2_ppm", 0), color="C1")
        ax[2].set_title(f"per-event depth (red = gap edges); gap_frac<1d={res.get('gap_frac_1d', 0):.2f}")
        ax[2].set_xlabel("BKJD")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fold_events.png"), dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--config", default="pipeline/prereg_poolA.json")
    ap.add_argument("--targets", default="targets/poolA_kbonus_bg_dwarfs_ranked.csv")
    ap.add_argument("--fnames", nargs="*")
    ap.add_argument("--top", type=int, default=None)
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args()
    cfg = load_config(a.config)
    run_dir = os.path.join("runs", a.run_id)
    sv = pd.read_csv(os.path.join(run_dir, "survivors.csv"))
    if a.fnames:
        sv = sv[sv.fname.isin(a.fnames)]
    if a.top:
        sv = sv.head(a.top)
    tg = pd.read_csv(a.targets).set_index("fname")
    rows = []
    for _, s in sv.iterrows():
        star = tg.loc[s.fname].to_dict() if s.fname in tg.index else {}
        star["_both_polarity"] = bool(s.get("both_polarity", False))
        out_dir = os.path.join(run_dir, "vet", s.fname)
        print(f"vetting {s.fname} P={s.real_period_d:.4f} ...", flush=True)
        r = vet_one(s.fname, float(s.real_period_d), float(s.real_duration_h), float(s.real_t0_bjd) - BKJD_OFFSET, star, cfg, out_dir, make_plots=not a.no_plots)
        flat = {k: v for k, v in r.items() if not isinstance(v, (list, dict))}
        if isinstance(r.get("tls"), dict):
            flat.update({f"tls_{k}": v for k, v in r["tls"].items()})
        flat["n_known_matches"] = len(r.get("known_matches", []))
        flat["stage2_fails"] = ",".join(r.get("stage2_fails", []))
        rows.append(flat)
        print(f"  -> depth2={r.get('depth_pass2_ppm')} snr2={r.get('snr_emp_pass2')} gap1d={r.get('gap_frac_1d')} oe={r.get('odd_even_sigma')} sec={r.get('secondary_sigma')} ls_amp/depth={r.get('ls_amp_over_depth')} durRatio={r.get('duration_ratio')} tlsSDE={r.get('tls', {}).get('SDE') if isinstance(r.get('tls'), dict) else None} rt={r.get('runtime_s')}s err={r.get('error')} | STAGE2 {'PASS' if r.get('stage2_pass') else 'fail: ' + ','.join(r.get('stage2_fails', []))}", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(run_dir, "vet_summary.csv"), index=False)


if __name__ == "__main__":
    main()
