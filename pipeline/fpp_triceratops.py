"""False-positive probability for one stage-2 survivor with TRICERATOPS (Giacalone et al. 2021).

    python -m pipeline.fpp_triceratops --run-id poolA_main --fname kic-006699887 --seeds 4 --N 1000000

What it does (every number lands in runs/<run-id>/fpp/<fname>/):
  1. Reads the stage-2 ephemeris from runs/<run-id>/vet/<fname>/vet.json (P, duration, t0, on-target depth).
  2. Rebuilds the stage-2 light curve (quality == 0, per-quarter normalise, biweight re-detrend with the
     transit masked), folds it, keeps +-3 durations around mid-transit and bins it to half a long cadence.
  3. Builds a TRICERATOPS target in Kepler mode: the TIC field around the target's Gaia position, the
     originating Kepler TPF (the target itself was never a Kepler target, so the host's TPF supplies the
     WCS and image), and a 3x3-pixel aperture centred on the target as the stand-in for the KBONUS PSF fit.
  4. Runs calc_probs with `seeds` different random seeds and reports FPP and NFPP (mean, std, per seed) with
     the per-scenario probabilities, plus the star table TRICERATOPS used.

Approximations that travel with the result (also written to fpp_summary.json):
  - TRICERATOPS models every PSF as a circular Gaussian of 0.75 px; the Kepler PRF is not used.
  - KBONUS FLUX is PSF photometry, not an aperture sum. A 3x3 aperture on the target is used only to set
    the flux ratio of each nearby star, i.e. the depth each would need to host; the light curve itself is
    treated as dilution-corrected (on-target depth).
  - Stellar parameters come from the TIC (Stassun et al. 2019), not from our colour-based estimate; both are
    recorded. No contrast curve and no MOLUSC companion file exist for these G 16-17 stars, so the unresolved
    bound-companion scenarios carry their prior weight.
  - "Validated" needs FPP < 0.015 and NFPP < 0.001 in a refereed paper (docs/05); this script only produces
    the numbers.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

from pipeline.batch_search import fetch_kbonus, read_kbonus, _phase_mask, BKJD_OFFSET
from pipeline.vet_batch import _normalise, _flatten

CADENCE_D = 29.4244 / 1440.0


def build_lc(fname: str, cfg: dict, P: float, D: float, T0: float, window_durations: float = 3.0, bin_d: float = CADENCE_D / 2):
    path = fetch_kbonus(fname, cache_dir="data/kbonus_batch", keep=True)
    t, f, fe, q, quarter, meta = read_kbonus(path, flux_col=cfg.get("flux_col", "FLUX"))
    t, y, e, quarter = _normalise(t, f, fe, q, quarter)
    mask = _phase_mask(t, P, 2.0 * D, T0)
    y2, _ = _flatten(t, y, cfg.get("window_length_d", 1.0), mask=mask)
    g = np.isfinite(y2)
    t, y2 = t[g], y2[g]
    dt = ((t - T0 + 0.5 * P) % P) - 0.5 * P
    w = np.abs(dt) < window_durations * D
    dt, yy = dt[w], y2[w]
    edges = np.arange(-window_durations * D, window_durations * D + bin_d, bin_d)
    idx = np.digitize(dt, edges) - 1
    rows = []
    for i in range(len(edges) - 1):
        m = idx == i
        if m.sum() < 3:
            continue
        rows.append((0.5 * (edges[i] + edges[i + 1]), float(np.mean(yy[m])), float(np.std(yy[m]) / np.sqrt(m.sum())), int(m.sum())))
    b = pd.DataFrame(rows, columns=["time", "flux", "flux_err", "n"])
    oot = np.abs(b.time) > 0.75 * D
    err = float(np.std(b.flux[oot])) if oot.sum() > 5 else float(np.median(b.flux_err))
    return b, err, {"n_unbinned": int(w.sum()), "n_total_cadences": int(len(t)), "bin_d": bin_d, "window_durations": window_durations,
                    "flux_err_used": err, "flux_err_source": "std of binned out-of-transit points (|dt| > 0.75 D)"}


def neighbour_depth(fname: str, cfg: dict, P: float, D: float, T0: float):
    """Box depth (fraction) and empirical error of another KBONUS source at the target's ephemeris; None if no light curve."""
    from pipeline.batch_search import _empirical_box_noise
    try:
        path = fetch_kbonus(fname, cache_dir="data/kbonus_batch", keep=True)
    except Exception as ex:  # noqa: BLE001  (404: not a KBONUS source)
        return None, None, f"no KBONUS light curve ({type(ex).__name__})"
    t, f, fe, q, quarter, _ = read_kbonus(path, flux_col=cfg.get("flux_col", "FLUX"))
    t, y, e, quarter = _normalise(t, f, fe, q, quarter)
    mask = _phase_mask(t, P, 2.0 * D, T0)
    y2, _ = _flatten(t, y, cfg.get("window_length_d", 1.0), mask=mask)
    g = np.isfinite(y2)
    t, y2 = t[g], y2[g]
    inm = _phase_mask(t, P, D, T0)
    if inm.sum() < 10:
        return None, None, "fewer than 10 in-transit cadences"
    depth = float(1.0 - np.mean(y2[inm]))
    err = float(_empirical_box_noise(t, y2, P, D, T0, n_boxes=80))
    return depth, err, f"{int(inm.sum())} in-transit cadences"


def tic_identifiers(tic_ids):
    """KIC and Gaia DR3 ids for a list of TIC ids (one MAST query)."""
    from astroquery.mast import Catalogs
    out = {}
    try:
        tab = Catalogs.query_criteria(catalog="TIC", ID=[int(i) for i in tic_ids])
        for r in tab:
            kic = r["KIC"]; gaia = r["GAIA"]
            out[int(r["ID"])] = {"KIC": None if (kic is None or str(kic) in ("--", "", "nan")) else int(kic),
                                 "GAIA": None if (gaia is None or str(gaia) in ("--", "", "nan")) else int(gaia)}
    except Exception as ex:  # noqa: BLE001
        print("TIC identifier lookup failed:", ex)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--fname", required=True)
    ap.add_argument("--config", default="pipeline/prereg_poolA.json")
    ap.add_argument("--targets", default="targets/poolA_kbonus_bg_dwarfs_ranked.csv")
    ap.add_argument("--quarters", type=int, nargs="+", default=[4, 5, 6, 7], help="one Kepler quarter per season by default")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--seed-start", type=int, default=0, help="first seed index; existing runs in fpp_summary.json are kept and appended to")
    ap.add_argument("--N", type=int, default=1_000_000)
    ap.add_argument("--parallel", action="store_true")
    ap.add_argument("--aperture-half-width", type=int, default=1, help="1 -> 3x3 pixels around the target")
    ap.add_argument("--out-root", default=None)
    ap.add_argument("--exclude-sigma", type=float, default=5.0,
                    help="drop a nearby star when its own KBONUS light curve excludes the eclipse it would need at this many sigma")
    ap.add_argument("--exclude-tic", type=int, nargs="*", default=[], help="manual exclusions (give --exclude-reason)")
    ap.add_argument("--exclude-reason", default="")
    a = ap.parse_args()

    warnings.filterwarnings("ignore")
    cfg = json.load(open(a.config))
    vet = json.load(open(f"runs/{a.run_id}/vet/{a.fname}/vet.json"))
    tg = pd.read_csv(a.targets)
    row = tg[tg.fname == a.fname].iloc[0]
    out = a.out_root or f"runs/{a.run_id}/fpp/{a.fname}"
    os.makedirs(out, exist_ok=True)
    log = open(os.path.join(out, "fpp.log"), "a")

    def say(*x):
        s = " ".join(str(i) for i in x)
        print(s, flush=True)
        log.write(s + "\n"); log.flush()

    P, D, T0 = float(vet["P"]), float(vet["dur_h"]) / 24.0, float(vet["t0_bkjd"])
    depth = float(vet["depth_pass2_ppm"]) * 1e-6
    say(f"[{_dt.datetime.utcnow().isoformat()}Z] {a.fname}: P={P:.6f} d, D={D*24:.3f} h, t0_bkjd={T0:.5f}, depth={depth*1e6:.0f} ppm (stage 2, {vet.get('stage2_rules')})")

    lc, err, lcmeta = build_lc(a.fname, cfg, P, D, T0)
    lc.to_csv(os.path.join(out, "lc_binned.csv"), index=False)
    say(f"light curve: {len(lc)} binned points from {lcmeta['n_unbinned']} cadences; flux_err {err*1e6:.0f} ppm")

    from triceratops import triceratops as tr
    t_build = time.time()
    tgt = tr.target(ID=int(row.tpforg), sectors=np.array(a.quarters), mission="Kepler", ra=float(row.ra), dec=float(row.dec),
                    lightkurve_cache_dir=os.path.abspath("data/lightkurve_cache"), search_radius=10)
    stars = tgt.stars.copy()
    stars.to_csv(os.path.join(out, "stars_tic.csv"), index=False)
    s0 = stars.iloc[0]
    say(f"TRICERATOPS target row 0: TIC {int(s0.ID)} Tmag {s0.Tmag:.2f} Teff {s0.Teff} R {s0.rad} M {s0['mass']} plx {s0.plx}; "
        f"{len(stars)} TIC stars within 10 px; build {time.time()-t_build:.0f}s")
    sep0 = 3600 * np.hypot((float(s0.ra) - float(row.ra)) * np.cos(np.radians(float(row.dec))), float(s0.dec) - float(row.dec))
    say(f"row-0 offset from Gaia position: {sep0:.2f} arcsec (must be < 1)")
    if sep0 > 1.0:
        say("ABORT: TIC row 0 is not the target"); sys.exit(2)

    # 3x3 aperture on the target in every quarter (stand-in for the KBONUS PSF fit)
    hw = a.aperture_half_width
    aps = []
    for k in range(len(tgt.pix_coords)):
        cx, cy = np.round(tgt.pix_coords[k][0])
        xs = np.repeat(np.arange(cx - hw, cx + hw + 1), 2 * hw + 1)
        ys = np.tile(np.arange(cy - hw, cy + hw + 1), 2 * hw + 1)
        aps.append(np.array([xs, ys]).T)
    tgt.calc_depths(tdepth=depth, all_ap_pixels=aps, dilution_corrected=True)
    tgt.stars.to_csv(os.path.join(out, "stars_tic.csv"), index=False)
    need = tgt.stars[tgt.stars.tdepth > 0][["ID", "Tmag", "fluxratio", "tdepth"]]
    say("stars that could host the dip (tdepth = eclipse depth they would need):")
    for _, r in need.iterrows():
        say(f"   TIC {int(r.ID)} Tmag {r.Tmag:.2f} fluxratio {r.fluxratio:.4f} tdepth {r.tdepth:.4f}")
    # data-driven exclusions: fold each possible host's own KBONUS light curve at the ephemeris
    exclusions = []
    others = need[need.ID != s0.ID]
    ids = tic_identifiers(others.ID.values) if len(others) else {}
    for _, r in others.iterrows():
        tid = int(r.ID)
        info = ids.get(tid, {})
        cand = []
        if info.get("KIC"):
            cand.append(f"kic-{info['KIC']:09d}")
        if info.get("GAIA"):
            cand.append(f"gaia-dr3-{info['GAIA']}")
        rec = {"TIC": tid, "Tmag": float(r.Tmag), "fluxratio": float(r.fluxratio), "required_depth": float(r.tdepth), "ids": info,
               "lc_tried": cand, "excluded": False}
        for fn in cand:
            d, e, note = neighbour_depth(fn, cfg, P, D, T0)
            rec["lc_used"], rec["note"] = fn, note
            if d is None:
                continue
            rec.update({"measured_depth_ppm": d * 1e6, "measured_err_ppm": e * 1e6,
                        "exclusion_sigma": (float(r.tdepth) - d) / e if e and e > 0 else None})
            break
        if tid in a.exclude_tic:
            rec.update({"excluded": True, "why": "manual: " + a.exclude_reason})
        elif rec.get("exclusion_sigma") is not None and rec["exclusion_sigma"] > a.exclude_sigma:
            rec.update({"excluded": True, "why": f"own light curve at the ephemeris: {rec['measured_depth_ppm']:.0f} +- {rec['measured_err_ppm']:.0f} ppm "
                                                 f"vs {rec['required_depth']*1e6:.0f} ppm required ({rec['exclusion_sigma']:.0f} sigma)"})
        exclusions.append(rec)
        say(f"   neighbour TIC {tid}: {json.dumps({k: v for k, v in rec.items() if k not in ('ids', 'lc_tried')}, default=str)}")
    for rec in exclusions:
        if rec["excluded"]:
            # zero the required depth instead of remove_star(): pix_coords are indexed by row position and
            # calc_probs only considers stars with tdepth > 0, so this drops the scenario without misaligning anything
            tgt.stars.loc[tgt.stars.ID.astype(str) == str(rec["TIC"]), "tdepth"] = 0.0
            say(f"   dropped TIC {rec['TIC']} from the scenario list: {rec['why']}")
    tgt.stars.to_csv(os.path.join(out, "stars_tic.csv"), index=False)
    need = tgt.stars[tgt.stars.tdepth > 0][["ID", "Tmag", "fluxratio", "tdepth"]]
    say("scenario hosts after exclusions: " + ", ".join(f"TIC {int(x.ID)} ({x.tdepth:.4f})" for _, x in need.iterrows()))
    try:
        tgt.plot_field(sector=a.quarters[0], ap_pixels=aps[0], save=True, fname=os.path.join(out, "field_q%d" % a.quarters[0]))
    except Exception as ex:  # noqa: BLE001
        say(f"plot_field failed: {type(ex).__name__}: {ex}")

    runs = []
    prev = os.path.join(out, "fpp_summary.json")
    if a.seed_start > 0 and os.path.exists(prev):
        runs = json.load(open(prev)).get("runs", [])
        say(f"appending to {len(runs)} existing seed runs")
    for s in range(a.seed_start, a.seed_start + a.seeds):
        np.random.seed(1000 + s)
        t1 = time.time()
        tgt.calc_probs(time=lc.time.values, flux_0=lc.flux.values, flux_err_0=err, P_orb=P, N=a.N, parallel=a.parallel,
                       verbose=0, exptime=CADENCE_D, nsamples=20, dilution_corrected=True)
        probs = tgt.probs.copy()
        probs.to_csv(os.path.join(out, f"probs_seed{s}.csv"), index=False)
        runs.append({"seed": 1000 + s, "FPP": float(tgt.FPP), "NFPP": float(tgt.NFPP), "runtime_s": round(time.time() - t1, 1),
                     "scenario_probs": {str(r.scenario) + ":" + str(int(r.ID)): float(r.prob) for _, r in probs.iterrows()}})
        say(f"seed {1000+s}: FPP {tgt.FPP:.4f}  NFPP {tgt.NFPP:.5f}  ({time.time()-t1:.0f}s)")
        if s == 0 and not os.path.exists(os.path.join(out, "fits_seed0.pdf")):
            try:
                tgt.plot_fits(time=lc.time.values, flux_0=lc.flux.values, flux_err_0=err, save=True, fname=os.path.join(out, "fits_seed0"))
            except Exception as ex:  # noqa: BLE001
                say(f"plot_fits failed: {type(ex).__name__}: {ex}")
        summary = {
            "fname": a.fname, "run_id": a.run_id, "utc": _dt.datetime.utcnow().isoformat() + "Z",
            "triceratops_version": getattr(__import__("triceratops"), "__version__", "1.1.0"),
            "ephemeris": {"P_d": P, "dur_h": D * 24, "t0_bkjd": T0, "t0_bjd": T0 + BKJD_OFFSET, "depth_ppm_on_target": depth * 1e6,
                          "source": f"runs/{a.run_id}/vet/{a.fname}/vet.json ({vet.get('stage2_rules')})"},
            "light_curve": lcmeta | {"n_binned": int(len(lc))},
            "tic_target": {k: (None if pd.isna(v) else (float(v) if isinstance(v, (int, float, np.floating, np.integer)) else str(v)))
                           for k, v in s0.items()},
            "our_host_estimate": {"Teff": float(row.Teff), "Rstar": float(row.Rstar), "dist_pc": float(row.dist_pc), "RUWE": float(row.RUWE), "sptype": str(row.sptype)},
            "quarters_used_for_field": a.quarters, "aperture": f"{2*hw+1}x{2*hw+1} px on target (stand-in for PSF photometry)",
            "N_draws": a.N, "seeds_done": len(runs),
            "FPP_mean": float(np.mean([r["FPP"] for r in runs])), "FPP_std": float(np.std([r["FPP"] for r in runs])),
            "NFPP_mean": float(np.mean([r["NFPP"] for r in runs])), "NFPP_std": float(np.std([r["NFPP"] for r in runs])),
            "runs": runs,
            "possible_hosts": need.to_dict(orient="records"),
            "neighbour_tests": exclusions, "exclude_sigma": a.exclude_sigma,
            "approximations": [
                "Gaussian 0.75-px PSF for every star (TRICERATOPS default), not the Kepler PRF",
                "KBONUS PSF flux treated as dilution-corrected; a 3x3 aperture on the target sets nearby-star flux ratios",
                "stellar parameters from the TIC; no contrast curve; no MOLUSC file",
                "folded light curve binned to half a long cadence; model supersampled 20x over 29.4 min",
            ],
            "thresholds": {"validated_FPP": 0.015, "validated_NFPP": 0.001, "note": "thresholds apply only in a refereed paper; this is a pipeline product"},
        }
        json.dump(summary, open(os.path.join(out, "fpp_summary.json"), "w"), indent=1)
    say(f"DONE {a.fname}: FPP {summary['FPP_mean']:.4f} +- {summary['FPP_std']:.4f}, NFPP {summary['NFPP_mean']:.5f} +- {summary['NFPP_std']:.5f} over {len(runs)} seeds")


if __name__ == "__main__":
    main()
