"""
Pixel-level localisation for a Kepler Bonus signal: which star actually dims?

The KBONUS light curve was extracted from the target-pixel file (TPF) of a *different*
star (tpforg).  For every quarter we build a per-pixel "missing light" map by local
in/out-of-transit differencing around each event, a per-pixel noise map (random-phase
ephemerides, or an analytic fallback when a quarter has too few events), a
significance-weighted deficit centroid with a bootstrap error over events, and a
model comparison (point source at the TARGET's Gaia position vs a scaled copy of the
stamp, i.e. a dip on whatever dominates the aperture).

Registered rules (pipeline/prereg_poolA.json -> pixel_localisation, post-hoc, v1.5):
  * a quarter is INCONCLUSIVE when the target lies within 0.5 px of the stamp edge or
    outside it, when < 3 events are usable, or when < 2 pixels exceed SNR 2;
  * verdict 'on-target' needs >= 5 conclusive quarters in which, for >= 70 %, the
    target-to-deficit offset is below 2 sigma AND below the offset to every other Gaia
    source within 2 px of the stamp; 'offset' needs the opposite in >= 70 %;
    otherwise 'inconclusive'.

    python -m pipeline.pixel_vet --fname kic-005689121 --period 13.769448 --t0-bjd 2455630.2454 --dur-h 2.3 \
        --ra 285.907827 --dec 40.941988 --tpforg 5689116 --out runs/controls/pixel/kic-005689121

Positive control (this module, 2026-10-17): KIC 5689121 -> deficit on the target in
13/13 quarters (an independent agent's own implementation found 12/13).
"""
from __future__ import annotations

import argparse
import json
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
BKJD_OFFSET = 2454833.0
PX_ARCSEC = 3.98


def _gaia_neighbours(ra, dec, radius_arcsec=30.0):
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "BP-RP", "Plx", "RUWE"], row_limit=300)
    r = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=radius_arcsec * u.arcsec, catalog="I/355/gaiadr3")
    if not r:
        return []
    return [{"source_id": str(x["Source"]), "ra": float(x["RA_ICRS"]), "dec": float(x["DE_ICRS"]),
             "G": float(x["Gmag"]) if np.isfinite(x["Gmag"]) else None} for x in r[0]]


def _wcentroid(wmap):
    tot = wmap.sum()
    if tot <= 0:
        return np.nan, np.nan
    yy, xx = np.indices(wmap.shape)
    return float((xx * wmap).sum() / tot), float((yy * wmap).sum() / tot)


def _load_rules(cfg_path):
    try:
        c = json.load(open(cfg_path)).get("pixel_localisation", {})
    except Exception:  # noqa: BLE001
        c = {}
    return {"edge_px": c.get("edge_px", 0.5), "min_events": c.get("min_events_per_quarter", 3), "min_sig": c.get("min_sig_pixels", 2),
            "min_q": c.get("min_conclusive_quarters", 5), "frac": c.get("on_target_fraction", 0.7),
            "nb_px": c.get("neighbour_search_px_outside_stamp", 2.0), "nboot": c.get("bootstrap_n", 200)}


def _load_rules_v16():
    try:
        c = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline", "prereg_poolA.json")))["pixel_localisation"]["v16_chi2_rule"]
    except Exception:  # noqa: BLE001
        c = {}
    return {"min_q": c.get("min_quarters", 5), "d_null": c.get("delta_chi2_vs_null_min", 25.0), "d_alt": c.get("delta_chi2_vs_best_alternative_min", 25.0)}


def pixel_vet(fname, period, t0_bjd, dur_h, ra, dec, tpforg, out_dir, make_plot=True, max_quarters=None,
              force_analytic=False, cfg_path=None, invert=False):
    sgn = -1.0 if invert else 1.0
    import lightkurve as lk
    rules = _load_rules(cfg_path or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline", "prereg_poolA.json"))
    os.makedirs(out_dir, exist_ok=True)
    res = {"fname": fname, "period": period, "t0_bjd": t0_bjd, "dur_h": dur_h, "tpforg": int(tpforg), "rules": rules,
           "force_analytic": bool(force_analytic), "quarters": []}
    sr = lk.search_targetpixelfile(f"KIC {int(tpforg)}", mission="Kepler", cadence="long")
    res["n_tpf_found"] = len(sr)
    neighbours = _gaia_neighbours(ra, dec)
    res["gaia_neighbours"] = neighbours
    t0_bkjd = t0_bjd - BKJD_OFFSET
    D = dur_h / 24.0
    rows, imgs = [], []
    rng = np.random.default_rng(1)
    for i in range(len(sr)):
        if max_quarters and i >= max_quarters:
            break
        try:
            tpf = sr[i].download(quality_bitmask="default")
        except Exception as ex:  # noqa: BLE001
            rows.append({"quarter": None, "error": str(ex)[:120], "conclusive": False})
            continue
        q = int(tpf.quarter) if hasattr(tpf, "quarter") else i
        t = tpf.time.value
        cube = np.array(tpf.flux.value, float)
        ny, nx = cube.shape[1:]
        valid_px = np.mean(np.isfinite(cube), axis=0) > 0.9
        cube[:, ~valid_px] = np.nan
        okc = np.isfinite(t) & (np.sum(np.isfinite(cube.reshape(len(t), -1)), axis=1) >= valid_px.sum() * 0.9)
        t, cube = t[okc], cube[okc]
        if len(t) < 200 or valid_px.sum() < 4:
            rows.append({"quarter": q, "error": f"too few cadences ({len(t)}) or pixels ({int(valid_px.sum())})", "conclusive": False})
            continue
        img_out = np.nanmean(cube, axis=0)
        w = tpf.wcs
        px_t = w.all_world2pix([[ra, dec]], 0)[0]
        e = rules["edge_px"]
        edge_flag = not (e <= px_t[0] <= nx - 1 - e and e <= px_t[1] <= ny - 1 - e)

        def event_diff_map(tc_list):
            acc = np.zeros((ny, nx)); nev = 0
            for tc in tc_list:
                inm = np.abs(t - tc) < 0.5 * D
                outm = (np.abs(t - tc) > 0.75 * D) & (np.abs(t - tc) < 3.0 * D)
                if inm.sum() < 2 or outm.sum() < 6:
                    continue
                # out - in = light missing during the event; invert=True mirrors it for the inverted-flux null
                acc += sgn * (np.nanmean(cube[outm], axis=0) - np.nanmean(cube[inm], axis=0))
                nev += 1
            return (acc / nev if nev else acc * np.nan), nev

        n0 = int(np.floor((t.min() - t0_bkjd) / period)) - 1
        n1 = int(np.ceil((t.max() - t0_bkjd) / period)) + 1
        tcs = [t0_bkjd + n * period for n in range(n0, n1 + 1) if t.min() - D < t0_bkjd + n * period < t.max() + D]
        usable = [tc for tc in tcs if (np.abs(t - tc) < 0.5 * D).sum() >= 2 and ((np.abs(t - tc) > 0.75 * D) & (np.abs(t - tc) < 3.0 * D)).sum() >= 6]
        diff, nev = event_diff_map(usable)
        if nev < 1:
            rows.append({"quarter": q, "error": "no usable events", "conclusive": False})
            continue
        rand_maps = []
        if not force_analytic:
            for _ in range(24):
                shift = rng.uniform(0.15, 0.85) * period
                m, ne = event_diff_map([tc + shift for tc in usable])
                if ne >= 2:
                    rand_maps.append(m)
        if len(rand_maps) >= 8:
            sig_map = np.nanstd(np.array(rand_maps), axis=0); noise_method = "random_phase"
        else:
            cad_std = np.nanstd(cube, axis=0)
            fac = [1.0 / (np.abs(t - tc) < 0.5 * D).sum() + 1.0 / ((np.abs(t - tc) > 0.75 * D) & (np.abs(t - tc) < 3.0 * D)).sum() for tc in usable]
            sig_map = cad_std * np.sqrt(np.mean(fac)) / np.sqrt(len(fac)); noise_method = "analytic"
        sig_map = np.where(valid_px & np.isfinite(sig_map) & (sig_map > 0), sig_map, np.nan)
        snr_map = diff / sig_map
        if not np.isfinite(snr_map).any():
            rows.append({"quarter": q, "error": "no usable noise map", "conclusive": False})
            continue
        sel = (snr_map > 2.0) & valid_px
        n_sig = int(sel.sum())
        wmap = np.where(sel, np.clip(diff, 0, None), 0.0)
        cx, cy = _wcentroid(wmap)
        bxs, bys = [], []
        if nev >= 2 and n_sig >= 1:
            for _ in range(rules["nboot"]):
                samp = list(rng.choice(usable, size=len(usable), replace=True))
                m, ne = event_diff_map(samp)
                if ne < 1:
                    continue
                x, y = _wcentroid(np.where(sel, np.clip(m, 0, None), 0.0))
                if np.isfinite(x):
                    bxs.append(x); bys.append(y)
        sig_cent_px = float(np.hypot(np.std(bxs), np.std(bys))) if len(bxs) >= 20 else np.nan
        yy, xx = np.indices((ny, nx))
        fx, fy = _wcentroid(np.where(valid_px, img_out, 0))
        by, bx = np.unravel_index(np.nanargmax(np.where(valid_px, snr_map, -np.inf)), (ny, nx))
        ok = valid_px & np.isfinite(sig_map) & (sig_map > 0)

        def chi2_fit(model):
            mm = model[ok]; dd = diff[ok]; ss = sig_map[ok]
            den = np.sum(mm * mm / ss**2)
            if den <= 1e-12:
                return float(np.sum((dd / ss) ** 2)), 0.0
            a = np.sum(mm * dd / ss**2) / den
            return float(np.sum(((dd - a * mm) / ss) ** 2)), float(a)

        chi2_null = float(np.sum((diff[ok] / sig_map[ok]) ** 2))
        chi2A, aA = chi2_fit(np.where(valid_px, img_out, 0.0))
        chi2B, aB = chi2_fit(np.exp(-0.5 * ((xx - px_t[0]) ** 2 + (yy - px_t[1]) ** 2) / 0.9**2))

        def dist(px, py, x, y):
            return float(np.hypot(px - x, py - y) * PX_ARCSEC) if np.isfinite(x) else None

        nb = []
        for n in neighbours:
            px_n = w.all_world2pix([[n["ra"], n["dec"]]], 0)[0]
            near = (-rules["nb_px"] <= px_n[0] < nx + rules["nb_px"]) and (-rules["nb_px"] <= px_n[1] < ny + rules["nb_px"])
            is_target = np.hypot(px_n[0] - px_t[0], px_n[1] - px_t[1]) < 0.3
            c2n = None
            if near and not is_target:
                c2n, _ = chi2_fit(np.exp(-0.5 * ((xx - px_n[0]) ** 2 + (yy - px_n[1]) ** 2) / 0.9**2))
            nb.append({"source_id": n["source_id"], "G": n["G"], "px": [float(px_n[0]), float(px_n[1])], "near_stamp": bool(near),
                       "is_target": bool(is_target), "deficit_dist_arcsec": dist(px_n[0], px_n[1], cx, cy), "chi2_gauss_at_source": c2n})
        d_t = dist(px_t[0], px_t[1], cx, cy)
        others = [x["deficit_dist_arcsec"] for x in nb if x["near_stamp"] and not x["is_target"] and x["deficit_dist_arcsec"] is not None]
        d_nb = min(others) if others else None
        sig_arcsec = sig_cent_px * PX_ARCSEC if np.isfinite(sig_cent_px) else None
        conclusive = (not edge_flag) and nev >= rules["min_events"] and n_sig >= rules["min_sig"] and d_t is not None
        on_q = bool(conclusive and sig_arcsec is not None and d_t < 2 * sig_arcsec and (d_nb is None or d_t < d_nb))
        off_q = bool(conclusive and sig_arcsec is not None and d_t > 2 * sig_arcsec and d_nb is not None and d_nb < d_t)
        aperture = tpf.pipeline_mask if tpf.pipeline_mask is not None else valid_px
        ap = aperture & valid_px
        alt = [x["chi2_gauss_at_source"] for x in nb if x["chi2_gauss_at_source"] is not None] + [chi2A]
        chi2_best_alt = float(min(alt))
        best_alt_name = "scaled_image" if chi2_best_alt == chi2A else [x["source_id"] for x in nb if x["chi2_gauss_at_source"] == chi2_best_alt][0]
        row = {"quarter": q, "n_events": nev, "n_cadences": int(len(t)), "stamp_shape": [int(ny), int(nx)], "n_valid_px": int(valid_px.sum()),
               "chi2_best_alternative": chi2_best_alt, "best_alternative": best_alt_name,
               "target_px": [float(px_t[0]), float(px_t[1])], "target_at_edge_or_outside": bool(edge_flag),
               "noise_method": noise_method, "n_random_maps": len(rand_maps),
               "deficit_centroid_px": [cx, cy], "centroid_bootstrap_sigma_px": sig_cent_px, "n_bootstrap": len(bxs), "n_sig_pixels": n_sig,
               "flux_centroid_px": [fx, fy], "deficit_peak_px": [int(bx), int(by)], "peak_snr": float(snr_map[by, bx]) if np.isfinite(snr_map[by, bx]) else None,
               "target_to_deficit_arcsec": d_t, "target_to_deficit_sigma": (d_t / sig_arcsec) if (sig_arcsec and d_t is not None) else None,
               "nearest_other_source_to_deficit_arcsec": d_nb,
               "fluxcentroid_to_deficit_arcsec": dist(fx, fy, cx, cy),
               "deficit_total_e_s": float(np.nansum(np.clip(np.where(valid_px, diff, 0), 0, None))),
               "depth_in_aperture_ppm": float(1e6 * np.nansum(diff[ap]) / np.nansum(img_out[ap])),
               "chi2_null": chi2_null, "chi2_scaled_image": chi2A, "chi2_gauss_at_target": chi2B, "n_pix_fit": int(ok.sum()),
               "conclusive": bool(conclusive), "on_target_quarter": on_q, "offset_quarter": off_q, "neighbours": nb}
        rows.append(row)
        imgs.append((q, np.where(valid_px, img_out, np.nan), np.where(valid_px, snr_map, np.nan), px_t, [n["px"] for n in nb if n["near_stamp"] and not n["is_target"]], conclusive))
    res["quarters"] = rows
    concl = [r for r in rows if r.get("conclusive")]
    allq = [r for r in rows if "target_to_deficit_arcsec" in r and r["target_to_deficit_arcsec"] is not None]
    summ = {"n_quarters_with_maps": len(allq), "n_conclusive": len(concl),
            "n_edge_or_outside": int(sum(1 for r in rows if r.get("target_at_edge_or_outside"))),
            "n_on_target_quarters": int(sum(1 for r in concl if r["on_target_quarter"])),
            "n_offset_quarters": int(sum(1 for r in concl if r["offset_quarter"]))}
    if concl:
        summ.update({
            "median_target_to_deficit_arcsec": float(np.median([r["target_to_deficit_arcsec"] for r in concl])),
            "median_target_to_deficit_sigma": float(np.nanmedian([r["target_to_deficit_sigma"] if r["target_to_deficit_sigma"] is not None else np.nan for r in concl])),
            "median_nearest_other_source_arcsec": float(np.nanmedian([r["nearest_other_source_to_deficit_arcsec"] if r["nearest_other_source_to_deficit_arcsec"] is not None else np.nan for r in concl])),
            "quarters_target_gaussian_beats_scaled_image": int(sum(1 for r in concl if r["chi2_gauss_at_target"] < r["chi2_scaled_image"])),
            "sum_chi2_null": float(sum(r["chi2_null"] for r in concl)), "sum_chi2_scaled_image": float(sum(r["chi2_scaled_image"] for r in concl)),
            "sum_chi2_gauss_at_target": float(sum(r["chi2_gauss_at_target"] for r in concl)),
            "median_peak_snr": float(np.median([r["peak_snr"] for r in concl if r["peak_snr"] is not None])) if any(r["peak_snr"] is not None for r in concl) else None,
            "median_depth_in_aperture_ppm": float(np.median([r["depth_in_aperture_ppm"] for r in concl]))})
    n = len(concl)
    if n >= rules["min_q"] and summ["n_on_target_quarters"] / n >= rules["frac"]:
        cverdict = "on-target"
    elif n >= rules["min_q"] and summ["n_offset_quarters"] / n >= rules["frac"]:
        cverdict = "offset"
    else:
        cverdict = "inconclusive"
    summ["centroid_verdict_v15"] = cverdict
    summ["centroid_verdict_reason"] = f"{summ['n_on_target_quarters']} on-target and {summ['n_offset_quarters']} offset of {n} conclusive quarters ({summ['n_edge_or_outside']} edge/outside)"
    # v1.6 chi-square model comparison (does not need a centroid; works at stamp edges)
    qq = [r for r in allq if r.get("n_events", 0) >= 2 and np.isfinite(r.get("chi2_gauss_at_target", np.nan))]
    d_null = float(sum(r["chi2_null"] - r["chi2_gauss_at_target"] for r in qq))
    d_alt = float(sum(r["chi2_best_alternative"] - r["chi2_gauss_at_target"] for r in qq))
    worst_alt = None
    if qq:
        names = {}
        for r in qq:
            names[r["best_alternative"]] = names.get(r["best_alternative"], 0.0) + (r["chi2_gauss_at_target"] - r["chi2_best_alternative"])
        worst_alt = max(names.items(), key=lambda kv: kv[1])
    r16 = _load_rules_v16()
    if len(qq) >= r16["min_q"] and d_null >= r16["d_null"] and d_alt >= r16["d_alt"]:
        verdict = "on-target"
    elif len(qq) >= r16["min_q"] and d_alt <= -r16["d_alt"]:
        verdict = "offset"
    else:
        verdict = "inconclusive"
    summ.update({"verdict": verdict, "verdict_rule": "v1.6 chi2 model comparison", "n_quarters_chi2": len(qq),
                 "sum_dchi2_null_minus_target": d_null, "sum_dchi2_bestalt_minus_target": d_alt,
                 "strongest_alternative": worst_alt[0] if worst_alt else None, "strongest_alternative_dchi2_over_target": worst_alt[1] if worst_alt else None,
                 "verdict_reason": f"{len(qq)} quarters: null-target dchi2 {d_null:.0f}, best-alternative-target dchi2 {d_alt:.0f} (rule: both >= {r16['d_alt']:.0f} over >= {r16['min_q']} quarters)"})
    res["summary"] = summ
    if make_plot and imgs:
        _plot(out_dir, imgs, fname, verdict)
    with open(os.path.join(out_dir, "pixel.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else (bool(o) if isinstance(o, np.bool_) else str(o)))
    return res


def _plot(out_dir, imgs, fname, verdict):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(imgs)
    cols = min(n, 6)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(2 * rows, cols, figsize=(2.6 * cols, 5.2 * rows))
    axes = np.array(axes).reshape(2 * rows, cols)
    for k, (q, img_out, snr, px_t, nbs, concl) in enumerate(imgs):
        r, c = divmod(k, cols)
        ax1, ax2 = axes[2 * r, c], axes[2 * r + 1, c]
        ax1.imshow(img_out, origin="lower", cmap="gray")
        ax1.set_title(f"Q{q} out-of-transit" + ("" if concl else " (inconclusive)"), fontsize=8)
        v = np.nanmax(np.abs(snr)) if np.isfinite(snr).any() else 1
        ax2.imshow(snr, origin="lower", cmap="RdBu_r", vmin=-v, vmax=v)
        ax2.set_title("missing-light SNR (red = lost in transit)", fontsize=8)
        for ax in (ax1, ax2):
            ax.plot(px_t[0], px_t[1], "o", mfc="none", mec="lime", ms=12, mew=1.5)
            for p in nbs:
                ax.plot(p[0], p[1], "x", color="orange", ms=6)
            ax.set_xlim(-0.5, img_out.shape[1] - 0.5); ax.set_ylim(-0.5, img_out.shape[0] - 0.5)
            ax.set_xticks([]); ax.set_yticks([])
    for k in range(n, rows * cols):
        r, c = divmod(k, cols)
        axes[2 * r, c].axis("off"); axes[2 * r + 1, c].axis("off")
    fig.suptitle(f"{fname}: verdict {verdict}. green circle = target (Gaia); orange x = other Gaia sources", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "diffimages.png"), dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fname", required=True)
    ap.add_argument("--period", type=float, required=True)
    ap.add_argument("--t0-bjd", type=float, required=True)
    ap.add_argument("--dur-h", type=float, required=True)
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--tpforg", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-quarters", type=int, default=None)
    ap.add_argument("--force-analytic", action="store_true", help="use the analytic noise map in every quarter (fallback test)")
    ap.add_argument("--invert", action="store_true", help="localise excess instead of missing light (inverted-flux null)")
    a = ap.parse_args()
    r = pixel_vet(a.fname, a.period, a.t0_bjd, a.dur_h, a.ra, a.dec, a.tpforg, a.out, max_quarters=a.max_quarters, force_analytic=a.force_analytic,
                  invert=a.invert)
    print(json.dumps(r.get("summary", {}), indent=1))
    for q in r["quarters"]:
        if "target_to_deficit_arcsec" in q:
            td, sg = q["target_to_deficit_arcsec"], q["target_to_deficit_sigma"]
            print(f"  Q{q['quarter']}: ev={q['n_events']} sig={q['n_sig_pixels']} edge={q['target_at_edge_or_outside']} noise={q['noise_method']} "
                  f"target->deficit {td and round(td, 1)}\" ({sg and round(sg, 1)} sigma) other->deficit {q['nearest_other_source_to_deficit_arcsec'] and round(q['nearest_other_source_to_deficit_arcsec'], 1)}\" "
                  f"chi2 null/img/gauss {q['chi2_null']:.0f}/{q['chi2_scaled_image']:.0f}/{q['chi2_gauss_at_target']:.0f} -> {'ON' if q['on_target_quarter'] else ('OFF' if q['offset_quarter'] else ('inconclusive' if not q['conclusive'] else 'undecided'))}")
        else:
            print(f"  Q{q.get('quarter')}: {q.get('error')}")


if __name__ == "__main__":
    main()
