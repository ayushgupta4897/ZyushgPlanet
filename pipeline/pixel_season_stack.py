"""Season-stacked pixel localisation (method extension of pipeline/pixel_vet.py v1.6; POST-HOC, prereg v1.8).

Why: the v1.6 chi-square test only uses quarters with >= 2 usable transits. Long-period signals (e.g. 47.69 d, about
two transits per 93-d quarter, often one) leave too few quarters. Kepler rolls every quarter and returns to the same
CCD module/output every fourth quarter (a "season"), so the event-differenced missing-light maps of all quarters in a
season can be placed on one CCD-pixel grid and averaged, giving one map per season with ~4x the events.

Per quarter: quality-filtered TPF cube of the originating star (tpforg); for every usable event (>= 2 in-transit and
>= 6 out-of-transit cadences within 3 durations) the map out - in is added to a per-season canvas at the stamp's CCD
(column, row) offset, with a per-pixel event count. The same is done for K random-phase ephemerides (same events,
shifted by a random fraction 0.15-0.85 of the period; one shift per (quarter, k)), which gives a per-pixel noise map
as the scatter of the K averaged random maps. Per season, the averaged missing-light map is fitted (amplitude only)
by: nothing (null); a Gaussian point source (sigma 0.9 px) at the target's Gaia position; the same at every other
Gaia source within 2 px of the canvas; and a scaled copy of the season's mean image. Positions are converted to the
canvas through each quarter's WCS and averaged within the season.

Rule (v1.8, registered before being applied to any signal whose verdict it could change): sum over seasons with a
map; 'on-target' needs >= 3 seasons AND sum(chi2_null - chi2_target) >= 25 AND sum(chi2_best_alt - chi2_target)
>= 25; 'offset' needs >= 3 seasons and some alternative to beat the target by >= 25 summed; otherwise inconclusive.
It must reproduce: KIC 5689121 on-target; kic-006699887 (15.11 d) on-target; kic-006367225 and kic-001577049 offset.

    python -m pipeline.pixel_season_stack --fname kic-006699887 --period 47.69006 --t0-bjd 2454968.71223 --dur-h 4.5 \
        --ra 294.14563 --dec 42.18799 --tpforg 6699878 --out runs/poolA_main/pixel_season/kic-006699887_P47
"""
from __future__ import annotations

import argparse
import json
import os
import warnings

import numpy as np

from pipeline.pixel_vet import BKJD_OFFSET, PX_ARCSEC, _gaia_neighbours

warnings.filterwarnings("ignore")
SIGMA_PX = 0.9


def _quarter(tpf, period, t0_bkjd, D, ra, dec, neighbours, K, rng):
    t = tpf.time.value
    cube = np.array(tpf.flux.value, float)
    ny, nx = cube.shape[1:]
    valid = np.mean(np.isfinite(cube), axis=0) > 0.9
    cube[:, ~valid] = np.nan
    okc = np.isfinite(t) & (np.sum(np.isfinite(cube.reshape(len(t), -1)), axis=1) >= valid.sum() * 0.9)
    t, cube = t[okc], cube[okc]
    if len(t) < 200 or valid.sum() < 4:
        return None
    hdr = tpf.hdu[0].header
    season = hdr.get("SEASON")
    channel = hdr.get("CHANNEL")
    col0, row0 = int(tpf.column), int(tpf.row)

    def usable(tc):
        return (np.abs(t - tc) < 0.5 * D).sum() >= 2 and ((np.abs(t - tc) > 0.75 * D) & (np.abs(t - tc) < 3.0 * D)).sum() >= 6

    def diff(tc):
        inm = np.abs(t - tc) < 0.5 * D
        outm = (np.abs(t - tc) > 0.75 * D) & (np.abs(t - tc) < 3.0 * D)
        return np.nanmean(cube[outm], axis=0) - np.nanmean(cube[inm], axis=0)

    n0 = int(np.floor((t.min() - t0_bkjd) / period)) - 1
    n1 = int(np.ceil((t.max() - t0_bkjd) / period)) + 1
    tcs = [t0_bkjd + n * period for n in range(n0, n1 + 1) if t.min() - D < t0_bkjd + n * period < t.max() + D]
    ev = [tc for tc in tcs if usable(tc)]
    real = np.zeros((ny, nx)); nreal = 0
    for tc in ev:
        real += diff(tc); nreal += 1
    rsum = np.zeros((K, ny, nx)); rn = np.zeros(K, int)
    if nreal:
        for k in range(K):
            s = rng.uniform(0.15, 0.85) * period
            for tc in ev:
                if usable(tc + s):
                    rsum[k] += diff(tc + s); rn[k] += 1
    w = tpf.wcs
    px_t = w.all_world2pix([[ra, dec]], 0)[0]
    nb = []
    for n in neighbours:
        p = w.all_world2pix([[n["ra"], n["dec"]]], 0)[0]
        nb.append({"source_id": n["source_id"], "G": n["G"], "ccd": [float(col0 + p[0]), float(row0 + p[1])]})
    return {"quarter": int(tpf.quarter), "season": season, "channel": channel, "col0": col0, "row0": row0, "ny": ny, "nx": nx,
            "valid": valid, "img": np.where(valid, np.nanmean(cube, axis=0), np.nan), "real": real, "nreal": nreal,
            "rsum": rsum, "rn": rn, "target_ccd": [float(col0 + px_t[0]), float(row0 + px_t[1])], "nb": nb, "n_tcs": len(tcs)}


def _stack(qs, K):
    cmin = min(q["col0"] for q in qs); rmin = min(q["row0"] for q in qs)
    cmax = max(q["col0"] + q["nx"] for q in qs); rmax = max(q["row0"] + q["ny"] for q in qs)
    NX, NY = cmax - cmin, rmax - rmin
    real = np.zeros((NY, NX)); cnt = np.zeros((NY, NX))
    rnd = np.zeros((K, NY, NX)); rcnt = np.zeros((K, NY, NX))
    img = np.zeros((NY, NX)); icnt = np.zeros((NY, NX))
    for q in qs:
        ys, xs = slice(q["row0"] - rmin, q["row0"] - rmin + q["ny"]), slice(q["col0"] - cmin, q["col0"] - cmin + q["nx"])
        v = q["valid"]
        real[ys, xs] += np.where(v, q["real"], 0); cnt[ys, xs] += np.where(v, q["nreal"], 0)
        for k in range(K):
            rnd[k, ys, xs] += np.where(v, q["rsum"][k], 0); rcnt[k, ys, xs] += np.where(v, q["rn"][k], 0)
        img[ys, xs] += np.where(v, q["img"], 0); icnt[ys, xs] += v
    with np.errstate(invalid="ignore", divide="ignore"):
        dmap = np.where(cnt > 0, real / cnt, np.nan)
        rmaps = np.where(rcnt > 0, rnd / rcnt, np.nan)
        imap = np.where(icnt > 0, img / icnt, np.nan)
    sig = np.nanstd(rmaps, axis=0)
    tx = np.mean([q["target_ccd"][0] for q in qs]) - cmin
    ty = np.mean([q["target_ccd"][1] for q in qs]) - rmin
    nbpos = {}
    for q in qs:
        for n in q["nb"]:
            nbpos.setdefault(n["source_id"], {"G": n["G"], "xs": [], "ys": []})
            nbpos[n["source_id"]]["xs"].append(n["ccd"][0] - cmin); nbpos[n["source_id"]]["ys"].append(n["ccd"][1] - rmin)
    return dmap, sig, imap, (tx, ty), {k: (np.mean(v["xs"]), np.mean(v["ys"]), v["G"]) for k, v in nbpos.items()}, cnt


def _fit(dmap, sig, model, ok):
    mm, dd, ss = model[ok], dmap[ok], sig[ok]
    den = np.sum(mm * mm / ss ** 2)
    if den <= 1e-12:
        return float(np.sum((dd / ss) ** 2)), 0.0
    a = np.sum(mm * dd / ss ** 2) / den
    return float(np.sum(((dd - a * mm) / ss) ** 2)), float(a)


def run(fname, period, t0_bjd, dur_h, ra, dec, tpforg, out, K=48, seed=1, max_quarters=None):
    import lightkurve as lk
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(seed)
    D = dur_h / 24.0
    t0 = t0_bjd - BKJD_OFFSET
    neighbours = _gaia_neighbours(ra, dec)
    sr = lk.search_targetpixelfile(f"KIC {int(tpforg)}", mission="Kepler", cadence="long")
    qs, errors = [], []
    for i in range(len(sr)):
        if max_quarters and i >= max_quarters:
            break
        try:
            tpf = sr[i].download(quality_bitmask="default")
            q = _quarter(tpf, period, t0, D, ra, dec, neighbours, K, rng)
            if q is not None and q["nreal"] > 0:
                qs.append(q)
        except Exception as ex:  # noqa: BLE001
            errors.append(str(ex)[:120])
    seasons = {}
    for q in qs:
        seasons.setdefault((q["season"], q["channel"]), []).append(q)
    rows, imgs = [], []
    for (season, channel), group in sorted(seasons.items(), key=lambda kv: (kv[0][0] is None, kv[0])):
        dmap, sig, imap, (tx, ty), nbpos, cnt = _stack(group, K)
        ok = np.isfinite(dmap) & np.isfinite(sig) & (sig > 0)
        if ok.sum() < 6:
            continue
        yy, xx = np.indices(dmap.shape)
        g = lambda x0, y0: np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / SIGMA_PX ** 2)  # noqa: E731
        c_null = float(np.sum((dmap[ok] / sig[ok]) ** 2))
        c_t, a_t = _fit(dmap, sig, g(tx, ty), ok)
        c_img, a_img = _fit(dmap, sig, np.where(np.isfinite(imap), imap, 0.0), ok)
        alts = [("scaled_image", c_img, None)]
        H, W = dmap.shape
        for sid, (x, y, G) in nbpos.items():
            if np.hypot(x - tx, y - ty) < 0.3 or not (-2 <= x < W + 2 and -2 <= y < H + 2):
                continue
            c, _ = _fit(dmap, sig, g(x, y), ok)
            alts.append((sid, c, G))
        best = min(alts, key=lambda z: z[1])
        rows.append({"season": season, "channel": channel, "quarters": [q["quarter"] for q in group],
                     "n_events": int(sum(q["nreal"] for q in group)), "n_pix_fit": int(ok.sum()), "canvas_shape": list(dmap.shape),
                     "target_canvas_px": [float(tx), float(ty)], "chi2_null": c_null, "chi2_target": c_t, "amp_target_e_s": a_t,
                     "chi2_scaled_image": c_img, "chi2_best_alternative": best[1], "best_alternative": best[0], "best_alternative_G": best[2],
                     "dchi2_null_minus_target": c_null - c_t, "dchi2_bestalt_minus_target": best[1] - c_t,
                     "alternatives": [{"name": n, "chi2": c, "G": G} for n, c, G in alts]})
        imgs.append((season, imap, dmap / sig, (tx, ty), [(x, y) for sid, (x, y, G) in nbpos.items()]))
    d_null = float(sum(r["dchi2_null_minus_target"] for r in rows))
    d_alt = float(sum(r["dchi2_bestalt_minus_target"] for r in rows))
    by_alt = {}
    for r in rows:
        for alt in r["alternatives"]:
            by_alt[alt["name"]] = by_alt.get(alt["name"], 0.0) + (r["chi2_target"] - alt["chi2"])
    worst = max(by_alt.items(), key=lambda kv: kv[1]) if by_alt else (None, None)
    if len(rows) >= 3 and d_null >= 25 and d_alt >= 25:
        verdict = "on-target"
    elif len(rows) >= 3 and worst[1] is not None and worst[1] >= 25:
        verdict = "offset"
    else:
        verdict = "inconclusive"
    res = {"fname": fname, "period": period, "t0_bjd": t0_bjd, "dur_h": dur_h, "tpforg": int(tpforg), "K_random": K,
           "rule": "v1.8 season-stacked chi2 (post-hoc method extension; >= 3 seasons, summed dchi2 >= 25 vs null and vs best alternative)",
           "n_quarters_used": len(qs), "errors": errors, "seasons": rows,
           "summary": {"verdict": verdict, "n_seasons": len(rows), "n_events_total": int(sum(r["n_events"] for r in rows)),
                       "sum_dchi2_null_minus_target": d_null, "sum_dchi2_bestalt_minus_target": d_alt,
                       "strongest_alternative": worst[0], "strongest_alternative_dchi2_over_target": worst[1],
                       "verdict_reason": f"{len(rows)} seasons, {sum(r['n_events'] for r in rows)} events: null-target {d_null:.0f}, best-alt-target {d_alt:.0f}"}}
    json.dump(res, open(os.path.join(out, "pixel_season.json"), "w"), indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    arrs = {}
    for j, (season, im, snr, (tx, ty), nbs) in enumerate(imgs):
        arrs[f"s{j}_image"] = im.astype(np.float32); arrs[f"s{j}_snr"] = snr.astype(np.float32)
        arrs[f"s{j}_target"] = np.array([tx, ty], np.float32); arrs[f"s{j}_neighbours"] = np.array(nbs, np.float32).reshape(-1, 2)
        arrs[f"s{j}_season"] = np.array([-1 if season is None else int(season)])
    np.savez_compressed(os.path.join(out, "season_maps.npz"), **arrs)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = len(imgs)
        fig, ax = plt.subplots(2, max(n, 1), figsize=(3.2 * max(n, 1), 6.4), squeeze=False)
        for j, (season, im, snr, (tx, ty), nbs) in enumerate(imgs):
            ax[0, j].imshow(im, origin="lower", cmap="gray"); ax[0, j].set_title(f"season {season}: mean image")
            ax[1, j].imshow(snr, origin="lower", cmap="RdBu_r", vmin=-5, vmax=5); ax[1, j].set_title("missing light / sigma")
            for a_ in ax[:, j]:
                a_.plot(tx, ty, "x", color="lime", ms=10, mew=2)
                for x, y in nbs:
                    a_.plot(x, y, "o", mfc="none", color="orange", ms=6)
        fig.suptitle(f"{fname} P={period:.4f} d: {verdict}"); fig.tight_layout(); fig.savefig(os.path.join(out, "season_maps.png"), dpi=110)
    except Exception as ex:  # noqa: BLE001
        print("plot failed:", ex)
    return res


def main():
    ap = argparse.ArgumentParser()
    for k in ("fname", "out"):
        ap.add_argument(f"--{k}", required=True)
    for k in ("period", "t0-bjd", "dur-h", "ra", "dec"):
        ap.add_argument(f"--{k}", type=float, required=True)
    ap.add_argument("--tpforg", type=int, required=True)
    ap.add_argument("--K", type=int, default=48)
    ap.add_argument("--max-quarters", type=int, default=None)
    a = ap.parse_args()
    r = run(a.fname, a.period, a.t0_bjd, a.dur_h, a.ra, a.dec, a.tpforg, a.out, K=a.K, max_quarters=a.max_quarters)
    for s in r["seasons"]:
        print(f"  season {s['season']} ch {s['channel']} q{s['quarters']}: {s['n_events']} events, dchi2 null-target {s['dchi2_null_minus_target']:.0f}, "
              f"best-alt({s['best_alternative']})-target {s['dchi2_bestalt_minus_target']:.0f}")
    print(json.dumps(r["summary"]))


if __name__ == "__main__":
    main()
