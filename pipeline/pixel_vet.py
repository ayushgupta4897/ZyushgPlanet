"""
Pixel-level localisation for a Kepler Bonus signal: which star actually dims?

The KBONUS light curve was extracted from the target-pixel file (TPF) of a *different*
star (tpforg).  We download that TPF for every quarter, build in-transit and
out-of-transit mean images at the signal's ephemeris, form the difference image
(out - in = the light that went missing), and measure where the deficit sits relative
to the Gaia positions of the target and of every neighbour in the stamp.  A real
transit on the target puts the deficit on the target's pixel; a blended eclipsing
binary puts it on the neighbour.

    python -m pipeline.pixel_vet --fname kic-005689121 --period 13.769448 --t0-bjd 2455630.2454 --dur-h 2.3 \
        --ra 285.907827 --dec 40.941988 --tpforg 5689116 --out runs/controls/pixel/kic-005689121

Outputs: <out>/pixel.json (per-quarter offsets and the summary), <out>/diffimages.png.
Positive control (2026-10-17): KIC 5689121 -> deficit closer to the target than to the
7.1-arcsec neighbour in 12 of 13 quarters (independent reproduction agent); this module
must reproduce that before being trusted on pool stars.
"""
from __future__ import annotations

import argparse
import json
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
BKJD_OFFSET = 2454833.0


def _gaia_neighbours(ra, dec, radius_arcsec=25.0):
    """Gaia DR3 sources within the stamp area via CDS XMatch-free Vizier cone (fast, robust)."""
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "BP-RP", "Plx", "RUWE"], row_limit=200)
    r = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=radius_arcsec * u.arcsec, catalog="I/355/gaiadr3")
    if not r:
        return []
    t = r[0]
    return [{"source_id": str(x["Source"]), "ra": float(x["RA_ICRS"]), "dec": float(x["DE_ICRS"]),
             "G": float(x["Gmag"]) if np.isfinite(x["Gmag"]) else None} for x in t]


def _centroid(img, mask=None):
    img = np.array(img, float)
    if mask is not None:
        img = np.where(mask, img, 0.0)
    img = np.clip(img, 0, None)
    tot = img.sum()
    if tot <= 0:
        return np.nan, np.nan
    yy, xx = np.indices(img.shape)
    return float((xx * img).sum() / tot), float((yy * img).sum() / tot)


def pixel_vet(fname, period, t0_bjd, dur_h, ra, dec, tpforg, out_dir, make_plot=True, max_quarters=None):
    import lightkurve as lk
    from astropy.wcs import WCS
    os.makedirs(out_dir, exist_ok=True)
    res = {"fname": fname, "period": period, "t0_bjd": t0_bjd, "dur_h": dur_h, "tpforg": int(tpforg), "quarters": []}
    sr = lk.search_targetpixelfile(f"KIC {int(tpforg)}", mission="Kepler", cadence="long")
    res["n_tpf_found"] = len(sr)
    neighbours = _gaia_neighbours(ra, dec)
    res["gaia_neighbours"] = neighbours
    t0_bkjd = t0_bjd - BKJD_OFFSET
    D = dur_h / 24.0
    rows = []
    imgs = []
    for i in range(len(sr)):
        if max_quarters and i >= max_quarters:
            break
        try:
            tpf = sr[i].download(quality_bitmask="default")
        except Exception as ex:  # noqa: BLE001
            rows.append({"quarter": None, "error": str(ex)[:120]})
            continue
        q = int(tpf.quarter) if hasattr(tpf, "quarter") else i
        t = tpf.time.value  # BKJD
        cube = np.array(tpf.flux.value, float)  # (n, ny, nx), e/s
        ny, nx = cube.shape[1:]
        valid_px = np.mean(np.isfinite(cube), axis=0) > 0.9  # collected pixels
        cube[:, ~valid_px] = np.nan
        okc = np.isfinite(t) & (np.sum(np.isfinite(cube.reshape(len(t), -1)), axis=1) >= valid_px.sum() * 0.9)
        t, cube = t[okc], cube[okc]
        if len(t) < 200 or valid_px.sum() < 4:
            rows.append({"quarter": q, "error": f"too few cadences ({len(t)}) or pixels ({int(valid_px.sum())})"})
            continue
        img_out = np.nanmean(cube, axis=0)

        def event_diff_map(tc_list):
            """Per-pixel (local out - in) averaged over events; returns map and number of events used."""
            acc = np.zeros((ny, nx)); nev = 0
            for tc in tc_list:
                inm = np.abs(t - tc) < 0.5 * D
                outm = (np.abs(t - tc) > 0.75 * D) & (np.abs(t - tc) < 3.0 * D)
                if inm.sum() < 2 or outm.sum() < 6:
                    continue
                acc += np.nanmean(cube[outm], axis=0) - np.nanmean(cube[inm], axis=0)
                nev += 1
            return (acc / nev if nev else acc * np.nan), nev

        n0 = int(np.floor((t.min() - t0_bkjd) / period)) - 1
        n1 = int(np.ceil((t.max() - t0_bkjd) / period)) + 1
        tcs = [t0_bkjd + n * period for n in range(n0, n1 + 1) if t.min() - D < t0_bkjd + n * period < t.max() + D]
        diff, nev = event_diff_map(tcs)
        if nev < 2:
            rows.append({"quarter": q, "error": f"only {nev} usable events"})
            continue
        # per-pixel noise: the same statistic at random phases (shifted ephemerides)
        rng = np.random.default_rng(1)
        rand_maps = []
        for _ in range(24):
            shift = rng.uniform(0.15, 0.85) * period
            m, ne = event_diff_map([tc + shift for tc in tcs])
            if ne >= 2:
                rand_maps.append(m)
        if len(rand_maps) >= 8:
            sig_map = np.nanstd(np.array(rand_maps), axis=0)
        else:
            # analytic fallback (long periods -> few events per quarter): per-pixel cadence scatter
            # propagated through the local in/out differencing and averaged over events
            cad_std = np.nanstd(cube, axis=0)
            fac = []
            for tc in tcs:
                inm = np.abs(t - tc) < 0.5 * D
                outm = (np.abs(t - tc) > 0.75 * D) & (np.abs(t - tc) < 3.0 * D)
                if inm.sum() >= 2 and outm.sum() >= 6:
                    fac.append(1.0 / inm.sum() + 1.0 / outm.sum())
            sig_map = cad_std * np.sqrt(np.mean(fac)) / np.sqrt(len(fac)) if fac else np.full((ny, nx), np.nan)
        sig_map = np.where(valid_px & np.isfinite(sig_map) & (sig_map > 0), sig_map, np.nan)
        snr_map = diff / sig_map
        if not np.isfinite(snr_map).any():
            rows.append({"quarter": q, "error": "no usable noise map"})
            continue
        # pixel coordinates of target and neighbours through the TPF WCS
        w = tpf.wcs
        px_t = w.all_world2pix([[ra, dec]], 0)[0]
        yy, xx = np.indices((ny, nx))
        # significance-weighted deficit centroid (pixels with SNR > 2 only)
        wmap = np.where((snr_map > 2.0) & valid_px, np.clip(diff, 0, None), 0.0)
        n_sig = int(((snr_map > 2.0) & valid_px).sum())
        if wmap.sum() > 0:
            cx, cy = float((xx * wmap).sum() / wmap.sum()), float((yy * wmap).sum() / wmap.sum())
        else:
            cx = cy = np.nan
        fx, fy = _centroid(np.where(valid_px, img_out, 0))
        by, bx = np.unravel_index(np.nanargmax(np.where(valid_px, snr_map, -np.inf)), (ny, nx))
        # model comparison on the deficit map: (A) scaled out-of-transit image (dip on whatever dominates the
        # stamp, i.e. the bright neighbour) vs (B) Gaussian PSF (sigma 0.9 px) at the TARGET position
        ok = valid_px & np.isfinite(sig_map) & (sig_map > 0)
        def chi2_fit(model):
            mm = model[ok]; dd = diff[ok]; ss = sig_map[ok]
            a = np.sum(mm * dd / ss**2) / np.sum(mm * mm / ss**2)
            return float(np.sum(((dd - a * mm) / ss) ** 2)), float(a)
        modelA = np.where(valid_px, img_out, 0.0)
        modelB = np.exp(-0.5 * ((xx - px_t[0]) ** 2 + (yy - px_t[1]) ** 2) / 0.9**2)
        chi2A, aA = chi2_fit(modelA)
        chi2B, aB = chi2_fit(modelB)
        chi2_null = float(np.sum((diff[ok] / sig_map[ok]) ** 2))
        scale = 3.98  # arcsec per Kepler pixel
        def dist(px, py, x, y):
            return float(np.hypot(px - x, py - y) * scale) if np.isfinite(x) else None
        nb = []
        for n in neighbours:
            px_n = w.all_world2pix([[n["ra"], n["dec"]]], 0)[0]
            inside = (0 <= px_n[0] < nx) and (0 <= px_n[1] < ny)
            modelN = np.exp(-0.5 * ((xx - px_n[0]) ** 2 + (yy - px_n[1]) ** 2) / 0.9**2)
            c2n, _ = chi2_fit(modelN) if inside else (None, None)
            nb.append({"source_id": n["source_id"], "G": n["G"], "px": [float(px_n[0]), float(px_n[1])], "in_stamp": bool(inside),
                       "deficit_dist_arcsec": dist(px_n[0], px_n[1], cx, cy), "chi2_gauss_at_source": c2n})
        aperture = tpf.pipeline_mask if tpf.pipeline_mask is not None else valid_px
        ap = aperture & valid_px
        row = {"quarter": q, "n_events": nev, "n_cadences": int(len(t)), "stamp_shape": [int(ny), int(nx)], "n_valid_px": int(valid_px.sum()),
               "target_px": [float(px_t[0]), float(px_t[1])], "deficit_centroid_px": [cx, cy], "n_sig_pixels": n_sig,
               "flux_centroid_px": [fx, fy], "deficit_peak_px": [int(bx), int(by)], "peak_snr": float(snr_map[by, bx]) if np.isfinite(snr_map[by, bx]) else None,
               "target_to_deficit_arcsec": dist(px_t[0], px_t[1], cx, cy),
               "fluxcentroid_to_deficit_arcsec": dist(fx, fy, cx, cy),
               "deficit_total_e_s": float(np.nansum(np.clip(np.where(valid_px, diff, 0), 0, None))),
               "depth_in_aperture_ppm": float(1e6 * np.nansum(diff[ap]) / np.nansum(img_out[ap])),
               "chi2_null": chi2_null, "chi2_scaled_image": chi2A, "chi2_gauss_at_target": chi2B, "n_pix_fit": int(ok.sum()),
               "amp_scaled_image_ppm": aA * 1e6, "amp_gauss_at_target_e_s": aB,
               "neighbours": nb}
        rows.append(row)
        imgs.append((q, np.where(valid_px, img_out, np.nan), np.where(valid_px, snr_map, np.nan), px_t, [n["px"] for n in nb if n["in_stamp"]]))
    res["quarters"] = rows
    good = [r for r in rows if "target_to_deficit_arcsec" in r and r["target_to_deficit_arcsec"] is not None]
    if good:
        d_t = np.array([r["target_to_deficit_arcsec"] for r in good])
        d_f = np.array([r["fluxcentroid_to_deficit_arcsec"] for r in good])
        # nearest neighbour (other than the target itself) distance per quarter
        d_nb = []
        for r in good:
            others = [n["deficit_dist_arcsec"] for n in r["neighbours"] if n["in_stamp"] and n["deficit_dist_arcsec"] is not None
                      and np.hypot(n["px"][0] - r["target_px"][0], n["px"][1] - r["target_px"][1]) > 0.5]
            d_nb.append(min(others) if others else np.nan)
        d_nb = np.array(d_nb, float)
        res["summary"] = {
            "n_quarters_used": len(good),
            "median_target_to_deficit_arcsec": float(np.median(d_t)),
            "median_fluxcentroid_to_deficit_arcsec": float(np.median(d_f)),
            "median_nearest_other_source_to_deficit_arcsec": float(np.nanmedian(d_nb)) if np.isfinite(d_nb).any() else None,
            "quarters_deficit_closer_to_target_than_any_other_source": int(np.nansum(d_t < d_nb)) if np.isfinite(d_nb).any() else None,
            "quarters_deficit_closer_to_target_than_flux_centroid": int(np.sum(d_t < d_f)),
            "median_depth_in_aperture_ppm": float(np.median([r["depth_in_aperture_ppm"] for r in good])),
            "quarters_target_gaussian_beats_scaled_image": int(sum(1 for r in good if r["chi2_gauss_at_target"] < r["chi2_scaled_image"])),
            "sum_chi2_scaled_image": float(sum(r["chi2_scaled_image"] for r in good)),
            "sum_chi2_gauss_at_target": float(sum(r["chi2_gauss_at_target"] for r in good)),
            "sum_chi2_null": float(sum(r["chi2_null"] for r in good)),
            "median_peak_snr": float(np.median([r["peak_snr"] for r in good if r["peak_snr"] is not None])) if any(r["peak_snr"] is not None for r in good) else None,
        }
    if make_plot and imgs:
        _plot(out_dir, imgs, fname)
    with open(os.path.join(out_dir, "pixel.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    return res


def _plot(out_dir, imgs, fname):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(imgs)
    cols = min(n, 6)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(2 * rows, cols, figsize=(2.6 * cols, 5.2 * rows))
    axes = np.array(axes).reshape(2 * rows, cols)
    for k, (q, img_out, diff, px_t, nbs) in enumerate(imgs):
        r, c = divmod(k, cols)
        ax1, ax2 = axes[2 * r, c], axes[2 * r + 1, c]
        ax1.imshow(img_out, origin="lower", cmap="gray")
        ax1.set_title(f"Q{q} out-of-transit", fontsize=8)
        v = np.nanmax(np.abs(diff)) if np.isfinite(diff).any() else 1
        ax2.imshow(diff, origin="lower", cmap="RdBu_r", vmin=-v, vmax=v)
        ax2.set_title("missing-light SNR map (red = light lost in transit)", fontsize=8)
        for ax in (ax1, ax2):
            ax.plot(px_t[0], px_t[1], "o", mfc="none", mec="lime", ms=12, mew=1.5)
            for p in nbs:
                ax.plot(p[0], p[1], "x", color="orange", ms=6)
            ax.set_xticks([])
            ax.set_yticks([])
    for k in range(n, rows * cols):
        r, c = divmod(k, cols)
        axes[2 * r, c].axis("off")
        axes[2 * r + 1, c].axis("off")
    fig.suptitle(f"{fname}: green circle = target (Gaia position); orange x = other Gaia sources", fontsize=9)
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
    a = ap.parse_args()
    r = pixel_vet(a.fname, a.period, a.t0_bjd, a.dur_h, a.ra, a.dec, a.tpforg, a.out, max_quarters=a.max_quarters)
    print(json.dumps(r.get("summary", {}), indent=1))
    for q in r["quarters"]:
        if "target_to_deficit_arcsec" in q:
            td = q['target_to_deficit_arcsec']; fd = q['fluxcentroid_to_deficit_arcsec']
            print(f"  Q{q['quarter']}: events={q['n_events']} sigpix={q['n_sig_pixels']} peakSNR={q['peak_snr'] and round(q['peak_snr'],1)} target->deficit {td and round(td,1)}\"  fluxcentroid->deficit {fd and round(fd,1)}\"  depth_in_ap {q['depth_in_aperture_ppm']:.0f} ppm  chi2 null/scaledimg/gaussTarget = {q['chi2_null']:.0f}/{q['chi2_scaled_image']:.0f}/{q['chi2_gauss_at_target']:.0f}")
        else:
            print(f"  Q{q.get('quarter')}: {q.get('error')}")


if __name__ == "__main__":
    main()
