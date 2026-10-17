"""Red-team pixel test for kic-006699887 (target sits in the edge column of the
KIC 6699878 postage stamp in every season).

Own implementation, no pipeline/ code.  For every quarter:
  * per-pixel event-differenced deficit map (out-of-transit minus in-transit,
    e/s), using a local window around each predicted transit;
  * a Gaussian PSF width fitted on the out-of-transit image (sum of Gaussians
    at the Gaia positions, common sigma, free amplitudes, constant sky);
  * a free-position Gaussian fit to the deficit map, with the position allowed
    to go OUTSIDE the stamp (beyond the edge the target sits on);
  * chi2 profile along the direction perpendicular to that edge, so we can say
    how far beyond the edge the source could sit;
  * implied fractional depth on the target = fitted deficit amplitude / target
    flux (from the image fit), and on the host.
Season stacks (same geometry) are also fitted.
"""
import glob, json, sys, warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares, minimize_scalar

warnings.filterwarnings("ignore")
OUT = sys.argv[1] if len(sys.argv) > 1 else "."
P = float(sys.argv[2]); T0 = float(sys.argv[3]); T14 = float(sys.argv[4]) / 24.0; TAG = sys.argv[5]
PIX_ARCSEC = 3.98
TPFS = sorted(glob.glob("/root/.cache/lightkurve/mastDownload/Kepler/kplr006699878_lc_Q011111111111111111/*targ.fits.gz"))
GAIA = json.load(open("/home/user/ZyushgPlanet/runs/poolA_main/pixel/kic-006699887/pixel.json"))["gaia_neighbours"]
TARGET_ID = "2077554781499350144"; HOST_ID = "2077554781499349120"


def gauss2d(x, y, x0, y0, sig):
    return np.exp(-0.5 * ((x - x0) ** 2 + (y - y0) ** 2) / sig ** 2) / (2 * np.pi * sig ** 2)


def fit_image(img, ivar, srcs):
    """img: (ny,nx); srcs: list of (x,y) pixel positions. Common sigma, free amps >=0, sky."""
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    good = np.isfinite(img) & (ivar > 0)

    def model(p):
        sig = p[0]; sky = p[1]; amps = p[2:]
        m = np.full(img.shape, sky, float)
        for a, (sx, sy) in zip(amps, srcs):
            m += a * gauss2d(xx, yy, sx, sy, sig)
        return m

    def resid(p):
        return ((img - model(p)) * np.sqrt(ivar))[good]

    p0 = np.r_[0.7, np.nanmedian(img), [max(np.nanmax(img), 1.0) for _ in srcs]]
    lb = np.r_[0.35, -np.inf, [0.0] * len(srcs)]
    ub = np.r_[2.5, np.inf, [np.inf] * len(srcs)]
    r = least_squares(resid, p0, bounds=(lb, ub))
    return r.x, model(r.x), np.sum(r.fun ** 2), good.sum()


def fit_deficit(dmap, dvar, sig, x0=None, y0=None, free_pos=True):
    ny, nx = dmap.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    good = np.isfinite(dmap) & np.isfinite(dvar) & (dvar > 0)
    w = np.sqrt(1.0 / dvar[good])

    def amp_for(px, py):
        g = gauss2d(xx, yy, px, py, sig)[good]
        d = dmap[good]
        a = np.sum(g * d / dvar[good]) / np.sum(g * g / dvar[good])
        chi2 = np.sum(((d - a * g) * w) ** 2)
        return a, chi2

    if not free_pos:
        a, c = amp_for(x0, y0)
        return dict(x=x0, y=y0, amp=a, chi2=c, n=int(good.sum()))
    # grid search over positions (allow outside the stamp by up to 3 px)
    best = None
    for px in np.arange(-3.0, nx + 2.0, 0.1):
        for py in np.arange(-2.0, ny + 1.0, 0.1):
            a, c = amp_for(px, py)
            if best is None or c < best[2]:
                best = (px, py, c, a)
    return dict(x=best[0], y=best[1], amp=best[3], chi2=best[2], n=int(good.sum()))


results = []
season_stacks = {}
for f in TPFS:
    hdr0 = fits.getheader(f, 0); q = hdr0["QUARTER"]; season = hdr0["SEASON"]
    d = fits.getdata(f, 1); ah = fits.getheader(f, 2)
    t = d["TIME"]; flux = d["FLUX"]; ferr = d["FLUX_ERR"]; qual = d["QUALITY"]
    ok = np.isfinite(t) & (qual == 0) & np.all(np.isfinite(flux.reshape(len(t), -1)) | ~np.isfinite(flux.reshape(len(t), -1)), axis=1)
    ok &= np.nansum(np.isfinite(flux.reshape(len(t), -1)), axis=1) > 0
    t = t[ok]; flux = flux[ok]; ferr = ferr[ok]
    ny, nx = flux.shape[1:]
    wcs = WCS(ah)
    srcs = {}
    for g in GAIA:
        px, py = wcs.all_world2pix(g["ra"], g["dec"], 0)
        srcs[g["source_id"]] = (float(px), float(py), g["G"])
    tx, ty, _ = srcs[TARGET_ID]; hx, hy, _ = srcs[HOST_ID]
    # which sources to include in the image fit: inside or within 1.5 px of the stamp
    use = [(sid, v) for sid, v in srcs.items() if -1.5 <= v[0] <= nx + 0.5 and -1.5 <= v[1] <= ny + 0.5]
    img = np.nanmedian(flux, axis=0)
    ivar = 1.0 / (np.nanmedian(ferr, axis=0) ** 2 / max(len(t), 1) + (0.01 * img) ** 2)
    pimg, mimg, chi2img, nimg = fit_image(img, ivar, [(v[0], v[1]) for _, v in use])
    sig = pimg[0]
    amps = {sid: a for (sid, _), a in zip(use, pimg[2:])}
    target_flux_fit = amps.get(TARGET_ID, np.nan); host_flux_fit = amps.get(HOST_ID, np.nan)
    # event-differenced deficit map
    n0 = int(np.floor((t.min() - T0) / P)) - 1; n1 = int(np.ceil((t.max() - T0) / P)) + 1
    num = np.zeros((ny, nx)); den = np.zeros((ny, nx)); nev = 0
    ev_list = []
    for n in range(n0, n1 + 1):
        tc = T0 + n * P
        intr = np.abs(t - tc) < 0.4 * T14          # core of the transit (avoid ingress/egress)
        oot = (np.abs(t - tc) > 0.8 * T14) & (np.abs(t - tc) < 3.0 * T14)
        if intr.sum() < 3 or oot.sum() < 10:
            continue
        # per-pixel linear trend from oot, evaluated at in-transit times
        X = np.vstack([np.ones(oot.sum()), t[oot] - tc]).T
        F = flux[oot].reshape(oot.sum(), -1)
        goodpx = np.all(np.isfinite(F), axis=0)
        beta = np.linalg.lstsq(X, np.where(goodpx, F, 0.0), rcond=None)[0]
        pred_in = (np.vstack([np.ones(intr.sum()), t[intr] - tc]).T @ beta)
        Fin = flux[intr].reshape(intr.sum(), -1)
        diff = np.nanmean(pred_in - Fin, axis=0)           # e/s lost per pixel
        resid_oot = F - X @ beta
        var = np.nanvar(resid_oot, axis=0, ddof=2) / intr.sum() + np.nanvar(resid_oot, axis=0, ddof=2) / oot.sum()
        var = np.where(goodpx & (var > 0), var, np.nan)
        diff = np.where(goodpx, diff, np.nan)
        wv = 1.0 / var
        num += np.where(np.isfinite(wv * diff), wv * diff, 0).reshape(ny, nx)
        den += np.where(np.isfinite(wv), wv, 0).reshape(ny, nx)
        nev += 1
        ev_list.append(dict(n=n, in_total=float(np.nansum(diff)), ))
    if nev < 2:
        results.append(dict(quarter=q, season=season, skipped="too few events")); continue
    dmap = num / den; dvar = 1.0 / den
    dmap[den == 0] = np.nan; dvar[den == 0] = np.nan
    # significance map
    snr = dmap / np.sqrt(dvar)
    # fits
    at_t = fit_deficit(dmap, dvar, sig, tx, ty, free_pos=False)
    at_h = fit_deficit(dmap, dvar, sig, hx, hy, free_pos=False)
    free = fit_deficit(dmap, dvar, sig)
    chi2_null = float(np.nansum(dmap ** 2 / dvar))
    # 1-D chi2 profile along x (the edge direction), y fixed at target y
    prof = []
    for px in np.arange(-3.0, nx + 1.0, 0.25):
        r = fit_deficit(dmap, dvar, sig, px, ty, free_pos=False)
        prof.append((float(px), float(r["chi2"]), float(r["amp"])))
    chi2min_x = min(p[1] for p in prof)
    allowed = [p[0] for p in prof if p[1] - chi2min_x < 4.0]  # 2-sigma range along x
    # target-pixel fraction check: ratio of deficit col0 vs col1 at target row
    r0 = int(round(ty)); r0 = min(max(r0, 0), ny - 1)
    results.append(dict(
        quarter=int(q), season=int(season), stamp=[int(ny), int(nx)], n_events=nev,
        psf_sigma_px=float(sig), image_chi2=float(chi2img), image_npix=int(nimg),
        target_px=[tx, ty], host_px=[hx, hy],
        target_flux_fit_e_s=float(target_flux_fit), host_flux_fit_e_s=float(host_flux_fit),
        deficit_at_target_e_s=float(at_t["amp"]), chi2_at_target=float(at_t["chi2"]),
        deficit_at_host_e_s=float(at_h["amp"]), chi2_at_host=float(at_h["chi2"]),
        free_fit=dict(x=float(free["x"]), y=float(free["y"]), amp=float(free["amp"]), chi2=float(free["chi2"])),
        chi2_null=chi2_null, npix=int(at_t["n"]),
        implied_depth_on_target_ppm=float(1e6 * at_t["amp"] / target_flux_fit) if target_flux_fit > 0 else None,
        implied_depth_on_host_ppm=float(1e6 * at_h["amp"] / host_flux_fit) if host_flux_fit > 0 else None,
        x_profile=prof, x_allowed_2sig=[min(allowed), max(allowed)] if allowed else None,
        dist_free_to_target_arcsec=float(np.hypot(free["x"] - tx, free["y"] - ty) * PIX_ARCSEC),
        dist_free_to_host_arcsec=float(np.hypot(free["x"] - hx, free["y"] - hy) * PIX_ARCSEC),
        deficit_col0_row_t=float(dmap[r0, 0]), deficit_col1_row_t=float(dmap[r0, 1]) if nx > 1 else None,
        snr_col0_row_t=float(snr[r0, 0]), snr_col1_row_t=float(snr[r0, 1]) if nx > 1 else None,
        model_col0_col1_ratio_target=float(gauss2d(0, r0, tx, ty, sig) / gauss2d(1, r0, tx, ty, sig)),
        sum_deficit_e_s=float(np.nansum(dmap)),
    ))
    # season stack (same stamp geometry within a season)
    key = int(season)
    if key not in season_stacks:
        season_stacks[key] = dict(num=np.zeros_like(num), den=np.zeros_like(den), sig=[], tx=tx, ty=ty, hx=hx, hy=hy, tf=[], hf=[], q=[])
    S = season_stacks[key]
    if S["num"].shape == num.shape:
        S["num"] += num; S["den"] += den; S["sig"].append(sig); S["tf"].append(target_flux_fit); S["hf"].append(host_flux_fit); S["q"].append(int(q))
    print(f"Q{q:2d} s{season} sig={sig:.2f}px Ftarget={target_flux_fit:7.0f} Fhost={host_flux_fit:7.0f} nev={nev} "
          f"def@target={at_t['amp']:6.2f} e/s (chi2 {at_t['chi2']:.1f}) def@host={at_h['amp']:6.2f} (chi2 {at_h['chi2']:.1f}) "
          f"null {chi2_null:.1f} free=({free['x']:.2f},{free['y']:.2f}) amp {free['amp']:.2f} chi2 {free['chi2']:.1f} "
          f"| implied depth target {1e6*at_t['amp']/target_flux_fit:6.0f} ppm; x allowed(2sig) {results[-1]['x_allowed_2sig']} "
          f"col0/col1 obs {dmap[r0,0]:.2f}/{dmap[r0,1]:.2f} model ratio {results[-1]['model_col0_col1_ratio_target']:.2f}", flush=True)

stack_res = []
for key, S in sorted(season_stacks.items()):
    dmap = S["num"] / S["den"]; dvar = 1.0 / S["den"]; dmap[S["den"] == 0] = np.nan; dvar[S["den"] == 0] = np.nan
    sig = float(np.median(S["sig"])); tx, ty, hx, hy = S["tx"], S["ty"], S["hx"], S["hy"]
    at_t = fit_deficit(dmap, dvar, sig, tx, ty, free_pos=False); at_h = fit_deficit(dmap, dvar, sig, hx, hy, free_pos=False)
    free = fit_deficit(dmap, dvar, sig)
    ny, nx = dmap.shape
    prof = []
    for px in np.arange(-3.0, nx + 1.0, 0.1):
        # at each x, optimise y too
        best = None
        for py in np.arange(-1.5, ny + 0.5, 0.1):
            r = fit_deficit(dmap, dvar, sig, px, py, free_pos=False)
            if best is None or r["chi2"] < best[0]:
                best = (r["chi2"], py, r["amp"])
        prof.append((float(px), float(best[0]), float(best[1]), float(best[2])))
    cmin = min(p[1] for p in prof)
    allowed1 = [p[0] for p in prof if p[1] - cmin < 1.0]
    allowed2 = [p[0] for p in prof if p[1] - cmin < 4.0]
    allowed3 = [p[0] for p in prof if p[1] - cmin < 9.0]
    tf = float(np.median(S["tf"])); hf = float(np.median(S["hf"]))
    chi2_null = float(np.nansum(dmap ** 2 / dvar))
    sr = dict(season=key, quarters=S["q"], psf_sigma_px=sig, target_px=[tx, ty], host_px=[hx, hy], stamp=[int(ny), int(nx)],
              deficit_at_target_e_s=float(at_t["amp"]), chi2_at_target=float(at_t["chi2"]),
              deficit_at_host_e_s=float(at_h["amp"]), chi2_at_host=float(at_h["chi2"]), chi2_null=chi2_null,
              free_fit=dict(x=float(free["x"]), y=float(free["y"]), amp=float(free["amp"]), chi2=float(free["chi2"])),
              dist_free_to_target_arcsec=float(np.hypot(free["x"] - tx, free["y"] - ty) * PIX_ARCSEC),
              implied_depth_on_target_ppm=float(1e6 * at_t["amp"] / tf), implied_depth_on_host_ppm=float(1e6 * at_h["amp"] / hf),
              x_allowed_1sig=[min(allowed1), max(allowed1)], x_allowed_2sig=[min(allowed2), max(allowed2)], x_allowed_3sig=[min(allowed3), max(allowed3)],
              x_profile=prof, dmap=np.where(np.isfinite(dmap), dmap, None).tolist(), snr=np.where(np.isfinite(dmap), dmap / np.sqrt(dvar), None).tolist())
    stack_res.append(sr)
    print(f"SEASON {key} Q{S['q']} sig={sig:.2f} target x={tx:.2f} y={ty:.2f}: def@target {at_t['amp']:.2f} e/s chi2 {at_t['chi2']:.1f} | "
          f"def@host {at_h['amp']:.2f} chi2 {at_h['chi2']:.1f} | null {chi2_null:.1f} | free ({free['x']:.2f},{free['y']:.2f}) amp {free['amp']:.2f} chi2 {free['chi2']:.1f} "
          f"| implied depth on target {1e6*at_t['amp']/tf:.0f} ppm | x allowed 1/2/3 sig: {sr['x_allowed_1sig']} {sr['x_allowed_2sig']} {sr['x_allowed_3sig']}", flush=True)

json.dump(dict(per_quarter=results, season_stacks=stack_res, ephemeris=dict(P=P, t0_bkjd=T0, T14_d=T14)),
          open(f"{OUT}/pixel_edge_test_{TAG}.json", "w"), indent=1, default=float)
print("written", f"{OUT}/pixel_edge_test_{TAG}.json")
