"""Red-team per-pixel check for kic-011138423 on the cached TPFs of KIC 11138418.
Season-stacked in-transit deficit maps; empirical-PRF template (from the bright host) placed at the
target position and at positions further out along the host->target line (i.e. beyond the stamp edge)."""
import glob, json, numpy as np
from astropy.io import fits
from astropy.wcs import WCS

P = 4.894451164394656; T0 = 125.33662554571406; DUR = 1.95 / 24.0
RA_T, DEC_T = 292.3076266151539, 48.75114527313966
RA_H, DEC_H = 292.30547569164, 48.74829473264     # KIC 11138418 (Gaia)
OUT = '/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/redteam-kic011138423/out/pixel_checks.json'
files = sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr011138418_lc_Q*/kplr011138418-*_lpd-targ.fits.gz'))
res = {'quarters': [], 'seasons': {}}
rng = np.random.default_rng(7)

def pix_of(hdr_ap, ra, dec):
    w = WCS(hdr_ap)
    x, y = w.all_world2pix(ra, dec, 0)
    return float(x), float(y)

stack = {}
for fn in files:
    h = fits.open(fn)
    hd = h[0].header; q = hd['QUARTER']; season = hd.get('SEASON'); ch = hd['CHANNEL']
    d = h[1].data
    t = d['TIME']; fl = d['FLUX']; ql = d['QUALITY']
    ok = np.isfinite(t) & (ql == 0)
    t = t[ok]; fl = fl[ok]
    ap = h[2]
    xt, yt = pix_of(ap.header, RA_T, DEC_T); xh, yh = pix_of(ap.header, RA_H, DEC_H)
    ny, nx = fl.shape[1:]
    # per-pixel local detrend: for each in-transit cadence, compare to median of out-of-transit cadences
    # within +-1 d (event-differenced), excluding +-1.5 dur around transit
    ph = ((t - T0 + P / 2) % P) - P / 2
    intr = np.abs(ph) < DUR / 2
    oot = (np.abs(ph) > 1.0 * DUR)
    deficit = np.zeros((ny, nx)); var = np.zeros((ny, nx)); nev = 0
    ncad = 0
    epochs = np.unique(np.round((t[intr] - T0) / P).astype(int))
    for n in epochs:
        tc = T0 + n * P
        a = intr & (np.abs(t - tc) < DUR)
        b = oot & (np.abs(t - tc) < 1.0) & (np.abs(t - tc) > 0.75 * DUR)
        if a.sum() < 2 or b.sum() < 15:
            continue
        img_in = np.nanmean(fl[a], axis=0); img_out = np.nanmedian(fl[b], axis=0)
        sig = np.nanstd(fl[b], axis=0) / np.sqrt(a.sum())
        dmap = img_out - img_in
        good = np.isfinite(dmap) & np.isfinite(sig) & (sig > 0)
        deficit[good] += dmap[good]; var[good] += sig[good]**2; nev += 1; ncad += a.sum()
    if nev < 3:
        continue
    img_oot = np.nanmedian(fl[oot], axis=0)
    key = int(season)
    if key not in stack:
        stack[key] = {'def': np.zeros((ny, nx)), 'var': np.zeros((ny, nx)), 'oot': np.zeros((ny, nx)), 'n': 0, 'xt': [], 'yt': [], 'xh': [], 'yh': [], 'shape': (ny, nx), 'quarters': []}
    s = stack[key]
    if s['shape'] != (ny, nx):
        # shape differs within season: skip stacking but still record
        res['quarters'].append({'quarter': int(q), 'season': key, 'note': 'shape differs, not stacked', 'shape': [ny, nx]})
        continue
    s['def'] += deficit; s['var'] += var; s['oot'] += img_oot * nev; s['n'] += nev
    s['xt'].append(xt); s['yt'].append(yt); s['xh'].append(xh); s['yh'].append(yh); s['quarters'].append(int(q))
    snr = deficit / np.sqrt(var + 1e-30)
    res['quarters'].append({'quarter': int(q), 'season': key, 'channel': int(ch), 'n_events': nev, 'shape': [ny, nx],
                            'target_px': [xt, yt], 'host_px': [xh, yh],
                            'target_pixel_snr': float(snr[int(round(yt)), int(round(xt))]) if 0 <= round(yt) < ny and 0 <= round(xt) < nx else None,
                            'target_pixel_deficit_e_s': float(deficit[int(round(yt)), int(round(xt))] / nev) if 0 <= round(yt) < ny and 0 <= round(xt) < nx else None,
                            'max_snr': float(np.nanmax(snr)), 'max_snr_px': [int(np.argmax(np.nanmax(snr, axis=0))), int(np.argmax(np.nanmax(snr, axis=1)))]})

# ---- season stacks: empirical PRF template from the host, placed at trial positions ----
for key, s in stack.items():
    ny, nx = s['shape']; n = s['n']
    dmap = s['def'] / n; sig = np.sqrt(s['var']) / n; oot = s['oot'] / n
    snr = dmap / sig
    xt, yt = np.mean(s['xt']), np.mean(s['yt']); xh, yh = np.mean(s['xh']), np.mean(s['yh'])
    yy, xx = np.mgrid[0:ny, 0:nx]
    # empirical PRF: oot image, host-centred; measure radial profile of the host (fraction of host total in a pixel vs distance)
    rh = np.hypot(xx - xh, yy - yh)
    tot = np.nansum(oot)
    prof = [(float(r), float(v / tot)) for r, v in zip(rh.ravel(), oot.ravel()) if np.isfinite(v)]
    # template: shift the oot image (dominated by the host) so that the host lands at (x0,y0); use bilinear shift via scipy
    from scipy.ndimage import shift as ndshift
    def template(x0, y0):
        return ndshift(np.nan_to_num(oot), (y0 - yh, x0 - xh), order=1, mode='constant', cval=0.0)
    good = np.isfinite(dmap) & np.isfinite(sig) & (sig > 0)
    def chi2_at(x0, y0):
        m = template(x0, y0)
        A = np.sum(m[good] * dmap[good] / sig[good]**2) / np.sum(m[good]**2 / sig[good]**2)
        c2 = np.sum(((dmap[good] - A * m[good]) / sig[good])**2)
        Aerr = 1.0 / np.sqrt(np.sum(m[good]**2 / sig[good]**2))
        return float(c2), float(A * tot), float(Aerr * tot)
    chi2_null = float(np.sum((dmap[good] / sig[good])**2))
    c2_t, amp_t, amperr_t = chi2_at(xt, yt)
    c2_h, amp_h, _ = chi2_at(xh, yh)
    # along the host->target line, beyond the target (outside the stamp)
    ux, uy = (xt - xh), (yt - yh); L = np.hypot(ux, uy); ux /= L; uy /= L
    line = []
    for dpx in [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0]:
        c2, amp, ae = chi2_at(xt + dpx * ux, yt + dpx * uy)
        line.append([dpx, round(c2, 2), round(amp, 2), round(ae, 2)])
    # grid search for best position over the stamp and 2 px beyond it
    best = (1e30, None, None)
    for gx in np.arange(-2, nx + 2, 0.25):
        for gy in np.arange(-2, ny + 2, 0.25):
            c2, amp, _ = chi2_at(gx, gy)
            if amp > 0 and c2 < best[0]:
                best = (c2, gx, gy)
    # 1-sigma region: delta chi2 < 2.3 (2 params)
    region = []
    for gx in np.arange(-2, nx + 2, 0.25):
        for gy in np.arange(-2, ny + 2, 0.25):
            c2, amp, _ = chi2_at(gx, gy)
            if amp > 0 and c2 - best[0] < 2.3:
                region.append((gx, gy))
    region = np.array(region)
    dist_best_target = float(np.hypot(best[1] - xt, best[2] - yt))
    # host-wing fraction at the target position: profile value at r = |target-host|
    r_t = float(np.hypot(xt - xh, yt - yh))
    near = [v for r, v in prof if abs(r - r_t) < 0.5]
    # flux in the target pixel and its contrast
    res['seasons'][str(key)] = {
        'quarters': s['quarters'], 'n_events_stacked': int(n), 'shape': [ny, nx], 'target_px': [xt, yt], 'host_px': [xh, yh],
        'target_to_host_px': r_t, 'host_wing_fraction_per_pixel_at_target_distance_median': float(np.median(near)) if near else None,
        'snr_map': np.round(snr, 1).tolist(), 'deficit_map_e_s': np.round(dmap, 3).tolist(),
        'chi2_null': chi2_null, 'chi2_prf_at_target': c2_t, 'amp_at_target_e_s': amp_t, 'amp_err_at_target_e_s': amperr_t,
        'chi2_prf_at_host': c2_h, 'amp_at_host_e_s': amp_h,
        'line_host_to_target_dpx_chi2_amp': line,
        'best_fit_px': [best[1], best[2]], 'chi2_best': best[0], 'best_to_target_px': dist_best_target,
        'best_to_target_arcsec': dist_best_target * 3.98,
        'onesigma_region_extent_px': [float(region[:, 0].min()), float(region[:, 0].max()), float(region[:, 1].min()), float(region[:, 1].max())] if len(region) else None,
        'n_good_px': int(good.sum()),
    }
    print(f"season {key} Q{s['quarters']} n_ev={n} target=({xt:.2f},{yt:.2f}) host=({xh:.2f},{yh:.2f}) chi2 null {chi2_null:.1f} target {c2_t:.1f} (amp {amp_t:.2f}+-{amperr_t:.2f} e/s) host {c2_h:.1f}; best ({best[1]:.2f},{best[2]:.2f}) chi2 {best[0]:.1f}, {dist_best_target:.2f} px from target; host wing/pixel at target distance {np.median(near) if near else None}")
    print('   line', line)
    print('   snr map:\n', np.round(snr, 1))

json.dump(res, open(OUT, 'w'), indent=1, default=float)
print('written', OUT)
