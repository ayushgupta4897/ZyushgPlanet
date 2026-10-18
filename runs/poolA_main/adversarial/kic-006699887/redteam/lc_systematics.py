"""Red-team light-curve checks for kic-006699887 (own code, no pipeline/ imports).
 - event times vs data gaps / quarter edges / thermal ramps; depth near vs far from gaps
 - depth per quarter and per season (stamp geometry changes each season)
 - depth with SAP vs PSF per season
 - Lomb-Scargle on out-of-transit flux (variability / rotation aliases)
 - secondary at every phase; box depth at P/2, 2P phases
 - inverted-flux BLS in a narrow window around P (masked detrend)
 - host star KIC 6699878 KBONUS light curve at the same ephemeris (downloaded)
 - stellar density / duration consistency
"""
import json, os, sys, urllib.request, warnings
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares, LombScargle
warnings.filterwarnings("ignore")
OUT = sys.argv[1]
P = 15.1136115; T0 = 134.72786; T14 = 3.75 / 24.0
FN = "/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-006699887_kepler_v1.0_lc.fits"
HOST_URL = "https://archive.stsci.edu/hlsps/kbonus-bkg/lcs/0066/006699878/hlsp_kbonus-bkg_kepler_kepler_kic-006699878_kepler_v1.0_lc.fits"
HOST_FN = f"{OUT}/dl/hlsp_kbonus-bkg_kepler_kepler_kic-006699878_kepler_v1.0_lc.fits"


def load(fn, col="FLUX"):
    d = fits.getdata(fn, "LIGHTCURVE_STITCHED")
    t = d["TIME"]; f = d[col]; q = d["QUARTER"]; qual = d["SAP_QUALITY"]
    ok = np.isfinite(t) & np.isfinite(f) & (qual == 0) & (f > 0)
    return t[ok], f[ok], q[ok]


def phase_h(t):
    ph = ((t - T0) / P + 0.5) % 1.0 - 0.5
    return ph * P * 24.0


def detrend(t, f, q, win=1.0, mask_transit=True):
    """per-quarter median-normalise, then divide by a running median (time window `win` d)
    computed from out-of-transit points only."""
    out = np.empty_like(f)
    intr = np.abs(phase_h(t)) < 0.75 * T14 * 24 if mask_transit else np.zeros(len(t), bool)
    for qq in np.unique(q):
        m = q == qq
        tt = t[m]; ff = f[m] / np.nanmedian(f[m]); it = intr[m]
        tr = np.empty_like(ff)
        to = tt[~it]; fo = ff[~it]
        lo = np.searchsorted(to, tt - win / 2); hi = np.searchsorted(to, tt + win / 2)
        for i in range(len(tt)):
            seg = fo[lo[i]:hi[i]]
            tr[i] = np.median(seg) if len(seg) >= 5 else np.nan
        out[m] = ff / tr
    ok = np.isfinite(out)
    # clip positive outliers only
    s = 1.4826 * np.nanmedian(np.abs(out[ok] - 1))
    ok &= out < 1 + 4 * s
    return out, ok


def box_depth(t, f, tc_list=None, phase_centre_h=0.0, width_h=3.75):
    ph = phase_h(t) - phase_centre_h
    ph = (ph + P * 12) % (P * 24) - P * 12
    inb = np.abs(ph) < width_h / 2
    oob = (np.abs(ph) > width_h) & (np.abs(ph) < 3 * width_h)
    if inb.sum() < 5:
        return np.nan, np.nan, 0
    d = np.mean(f[oob]) - np.mean(f[inb])
    e = np.sqrt(np.var(f[oob]) / inb.sum() + np.var(f[oob]) / oob.sum())
    return 1e6 * d, 1e6 * e, int(inb.sum())


res = {}
t, f, q = load(FN, "FLUX")
fd, ok = detrend(t, f, q)
t, fd, q, fraw = t[ok], fd[ok], q[ok], f[ok]
tS, fS, qS = load(FN, "SAP_FLUX")
fSd, okS = detrend(tS, fS, qS); tS, fSd, qS = tS[okS], fSd[okS], qS[okS]
print("N psf", len(t), "N sap", len(tS))
res["depth_box_psf"] = box_depth(t, fd)
res["depth_box_sap"] = box_depth(tS, fSd)
print("PSF box depth", res["depth_box_psf"], "SAP", res["depth_box_sap"])

# ---------- gaps / quarter edges / ramps ----------
dt = np.diff(t); gaps = np.where(dt > 0.5)[0]
gap_edges = []  # (end_of_segment, start_of_next)
for i in gaps:
    gap_edges.append((float(t[i]), float(t[i + 1]), float(dt[i])))
res["n_gaps_gt_0.5d"] = len(gap_edges)
qstart = {int(qq): float(t[q == qq].min()) for qq in np.unique(q)}
qend = {int(qq): float(t[q == qq].max()) for qq in np.unique(q)}
events = []
n0 = int(np.floor((t.min() - T0) / P)); n1 = int(np.ceil((t.max() - T0) / P))
for n in range(n0, n1 + 1):
    tc = T0 + n * P
    w = np.abs(t - tc) < 1.5 * T14
    inb = np.abs(t - tc) < 0.5 * T14
    oob = (np.abs(t - tc) > 0.75 * T14) & (np.abs(t - tc) < 3 * T14)
    if inb.sum() < 4 or oob.sum() < 10:
        continue
    d = 1e6 * (np.mean(fd[oob]) - np.mean(fd[inb])); e = 1e6 * np.sqrt(np.var(fd[oob]) * (1 / inb.sum() + 1 / oob.sum()))
    # distance to nearest gap edge and to nearest quarter edge; time since last gap start (ramp)
    dg = min([min(abs(tc - a), abs(tc - b)) for a, b, _ in gap_edges] + [99])
    after_gap = min([tc - b for a, b, _ in gap_edges if tc >= b] + [99])
    dq = min([min(abs(tc - qstart[k]), abs(tc - qend[k])) for k in qstart])
    qq = int(np.median(q[inb]))
    events.append(dict(n=n, tc=tc, depth=d, err=e, n_in=int(inb.sum()), gap_dist=dg, days_after_gap=after_gap, q_edge_dist=dq, quarter=qq, season=(qq - 2) % 4 if qq >= 2 else 3))
res["events"] = events
ev = np.array([(e["depth"], e["err"], e["gap_dist"], e["days_after_gap"], e["q_edge_dist"]) for e in events])
near = ev[:, 2] < 1.0; ramp = ev[:, 3] < 2.0
def wmean(x, e):
    w = 1 / e ** 2; return float(np.sum(w * x) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))
res["depth_near_gap_1d"] = wmean(ev[near, 0], ev[near, 1]) + (int(near.sum()),)
res["depth_far_gap_1d"] = wmean(ev[~near, 0], ev[~near, 1]) + (int((~near).sum()),)
res["depth_within_2d_after_gap"] = wmean(ev[ramp, 0], ev[ramp, 1]) + (int(ramp.sum()),)
res["depth_not_after_gap"] = wmean(ev[~ramp, 0], ev[~ramp, 1]) + (int((~ramp).sum()),)
res["n_events_within_0.5d_of_quarter_edge"] = int((ev[:, 4] < 0.5).sum())
print("events", len(events), "near gap(1d)", res["depth_near_gap_1d"], "far", res["depth_far_gap_1d"], "ramp", res["depth_within_2d_after_gap"], "not", res["depth_not_after_gap"])
# fraction of events in first 2 days after a gap vs expectation
tot_span = t.max() - t.min()
res["expected_frac_within_2d_after_gap"] = float(2.0 * len(gap_edges) / tot_span)
res["observed_frac_within_2d_after_gap"] = float(ramp.mean())
# per-event chi2 vs constant and against a model where depth = 0
w = 1 / ev[:, 1] ** 2; mu = np.sum(w * ev[:, 0]) / np.sum(w)
res["event_chi2_const"] = float(np.sum(w * (ev[:, 0] - mu) ** 2)); res["event_dof"] = len(ev) - 1
res["events_individually_gt_2sig"] = int((ev[:, 0] / ev[:, 1] > 2).sum())
res["events_negative"] = int((ev[:, 0] < 0).sum())
# leave-one-out: largest single-event contribution to the stacked S/N
snr_all = mu / (1 / np.sqrt(np.sum(w)))
loo = []
for i in range(len(ev)):
    m = np.ones(len(ev), bool); m[i] = False
    loo.append(float(np.sum(w[m] * ev[m, 0]) / np.sum(w[m])))
res["loo_depth_min_max"] = [min(loo), max(loo)]

# ---------- per quarter / per season ----------
pq = {}
for qq in np.unique(q):
    m = q == qq
    pq[int(qq)] = box_depth(t[m], fd[m])
res["depth_per_quarter_psf"] = pq
ps = {}
for s in range(4):
    m = np.isin(q, [qq for qq in np.unique(q) if qq >= 2 and (qq - 2) % 4 == s] + ([1, 5, 9, 13, 17] if s == 3 else []))
    # season definition from TPF headers: Q2,6,10,14 -> season 0; Q3,7,11,15 -> 1; Q4,8,12,16 -> 2; Q1,5,9,13,17 -> 3
    ps[s] = dict(psf=box_depth(t[m], fd[m]), sap=box_depth(tS[np.isin(qS, np.unique(q[m]))], fSd[np.isin(qS, np.unique(q[m]))]), quarters=[int(x) for x in np.unique(q[m])])
res["depth_per_season"] = ps
print("per season", json.dumps(ps, default=float))
# season chi2
sd = np.array([ps[s]["psf"][0] for s in range(4)]); se = np.array([ps[s]["psf"][1] for s in range(4)])
res["season_chi2_const"] = float(np.sum(((sd - np.average(sd, weights=1 / se ** 2)) / se) ** 2))
qd = np.array([pq[k][0] for k in pq]); qe = np.array([pq[k][1] for k in pq]); mq = np.isfinite(qd) & (qe > 0)
res["quarter_chi2_const"] = float(np.sum(((qd[mq] - np.average(qd[mq], weights=1 / qe[mq] ** 2)) / qe[mq]) ** 2)); res["quarter_dof"] = int(mq.sum() - 1)

# ---------- secondary / phase scan ----------
scan = []
for ph in np.arange(-P * 12, P * 12, 0.5):
    d, e, n = box_depth(t, fd, phase_centre_h=ph)
    scan.append((float(ph), float(d), float(e)))
scan = np.array(scan)
outside = np.abs(scan[:, 0]) > 1.5 * T14 * 24
res["phase_scan_max_depth_outside_transit"] = [float(x) for x in scan[outside][np.argmax(scan[outside, 1])]]
res["phase_scan_min_depth_outside_transit"] = [float(x) for x in scan[outside][np.argmin(scan[outside, 1])]]
res["phase_scan_sigma_outside"] = float(np.std(scan[outside, 1] / scan[outside, 2]))
res["secondary_phase0.5"] = box_depth(t, fd, phase_centre_h=P * 12)
res["secondary_phase0.5_wide_12h"] = box_depth(t, fd, phase_centre_h=P * 12, width_h=12)
print("secondary", res["secondary_phase0.5"], "max outside", res["phase_scan_max_depth_outside_transit"])
# depth at 2P odd/even and at P/3 phases
res["depth_phase_third"] = [box_depth(t, fd, phase_centre_h=P * 24 / 3), box_depth(t, fd, phase_centre_h=-P * 24 / 3)]

# ---------- Lomb-Scargle (out-of-transit, quarter-normalised, no detrend) ----------
intr = np.abs(phase_h(t)) < T14 * 24
fn = np.empty_like(fraw)
for qq in np.unique(q):
    m = q == qq; fn[m] = fraw[m] / np.nanmedian(fraw[m])
freq = np.linspace(1 / 60.0, 1 / 0.1, 200000)
ls = LombScargle(t[~intr], fn[~intr])
pw = ls.power(freq)
i = np.argsort(pw)[::-1]
top = []
for j in i[:3000]:
    if all(abs(1 / freq[j] - p) > 0.02 * p for p, _, _ in top):
        top.append((float(1 / freq[j]), float(pw[j]), float(np.sqrt(pw[j]) * np.std(fn[~intr]) * 1e6)))
    if len(top) >= 8: break
res["ls_top_periods"] = top
res["ls_power_at_P"] = float(ls.power(np.array([1 / P]))[0]); res["ls_power_at_2P"] = float(ls.power(np.array([1 / (2 * P)]))[0])
print("LS top", top)

# ---------- inverted BLS near P (masked detrend) ----------
bls = BoxLeastSquares(t, 2 - fd)
periods = np.linspace(14.5, 15.8, 20000)
r = bls.power(periods, [2.5 / 24, 3.75 / 24, 5 / 24], objective="snr")
res["inverted_bls_14.5_15.8"] = dict(best_period=float(periods[np.argmax(r.power)]), best_snr=float(np.max(r.power)), depth_ppm=float(1e6 * r.depth[np.argmax(r.power)]), snr_at_P=float(r.power[np.argmin(np.abs(periods - P))]))
blsr = BoxLeastSquares(t, fd)
rr = blsr.power(periods, [2.5 / 24, 3.75 / 24, 5 / 24], objective="snr")
res["real_bls_14.5_15.8"] = dict(best_period=float(periods[np.argmax(rr.power)]), best_snr=float(np.max(rr.power)), depth_ppm=float(1e6 * rr.depth[np.argmax(rr.power)]))
print("inverted near P", res["inverted_bls_14.5_15.8"], "real", res["real_bls_14.5_15.8"])

# ---------- host star light curve ----------
os.makedirs(f"{OUT}/dl", exist_ok=True)
if not os.path.exists(HOST_FN):
    try:
        urllib.request.urlretrieve(HOST_URL, HOST_FN)
    except Exception as e:
        print("host download failed", e)
if os.path.exists(HOST_FN):
    th, fh, qh = load(HOST_FN, "FLUX"); fhd, okh = detrend(th, fh, qh); th, fhd, qh = th[okh], fhd[okh], qh[okh]
    thS, fhS, qhS = load(HOST_FN, "SAP_FLUX"); fhSd, okhS = detrend(thS, fhS, qhS); thS, fhSd = thS[okhS], fhSd[okhS]
    res["host_psf_depth_at_ephemeris"] = box_depth(th, fhd)
    res["host_sap_depth_at_ephemeris"] = box_depth(thS, fhSd)
    res["host_psf_scatter_ppm"] = float(1e6 * 1.4826 * np.median(np.abs(fhd - 1)))
    hdr = fits.getheader(HOST_FN, 0)
    res["host_header"] = {k: hdr.get(k) for k in ["KEPMAG", "GMAG", "GAIAID", "TPFORG", "PARALLAX", "PMRA", "PMDEC"]}
    cat = fits.getdata("/home/user/ZyushgPlanet/data/kbonus-catalogs/hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits")
    print("host PSF depth", res["host_psf_depth_at_ephemeris"], "host SAP", res["host_sap_depth_at_ephemeris"], "scatter", res["host_psf_scatter_ppm"])
    # expected host depth if the signal were on the HOST and leaked into the target PSF flux:
    # target depth 1086 ppm * Ftarget / Fhost = 1086*4629/12632 = 398 ppm would then be the *leak*; i.e. host would show >= ~1000 ppm.
    # expected host depth if the signal is on the target and the PSF fit leaks ~x% of the target's dip into the host: small.
    res["note_host"] = "If the dip were on the host, the host's own PSF/SAP flux would show the full depth (>= the target's 1086 ppm x F_t/F_h would be the leak INTO the target, so host depth >= 1086*4629/12632*... ). See redteam.md."

# ---------- density / duration ----------
G = 6.674e-8; Psec = P * 86400
for T14h, lab in [(3.5, "T14=3.5h"), (3.75, "T14=3.75h"), (4.0, "T14=4.0h")]:
    k = 0.034
    aR = Psec / (np.pi * T14h * 3600) * (1 + k)  # b = 0
    rho = 3 * np.pi / (G * Psec ** 2) * aR ** 3
    res[f"rho_star_b0_{lab}_gcc"] = float(rho)
res["rho_catalogue_target_gcc"] = float(0.7328 / 0.6599 ** 3 * 1.41)
res["rho_host_6699878_gcc_R0.60_0.64"] = [float(0.7 / 0.60 ** 3 * 1.41), float(0.7 / 0.64 ** 3 * 1.41)]
# planet radius with dilution scenarios: PSF flux is deblended, but if the target PSF flux itself contains x% of host light
res["rp_earth_depth1086_R0.66"] = float(np.sqrt(1086e-6) * 0.66 * 109.1)
res["rp_earth_depth1141_R0.66"] = float(np.sqrt(1141e-6) * 0.66 * 109.1)

json.dump(res, open(f"{OUT}/lc_systematics.json", "w"), indent=1, default=float)
print("written")
