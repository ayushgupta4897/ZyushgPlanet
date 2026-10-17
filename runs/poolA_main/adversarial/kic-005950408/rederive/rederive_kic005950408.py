#!/usr/bin/env python
"""
Independent re-derivation for kic-005950408 (KBONUS-BKG PSF light curve).

Written from scratch with numpy / scipy / astropy only.  Does NOT import or read
anything under pipeline/.  Nothing here is a planet; this re-derives a first-screen
signal's period, depth, duration and consistency statistics.

Data sets
  A : all SAP_QUALITY==0 cadences (what a blind pipeline would use)
  B : A minus quarters whose PSF-fit median RED_CHI2 > 1e5 (data-driven, decided
      AFTER looking at the light curve; labelled post-hoc).  For this star those
      are Q0,1,5,9,13,17 = every season-3 quarter, all on CCD channel 25.

Steps
 1. Load LIGHTCURVE_STITCHED, keep SAP_QUALITY==0 and finite values, normalise
    per quarter, detrend with a 1-d time-windowed running median, two-pass,
    upward-only 3-MAD clip (own code, not wotan).
 2. Blind BLS 1-50 d (astropy BoxLeastSquares, own log-uniform frequency grid),
    real flux and inverted flux (2-f) null, on A and on B.
 3. Characterise (box depth with random-phase empirical error, odd/even,
    secondary at phase 0.5, per-year depths, N transits, duration) at
    (a) the blind top peak and (b) the pipeline period 38.60670 d with my own
    t0/duration, on A and on B.
 4. Same depth on SAP_FLUX (dilution check).
Outputs: rederive_results.json, rederive_fold.png, bls_*.npz (cached spectra).
"""
import json, os, sys, time
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FNAME = "kic-005950408"
FITS = "/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-005950408_kepler_v1.0_lc.fits"
OUTDIR = os.path.dirname(os.path.abspath(__file__))
RNG = np.random.default_rng(20261017)
PIPE_P, PIPE_DEPTH_PPM, PIPE_DUR_H = 38.60670, 4104.0, 0.9     # task text, comparison only
BAD_RCHI2 = 1e5

# ----------------------------------------------------------------------------- loading
def load(col="FLUX"):
    with fits.open(FITS) as F:
        d = F["LIGHTCURVE_STITCHED"].data
        t = np.asarray(d["TIME"], float); y = np.asarray(d[col], float)
        e = np.asarray(d[col + "_ERR"], float); q = np.asarray(d["SAP_QUALITY"], int)
        qu = np.asarray(d["QUARTER"], int)
    m = (q == 0) & np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (y > 0)
    t, y, e, qu = t[m], y[m], e[m], qu[m]
    o = np.argsort(t)
    return t[o], y[o], e[o], qu[o]

def quarter_table():
    rows = {}
    with fits.open(FITS) as F:
        for q in range(18):
            h = F[f"LIGHTCURVE_Q{q}"].header; d = F[f"LIGHTCURVE_Q{q}"].data
            g = d["SAP_QUALITY"] == 0
            def mad(x):
                x = x[g]; x = x[np.isfinite(x)]
                if len(x) < 10: return None
                x = x / np.median(x); return float(1.4826 * np.median(np.abs(x - np.median(x))) * 1e6)
            rows[q] = dict(season=h.get("SEASON"), module=h.get("MODULE"), output=h.get("OUTPUT"),
                           channel=h.get("CHANNEL"), psffrac=h.get("PSFFRAC"), flfrcsap=h.get("FLFRCSAP"),
                           npixsap=h.get("NPIXSAP"), median_red_chi2=float(np.nanmedian(d["RED_CHI2"][g])),
                           raw_mad_psf_ppm=mad(d["FLUX"]), raw_mad_sap_ppm=mad(d["SAP_FLUX"]),
                           median_psf_flux=float(np.nanmedian(d["FLUX"][g])),
                           median_sap_flux=float(np.nanmedian(d["SAP_FLUX"][g])), n_quality0=int(g.sum()))
    return rows

def running_median(t, y, window=1.0):
    half = window / 2.0
    lo = np.searchsorted(t, t - half, side="left"); hi = np.searchsorted(t, t + half, side="right")
    out = np.empty_like(y)
    for i in range(len(t)):
        out[i] = np.median(y[lo[i]:hi[i]])
    return out

def detrend(t, y, e, qu, window=1.0, nclip=3.0):
    yn, en = np.empty_like(y), np.empty_like(e)
    for Q in np.unique(qu):
        s = qu == Q; med = np.median(y[s]); yn[s] = y[s] / med; en[s] = e[s] / med
    trend = np.ones_like(yn)
    for Q in np.unique(qu):
        s = qu == Q; trend[s] = running_median(t[s], yn[s], window)
    r = yn / trend - 1.0
    mad = 1.4826 * np.median(np.abs(r - np.median(r)))
    good = r < nclip * mad                       # upward-only clip
    use = good & (r > -5 * mad)                  # deep dips excluded from the trend, kept in data
    trend2 = np.ones_like(yn)
    for Q in np.unique(qu):
        s = qu == Q; su = s & use
        trend2[s] = np.interp(t[s], t[su], running_median(t[su], yn[su], window))
    f = yn / trend2; ef = en / trend2
    return t[good], f[good], ef[good], qu[good], float(mad)

# ----------------------------------------------------------------------------- BLS
DURATIONS = np.array([0.03, 0.045, 0.06, 0.08, 0.11, 0.15, 0.20, 0.27])

def make_grid(t, pmin=1.0, pmax=50.0, dur_min=0.04, oversample=2.0):
    T = t.max() - t.min()
    dlnf = dur_min / (oversample * T)
    lnf = np.arange(np.log(1.0 / pmax), np.log(1.0 / pmin), dlnf)
    return 1.0 / np.exp(lnf)[::-1]

def run_bls_cached(tag, t, f, ef, periods):
    fn = os.path.join(OUTDIR, f"bls_{tag}.npz")
    key = np.array([len(t), float(t.sum()), float(f.sum()), len(periods)])
    if os.path.exists(fn):
        z = np.load(fn)
        if np.allclose(z["key"], key):
            return dict(power=z["power"], depth=z["depth"], duration=z["duration"], t0=z["t0"])
    bls = BoxLeastSquares(t, f, dy=ef)
    res = bls.power(periods, DURATIONS, objective="likelihood")
    out = dict(power=np.asarray(res.power), depth=np.asarray(res.depth),
               duration=np.asarray(res.duration), t0=np.asarray(res.transit_time))
    np.savez(fn, key=key, periods=periods, **out)
    return out

def sde_of(power):
    return float((power.max() - power.mean()) / power.std())

def sde_at(power, i):
    return float((power[i] - power.mean()) / power.std())

def top_peaks(periods, power, n=3, exclude_frac=0.02):
    idx = np.argsort(power)[::-1]; peaks = []
    for i in idx:
        P = periods[i]
        if all(abs(P - q) / q > exclude_frac for q, _ in peaks):
            peaks.append((float(P), int(i)))
        if len(peaks) == n: break
    return peaks

def coherent_with_pipeline(P, T, dur=0.0375):
    """True if P is within one transit duration of phase drift over the baseline of the pipeline P
    (or of an integer alias of it)."""
    for k in (1, 2, 3, 1 / 2, 1 / 3):
        P0 = PIPE_P * k
        if abs(P - P0) < dur * P0 / T: return k
    return None

def near_pipeline(P, tol=0.01):
    for k in (1, 2, 3, 4, 1/2, 1/3, 1/4, 3/2, 2/3):
        if abs(P / PIPE_P - k) < tol: return k
    return None

# ----------------------------------------------------------------------------- characterisation
def phase_fold(t, P, t0):
    return ((t - t0) / P + 0.5) % 1.0 - 0.5

def box_depth(t, f, P, t0, dur, phase_center=0.0, n_random=2000, rng=RNG):
    ph = phase_fold(t, P, t0); hw = 0.5 * dur / P
    inb = np.abs(((ph - phase_center + 0.5) % 1.0) - 0.5) < hw
    oot = np.abs(ph) > 2 * hw
    if inb.sum() == 0: return np.nan, np.nan, 0
    base = np.median(f[oot]); depth = base - np.mean(f[inb])
    rs = []
    for c in rng.uniform(-0.5, 0.5, n_random * 3):
        if abs(c) < 3 * hw or abs(((c - phase_center + 0.5) % 1.0) - 0.5) < 3 * hw: continue
        m = np.abs(((ph - c + 0.5) % 1.0) - 0.5) < hw
        if m.sum() < max(3, 0.3 * inb.sum()): continue
        rs.append(base - np.mean(f[m]))
        if len(rs) >= n_random: break
    rs = np.array(rs)
    return float(depth), float(np.std(rs)) if len(rs) > 10 else np.nan, int(inb.sum())

def per_event(t, f, P, t0, dur, min_loc=3):
    n = np.round((t - t0) / P); ev = []
    for k in np.unique(n):
        tc = t0 + k * P
        inb = np.abs(t - tc) < dur / 2
        loc = (np.abs(t - tc) > dur) & (np.abs(t - tc) < 3 * dur)
        if inb.sum() >= 1 and loc.sum() >= min_loc:
            d = np.median(f[loc]) - np.mean(f[inb]); scat = np.std(f[loc]) / np.sqrt(inb.sum())
            ev.append(dict(epoch=int(k), tc_bkjd=float(tc), n_in=int(inb.sum()), n_loc=int(loc.sum()),
                           depth_ppm=float(d * 1e6), err_ppm=float(scat * 1e6),
                           in_transit_flux_ppm=[float(x) for x in ((f[inb] - 1) * 1e6).round(0)]))
    return ev

def refine(t, f, ef, P0, dur0, halfwidth=0.005):
    bls = BoxLeastSquares(t, f, dy=ef); T = t.max() - t.min()
    dP = 0.02 * dur0 * P0 / T
    periods = np.arange(P0 * (1 - halfwidth), P0 * (1 + halfwidth), dP)
    durs = np.array([0.02, 0.025, 0.03, 0.035, 0.0375, 0.04, 0.045, 0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.15, 0.2])
    res = bls.power(periods, durs, objective="likelihood", oversample=20)
    i = int(np.argmax(res.power))
    return dict(period_d=float(res.period[i]), t0_bkjd=float(res.transit_time[i]), duration_d=float(res.duration[i]),
                bls_depth_ppm=float(res.depth[i] * 1e6), bls_depth_err_ppm=float(res.depth_err[i] * 1e6),
                bls_depth_snr=float(res.depth_snr[i]), power=float(res.power[i]))

def at_fixed_period(t, f, ef, P, durs=(0.0375,)):
    bls = BoxLeastSquares(t, f, dy=ef)
    res = bls.power(np.array([P]), np.array(durs), objective="likelihood", oversample=40)
    i = int(np.argmax(res.power))
    return dict(period_d=P, t0_bkjd=float(res.transit_time[i]), duration_d=float(res.duration[i]),
                bls_depth_ppm=float(res.depth[i] * 1e6), bls_depth_err_ppm=float(res.depth_err[i] * 1e6),
                bls_depth_snr=float(res.depth_snr[i]), power=float(res.power[i]))

def duration_scan(t, f, P, t0):
    durs = np.arange(0.01, 0.25, 0.0025); ph = phase_fold(t, P, t0) * P
    oot = np.abs(ph) > 0.4; base, sig = np.median(f[oot]), np.std(f[oot]); snr = []
    for d in durs:
        m = np.abs(ph) < d / 2
        snr.append(0 if m.sum() < 2 else (base - f[m].mean()) * np.sqrt(m.sum()) / sig)
    snr = np.array(snr)
    return float(durs[np.argmax(snr)]), float(snr.max())

def year_of_quarter(q):
    return 1 if q <= 4 else 2 if q <= 8 else 3 if q <= 12 else 4

def characterise(label, t, f, ef, qu, eph, rng=RNG):
    P, t0, dur = eph["period_d"], eph["t0_bkjd"], eph["duration_d"]
    out = dict(label=label, ephemeris=dict(eph, t0_bjd=t0 + 2454833.0, duration_h=dur * 24))
    dsc, _ = duration_scan(t, f, P, t0); out["ephemeris"]["duration_h_snr_scan"] = dsc * 24
    depth, derr, nin = box_depth(t, f, P, t0, dur)
    out["depth"] = dict(depth_ppm=depth * 1e6, err_ppm_random_phase=derr * 1e6,
                        empirical_snr=depth / derr if derr > 0 else None, n_in=nin)
    ph = phase_fold(t, P, t0); n = np.round((t - t0) / P).astype(int); hw = 0.5 * dur / P
    inb = np.abs(ph) < hw; oot = np.abs(ph) > 2 * hw
    base = np.median(f[oot]); sig = np.std(f[oot])
    odd = inb & (n % 2 == 1); even = inb & (n % 2 == 0)
    if odd.sum() and even.sum():
        d_odd = base - f[odd].mean(); d_even = base - f[even].mean()
        e_odd = sig / np.sqrt(odd.sum()); e_even = sig / np.sqrt(even.sum())
        out["odd_even"] = dict(depth_odd_ppm=d_odd * 1e6, err_odd_ppm=e_odd * 1e6, n_odd=int(odd.sum()),
                               depth_even_ppm=d_even * 1e6, err_even_ppm=e_even * 1e6, n_even=int(even.sum()),
                               difference_sigma=float(abs(d_odd - d_even) / np.hypot(e_odd, e_even)))
    sec, sec_err, nsec = box_depth(t, f, P, t0, dur, phase_center=0.5)
    sec_w, sec_w_err, nsec_w = box_depth(t, f, P, t0, 3 * dur, phase_center=0.5)
    out["secondary"] = dict(depth_ppm=sec * 1e6, err_ppm=sec_err * 1e6, sigma=sec / sec_err if sec_err > 0 else None,
                            n_in=nsec, wide3x_depth_ppm=sec_w * 1e6, wide3x_err_ppm=sec_w_err * 1e6,
                            wide3x_sigma=sec_w / sec_w_err if sec_w_err > 0 else None, wide3x_n_in=nsec_w)
    ev = per_event(t, f, P, t0, dur)
    out["events"] = ev; out["n_transits_with_data"] = len(ev)
    out["n_epochs_in_baseline"] = int(np.floor((t.max() - t0) / P) - np.ceil((t.min() - t0) / P) + 1)
    yrs = {}; yr = np.array([year_of_quarter(q) for q in qu])
    for Y in (1, 2, 3, 4):
        s = yr == Y
        if inb[s].sum() == 0 or (s & oot).sum() < 10:
            yrs[f"year{Y}"] = dict(depth_ppm=None, err_ppm=None, n_in=int(inb[s].sum())); continue
        yrs[f"year{Y}"] = dict(depth_ppm=float((np.median(f[s & oot]) - f[s & inb].mean()) * 1e6),
                               err_ppm=float(np.std(f[s & oot]) / np.sqrt(inb[s].sum()) * 1e6),
                               n_in=int(inb[s].sum()), quarters=sorted(int(x) for x in np.unique(qu[s])))
    out["per_year"] = yrs
    if ev:
        dd = np.array([x["depth_ppm"] for x in ev]); ee = np.array([x["err_ppm"] for x in ev])
        snr_ev = dd / np.maximum(ee, 1.0); kmax = int(np.argmax(snr_ev))
        drop = np.abs(t - ev[kmax]["tc_bkjd"]) > dur
        d_wo, e_wo, _ = box_depth(t[drop], f[drop], P, t0, dur, n_random=800)
        out["single_event_check"] = dict(strongest_event_epoch=ev[kmax]["epoch"], strongest_event_depth_ppm=float(dd[kmax]),
                                         depth_without_it_ppm=d_wo * 1e6, err_without_it_ppm=e_wo * 1e6,
                                         snr_without_it=d_wo / e_wo if e_wo > 0 else None,
                                         n_events_depth_over_2sigma=int(np.sum(snr_ev > 2)),
                                         n_events_total=len(ev),
                                         per_event_snr_sorted=sorted(snr_ev.round(2).tolist(), reverse=True))
    core = np.abs(ph) < hw / 3
    if core.sum():
        out["shape"] = dict(depth_core_third_ppm=float((base - f[core].mean()) * 1e6), n_core=int(core.sum()),
                            core_over_full=float((base - f[core].mean()) / depth) if depth > 0 else None)
    return out, (ph, inb, odd, even, base)

# ----------------------------------------------------------------------------- main
def main():
    T0 = time.time()
    out = dict(fname=FNAME, fits=FITS, script=os.path.abspath(__file__), date="2026-10-17",
               pipeline_numbers_from_task=dict(period_d=PIPE_P, depth_ppm=PIPE_DEPTH_PPM, duration_h=PIPE_DUR_H),
               method=dict(detrend="per-quarter median normalisation; 1.0 d time-windowed running median, two-pass "
                                   "(upward 3-MAD clip removed; >5-MAD dips excluded from trend estimate)",
                           bls="astropy BoxLeastSquares, objective=likelihood, durations " + str(DURATIONS.tolist()) +
                               " d, log-uniform frequency grid dlnf = 0.04 d / (2 * baseline), 1-50 d",
                           sde="(max power - mean power) / std power over the 1-50 d grid",
                           depth="median(out-of-box) - mean(in-box); error = std of the same statistic over "
                                 "random-phase boxes of the same width (empirical)",
                           dataset_B="post-hoc: drop quarters whose median PSF-fit RED_CHI2 > %g" % BAD_RCHI2))
    qt = quarter_table(); out["quarter_table"] = qt
    bad_q = sorted(q for q, r in qt.items() if r["median_red_chi2"] > BAD_RCHI2)
    out["dataset_B_dropped_quarters"] = bad_q
    out["dataset_B_dropped_channels"] = sorted(set(qt[q]["channel"] for q in bad_q))
    print("bad quarters (RED_CHI2>1e5):", bad_q, "channels", out["dataset_B_dropped_channels"])

    t, y, e, qu = load("FLUX")
    out["n_cadences_quality0"] = int(len(t))
    datasets = {}
    tA, fA, eA, qA, madA = detrend(t, y, e, qu); datasets["A"] = (tA, fA, eA, qA)
    mB = ~np.isin(qu, bad_q)
    tB, fB, eB, qB, madB = detrend(t[mB], y[mB], e[mB], qu[mB]); datasets["B"] = (tB, fB, eB, qB)
    out["datasets"] = dict(A=dict(n=int(len(tA)), mad_ppm=madA * 1e6, baseline_d=float(tA.max() - tA.min()),
                                  description="all SAP_QUALITY==0 cadences, after upward clip"),
                           B=dict(n=int(len(tB)), mad_ppm=madB * 1e6, baseline_d=float(tB.max() - tB.min()),
                                  quarters=sorted(int(x) for x in np.unique(qB)),
                                  description="A minus quarters %s" % bad_q))
    print(f"A: {len(tA)} pts MAD {madA*1e6:.0f} ppm;  B: {len(tB)} pts MAD {madB*1e6:.0f} ppm")

    # ---- blind BLS on A and B, real + inverted
    out["blind_bls"] = {}
    spectra = {}
    for tag in ("A", "B"):
        tt, ff, ee, qq = datasets[tag]
        periods = make_grid(tt); T = tt.max() - tt.min()
        t1 = time.time(); r = run_bls_cached(tag + "_real", tt, ff, ee, periods)
        ri = run_bls_cached(tag + "_inv", tt, 2.0 - ff, ee, periods)
        print(f"[{tag}] BLS real+inverted on {len(periods)} periods in {time.time()-t1:.0f}s")
        pw, pwi = r["power"], ri["power"]
        pk = top_peaks(periods, pw); pki = top_peaks(periods, pwi)
        j = int(np.argmin(np.abs(periods - PIPE_P)))
        # top peak within the pipeline's own range (0.5-40 d, here 1-40 d)
        m40 = periods <= 40.0; k40 = int(np.argmax(pw[m40])); P40 = float(periods[m40][k40])
        out["blind_bls"][tag] = dict(
            n_periods=int(len(periods)), sde_real=sde_of(pw), sde_inverted=sde_of(pwi),
            top3_real=[dict(period_d=P, power=float(pw[i]), sde=sde_at(pw, i), depth_ppm=float(r["depth"][i] * 1e6),
                            duration_h=float(r["duration"][i] * 24), t0_bkjd=float(r["t0"][i]),
                            within_1pct_of_pipeline_P_or_alias=near_pipeline(P),
                            phase_coherent_with_pipeline_P=coherent_with_pipeline(P, T)) for P, i in pk],
            top3_inverted=[dict(period_d=P, power=float(pwi[i]), sde=sde_at(pwi, i), depth_ppm=float(ri["depth"][i] * 1e6),
                                duration_h=float(ri["duration"][i] * 24)) for P, i in pki],
            at_pipeline_P=dict(period_d=float(periods[j]), power=float(pw[j]), sde=sde_at(pw, j),
                               duration_h=float(r["duration"][j] * 24), depth_ppm=float(r["depth"][j] * 1e6),
                               rank_among_grid=int((pw > pw[j]).sum() + 1)),
            top_peak_1_to_40d=dict(period_d=P40, power=float(pw[m40][k40]), sde=sde_at(pw, int(np.where(m40)[0][k40])),
                                   phase_coherent_with_pipeline_P=coherent_with_pipeline(P40, T)),
            recovered_blind=bool(coherent_with_pipeline(pk[0][0], T) == 1))
        spectra[tag] = (periods, pw, pwi)
        print(f"[{tag}] SDE real {sde_of(pw):.2f} inverted {sde_of(pwi):.2f}; top3 real:",
              [(round(P, 4), round(float(pw[i]), 1)) for P, i in pk], "; at pipeline P: pow %.1f rank %d" % (pw[j], (pw > pw[j]).sum() + 1))

    # ---- characterise
    out["characterisation"] = {}
    folds = {}
    for tag in ("A", "B"):
        tt, ff, ee, qq = datasets[tag]
        periods, pw, _ = spectra[tag]
        Ptop, itop = top_peaks(periods, pw, 1)[0]
        r = run_bls_cached(tag + "_real", tt, ff, ee, periods)
        eph_top = refine(tt, ff, ee, Ptop, float(r["duration"][itop]))
        c_top, fold_top = characterise(f"{tag}: blind top peak (refined)", tt, ff, ee, qq, eph_top)
        eph_pipe = at_fixed_period(tt, ff, ee, PIPE_P, durs=(0.025, 0.03, 0.0375, 0.045, 0.06, 0.08, 0.11, 0.15, 0.2))
        c_pipe, fold_pipe = characterise(f"{tag}: pipeline period fixed (own t0/duration)", tt, ff, ee, qq, eph_pipe)
        eph_pipe09 = at_fixed_period(tt, ff, ee, PIPE_P, durs=(0.0375,))
        c_pipe09, fold_pipe09 = characterise(f"{tag}: pipeline period and 0.9 h duration fixed (own t0)", tt, ff, ee, qq, eph_pipe09)
        eph_loc = refine(tt, ff, ee, PIPE_P, 0.0375, halfwidth=0.003)
        c_loc, _ = characterise(f"{tag}: local refine within 0.3% of pipeline P", tt, ff, ee, qq, eph_loc)
        out["characterisation"][tag] = dict(blind_top=c_top, pipeline_P_fixed=c_pipe, pipeline_P_and_dur_fixed=c_pipe09,
                                            local_refine_near_pipeline_P=c_loc)
        folds[tag] = dict(top=(eph_top, fold_top), pipe=(eph_pipe09, fold_pipe09))
        for c in (c_top, c_pipe, c_pipe09, c_loc):
            d = c["depth"]; e_ = c["ephemeris"]
            print(f"  {c['label']}: P={e_['period_d']:.5f} t0={e_['t0_bkjd']:.4f} dur={e_['duration_h']:.2f}h "
                  f"depth {d['depth_ppm']:.0f}+-{d['err_ppm_random_phase']:.0f} ppm (SNR {d['empirical_snr']:.1f}, n_in {d['n_in']}) "
                  f"oe {c.get('odd_even',{}).get('difference_sigma', float('nan')):.2f}s sec {c['secondary']['depth_ppm']:.0f}+-{c['secondary']['err_ppm']:.0f} "
                  f"N_ev {c['n_transits_with_data']} yrs " + str([round(v['depth_ppm']) if v['depth_ppm'] is not None else None for v in c['per_year'].values()]))

    # ---- SAP_FLUX dilution check at the pipeline ephemeris (dataset A and B, PSF vs SAP)
    ts, ys, es, qs = load("SAP_FLUX")
    out["sap_check"] = {}
    for tag in ("A", "B"):
        m = np.ones(len(ts), bool) if tag == "A" else ~np.isin(qs, bad_q)
        tS, fS, eS, qS, madS = detrend(ts[m], ys[m], es[m], qs[m])
        tt, ff, ee, qq = datasets[tag]
        res = {}
        for which, (eph, _) in folds[tag].items():
            P, t0, dur = eph["period_d"], eph["t0_bkjd"], eph["duration_d"]
            dS, eSd, nS = box_depth(tS, fS, P, t0, dur)
            dP, ePd, nP = box_depth(tt, ff, P, t0, dur)
            # expected SAP/PSF depth ratio if the dip is on this star: flux-weighted over quarters used
            qs_used = np.unique(qq)
            w = np.array([qt[int(q)]["n_quality0"] for q in qs_used], float)
            exp_ratio = float(np.sum(w * np.array([qt[int(q)]["flfrcsap"] * qt[int(q)]["median_psf_flux"] / qt[int(q)]["median_sap_flux"] for q in qs_used])) / w.sum())
            res[which] = dict(period_d=P, t0_bkjd=t0, duration_h=dur * 24,
                              psf_depth_ppm=dP * 1e6, psf_err_ppm=ePd * 1e6,
                              sap_depth_ppm=dS * 1e6, sap_err_ppm=eSd * 1e6, sap_n_in=nS,
                              ratio_sap_over_psf=float(dS / dP) if dP != 0 else None,
                              expected_ratio_if_on_target=exp_ratio, sap_mad_ppm=madS * 1e6)
            print(f"  SAP[{tag}/{which}] P={P:.4f}: PSF {dP*1e6:.0f}+-{ePd*1e6:.0f}  SAP {dS*1e6:.0f}+-{eSd*1e6:.0f}  ratio {dS/dP:.2f} (expected if on target {exp_ratio:.2f})")
        out["sap_check"][tag] = res
    out["sap_check"]["note"] = ("SAP_FLUX is a 2-4 pixel aperture sum that also contains neighbours' light; if the dip is on this "
                                "star, SAP depth / PSF depth ~ flfrcsap * F_psf / F_sap (flux-weighted over the quarters used). "
                                "NB: SAP_FLUX is very noisy in season-2 quarters (Q4,8,12,16: raw MAD 60-77k ppm) and season-3 quarters.")

    # ---- summary / comparison
    cA = out["characterisation"]["A"]["pipeline_P_and_dur_fixed"]; cB = out["characterisation"]["B"]["pipeline_P_and_dur_fixed"]
    out["summary"] = dict(
        blind_top_peak_A_d=out["blind_bls"]["A"]["top3_real"][0]["period_d"],
        blind_top_peak_B_d=out["blind_bls"]["B"]["top3_real"][0]["period_d"],
        recovered_blind_A=out["blind_bls"]["A"]["recovered_blind"], recovered_blind_B=out["blind_bls"]["B"]["recovered_blind"],
        depth_at_pipeline_ephemeris_A_ppm=cA["depth"]["depth_ppm"], depth_at_pipeline_ephemeris_B_ppm=cB["depth"]["depth_ppm"],
        pipeline_depth_ppm=PIPE_DEPTH_PPM)
    out["runtime_s"] = time.time() - T0
    def _clean(o):
        if isinstance(o, dict): return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [_clean(v) for v in o]
        if isinstance(o, (np.floating, float)): return None if not np.isfinite(o) else float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.bool_): return bool(o)
        return o
    with open(os.path.join(OUTDIR, "rederive_results.json"), "w") as fh:
        json.dump(_clean(out), fh, indent=2)

    # ---- figure
    fig, ax = plt.subplots(3, 2, figsize=(14, 13))
    for col, tag in enumerate(("A", "B")):
        periods, pw, pwi = spectra[tag]; tt, ff, ee, qq = datasets[tag]
        a = ax[0, col]
        a.plot(periods, pw, lw=0.4, color="0.3", label="real")
        a.plot(periods, pwi, lw=0.4, color="tab:orange", alpha=0.6, label="inverted (2-f)")
        a.axvline(PIPE_P, color="tab:red", ls="--", lw=1, label=f"pipeline P={PIPE_P}")
        a.set_xscale("log"); a.set_xlabel("period [d]"); a.set_ylabel("BLS power (dlogL)")
        a.set_title(f"dataset {tag} ({out['datasets'][tag]['n']} pts): blind BLS 1-50 d, SDE real {sde_of(pw):.1f} / inverted {sde_of(pwi):.1f}", fontsize=9)
        a.legend(fontsize=7)
        for row, which in ((1, "top"), (2, "pipe")):
            eph, (ph, inb, odd, even, base) = folds[tag][which]
            P, t0, dur = eph["period_d"], eph["t0_bkjd"], eph["duration_d"]
            c = out["characterisation"][tag]["blind_top" if which == "top" else "pipeline_P_and_dur_fixed"]
            a = ax[row, col]; phd = ph * P * 24; s = np.abs(phd) < 10
            a.plot(phd[s], (ff[s] - 1) * 1e6, ".", ms=2, color="0.65", alpha=0.6)
            bins = np.arange(-10, 10.01, 0.25); bc = 0.5 * (bins[1:] + bins[:-1])
            bm = [np.mean(ff[s][(phd[s] >= bins[i]) & (phd[s] < bins[i + 1])]) if np.any((phd[s] >= bins[i]) & (phd[s] < bins[i + 1])) else np.nan for i in range(len(bc))]
            a.plot(bc, (np.array(bm) - 1) * 1e6, "o-", ms=3, color="tab:blue", label="15-min bins")
            a.plot(phd[odd], (ff[odd] - 1) * 1e6, "^", ms=5, color="tab:purple", label=f"odd in-box (n={odd.sum()})")
            a.plot(phd[even], (ff[even] - 1) * 1e6, "v", ms=5, color="tab:green", label=f"even in-box (n={even.sum()})")
            a.axvspan(-dur * 12, dur * 12, color="tab:orange", alpha=0.15)
            d = c["depth"]
            a.axhline(-d["depth_ppm"], color="tab:red", lw=1, label=f"box depth {d['depth_ppm']:.0f}+-{d['err_ppm_random_phase']:.0f} ppm")
            a.set_ylim(-15000, 8000); a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm")
            ttl = "blind top peak" if which == "top" else "pipeline P (+0.9 h box), own t0"
            a.set_title(f"{tag}: {ttl}: P={P:.5f} d t0={t0:.3f} BKJD dur={dur*24:.2f} h; o/e {c.get('odd_even',{}).get('difference_sigma',float('nan')):.1f}s", fontsize=9)
            a.legend(fontsize=7, loc="lower right")
    fig.suptitle(f"{FNAME} independent re-derivation. A = all quality-0; B = without channel-25 quarters {bad_q}. First-screen signal, not a planet.", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "rederive_fold.png"), dpi=120)
    print("done in %.0f s" % out["runtime_s"])

if __name__ == "__main__":
    main()
