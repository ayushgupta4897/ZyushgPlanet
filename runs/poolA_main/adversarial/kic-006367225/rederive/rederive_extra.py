#!/usr/bin/env python
"""
Supplementary measurements for the independent re-derivation of kic-006367225.

Re-uses only the loader / detrend / box helpers from rederive_kic006367225.py
(this directory; written from scratch, nothing from ZyushgPlanet/pipeline/).

Adds, at BOTH the pipeline's published ephemeris (P 12.971528 d, t0 BKJD
142.7530, 2.25 h) and at this re-derivation's own ephemeris:
  1. depth under several definitions (box mean, box median, 1/err^2-weighted
     mean, 10 % trimmed mean, central-50 % "core") and as a function of box
     width 0.5-6 h, to see which definition could give ~3000 ppm;
  2. per-season (same module / aperture every 4th quarter) PSF and SAP depths,
     the SAP/PSF ratio per season, and the header CROWDSAP / PSFFRAC the
     ratio should equal if the signal were on the target and undiluted;
  3. absolute flux deficits in e-/s in PSF and SAP flux per season;
  4. the same box depth in the KBONUS-BKG light curves of the nearest Gaia
     neighbours (downloaded into an isolated directory, read as data only):
     the G=13.8 TPF-origin star 8" away and the G=20.5 / 20.9 stars 3.6" / 9.6"
     away that the pipeline's pixel deficit sits closest to in several quarters.

Run:  python -I rederive_extra.py
Writes rederive_extra.json next to this script and prints a summary.
"""
import json
import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rederive_kic006367225 as R  # noqa: E402  (own code, same directory)

OUT = R.OUT
DL = os.path.join(R.SCR, "downloads")
BKJD0 = R.BKJD0

PIPE = dict(P=12.971527521795288, t0=142.75300012519583, dur=2.25 / 24)   # vet.json
MINE = dict(P=12.971955910945276, t0=856.181950227186, dur=2.4 / 24)        # pass2 top peak
SEASON_Q = {0: [4, 8, 12, 16], 1: [1, 5, 9, 13, 17], 2: [2, 6, 10, 14], 3: [3, 7, 11, 15]}

NEIGHBOURS = {  # Gaia DR3 id -> (fname, G, separation from target in arcsec, note)
    "2077645628644904576": ("kic-006367220", 13.82, 8.1, "TPF-origin star (KIC 6367220), 23x brighter than target"),
    "2077645628637721728": ("gaia-dr3-2077645628637721728", 20.49, 3.6, "nearest to pixel deficit in Q2,4,7,8"),
    "2077645628637706112": ("gaia-dr3-2077645628637706112", 20.87, 9.6, "2nd nearest to pixel deficit in Q1,3"),
}


def load_file(path, col):
    h = fits.open(path)
    d = h["LIGHTCURVE_STITCHED"].data
    t = np.asarray(d["TIME"], float)
    f = np.asarray(d[col], float)
    e = np.asarray(d[col + "_ERR"], float)
    q = np.asarray(d["QUARTER"], int)
    qual = np.asarray(d["SAP_QUALITY"], int)
    m = (qual == 0) & np.isfinite(t) & np.isfinite(f) & (f > 0) & np.isfinite(e) & (e > 0)
    o = np.argsort(t[m])
    hdr = {}
    for x in h:
        if x.name.startswith("LIGHTCURVE_Q"):
            hh = x.header
            hdr[int(hh["QUARTER"])] = dict(CROWDSAP=float(hh.get("CROWDSAP", np.nan)),
                                           PSFFRAC=float(hh.get("PSFFRAC", np.nan)),
                                           FLFRCSAP=float(hh.get("FLFRCSAP", np.nan)),
                                           NPIXSAP=int(hh.get("NPIXSAP", -1)),
                                           MODULE=int(hh.get("MODULE", -1)))
    prim = h[0].header
    meta = dict(OBJECT=str(prim.get("OBJECT")), GAIAID=str(prim.get("GAIAID")), GMAG=float(prim.get("GMAG", np.nan)),
                KEPMAG=float(prim.get("KEPMAG", np.nan)), RA=float(prim.get("RA_OBJ", np.nan)),
                DEC=float(prim.get("DEC_OBJ", np.nan)), TPFORG=int(prim.get("TPFORG", -1)))
    return t[m][o], f[m][o], e[m][o], q[m][o], hdr, meta


def masked_detrend(t, f, q, eph, robust=True):
    ph = R.phase_of(t, eph["P"], eph["t0"])
    tmask = ~(np.abs(ph) * eph["P"] < 1.5 * eph["dur"])
    keep_idx = np.arange(len(t))
    # detrend() drops points; we need the error column aligned, so re-derive the
    # kept index set by running detrend on an index-carrying copy.
    t2, y2, q2, dy2, info = R.detrend(t, f, q, win=1.0, transit_mask=tmask, robust=robust)
    # align by time (times are unique and sorted)
    pos = np.searchsorted(t, t2)
    assert np.allclose(t[pos], t2)
    return t2, y2, q2, dy2, pos, info


def depth_defs(t, y, e, P, t0, dur, oot_excl=1.5):
    ph = R.phase_of(t, P, t0)
    half = 0.5 * dur / P
    inb = np.abs(ph) < half
    core = np.abs(ph) < 0.5 * half
    oot = np.abs(ph) > oot_excl * half
    yo = y[oot]
    base_mean, base_med = np.mean(yo), np.median(yo)
    w = 1.0 / e**2
    yi = y[inb]
    k = max(int(0.1 * len(yi)), 1)
    ys = np.sort(yi)
    return dict(
        n_in=int(inb.sum()),
        box_mean_ppm=float((base_mean - np.mean(yi)) * 1e6),
        box_median_ppm=float((base_med - np.median(yi)) * 1e6),
        box_wmean_fluxerr_ppm=float((np.sum(w[oot] * yo) / np.sum(w[oot]) - np.sum(w[inb] * yi) / np.sum(w[inb])) * 1e6),
        box_trimmed10_ppm=float((base_mean - np.mean(ys[k:-k])) * 1e6),
        core50_mean_ppm=float((base_mean - np.mean(y[core])) * 1e6),
        n_core=int(core.sum()),
    )


def depth_vs_width(t, y, P, t0, widths_h):
    out = []
    for wh in widths_h:
        d, nin = R.box_depth(t, y, P, t0, wh / 24)
        out.append(dict(width_h=float(wh), depth_ppm=float(d * 1e6), n_in=int(nin)))
    return out


def season_of(q):
    s = np.full(len(q), -1)
    for k, qs in SEASON_Q.items():
        s[np.isin(q, qs)] = k
    return s


def per_season(t, y, q, P, t0, dur, nrand=1500):
    s = season_of(q)
    res = {}
    for k, qs in SEASON_Q.items():
        sel = s == k
        r = R.subset_depth(t, y, P, t0, dur, sel, nrand=nrand)
        r["quarters"] = qs
        res[f"season_{k}"] = r
    return res


def abs_deficit(t, fnorm_y, f_raw_aligned, q, P, t0, dur):
    """Absolute in-transit deficit in e-/s = fractional box depth x median raw flux, per season."""
    s = season_of(q)
    out = {}
    for k in SEASON_Q:
        sel = s == k
        if sel.sum() < 50:
            continue
        d, nin = R.box_depth(t[sel], fnorm_y[sel], P, t0, dur)
        med = float(np.median(f_raw_aligned[sel]))
        out[f"season_{k}"] = dict(median_flux_e_s=med, frac_depth_ppm=float(d * 1e6),
                                  abs_deficit_e_s=float(d * med), n_in=int(nin))
    return out


def analyse_target(eph, label):
    t, f, e, q, hdr, meta = load_file(R.FITS, "FLUX")
    ts, fs, es, qs, _, _ = load_file(R.FITS, "SAP_FLUX")
    t2, y2, q2, dy2, pos, info = masked_detrend(t, f, q, eph)
    e2 = e[pos] / f[pos]                      # fractional errors aligned to kept points
    f2 = f[pos]
    ts2, ys2, qs2, dys2, poss, sinfo = masked_detrend(ts, fs, qs, eph)
    es2 = es[poss] / fs[poss]
    fs2 = fs[poss]
    P, t0, dur = eph["P"], eph["t0"], eph["dur"]
    out = dict(ephemeris=dict(P=P, t0_bkjd=t0, dur_h=dur * 24), n_used_psf=int(len(t2)), n_used_sap=int(len(ts2)))
    out["psf_depth_definitions"] = depth_defs(t2, y2, e2, P, t0, dur)
    out["sap_depth_definitions"] = depth_defs(ts2, ys2, es2, P, t0, dur)
    widths = [0.5, 1.0, 1.5, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]
    out["psf_depth_vs_width"] = depth_vs_width(t2, y2, P, t0, widths)
    out["sap_depth_vs_width"] = depth_vs_width(ts2, ys2, P, t0, widths)
    ps = per_season(t2, y2, q2, P, t0, dur)
    ss = per_season(ts2, ys2, qs2, P, t0, dur)
    seas = {}
    for k in SEASON_Q:
        key = f"season_{k}"
        cs = [hdr[qq]["CROWDSAP"] for qq in SEASON_Q[k] if qq in hdr]
        pf = [hdr[qq]["PSFFRAC"] for qq in SEASON_Q[k] if qq in hdr]
        npx = [hdr[qq]["NPIXSAP"] for qq in SEASON_Q[k] if qq in hdr]
        mod = sorted({hdr[qq]["MODULE"] for qq in SEASON_Q[k] if qq in hdr})
        dp, ep_ = ps[key]["depth_ppm"], ps[key]["err_ppm"]
        ds, es_ = ss[key]["depth_ppm"], ss[key]["err_ppm"]
        ratio = ds / dp if dp else np.nan
        rerr = abs(ratio) * np.hypot(es_ / ds if ds else np.inf, ep_ / dp if dp else np.inf)
        seas[key] = dict(quarters=SEASON_Q[k], module=mod, npixsap=npx,
                         crowdsap_header=cs, crowdsap_mean=float(np.mean(cs)) if cs else None,
                         psffrac_header=pf, psffrac_mean=float(np.mean(pf)) if pf else None,
                         psf_depth_ppm=dp, psf_err_ppm=ep_, psf_n_in=ps[key]["n_in"],
                         sap_depth_ppm=ds, sap_err_ppm=es_, sap_n_in=ss[key]["n_in"],
                         ratio_sap_over_psf=float(ratio), ratio_err=float(rerr),
                         expected_ratio_if_on_target=float(np.mean(cs)) if cs else None)
    out["per_season"] = seas
    # chi2 of per-season PSF and SAP depths against a constant
    for name, src in (("psf", ps), ("sap", ss)):
        d = np.array([src[f"season_{k}"]["depth_ppm"] for k in SEASON_Q])
        er = np.array([src[f"season_{k}"]["err_ppm"] for k in SEASON_Q])
        ok = np.isfinite(d) & np.isfinite(er) & (er > 0)
        wm = np.sum(d[ok] / er[ok]**2) / np.sum(1 / er[ok]**2)
        out[f"per_season_{name}_chi2_const"] = float(np.sum((d[ok] - wm)**2 / er[ok]**2))
        out[f"per_season_{name}_dof"] = int(ok.sum() - 1)
        out[f"per_season_{name}_weighted_mean_ppm"] = float(wm)
    out["abs_deficit_psf"] = abs_deficit(t2, y2, f2, q2, P, t0, dur)
    out["abs_deficit_sap"] = abs_deficit(ts2, ys2, fs2, qs2, P, t0, dur)
    out["header_per_quarter"] = hdr
    out["meta"] = meta
    out["point_scatter_ppm"] = dict(psf=float(R.robust_std(y2) * 1e6), sap=float(R.robust_std(ys2) * 1e6))
    out["median_flux_e_s"] = dict(psf=float(np.median(f2)), sap=float(np.median(fs2)))
    return out


def analyse_neighbour(gid, eph):
    fname, G, sep, note = NEIGHBOURS[gid]
    path = os.path.join(DL, f"{fname}.fits")
    res = dict(gaia_id=gid, fname=fname, G=G, sep_arcsec=sep, note=note, file=path, available=os.path.exists(path))
    if not res["available"] or os.path.getsize(path) < 10000:
        res["available"] = False
        return res
    try:
        P, t0, dur = eph["P"], eph["t0"], eph["dur"]
        for col in ("FLUX", "SAP_FLUX"):
            t, f, e, q, hdr, meta = load_file(path, col)
            res["meta"] = meta
            t2, y2, q2, dy2, pos, info = masked_detrend(t, f, q, eph)
            d, nin = R.box_depth(t2, y2, P, t0, dur)
            rd = R.random_phase_depths(t2, y2, P, t0, dur, n=2000)
            err = R.robust_std(rd)
            ds, _ = R.box_depth(t2, y2, P, t0, dur, phase_c=0.5)
            # also the deepest box anywhere in phase, same width, to see if a
            # shifted/related event exists
            scan = []
            for pc in np.linspace(-0.5, 0.5, 400, endpoint=False):
                dd, _ = R.box_depth(t2, y2, P, t0, dur, phase_c=pc)
                scan.append((pc, dd))
            scan = np.array(scan)
            k = int(np.nanargmax(scan[:, 1]))
            # odd/even on the neighbour too (an EB at 2P would show it)
            epoch = np.round((t2 - t0) / P).astype(int)
            do = R.subset_depth(t2, y2, P, t0, dur, (epoch % 2) != 0, nrand=800)
            de = R.subset_depth(t2, y2, P, t0, dur, (epoch % 2) == 0, nrand=800)
            res[col] = dict(n_used=int(len(t2)), point_scatter_ppm=float(R.robust_std(y2) * 1e6),
                            median_flux_e_s=float(np.median(f[pos])),
                            box_depth_ppm=float(d * 1e6), empirical_err_ppm=float(err * 1e6),
                            sigma=float(d / err), n_in=int(nin),
                            abs_deficit_e_s=float(d * np.median(f[pos])),
                            secondary_ppm=float(ds * 1e6), secondary_sigma=float(ds / err),
                            deepest_box_any_phase_ppm=float(scan[k, 1] * 1e6), deepest_box_phase=float(scan[k, 0]),
                            odd_ppm=do["depth_ppm"], odd_err_ppm=do["err_ppm"],
                            even_ppm=de["depth_ppm"], even_err_ppm=de["err_ppm"],
                            crowdsap_header=[hdr[qq]["CROWDSAP"] for qq in sorted(hdr)],
                            n_steps_split=len(info["steps"]), n_clipped=info["n_clipped"])
    except Exception as ex:  # file is untrusted data: never let it abort the run
        res["error"] = f"{type(ex).__name__}: {ex}"
    return res


# ----------------------------------------------------------------------------
# raw (un-stitched) per-quarter flux levels: the stitched columns are rescaled
# per quarter to a common median, which hides how much flux the PSF fit / the
# aperture actually assigned to this star in each quarter.
# ----------------------------------------------------------------------------
KP_ZP_E_S = 1.74e5   # approximate Kepler zero point: flux(e/s) ~ 1.74e5 * 10^(-0.4 (Kp - 12))


def raw_quarter_levels(path):
    h = fits.open(path)
    kp = float(h[0].header.get("KEPMAG", np.nan))
    expected = KP_ZP_E_S * 10 ** (-0.4 * (kp - 12.0))
    rows = {}
    for x in h:
        if not x.name.startswith("LIGHTCURVE_Q"):
            continue
        d, hh = x.data, x.header
        m = (np.asarray(d["SAP_QUALITY"]) == 0) & np.isfinite(d["FLUX"]) & np.isfinite(d["SAP_FLUX"])
        Q = int(hh["QUARTER"])
        rows[Q] = dict(module=int(hh.get("MODULE", -1)), npixsap=int(hh.get("NPIXSAP", -1)),
                       crowdsap=float(hh.get("CROWDSAP", np.nan)), psffrac=float(hh.get("PSFFRAC", np.nan)),
                       flfrcsap=float(hh.get("FLFRCSAP", np.nan)),
                       raw_psf_flux_e_s=float(np.nanmedian(d["FLUX"][m])),
                       raw_sap_flux_e_s=float(np.nanmedian(d["SAP_FLUX"][m])),
                       n_good=int(m.sum()))
        rows[Q]["psf_over_expected"] = rows[Q]["raw_psf_flux_e_s"] / expected
        rows[Q]["sap_over_expected"] = rows[Q]["raw_sap_flux_e_s"] / expected
        rows[Q]["sap_target_flux_if_crowdsap_right_e_s"] = rows[Q]["raw_sap_flux_e_s"] * rows[Q]["crowdsap"]
    return dict(kepmag=kp, expected_total_flux_e_s=float(expected), quarters=rows)


def raw_abs_deficits(levels, per_season):
    """Fractional per-season depth (from the stitched, rescaled flux) x median
    RAW flux of that season = absolute deficit in e-/s actually measured."""
    out = {}
    for key, s in per_season.items():
        qs = s["quarters"]
        raw_psf = float(np.median([levels["quarters"][q]["raw_psf_flux_e_s"] for q in qs if q in levels["quarters"]]))
        raw_sap = float(np.median([levels["quarters"][q]["raw_sap_flux_e_s"] for q in qs if q in levels["quarters"]]))
        exp_ = levels["expected_total_flux_e_s"]
        psffrac = s["psffrac_mean"]
        out[key] = dict(
            quarters=qs, raw_psf_flux_e_s=raw_psf, raw_sap_flux_e_s=raw_sap,
            psf_abs_deficit_e_s=s["psf_depth_ppm"] * 1e-6 * raw_psf,
            psf_abs_deficit_err_e_s=s["psf_err_ppm"] * 1e-6 * raw_psf,
            sap_abs_deficit_e_s=s["sap_depth_ppm"] * 1e-6 * raw_sap,
            sap_abs_deficit_err_e_s=s["sap_err_ppm"] * 1e-6 * raw_sap,
            # if the signal were on the target and the whole PSF-fit deficit were the target's:
            depth_if_on_target_from_psf_deficit_ppm=s["psf_depth_ppm"] * raw_psf / exp_,
            # if on target and the aperture captured PSFFRAC of the target's light:
            depth_if_on_target_from_sap_deficit_ppm=s["sap_depth_ppm"] * raw_sap / (exp_ * psffrac) if psffrac else None,
        )
    return out


if __name__ == "__main__" and "--raw-only" in sys.argv:
    jp = os.path.join(OUT, "rederive_extra.json")
    out = json.load(open(jp))
    lev = raw_quarter_levels(R.FITS)
    out["raw_flux_levels_target"] = lev
    nb = os.path.join(DL, "kic-006367220.fits")
    if os.path.exists(nb) and os.path.getsize(nb) > 10000:
        out["raw_flux_levels_kic006367220"] = raw_quarter_levels(nb)
    for lab in ("at_pipeline_ephemeris", "at_own_ephemeris"):
        out[lab]["abs_deficit_vs_raw_flux"] = raw_abs_deficits(lev, out[lab]["per_season"])
    with open(jp, "w") as fh:
        json.dump(R.to_jsonable(out), fh, indent=2)
    print("target Kp", lev["kepmag"], "expected total flux ~", round(lev["expected_total_flux_e_s"]), "e/s")
    for Q, r in sorted(lev["quarters"].items(), key=lambda kv: int(kv[0])):
        print(f" Q{int(Q):2d} mod{r['module']:2d} npix{r['npixsap']} CROWD {r['crowdsap']:.2f} PSFFRAC {r['psffrac']:.2f} | raw PSF {r['raw_psf_flux_e_s']:7.1f} ({r['psf_over_expected']:.2f}x exp) | raw SAP {r['raw_sap_flux_e_s']:7.1f} ({r['sap_over_expected']:.2f}x exp)")
    if "raw_flux_levels_kic006367220" in out:
        l2 = out["raw_flux_levels_kic006367220"]
        print("neighbour KIC 6367220 Kp", l2["kepmag"], "expected ~", round(l2["expected_total_flux_e_s"]))
        for Q, r in sorted(l2["quarters"].items(), key=lambda kv: int(kv[0])):
            print(f" Q{int(Q):2d} mod{r['module']:2d} npix{r['npixsap']} CROWD {r['crowdsap']:.2f} | raw PSF {r['raw_psf_flux_e_s']:8.1f} | raw SAP {r['raw_sap_flux_e_s']:8.1f}")
    for lab in ("at_pipeline_ephemeris", "at_own_ephemeris"):
        print("==", lab)
        for k, v in out[lab]["abs_deficit_vs_raw_flux"].items():
            print(f" {k} Q{v['quarters']}: raw PSF {v['raw_psf_flux_e_s']:.0f} e/s, PSF deficit {v['psf_abs_deficit_e_s']:.2f}+-{v['psf_abs_deficit_err_e_s']:.2f} e/s -> {v['depth_if_on_target_from_psf_deficit_ppm']:.0f} ppm if on target; "
                  f"raw SAP {v['raw_sap_flux_e_s']:.0f} e/s, SAP deficit {v['sap_abs_deficit_e_s']:.2f}+-{v['sap_abs_deficit_err_e_s']:.2f} e/s -> {v['depth_if_on_target_from_sap_deficit_ppm']:.0f} ppm if on target")
    sys.exit(0)


if __name__ == "__main__":
    out = dict(fname=R.FNAME, date="2026-10-17",
               note="supplementary to rederive_results.json; same own code; neighbour files downloaded from MAST "
                    "into an isolated scratch directory and read as data only")
    out["at_pipeline_ephemeris"] = analyse_target(PIPE, "pipe")
    out["at_own_ephemeris"] = analyse_target(MINE, "mine")
    out["neighbours_at_own_ephemeris"] = {gid: analyse_neighbour(gid, MINE) for gid in NEIGHBOURS}
    jp = os.path.join(OUT, "rederive_extra.json")
    with open(jp, "w") as fh:
        json.dump(R.to_jsonable(out), fh, indent=2)
    print("wrote", jp)
    for lab in ("at_pipeline_ephemeris", "at_own_ephemeris"):
        o = out[lab]
        print(f"\n=== {lab} {o['ephemeris']}")
        print(" PSF depth defs:", {k: round(v, 1) if isinstance(v, float) else v for k, v in o["psf_depth_definitions"].items()})
        print(" SAP depth defs:", {k: round(v, 1) if isinstance(v, float) else v for k, v in o["sap_depth_definitions"].items()})
        print(" PSF depth vs width:", [(w["width_h"], round(w["depth_ppm"])) for w in o["psf_depth_vs_width"]])
        print(" SAP depth vs width:", [(w["width_h"], round(w["depth_ppm"])) for w in o["sap_depth_vs_width"]])
        for k, s in o["per_season"].items():
            print(f" {k} Q{s['quarters']} mod{s['module']} npix{s['npixsap']} CROWDSAP~{s['crowdsap_mean']:.2f} PSFFRAC~{s['psffrac_mean']:.2f}: "
                  f"PSF {s['psf_depth_ppm']:.0f}+-{s['psf_err_ppm']:.0f}  SAP {s['sap_depth_ppm']:.0f}+-{s['sap_err_ppm']:.0f}  "
                  f"ratio {s['ratio_sap_over_psf']:.2f}+-{s['ratio_err']:.2f}")
        print(f" per-season chi2(const): PSF {o['per_season_psf_chi2_const']:.1f}/{o['per_season_psf_dof']}  "
              f"SAP {o['per_season_sap_chi2_const']:.1f}/{o['per_season_sap_dof']}; weighted means PSF {o['per_season_psf_weighted_mean_ppm']:.0f} SAP {o['per_season_sap_weighted_mean_ppm']:.0f}")
        print(" abs deficit PSF (e/s):", {k: (round(v["median_flux_e_s"], 1), round(v["abs_deficit_e_s"], 3)) for k, v in o["abs_deficit_psf"].items()})
        print(" abs deficit SAP (e/s):", {k: (round(v["median_flux_e_s"], 1), round(v["abs_deficit_e_s"], 3)) for k, v in o["abs_deficit_sap"].items()})
        print(" scatter ppm:", o["point_scatter_ppm"], " median flux e/s:", o["median_flux_e_s"])
    print("\n=== neighbours at own ephemeris")
    for gid, n in out["neighbours_at_own_ephemeris"].items():
        print(gid, n["fname"], "G=", n["G"], "sep=", n["sep_arcsec"], "available=", n["available"], n.get("error", ""))
        if n.get("meta"):
            print("   meta:", n["meta"])
        for col in ("FLUX", "SAP_FLUX"):
            if col in n:
                c = n[col]
                print(f"   {col}: depth {c['box_depth_ppm']:.0f} +- {c['empirical_err_ppm']:.0f} ppm ({c['sigma']:.1f} sig), "
                      f"abs {c['abs_deficit_e_s']:.2f} e/s of {c['median_flux_e_s']:.0f}; sec {c['secondary_ppm']:.0f} ({c['secondary_sigma']:.1f} sig); "
                      f"deepest any phase {c['deepest_box_any_phase_ppm']:.0f} at {c['deepest_box_phase']:.3f}; "
                      f"odd {c['odd_ppm']:.0f}+-{c['odd_err_ppm']:.0f} even {c['even_ppm']:.0f}+-{c['even_err_ppm']:.0f}; scatter {c['point_scatter_ppm']:.0f} ppm")
