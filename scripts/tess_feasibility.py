"""Feasibility of a TESS search of faint, nearby northern M dwarfs (docs/02 section 5, item 2). Nothing here is a search.

    python scripts/tess_feasibility.py pool        # Gaia DR3 pool + TOI/CTOI hosts      -> runs/tess_feasibility/pool_summary.json
    python scripts/tess_feasibility.py coverage    # light-curve products for a sample   -> runs/tess_feasibility/coverage_sample.csv
    python scripts/tess_feasibility.py controls    # known-TOI recovery (TARS, TGLC)     -> runs/tess_feasibility/controls.json
    python scripts/tess_feasibility.py inject      # injection-recovery on random stars  -> runs/tess_feasibility/inject.json

Pool: Gaia DR3, parallax > 10 mas at > 10 sigma, 1.8 < BP-RP < 5, M_G > 8, ecliptic latitude > 0, G 13.8-17.5, then
T 13.5-16 from the Stassun et al. (2019) G -> T colour relation (BP-RP clipped at 3, so the reddest stars are approximate).
Light curves: TARS (Boyle, Bouma & Mann 2026, DOI 10.17909/2sff-fn29; unpopular CPM, S1-96, T < 16) and TGLC
(Han & Brandt 2023, DOI 10.17909/610m-9474; Gaia-prior PSF photometry, S1-55). Downloads are cached in data/tess/.
"""
import argparse
import datetime as dt
import glob
import io
import json
import os
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

import numpy as np
import pandas as pd
import requests

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.search import bls_search  # noqa: E402

CACHE, OUT = "data/tess", "runs/tess_feasibility"
POOL = f"{CACHE}/pool_T13p5_16.csv"
GAIA_Q = """
SELECT g.source_id, g.ra, g.dec, g.l, g.b, g.ecl_lat, g.parallax, g.parallax_over_error,
       g.phot_g_mean_mag AS gmag, g.bp_rp, g.ruwe, g.phot_g_mean_mag + 5*LOG10(g.parallax/100.0) AS abs_g
FROM gaiadr3.gaia_source AS g
WHERE g.parallax > 10 AND g.parallax_over_error > 10 AND g.bp_rp > 1.8 AND g.bp_rp < 5.0
  AND g.phot_g_mean_mag BETWEEN 13.8 AND 17.5 AND g.phot_g_mean_mag + 5*LOG10(g.parallax/100.0) > 8.0 AND g.ecl_lat > 0
"""
EXOFOP = {"toi": "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
          "ctoi": "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"}
# Pecaut & Mamajek (2013, v2022) dwarf sequence, M0V-M8V: absolute G -> radius (feasibility-grade stellar radii)
MG = np.array([8.19, 8.66, 9.33, 10.27, 11.68, 13.13, 14.85, 16.3])
RS = np.array([0.588, 0.501, 0.446, 0.361, 0.274, 0.196, 0.137, 0.114])
DURS_H = (0.5, 0.75, 1, 1.5, 2, 3)
S = requests.Session()


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def tmag_from_gaia(g, bp_rp):
    c = np.clip(bp_rp, None, 3.0)
    return g - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473


def exofop(name):
    path = f"{CACHE}/{name}.csv"
    if not os.path.exists(path):
        r = S.get(EXOFOP[name], timeout=120); r.raise_for_status(); open(path, "wb").write(r.content)
    return pd.read_csv(path), path


def match_hosts(pool, cat, radius_arcsec=10.0):
    """Index into pool and separation for each catalogue row (loose radius: 100-pc dwarfs move several arcsec since J2000)."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    cp = SkyCoord(pool.ra.values * u.deg, pool.dec.values * u.deg)
    cc = SkyCoord(cat.RA.values, cat.Dec.values, unit=(u.hourangle, u.deg)) if cat.RA.dtype == object else \
        SkyCoord(cat.RA.values * u.deg, cat.Dec.values * u.deg)
    idx, sep, _ = cc.match_to_catalog_sky(cp)
    return idx, sep.arcsec < radius_arcsec


def cmd_pool(a):
    from astroquery.gaia import Gaia
    os.makedirs(CACHE, exist_ok=True); os.makedirs(OUT, exist_ok=True)
    raw = f"{CACHE}/gaia_md_north_100pc.csv"
    if not os.path.exists(raw):
        Gaia.ROW_LIMIT = -1
        Gaia.launch_job_async(GAIA_Q, dump_to_file=False).get_results().write(raw, overwrite=True)
    g = pd.read_csv(raw)
    g["tmag"] = tmag_from_gaia(g.gmag, g.bp_rp); g["dist_pc"] = 1000 / g.parallax
    summ = {"computed_utc": now(), "gaia_query": GAIA_Q.strip(), "n_gaia": len(g), "bins": {}}
    for name in EXOFOP:
        cat, path = exofop(name)
        idx, ok = match_hosts(g, cat)
        g[f"{name}_host"] = False; g.loc[g.index[idx[ok]], f"{name}_host"] = True
        summ[f"{name}_file"] = {"path": path, "rows": len(cat), "downloaded_utc": dt.datetime.fromtimestamp(os.path.getmtime(path), dt.timezone.utc).isoformat(timespec="seconds")}
    for lo, hi in [(10, 13.5), (13.5, 16), (16, 17.5)]:
        s = g[(g.tmag >= lo) & (g.tmag < hi)]
        summ["bins"][f"T{lo}-{hi}"] = {"n": len(s), "within_60pc": int((s.dist_pc < 60).sum()), "toi_hosts": int(s.toi_host.sum()),
                                       "ctoi_hosts": int(s.ctoi_host.sum())}
    p = g[(g.tmag >= 13.5) & (g.tmag < 16)].copy(); p.to_csv(POOL, index=False)
    summ["pool"] = {"n": len(p), "within_60pc": int((p.dist_pc < 60).sum()), "abs_b_lt_10": int((p.b.abs() < 10).sum()),
                    "T_gt_15p5": int((p.tmag > 15.5).sum()), "ruwe_lt_1p4": int((p.ruwe < 1.4).sum())}
    # the TOIs that sit on (or within 10" of) a pool star: these are the positive controls
    toi, _ = exofop("toi"); idx, ok = match_hosts(p, toi)
    m = toi[ok].copy(); m["source_id"] = p.source_id.values[idx[ok]]; m["pool_tmag"] = p.tmag.values[idx[ok]]; m["dist_pc"] = p.dist_pc.values[idx[ok]]
    m.to_csv(f"{OUT}/pool_toi_matches.csv", index=False)
    summ["pool"]["toi_matches"] = len(m)
    json.dump(summ, open(f"{OUT}/pool_summary.json", "w"), indent=1); print(json.dumps(summ["bins"]), json.dumps(summ["pool"]))


def tglc_url(gid, sec, cam, ccd):
    z = str(int(gid)).zfill(21)[:16]; d = "/".join(z[i:i + 4] for i in range(0, 16, 4))
    return (f"https://archive.stsci.edu/hlsps/tglc/s{sec:04d}/cam{cam}-ccd{ccd}/{d}/"
            f"hlsp_tglc_tess_ffi_gaiaid-{int(gid)}-s{sec:04d}-cam{cam}-ccd{ccd}_tess_v1_llc.fits")


def sectors_le(ra, dec, smax):
    from tess_stars2px import tess_stars2px_function_entry as t2p
    _, _, _, sec, cam, ccd, *_ = t2p(np.array([1]), np.array([ra]), np.array([dec]))
    return sorted({(int(s), int(c), int(d)) for s, c, d in zip(sec, cam, ccd) if s <= smax})


def products_near(ra, dec, radius="4 arcsec"):
    from astroquery.mast import Observations
    o = Observations.query_region(f"{ra} {dec}", radius=radius)
    return o[np.char.startswith(np.asarray(o["dataproduct_type"], dtype=str), "timeseries")]


def cmd_coverage(a):
    p = pd.read_csv(POOL).sample(n=a.n, random_state=a.seed)

    def one(r):
        o = products_near(r.ra, r.dec)
        nsec = lambda prov: len({int(s) for s, pr in zip(o["sequence_number"], o["provenance_name"]) if str(pr) == prov})
        secs = sectors_le(r.ra, r.dec, 55)
        tglc = [s for s, c, d in secs if S.head(tglc_url(r.source_id, s, c, d), timeout=30).status_code == 200]
        tic = next((str(t) for t, pr in zip(o["target_name"], o["provenance_name"]) if str(pr) == "TARS"), "")
        return dict(source_id=int(r.source_id), tic=tic, tmag=r.tmag, dist_pc=r.dist_pc, b=r.b, abs_g=r.abs_g, n_s1_55=len(secs), n_tglc=len(tglc),
                    **{f"n_{k.lower().replace('-', '')}": nsec(k) for k in ["TARS", "TESS-SPOC", "SPOC", "QLP"]})
    with ThreadPoolExecutor(6) as ex:
        d = pd.DataFrame(list(ex.map(one, [r for _, r in p.iterrows()])))
    d.to_csv(f"{OUT}/coverage_sample.csv", index=False)
    # the closest public search (zacharynapier/tess-mdwarf-search README, read 2026-10-19): T 13-15.5, |b| > 10,
    # >= 4 TESS-SPOC FFI sectors, no 2-min data. Teff/logg/contamination cuts not reproduced, so this over-counts its overlap.
    other = (d.tmag >= 13) & (d.tmag <= 15.5) & (d.b.abs() > 10) & (d.n_tessspoc >= 4) & (d.n_spoc == 0)
    comp = ~other & (d.n_spoc == 0) & (d.n_tars >= 2)
    s = {"computed_utc": now(), "n": len(d), "seed": a.seed,
         "frac_with": {k: round(float((d[f"n_{k}"] > 0).mean()), 3) for k in ["tars", "tessspoc", "spoc", "qlp", "tglc"]},
         "median_sectors": {k: float(d[f"n_{k}"].median()) for k in ["tars", "tessspoc", "spoc", "qlp", "tglc"]},
         "tglc_files_found_over_expected_s1_55": [int(d.n_tglc.sum()), int(d.n_s1_55.sum())],
         "in_napier_style_selection": int(other.sum()), "complement_no2min_ge2tars": int(comp.sum()),
         "complement_median_tars_sectors": float(d[comp].n_tars.median()), "complement_mean_tmag": round(float(d[comp].tmag.mean()), 2)}
    json.dump(s, open(f"{OUT}/coverage_summary.json", "w"), indent=1); print(json.dumps(s, indent=1))


def fetch_tars(tic):
    from astroquery.mast import Observations
    d = f"{CACHE}/tars/{tic}"
    if glob.glob(f"{d}/*_lc.fits"):
        return True
    obs = Observations.query_criteria(provenance_name="tars", target_name=str(tic))
    if len(obs) == 0:
        return False
    pl = Observations.get_product_list(obs); pl = pl[[str(f).endswith("_lc.fits") for f in pl["productFilename"]]]
    os.makedirs(d, exist_ok=True); Observations.download_products(pl, download_dir=d, flat=True)
    return True


def load_tars(tic):
    from astropy.io import fits
    T, F = [], []
    for f in sorted(glob.glob(f"{CACHE}/tars/{tic}/*_lc.fits")):
        d = fits.getdata(f, 1); t, x = np.asarray(d["TIME"], float), np.asarray(d["FLUX"], float)
        m = np.isfinite(t) & np.isfinite(x); t, x = t[m], x[m]
        if len(t) < 100:
            continue
        med = np.nanmedian(x); T.append(t); F.append(x / med if abs(med) > 1e-3 else 1 + x)
    return (np.concatenate(T), np.concatenate(F)) if T else (None, None)


def load_tglc(gid, ra, dec, col="cal_psf_flux"):
    from astropy.io import fits
    T, F = [], []
    for s, c, d in sectors_le(ra, dec, 55):
        r = S.get(tglc_url(gid, s, c, d), timeout=60)
        if r.status_code != 200:
            continue
        x = fits.open(io.BytesIO(r.content))[1].data
        m = (x["TGLC_flags"] == 0) & (x["TESS_flags"] == 0) & np.isfinite(x[col]) & np.isfinite(x["time"])
        if m.sum() < 100:
            continue
        f = np.asarray(x[col][m], float); T.append(np.asarray(x["time"][m], float)); F.append(f / np.nanmedian(f))
    return (np.concatenate(T), np.concatenate(F)) if T else (None, None)


def clean(t, f, bin_min=10.0):
    """biweight 0.5 d; trim 0.3 d after / 0.1 d before each gap > 0.5 d; 10-min bins; per-segment MAD noise; upward 4-sigma clip."""
    from wotan import flatten
    o = np.argsort(t); t, f = t[o], f[o]
    flat = flatten(t, f, method="biweight", window_length=0.5, break_tolerance=0.3)
    m = np.isfinite(flat); t, flat = t[m], flat[m]
    seg = np.concatenate([[0], np.cumsum(np.diff(t) > 0.5)])
    keep = np.ones(len(t), bool)
    for k in np.unique(seg):
        j = seg == k; keep[j] = (t[j] > t[j].min() + 0.3) & (t[j] < t[j].max() - 0.1)
    t, flat, seg = t[keep], flat[keep], seg[keep]
    tb, fb, eb, w = [], [], [], bin_min / 1440
    for k in np.unique(seg):
        j = seg == k
        if j.sum() < 50:
            continue
        b = np.floor((t[j] - t[j].min()) / w).astype(int); nb = np.bincount(b); ok = nb > 0
        tt = np.bincount(b, t[j])[ok] / nb[ok]; ff = np.bincount(b, flat[j])[ok] / nb[ok]
        sig = 1.4826 * np.median(np.abs(ff - np.median(ff))); m = ff < np.median(ff) + 4 * sig
        tb.append(tt[m]); fb.append(ff[m]); eb.append(np.full(m.sum(), sig))
    return np.concatenate(tb), np.concatenate(fb), np.concatenate(eb)


def search(t, f, refine=True):
    tc, fc, ec = clean(t, f)
    return bls_search(tc, fc, ec, pmin=0.5, pmax=15, durations_h=DURS_H, refine=refine, verbose=False), float(np.median(ec)), len(tc)


def depth_at(t, f, P, t0, dur_h):
    tc, fc, _ = clean(t, f)
    ph = ((tc - t0 + 0.5 * P) % P) - 0.5 * P
    inn = np.abs(ph) < 0.4 * dur_h / 24; out = (np.abs(ph) > dur_h / 24) & (np.abs(ph) < 3 * dur_h / 24)
    return float(1e6 * (np.median(fc[out]) - np.median(fc[inn])))


def same_period(p, q, tol=2e-3):
    return min(abs(p / q - h) for h in (1, 2, 0.5)) < tol


def cmd_controls(a):
    m = pd.read_csv(f"{OUT}/pool_toi_matches.csv"); p = pd.read_csv(POOL).set_index("source_id")
    # genuine pool members: the TOI's own TESS magnitude must agree with the pool estimate (a mismatch is a bright neighbour)
    m = m[(m["TESS Mag"] - m.pool_tmag).abs() < 0.6]
    rows = []
    for tic, grp in m.groupby("TIC ID"):
        fetch_tars(int(tic)); sid = int(grp.source_id.iloc[0]); g = p.loc[sid]
        lcs = {"tars": load_tars(int(tic)), "tglc": load_tglc(sid, g.ra, g.dec)}
        res = {k: (search(*v) if v[0] is not None else None) for k, v in lcs.items()}
        for _, k in grp.iterrows():
            P, t0, D = k["Period (days)"], k["Epoch (BJD)"] - 2457000, k["Duration (hours)"]
            row = dict(toi=k.TOI, tic=int(tic), tmag=k["TESS Mag"], disposition=k["TFOPWG Disposition"], period_d=P, depth_ppm_toi=k["Depth (ppm)"],
                       duration_h=D, source=k["Pipeline Signal ID"])
            for name, r in res.items():
                if r is None:
                    row[name] = None; continue
                bls, rms, n = r; t, f = lcs[name]
                row[name] = dict(n_points_10min=n, rms_10min_ppm=round(rms * 1e6), best_period_d=bls.period, sde=round(bls.sde, 2),
                                 depth_snr=round(bls.depth_snr, 2), recovered_as_top_peak=bool(same_period(bls.period, P)),
                                 depth_ppm_at_toi_ephemeris=round(depth_at(t, f, P, t0, D)))
            rows.append(row); print(json.dumps(row))
    json.dump({"computed_utc": now(), "durations_h": DURS_H, "period_range_d": [0.5, 15], "controls": rows}, open(f"{OUT}/controls.json", "w"), indent=1)


def _inject_one(arg):
    tic, absg, tmag, seed, radii = arg
    rng = np.random.default_rng(seed); t, f = load_tars(tic); R = float(np.interp(absg, MG, RS)); M = R; out = []
    for rp in [0.0] + list(radii):
        ff = f.copy(); P = float(np.exp(rng.uniform(np.log(0.5), np.log(10)))); t0 = rng.uniform(t.min(), t.min() + P)
        dep = T14 = 0.0
        if rp > 0:
            k = rp / 109.1 / R; a_rs = 215.03 * (M * (P / 365.25) ** 2) ** (1 / 3) / R; b = rng.uniform(0, 0.8)
            T14 = P / np.pi * np.arcsin(np.sqrt(max((1 + k) ** 2 - b * b, 1e-4)) / a_rs)
            dep = 0.9 * k * k  # box at ~the limb-darkened mean depth
            ph = ((t - t0 + 0.5 * P) % P) - 0.5 * P; ff = ff - dep * (np.abs(ph) < T14 / 2)
        bls, rms, n = search(t, ff, refine=False)
        out.append(dict(tic=tic, tmag=tmag, rstar=R, n_sectors=len(glob.glob(f"{CACHE}/tars/{tic}/*_lc.fits")), rms_10min_ppm=rms * 1e6,
                        rp_earth=rp, period_d=P, depth_ppm=dep * 1e6, T14_h=T14 * 24, best_period_d=bls.period, sde=bls.sde,
                        depth_snr=bls.depth_snr, period_match=bool(rp > 0 and same_period(bls.period, P))))
    return out


def cmd_inject(a):
    from astroquery.mast import Observations
    p = pd.read_csv(POOL).sample(n=4 * a.n, random_state=a.seed)

    def tic_of(r):
        try:
            o = products_near(r.ra, r.dec)
            tic = next((int(t) for t, pr in zip(o["target_name"], o["provenance_name"]) if str(pr) == "TARS"), None)
            return (tic, float(r.abs_g), float(r.tmag), int(r.source_id) % 100000) if tic and fetch_tars(tic) else None
        except Exception as ex:  # a failed lookup just drops the star
            print("skip", r.source_id, str(ex)[:80]); return None
    with ThreadPoolExecutor(6) as ex:
        stars = [s for s in ex.map(tic_of, [r for _, r in p.iterrows()]) if s][: a.n]
    print("stars:", len(stars), flush=True)
    radii = [1.0, 1.5, 2.0, 3.0, 4.0]
    with Pool(a.workers) as pool:
        res = [x for y in pool.map(_inject_one, [s + (radii,) for s in stars]) for x in y]
    json.dump({"computed_utc": now(), "seed": a.seed, "note": "box injected into the TARS light curve after CPM extraction: "
               "ignores CPM depth suppression, so recovery is an upper bound", "results": res}, open(f"{OUT}/inject.json", "w"), indent=1)
    d = pd.DataFrame(res); print(d.groupby("rp_earth")[["period_match", "sde", "depth_snr"]].mean())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["pool", "coverage", "controls", "inject"])
    ap.add_argument("--n", type=int, default=120); ap.add_argument("--seed", type=int, default=20261019); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(); os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)
    {"pool": cmd_pool, "coverage": cmd_coverage, "controls": cmd_controls, "inject": cmd_inject}[a.cmd](a)
