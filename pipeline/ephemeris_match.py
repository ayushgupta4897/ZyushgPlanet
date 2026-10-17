"""
Ephemeris / known-object matching for first-screen survivors.

Why: in Kepler, a periodic dip can be a copy of a neighbour's eclipse (PRF wings,
electronic cross-talk, column anomalies; Coughlin et al. 2014) or simply an object
somebody already catalogued.  Every survivor is checked against

  * KOI cumulative table (all dispositions, incl. false positives)      data/known/koi_cumulative.csv
  * Kepler DR25 threshold-crossing events (34k, incl. junk)              data/known/dr25_tce.csv
  * Kepler eclipsing-binary catalogue (Kirk et al. 2016, VizieR J/AJ/151/68)  data/known/J_AJ_151_68_catalog.csv
  * confirmed planets in the Kepler field (NASA Exoplanet Archive ps)   data/known/ps_keplerfield.csv

at period ratios 1, 2, 1/2, 3, 1/3 (tolerance `tol`, fractional) for every object
within `radius_arcmin`, and *always* against the originating target-pixel-file star
(tpforg) whatever its distance, because that star's photons are in the same pixels.

Usage
-----
    from pipeline.ephemeris_match import match_signal
    hits = match_signal(ra, dec, period_d, tpforg_kic=5689116, kic=5689121)
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KNOWN = os.path.join(ROOT, "data", "known")
_cache: dict[str, pd.DataFrame] = {}


def _load() -> pd.DataFrame:
    """One tidy table: catalog, name, kepid, ra, dec, period, t0_bkjd, disposition, depth."""
    if "all" in _cache:
        return _cache["all"]
    frames = []
    koi = pd.read_csv(os.path.join(KNOWN, "koi_cumulative.csv"))
    frames.append(pd.DataFrame({
        "catalog": "KOI", "name": koi.kepoi_name, "kepid": koi.kepid, "ra": koi.ra, "dec": koi.dec,
        "period": koi.koi_period, "t0_bkjd": koi.koi_time0bk, "disposition": koi.koi_disposition,
        "depth_ppm": koi.koi_depth}))
    tce = pd.read_csv(os.path.join(KNOWN, "dr25_tce.csv"))
    frames.append(pd.DataFrame({
        "catalog": "TCE", "name": tce.kepid.astype(str) + "-" + tce.tce_plnt_num.astype(str), "kepid": tce.kepid,
        "ra": tce.ra, "dec": tce.dec, "period": tce.tce_period, "t0_bkjd": tce.tce_time0bk,
        "disposition": "TCE", "depth_ppm": tce.tce_depth}))
    eb = pd.read_csv(os.path.join(KNOWN, "J_AJ_151_68_catalog.csv"))
    ra_col = "RAJ2000" if "RAJ2000" in eb.columns else None
    frames.append(pd.DataFrame({
        "catalog": "KeplerEB", "name": "KIC" + eb.KIC.astype(str), "kepid": eb.KIC,
        "ra": eb[ra_col] if ra_col else np.nan, "dec": eb["DEJ2000"] if ra_col else np.nan,
        "period": eb.Per, "t0_bkjd": eb.BJD0 if "BJD0" in eb.columns else np.nan,
        "disposition": "EB", "depth_ppm": np.nan}))
    ps = pd.read_csv(os.path.join(KNOWN, "ps_keplerfield.csv"))
    frames.append(pd.DataFrame({
        "catalog": "ConfirmedPlanet", "name": ps.pl_name, "kepid": np.nan, "ra": ps.ra, "dec": ps.dec,
        "period": ps.pl_orbper, "t0_bkjd": ps.pl_tranmid - 2454833.0, "disposition": "CONFIRMED",
        "depth_ppm": ps.pl_trandep * 1e4 if "pl_trandep" in ps.columns else np.nan}))
    allk = pd.concat(frames, ignore_index=True)
    allk = allk[np.isfinite(allk.period)]
    # fill EB coordinates from KOI/TCE kepid where missing
    kc = pd.concat([koi[["kepid", "ra", "dec"]], tce[["kepid", "ra", "dec"]]]).drop_duplicates("kepid").set_index("kepid")
    miss = allk.ra.isna() & allk.kepid.notna()
    allk.loc[miss, "ra"] = allk.loc[miss, "kepid"].map(kc.ra)
    allk.loc[miss, "dec"] = allk.loc[miss, "kepid"].map(kc.dec)
    _cache["all"] = allk
    return allk


def _sep_arcmin(ra1, dec1, ra2, dec2):
    d2r = np.pi / 180.0
    dra = (np.asarray(ra2) - ra1) * np.cos(dec1 * d2r)
    ddec = np.asarray(dec2) - dec1
    return np.sqrt(dra**2 + ddec**2) * 60.0


RATIOS = {"1:1": 1.0, "2:1": 2.0, "1:2": 0.5, "3:1": 3.0, "1:3": 1.0 / 3.0}


def match_signal(ra: float, dec: float, period_d: float, tpforg_kic: int | None = None, kic: int | None = None,
                 radius_arcmin: float = 3.0, tol: float = 0.002, t0_bkjd: float | None = None) -> pd.DataFrame:
    """Return every catalogue object that could explain this period.

    Columns: catalog, name, kepid, sep_arcmin, period, ratio, frac_diff, disposition,
    same_tpf (the originating TPF star), same_kic, epoch_offset_cycles (if t0 given).
    """
    allk = _load()
    sep = _sep_arcmin(ra, dec, allk.ra.values, allk.dec.values)
    near = (sep <= radius_arcmin) | (allk.kepid.values == (tpforg_kic or -1)) | (allk.kepid.values == (kic or -2))
    sub = allk[near].copy()
    sub["sep_arcmin"] = sep[near]
    rows = []
    for _, r in sub.iterrows():
        for lab, k in RATIOS.items():
            fd = abs(period_d - k * r.period) / (k * r.period)
            if fd <= tol:
                eo = np.nan
                if t0_bkjd is not None and np.isfinite(r.t0_bkjd):
                    eo = ((t0_bkjd - r.t0_bkjd) / period_d) % 1.0
                    eo = min(eo, 1 - eo)  # 0 = same epoch; 0.5 = anti-phase
                rows.append({**r.to_dict(), "ratio": lab, "frac_diff": fd,
                             "same_tpf": bool(r.kepid == tpforg_kic) if tpforg_kic else False,
                             "same_kic": bool(r.kepid == kic) if kic else False,
                             "epoch_offset_phase": eo})
    cols = ["catalog", "name", "kepid", "sep_arcmin", "period", "ratio", "frac_diff", "disposition",
            "depth_ppm", "same_tpf", "same_kic", "epoch_offset_phase"]
    return pd.DataFrame(rows, columns=cols).sort_values(["same_tpf", "same_kic", "sep_arcmin"], ascending=[False, False, True]) if rows else pd.DataFrame(columns=cols)


if __name__ == "__main__":
    # self-test on Kuber's candidate: expect Kepler-169 e (KOI-505.01) at 4.4 arcmin, 1:1, epoch offset ~0.21
    h = match_signal(285.907827, 40.941988, 13.769448, tpforg_kic=5689116, kic=5689121, radius_arcmin=5.0,
                     t0_bkjd=2455630.245431 - 2454833.0)
    print(h.to_string())
