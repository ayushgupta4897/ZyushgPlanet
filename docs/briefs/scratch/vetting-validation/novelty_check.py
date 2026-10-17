"""novelty_check.py -- catalogue-novelty + ephemeris-match check for a transit candidate.

Verified live 2026-10-17 from the ZyushgPlanet sandbox (see brief.md, section 4).
Usage:  python -I novelty_check.py RA_DEG DEC_DEG PERIOD_D T0_BJD [RADIUS_ARCMIN]
Notes:
  * NASA Exoplanet Archive TAP: the IPAC WAF returned HTTP 403 for ADQL
    `contains(point(...),circle(...))` from this sandbox's proxy, while plain
    `ra between .. and dec between ..` box queries succeeded; the legacy
    nstedAPI also works. The code uses boxes and falls back to the legacy API.
  * Kepler EB catalogue (Villanova) needs a CSRF-protected POST; it has GLon/GLat only.
  * All catalogue absences are provisional: snapshot date matters (Kuber, planethunt.md).
"""
import sys, io, re, time, requests, numpy as np, pandas as pd
from scipy.special import erfcinv
from astropy.coordinates import SkyCoord
import astropy.units as u

TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
API = "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI"
H = {"User-Agent": "Mozilla/5.0 ZyushgPlanet-novelty/0.1"}

def tap(sql, retries=3):
    for i in range(retries):
        r = requests.get(TAP, params=dict(query=sql, format="csv"), headers=H, timeout=180)
        if r.status_code == 200:
            return pd.read_csv(io.StringIO(r.text))
        time.sleep(10 * (i + 1))
    raise RuntimeError(f"TAP failed: {r.status_code} {r.text[:120]}")

def legacy(table, select, where):
    r = requests.get(API, params={"table": table, "select": select, "format": "csv", "where": where}, headers=H, timeout=180)
    r.raise_for_status(); return pd.read_csv(io.StringIO(r.text))

def box(ra, dec, rad_deg):
    dra = rad_deg / np.cos(np.radians(dec))
    return f"ra between {ra-dra} and {ra+dra} and dec between {dec-rad_deg} and {dec+rad_deg}"

def sigma_match(Pa, Ta, Pb, Tb):
    """Coughlin et al. 2014 (arXiv:1401.1240) ephemeris-match statistics.
    Returns (sigma_P, sigma_T); a match needs sigma_P >= 3.5 AND sigma_T >= 2.0.
    Pa, Ta = candidate period/epoch; Pb, Tb = comparison object (same time system!)."""
    dP = (Pa - Pb) / Pa; dT = (Ta - Tb) / Pa
    dPp = abs(dP - np.rint(dP)); dTp = abs(dT - np.rint(dT))
    return np.sqrt(2) * erfcinv(max(dPp, 1e-300)), np.sqrt(2) * erfcinv(max(dTp, 1e-300))

RATIOS = (1, 2, 0.5, 3, 1/3, 4, 0.25, 5, 0.2)

def ephem_table(P, T0, df, pcol, tcol, idcols, ratios=RATIOS):
    rows = []
    for _, r in df.iterrows():
        if not (np.isfinite(r[pcol]) and np.isfinite(r[tcol])): continue
        for q in ratios:
            sP, sT = sigma_match(P, T0, r[pcol] * q, r[tcol])
            rows.append({**{c: r[c] for c in idcols}, "P_other": r[pcol], "ratio": round(q, 3), "sigP": round(sP, 2), "sigT": round(sT, 2),
                         "match_C14": bool(sP >= 3.5 and sT >= 2.0)})
    return pd.DataFrame(rows).sort_values("sigP", ascending=False) if rows else pd.DataFrame()

def main(ra, dec, P, T0_bjd, rad_arcmin=10.0):
    rad = rad_arcmin / 60.0
    T0_bkjd = T0_bjd - 2454833.0          # Kepler tables use BKJD
    T0_btjd = T0_bjd - 2457000.0          # TESS tables use BTJD; ExoFOP/TOI table give full BJD
    out = {}
    # --- Kepler KOIs (cumulative) ---
    try:
        koi = tap("select kepid,kepoi_name,koi_disposition,koi_pdisposition,koi_score,koi_period,koi_time0bk,koi_depth,koi_duration,"
                  "koi_fpflag_nt,koi_fpflag_ss,koi_fpflag_co,koi_fpflag_ec,koi_comment,ra,dec from cumulative where " + box(ra, dec, rad))
    except Exception:
        koi = legacy("cumulative", "kepid,kepoi_name,koi_disposition,koi_pdisposition,koi_score,koi_period,koi_time0bk,koi_depth,koi_duration,koi_fpflag_nt,koi_fpflag_ss,koi_fpflag_co,koi_fpflag_ec,koi_comment,ra,dec", box(ra, dec, rad))
    out["koi"] = koi; out["koi_ephem"] = ephem_table(P, T0_bkjd, koi, "koi_period", "koi_time0bk", ["kepid", "kepoi_name", "koi_disposition"])
    time.sleep(2)
    # --- Kepler DR25 TCEs (includes FPs never promoted to KOI) ---
    try:
        tce = tap("select kepid,tce_plnt_num,tce_period,tce_time0bk,tce_depth,tce_duration,tce_max_mult_ev,ra,dec from q1_q17_dr25_tce where " + box(ra, dec, rad))
    except Exception:
        tce = legacy("q1_q17_dr25_tce", "kepid,tce_plnt_num,tce_period,tce_time0bk,tce_depth,tce_duration,tce_max_mult_ev,ra,dec", box(ra, dec, rad))
    out["tce"] = tce; out["tce_ephem"] = ephem_table(P, T0_bkjd, tce, "tce_period", "tce_time0bk", ["kepid", "tce_plnt_num"])
    time.sleep(2)
    # --- TOIs (TAP 'toi' table; pl_tranmid is full BJD) ---
    toi = tap("select tid,toi,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandep,pl_trandurh,ra,dec from toi where " + box(ra, dec, rad))
    out["toi"] = toi; out["toi_ephem"] = ephem_table(P, T0_bjd, toi, "pl_orbper", "pl_tranmid", ["tid", "toi", "tfopwg_disp"])
    time.sleep(2)
    # --- Confirmed planets (pscomppars) and K2 ---
    out["ps"] = tap("select pl_name,hostname,tic_id,gaia_dr3_id,pl_orbper,pl_tranmid,pl_rade,disc_facility,ra,dec from pscomppars where " + box(ra, dec, rad))
    time.sleep(2)
    out["k2"] = tap("select pl_name,epic_candname,disposition,pl_orbper,pl_tranmid,ra,dec from k2pandc where " + box(ra, dec, rad))
    # --- ExoFOP TOI + CTOI full lists (CSV) ---
    toi_csv = pd.read_csv(io.StringIO(requests.get("https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv", headers=H, timeout=300).text))
    ctoi_csv = pd.read_csv(io.StringIO(requests.get("https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv", headers=H, timeout=300).text))
    for name, df in (("exofop_toi", toi_csv), ("exofop_ctoi", ctoi_csv)):
        c = SkyCoord(df["RA"].values, df["Dec"].values, unit=(u.hourangle if ":" in str(df["RA"].iloc[0]) else u.deg, u.deg))
        sep = c.separation(SkyCoord(ra * u.deg, dec * u.deg)).arcmin
        out[name] = df[sep < rad_arcmin].assign(sep_arcmin=sep[sep < rad_arcmin])
    # --- Villanova Kepler EB catalogue (CSRF POST; catalogue has GLon/GLat, bjd0 = BJD-2400000) ---
    s = requests.Session(); home = s.get("https://keplerebs.villanova.edu/", headers=H, timeout=120).text
    tok = re.search(r"name='csrfmiddlewaretoken' value='([^']+)'", home).group(1)
    r = s.post("https://keplerebs.villanova.edu/", headers={**H, "Referer": "https://keplerebs.villanova.edu/"},
               files={"csrfmiddlewaretoken": (None, tok), "format": (None, "file_csv")}, timeout=300)
    lines = [l.lstrip("#") if l.startswith("#KIC") else l for l in r.text.splitlines() if not l.startswith("##")]
    eb = pd.read_csv(io.StringIO("\n".join(lines))); eb = eb.loc[:, ~eb.columns.str.startswith("Unnamed")]
    g = SkyCoord(ra * u.deg, dec * u.deg).galactic
    ebc = SkyCoord(l=eb.GLon.values * u.deg, b=eb.GLat.values * u.deg, frame="galactic")
    eb["sep_arcmin"] = ebc.separation(g).arcmin
    out["kepler_eb"] = eb[eb.sep_arcmin < max(rad_arcmin, 60)]
    out["kepler_eb_ephem"] = ephem_table(P, T0_bjd - 2400000.0, out["kepler_eb"], "period", "bjd0", ["KIC", "sep_arcmin"])
    # --- VizieR: Kirk+2016 Kepler EBs (J/AJ/151/68), Prsa+2022 TESS EBs (J/ApJS/258/16) ---
    from astroquery.vizier import Vizier
    v = Vizier(columns=["**"], row_limit=-1)
    c0 = SkyCoord(ra * u.deg, dec * u.deg)
    for key, cat in (("vizier_kepler_eb", "J/AJ/151/68"), ("vizier_tess_eb", "J/ApJS/258/16")):
        res = v.query_region(c0, radius=rad_arcmin * u.arcmin, catalog=cat)
        out[key] = res[0].to_pandas() if res else pd.DataFrame()
    # --- Gaia DR3 variability / NSS for the exact source (if Gaia ID known, do by source_id) ---
    from astroquery.gaia import Gaia
    job = Gaia.launch_job(f"select source_id,ruwe,astrometric_excess_noise,astrometric_excess_noise_sig,ipd_frac_multi_peak,non_single_star,phot_variable_flag,phot_g_mean_mag "
                          f"from gaiadr3.gaia_source where 1=contains(point('ICRS',ra,dec),circle('ICRS',{ra},{dec},{3/3600.}))")
    out["gaia_source"] = job.get_results().to_pandas()
    if len(out["gaia_source"]):
        sid = int(out["gaia_source"].source_id.iloc[0])
        for tbl in ("vari_summary", "vari_classifier_result", "vari_eclipsing_binary", "vari_planetary_transit", "nss_two_body_orbit", "nss_acceleration_astro", "nss_non_linear_spectro"):
            out["gaia_" + tbl] = Gaia.launch_job(f"select * from gaiadr3.{tbl} where source_id={sid}").get_results().to_pandas()
    # --- SIMBAD ---
    from astroquery.simbad import Simbad
    sb = Simbad(); sb.add_votable_fields("otype", "ids")
    out["simbad"] = sb.query_region(c0, radius=15 * u.arcsec)
    return out

if __name__ == "__main__":
    ra, dec, P, T0 = map(float, sys.argv[1:5]); rad = float(sys.argv[5]) if len(sys.argv) > 5 else 10.0
    res = main(ra, dec, P, T0, rad)
    for k, v in res.items():
        n = len(v) if v is not None else 0
        print(f"== {k}: {n} rows")
        if n and hasattr(v, "head"): print(v.head(12).to_string(index=False))
    for k in ("koi_ephem", "tce_ephem", "toi_ephem", "kepler_eb_ephem"):
        if len(res[k]): print(k, "any Coughlin-2014 match:", bool(res[k].match_C14.any()))
