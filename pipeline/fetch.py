"""Locate, download and describe Kepler Bonus background-source (KBONUS-BKG) light curves.

Data product facts (verified against the live MAST products on 2026-10-17):
* HLSP 'KBONUS-BKG' (Martinez-Palomera, Hedges & Dotson 2023, AJ 166, 265; arXiv:2310.17733;
  DOI 10.17909/7jbr-w430; https://archive.stsci.edu/hlsp/kbonus-bkg).
* MAST serves ONE multi-extension FITS file per source for all quarters, e.g.
  mast:HLSP/kbonus-bkg/lcs/0056/005689121/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits
  (sources with no KIC id are named gaia-dr3-<id>). The CAOM 'target_name' is the Gaia DR3 designation.
* HDUs: PRIMARY (no data; header carries KEPLERID, GAIAID, TPFORG (the TPF the source was extracted
  from), RA_OBJ/DEC_OBJ, GMAG/BPMAG/RPMAG/KEPMAG, PARALLAX, XPOSURE=1765.464 s, TIMESYS='TDB',
  QDETECT = per-quarter detection string Q0..Q17),
  LIGHTCURVE_STITCHED (all quarters; columns TIME, CADENCENO, QUARTER, FLUX, FLUX_ERR, SAP_FLUX,
  SAP_FLUX_ERR, PSF_FLAT_FLUX, PSF_FLAT_FLUX_ERR, SAP_QUALITY, FLATTEN_MASK),
  LIGHTCURVE_Q<n> one per detected quarter (columns CADENCENO, TIME, FLUX, FLUX_ERR, SAP_FLUX,
  SAP_FLUX_ERR, CENTROID_COLUMN, CENTROID_ROW, SAP_QUALITY, SAP_BKG, RED_CHI2, PSF_FLUX_NOVA,
  PSF_FLUX_ERR_NOVA; header keys QUARTER, SEASON, MODULE, OUTPUT, CHANNEL, FLFRCSAP, CROWDSAP,
  NPIXSAP, PSFFRAC, PERTRATI, PERTSTD, ITERNEG, SAP='optimal'),
  APERTURE_Q<n> (uint8 image = SAP aperture mask in the originating TPF frame).
* FLUX / FLUX_ERR = the corrected PSF (psfmachine 'Linearized Field Deblending') photometry.
  SAP_FLUX = simple aperture photometry in the APERTURE_Q<n> mask (diluted by neighbours).
  PSF_FLAT_FLUX = PSF flux after the authors' own flattening; PSF_FLUX_NOVA = PSF flux from the
  mean (un-perturbed) shape model. We search FLUX.
* SAP_QUALITY = the Kepler TPF QUALITY bitmask copied from the originating TPF (values seen for
  KIC 5689121: 0, 16, 128, 144, 1024, 2048, 2176, 8192, 262144, 393216, ...). FLATTEN_MASK is the
  authors' flattening mask (0/1), only in the stitched extension.
* TIME column unit is 'BJD - 2454833' (BKJD) on the TDB scale. BJD_TDB = TIME + 2454833.0.
  lightkurve.read() parses these files natively as KeplerLightCurve(time.format='bkjd', scale='tdb').
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.io import fits

BKJD_OFFSET = 2454833.0  # BJD_TDB = BKJD + BKJD_OFFSET  (Kepler Archive Manual, and TUNIT 'BJD - 2454833')
AUTHOR = "KBONUS-BKG"
REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT / "data" / "kbonus"

PRIMARY_KEYS = ["OBJECT", "KEPLERID", "GAIAID", "TPFORG", "RA_OBJ", "DEC_OBJ", "KEPMAG", "GMAG", "BPMAG",
                "RPMAG", "PARALLAX", "PMRA", "PMDEC", "XPOSURE", "TIMESYS", "MJD-BEG", "MJD-END", "QDETECT",
                "HLSPVER", "VERSION", "PROCVER", "DOI", "APERTURE", "DATE"]
QUARTER_KEYS = ["QUARTER", "SEASON", "MODULE", "OUTPUT", "CHANNEL", "FLFRCSAP", "CROWDSAP", "NPIXSAP",
                "PSFFRAC", "PERTRATI", "PERTSTD", "ITERNEG", "SAP", "ROW", "COLUMN", "TPF_ROW", "TPF_COL"]


@dataclass
class Quarter:
    quarter: int
    time: np.ndarray            # BKJD (BJD_TDB - 2454833)
    flux: np.ndarray            # PSF flux, e-/s
    flux_err: np.ndarray
    sap_flux: np.ndarray
    sap_flux_err: np.ndarray
    quality: np.ndarray         # SAP_QUALITY bitmask (Kepler TPF QUALITY)
    centroid_col: np.ndarray
    centroid_row: np.ndarray
    header: dict = field(default_factory=dict)


def find_products(gaia_id: Optional[str] = None, ra: Optional[float] = None, dec: Optional[float] = None,
                  radius_arcsec: float = 20.0):
    """Return an astropy Table of KBONUS-BKG CAOM observations for a Gaia DR3 id and/or a sky position.

    gaia_id -> astroquery.mast.Observations.query_criteria(provenance_name='KBONUS-BKG',
               target_name='Gaia DR3 <id>')  (exact target-name match, no coordinates needed)
    ra/dec  -> lightkurve.search_lightcurve(SkyCoord, author='KBONUS-BKG', radius=...) cone search.
    """
    from astroquery.mast import Observations
    tabs = []
    if gaia_id is not None:
        name = f"Gaia DR3 {str(gaia_id).split()[-1]}"
        t = Observations.query_criteria(provenance_name=AUTHOR, target_name=name)
        if len(t):
            tabs.append(t)
    if not tabs and ra is not None and dec is not None:
        import astropy.units as u
        import lightkurve as lk
        from astropy.coordinates import SkyCoord
        sr = lk.search_lightcurve(SkyCoord(ra, dec, unit="deg"), author=AUTHOR, radius=radius_arcsec * u.arcsec)
        if len(sr):
            tabs.append(sr.table)
    if not tabs:
        raise RuntimeError(f"No {AUTHOR} product found for gaia_id={gaia_id} ra={ra} dec={dec}")
    return tabs[0]


def pick_row(table, gaia_id: Optional[str] = None):
    """Pick the row for the requested Gaia source (or the nearest if only coordinates were given)."""
    names = [str(x) for x in table["target_name"]]
    if gaia_id is not None:
        want = f"Gaia DR3 {str(gaia_id).split()[-1]}"
        for i, n in enumerate(names):
            if n == want:
                return table[i]
        raise RuntimeError(f"{want} not among returned target_names: {names}")
    if "distance" in table.colnames:
        return table[int(np.argmin(np.asarray(table["distance"], float)))]
    return table[0]


def download(gaia_id: Optional[str] = None, ra: Optional[float] = None, dec: Optional[float] = None,
             data_root: Path = DATA_ROOT, verbose: bool = True) -> Path:
    """Download (or reuse cached) the single all-quarters KBONUS-BKG FITS file for one source.

    Cache location: <repo>/data/kbonus/<gaia_source_id>/<original MAST filename>
    """
    from astroquery.mast import Observations
    table = find_products(gaia_id, ra, dec)
    row = pick_row(table, gaia_id)
    target = str(row["target_name"])
    sid = target.split()[-1]
    uri = str(row["dataURL"])
    dest = Path(data_root) / sid / uri.split("/")[-1]
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        if verbose:
            print(f"[fetch] cached {dest}")
    else:
        t0 = time.time()
        status = Observations.download_file(uri, local_path=str(dest))
        if verbose:
            print(f"[fetch] downloaded {uri} -> {dest} ({dest.stat().st_size} bytes, {time.time()-t0:.1f}s, status={status[0]})")
        with open(dest, "rb") as fh:
            if fh.read(6) != b"SIMPLE":
                dest.unlink()
                raise RuntimeError("Downloaded file is not FITS")
    meta = {"target_name": target, "gaia_source_id": sid, "dataURL": uri, "file": str(dest),
            "s_ra": float(row["s_ra"]), "s_dec": float(row["s_dec"]),
            "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (dest.parent / "manifest.json").write_text(json.dumps(meta, indent=2))
    return dest


def describe_fits(path: Path, verbose: bool = True) -> dict:
    """Print and return the HDU / column / header structure of a KBONUS-BKG light-curve file."""
    out = {"file": str(path), "hdus": [], "primary": {}, "quarters": [], "stitched_columns": [],
           "quarter_columns": [], "quality_values_seen": []}
    with fits.open(path, memmap=False) as h:
        for i, hd in enumerate(h):
            shape = getattr(hd, "shape", None)
            n = len(hd.data) if isinstance(hd, fits.BinTableHDU) else (list(shape) if shape else None)
            out["hdus"].append({"index": i, "name": hd.name, "type": type(hd).__name__, "rows_or_shape": n})
        ph = h[0].header
        out["primary"] = {k: (ph[k] if not isinstance(ph[k], (np.floating, np.integer)) else float(ph[k]))
                          for k in PRIMARY_KEYS if k in ph}
        st = h["LIGHTCURVE_STITCHED"]
        out["stitched_columns"] = [{"name": c.name, "format": c.format, "unit": c.unit} for c in st.columns]
        out["time_unit"] = st.columns["TIME"].unit
        out["quality_values_seen"] = [int(v) for v in np.unique(st.data["SAP_QUALITY"])]
        out["n_rows_stitched"] = int(len(st.data))
        for hd in h:
            if hd.name.startswith("LIGHTCURVE_Q"):
                if not out["quarter_columns"]:
                    out["quarter_columns"] = [{"name": c.name, "format": c.format, "unit": c.unit} for c in hd.columns]
                q = {k: hd.header.get(k) for k in QUARTER_KEYS if k in hd.header}
                d = hd.data
                q.update(n=int(len(d)), n_quality0=int((d["SAP_QUALITY"] == 0).sum()),
                         med_psf_flux=float(np.nanmedian(d["FLUX"])), med_sap_flux=float(np.nanmedian(d["SAP_FLUX"])),
                         t_min_bkjd=float(np.nanmin(d["TIME"])), t_max_bkjd=float(np.nanmax(d["TIME"])))
                ap = h[f"APERTURE_Q{hd.header['QUARTER']}"]
                q["aperture_shape"] = list(ap.data.shape)
                q["aperture_npix"] = int(ap.data.sum())
                out["quarters"].append(q)
    if verbose:
        print(f"[fetch] {path}")
        for hd in out["hdus"]:
            print(f"   HDU {hd['index']:2d} {hd['name']:22s} {hd['type']:14s} {hd['rows_or_shape']}")
        print("   PRIMARY:", json.dumps(out["primary"], default=str))
        print("   TIME unit:", out["time_unit"], "| TIMESYS:", out["primary"].get("TIMESYS"),
              "| => BJD_TDB = TIME +", BKJD_OFFSET)
        print("   stitched columns:", [c["name"] for c in out["stitched_columns"]])
        print("   per-quarter columns:", [c["name"] for c in out["quarter_columns"]])
        print("   SAP_QUALITY values seen:", out["quality_values_seen"])
        print("   Q   n   q0   PSFFRAC CROWDSAP FLFRCSAP NPIX  medPSF  medSAP  tmin  tmax")
        for q in out["quarters"]:
            print(f"   {q['QUARTER']:2d} {q['n']:5d} {q['n_quality0']:5d} {q['PSFFRAC']:.3f}  {q['CROWDSAP']:.3f}   "
                  f"{q['FLFRCSAP']:.3f}  {q['NPIXSAP']:2d}  {q['med_psf_flux']:7.1f} {q['med_sap_flux']:7.1f} "
                  f"{q['t_min_bkjd']:7.1f} {q['t_max_bkjd']:7.1f}")
    return out


def load_quarters(path: Path, min_psffrac: float = 0.0) -> tuple[list[Quarter], dict]:
    """Read every LIGHTCURVE_Q<n> extension into Quarter objects (PSF flux = FLUX).

    min_psffrac: drop quarters whose header PSFFRAC (fraction of the PSF inside the data) is below this.
    Returns (quarters, primary_header_dict).
    """
    quarters = []
    with fits.open(path, memmap=False) as h:
        prim = {k: h[0].header.get(k) for k in PRIMARY_KEYS if k in h[0].header}
        for hd in h:
            if not hd.name.startswith("LIGHTCURVE_Q"):
                continue
            H = hd.header
            if float(H.get("PSFFRAC", 1.0)) < min_psffrac:
                continue
            d = hd.data
            quarters.append(Quarter(
                quarter=int(H["QUARTER"]),
                time=np.asarray(d["TIME"], float), flux=np.asarray(d["FLUX"], float),
                flux_err=np.asarray(d["FLUX_ERR"], float), sap_flux=np.asarray(d["SAP_FLUX"], float),
                sap_flux_err=np.asarray(d["SAP_FLUX_ERR"], float), quality=np.asarray(d["SAP_QUALITY"], int),
                centroid_col=np.asarray(d["CENTROID_COLUMN"], float), centroid_row=np.asarray(d["CENTROID_ROW"], float),
                header={k: H.get(k) for k in QUARTER_KEYS if k in H}))
    quarters.sort(key=lambda q: q.quarter)
    return quarters, prim


if __name__ == "__main__":  # python -m pipeline.fetch --gaia 2103792461826460288
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gaia", help="Gaia DR3 source id")
    ap.add_argument("--ra", type=float); ap.add_argument("--dec", type=float)
    a = ap.parse_args()
    p = download(a.gaia, a.ra, a.dec)
    describe_fits(p)
