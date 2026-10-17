"""Run MOLUSC (Wood et al. 2021, AJ 162, 128; github.com/woodml/MOLUSC) headless to constrain unseen bound
companions of a target with Gaia DR3/EDR3 RUWE and Gaia imaging completeness, and write a companion file that
TRICERATOPS accepts as `molusc_file`.

    python scripts/molusc_run.py <fname> --mass 0.73 [--age 5] [--n-trice 1000000] [--fmult 0.46] \
        [--molusc-dir <path to MOLUSC/code>] [--out runs/poolA_main/molusc/<fname>]

MOLUSC is not modified. Two shims make it run in a headless, current-SciPy environment:
  * a stub `tkinter` module (MOLUSC imports it at the top for its GUI; the command-line path never touches it);
  * `Application.star_jitter = -1`, an attribute MOLUSC's command-line path reads but never sets (-1 is its own
    "no stellar-jitter test" value);
  * `scipy.interpolate.interp2d` (removed in SciPy 1.14) replaced by a bilinear RegularGridInterpolator with
    nearest-neighbour extrapolation, which is what interp2d(kind='linear', fill_value=None) did on a rectangular grid.

Normalisation. TRICERATOPS draws N trial systems; with a MOLUSC file it takes the surviving companions (periastron
> 10 AU) as the companion draws and pads the rest with "no companion". MOLUSC generates companions assuming the star
HAS one. So we generate n = fmult x N companions, where fmult is the stellar multiplicity fraction (default 0.46,
Raghavan et al. 2010 for solar-type stars; conservative for K and M dwarfs, whose fractions are lower). The result
is a companion prior equal to fmult x (fraction of companions Gaia would not have revealed).
"""
import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import types

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("fname")
ap.add_argument("--mass", type=float, required=True)
ap.add_argument("--age", type=float, default=5.0)
ap.add_argument("--n-trice", type=int, default=1_000_000)
ap.add_argument("--fmult", type=float, default=0.46)
ap.add_argument("--molusc-dir", default="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/tools/MOLUSC/code")
ap.add_argument("--targets", default="targets/poolA_kbonus_bg_dwarfs_ranked.csv")
ap.add_argument("--no-ruwe", action="store_true")
ap.add_argument("--no-gaia", action="store_true")
ap.add_argument("--seed", type=int, default=7)
ap.add_argument("--gaia-limit", type=float, default=20.0, help="MOLUSC Gaia completeness limit: 18 (its conservative default) or 20 (its 'more complete' GUI option)")
ap.add_argument("--out", default=None)
a = ap.parse_args()

import pandas as pd  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
import astropy.units as u  # noqa: E402

repo = os.getcwd()
row = pd.read_csv(a.targets).set_index("fname").loc[a.fname]
c = SkyCoord(float(row.ra) * u.deg, float(row.dec) * u.deg)
ra_s = c.ra.to_string(unit=u.hourangle, sep="hms", precision=3)
dec_s = c.dec.to_string(unit=u.deg, sep="dms", precision=2, alwayssign=True)
out = os.path.abspath(a.out or f"runs/poolA_main/molusc/{a.fname}")
os.makedirs(out, exist_ok=True)
n_gen = int(round(a.fmult * a.n_trice))

# ---- shims (no change to MOLUSC's own code)
tk = types.ModuleType("tkinter")
tk.Frame = object
tk.messagebox = types.ModuleType("tkinter.messagebox")
sys.modules["tkinter"] = tk
sys.modules["tkinter.messagebox"] = tk.messagebox

import scipy.interpolate as _si  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402


class _One(np.ndarray):
    """1-element array that still converts with float(), as interp2d's output did under NumPy 1.x (NumPy 2 refuses
    float() on a shape-(1,) array; MOLUSC both concatenates these outputs and calls float() on them)."""

    def __float__(self):
        return float(np.asarray(self).reshape(-1)[0])


class _Interp2dShim:
    """interp2d(x, y, z, kind='linear') on a rectangular grid, z shaped (len(y), len(x)); nearest-neighbour outside."""

    def __init__(self, x, y, z, kind="linear", **kw):
        x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
        ox, oy = np.argsort(x), np.argsort(y)
        self.x, self.y = x[ox], y[oy]
        self.f = RegularGridInterpolator((self.y, self.x), z[np.ix_(oy, ox)], method="linear", bounds_error=False, fill_value=None)

    def __call__(self, xi, yi):
        xi = np.clip(np.atleast_1d(np.asarray(xi, float)), self.x[0], self.x[-1])
        yi = np.clip(np.atleast_1d(np.asarray(yi, float)), self.y[0], self.y[-1])
        if xi.size == 1 and yi.size == 1:
            return np.atleast_1d(self.f([[yi[0], xi[0]]])).view(_One)
        Y, X = np.meshgrid(yi, xi, indexing="ij")
        return self.f(np.c_[Y.ravel(), X.ravel()]).reshape(Y.shape)


_si.interp2d = _Interp2dShim

# ---- run MOLUSC's own command-line application from its directory (it reads its model tables by relative path)
np.random.seed(a.seed)
sys.path.insert(0, a.molusc_dir)
os.chdir(a.molusc_dir)
import BinaryStarGUI as M  # noqa: E402

prefix = os.path.join(out, "molusc")
argv = ["BinaryStarGUI.py", "-v", "--age", str(a.age)]
if not a.no_ruwe:
    argv.append("--ruwe")
if not a.no_gaia:
    argv.append("--gaia")
argv += [prefix, ra_s, "--", dec_s, str(n_gen), str(a.mass)] if dec_s.startswith("-") else [prefix, ra_s, dec_s, str(n_gen), str(a.mass)]
sys.argv = argv
app = M.Application(argv)
app.star_jitter = -1  # MOLUSC's CLI never sets this attribute; -1 is its own "no stellar-jitter test" value
app.gaia_limit = a.gaia_limit  # MOLUSC exposes 18 or 20 in its GUI only
app.start()
os.chdir(repo)

kept = pd.read_csv(prefix + "_kept.csv")
peri = kept["semi-major axis(AU)"] * (1 - kept["eccentricity"])
n_wide = int((peri > 10).sum())
commit = subprocess.run(["git", "-C", a.molusc_dir, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
summary = {
    "fname": a.fname, "utc": dt.datetime.utcnow().isoformat() + "Z", "molusc_repo": "https://github.com/woodml/MOLUSC", "molusc_commit": commit,
    "ra_hms": ra_s, "dec_dms": dec_s, "mass_msun": a.mass, "age_gyr": a.age, "gaia_limit_G": a.gaia_limit, "tests": [t for t in ("ruwe", "gaia") if not getattr(a, "no_" + t)],
    "fmult": a.fmult, "n_trice": a.n_trice, "n_generated": n_gen, "n_kept": int(len(kept)), "n_kept_periastron_gt_10au": n_wide,
    "fraction_kept": len(kept) / n_gen, "companion_prior_passed_to_triceratops": n_wide / a.n_trice,
    "target_ruwe": float(np.exp(getattr(app, "ruwe_ln", np.nan))) if hasattr(app, "ruwe_ln") else None,
    "shims": ["stub tkinter (GUI not used)", "Application.star_jitter = -1 (uninitialised in MOLUSC's CLI; -1 = no jitter test)", "scipy.interpolate.interp2d -> bilinear RegularGridInterpolator with nearest-neighbour extrapolation"],
    "kept_file": prefix + "_kept.csv",
}
json.dump(summary, open(os.path.join(out, "molusc_summary.json"), "w"), indent=1)
print(json.dumps({k: summary[k] for k in ("fname", "n_generated", "n_kept", "n_kept_periastron_gt_10au", "fraction_kept", "companion_prior_passed_to_triceratops")}))
