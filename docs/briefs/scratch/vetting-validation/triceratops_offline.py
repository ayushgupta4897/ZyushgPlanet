"""triceratops_offline.py -- run TRICERATOPS 1.1.0 WITHOUT target.__init__ (no MAST/Gaia/TRILEGAL at run time).

This is the pattern Kuber Mehta used for KIC 5689121 (analysis/validation-20261011/scripts/triceratops/run_kepler.py)
and the pattern you need whenever (a) the archive is unreachable, (b) the mission is Kepler/KBONUS and you want to
supply your own folded light curve, or (c) you want full control of the field-star population.
Verified 2026-10-17: triceratops 1.1.0 (PyPI 2026-09-02), pytransit 2.2.0, Python 3.13; N=20000 draws -> 17 s on 1 CPU
(default N=1e6 scales ~linearly, i.e. ~15 min). The synthetic field population below is a DEMO ONLY: build the real one
from a Gaia DR3 cone (triceratops.funcs.query_gaia_background does this for you when you CAN use target()), or from
TRILEGAL, resampled to the 0.1 deg^2 the priors assume.

Required attributes when bypassing __init__ (from inspecting calc_probs in 1.1.0):
   ID, mission ('TESS'|'Kepler'|'K2'), stars (DataFrame), trilegal_fname (CSV path or None), trilegal_url (None),
   sectors, search_radius, N_pix (cosmetic), background_population_source (string).
stars columns: ID, Tmag, Jmag, Hmag, Kmag, ra, dec, mass, rad, Teff, plx, sep, PA, fluxratio, tdepth
   (row 0 = target; add rows for resolved neighbours within ~10 px with sep/PA/fluxratio/tdepth from calc_depths
    if you want NFPP != 0 -- Kuber left them out, hence his NFPP=0 caveat).
background CSV columns: Mact, logg, logTe, [M/H], TESS, J, H, Ks   (TRILEGAL schema; Gaia path writes the same)
Light curve: time = days from mid-transit (folded), flux_0 normalised to 1 out of transit (dilution-corrected, PDCSAP-like,
   else pass dilution_corrected=False), flux_err_0 = scalar out-of-transit rms. exptime in days; nsamples supersampling.
Contrast curves: text file, col1 separation [arcsec], col2 delta-mag; filt in {'TESS','Vis','J','H','K'}; list allowed.
molusc_file: MOLUSC output CSV of surviving companions (used by PTP/PEB/STP/SEB only).
"""
import numpy as np, pandas as pd, json, sys, time as _t, warnings
warnings.filterwarnings("ignore")
import triceratops.triceratops as tr, triceratops.funcs as F

def run(fold_time, fold_flux, fold_err, P, stars, bg_csv, mission="Kepler", ID=0, N=1_000_000, exptime=1765.464/86400,
        nsamples=20, contrast_curve_file=None, filt="K", molusc_file=None, seed=0, dilution_corrected=True):
    T = object.__new__(tr.target)
    T.ID = int(ID); T.mission = mission; T.search_radius = 10; T.sectors = np.array([1]); T.N_pix = 22
    T.stars = stars; T.trilegal_fname = bg_csv; T.trilegal_url = None; T.background_population_source = "file"
    np.random.seed(seed)
    T.calc_probs(time=fold_time, flux_0=fold_flux, flux_err_0=float(fold_err), P_orb=P, N=N, parallel=False, verbose=0,
                 exptime=exptime, nsamples=nsamples, contrast_curve_file=contrast_curve_file,
                 filt=filt if contrast_curve_file else ("TESS"), molusc_file=molusc_file, dilution_corrected=dilution_corrected)
    return T

if __name__ == "__main__":
    import batman
    rng = np.random.default_rng(3); P = 13.769448
    ph = np.sort(rng.uniform(-0.25, 0.25, 1400))
    pm = batman.TransitParams(); pm.t0=0; pm.per=P; pm.rp=0.046; pm.a=39.7; pm.inc=89.3; pm.ecc=0; pm.w=90; pm.u=[0.4,0.25]; pm.limb_dark="quadratic"
    fl = batman.TransitModel(pm, ph, supersample_factor=7, exp_time=1765.464/86400).light_curve(pm) + rng.normal(0, 1500e-6, ph.size)
    err = float(np.std(fl[np.abs(ph) > 0.075])); depth = 1 - np.mean(fl[np.abs(ph) < 0.03])/np.mean(fl[np.abs(ph) > 0.075])
    stars = pd.DataFrame([dict(ID="120686465", Tmag=15.2114, Jmag=13.883, Hmag=13.213, Kmag=12.969, ra=285.907827, dec=40.941988,
                               mass=0.505, rad=0.507, Teff=3633., plx=3.40648, sep=0.0, PA=0.0, fluxratio=1.0, tdepth=depth)])
    # DEMO field population (replace with Gaia-cone or TRILEGAL population!)
    def seq(xc, yc):
        x, y = F.mamajek_column(xc), F.mamajek_column(yc); k = np.isfinite(x) & np.isfinite(y); x, y = x[k], y[k]; o = np.argsort(x); x, y = x[o], y[o]
        k = np.concatenate(([True], np.diff(x) > 1e-9)); return x[k], y[k]
    n = 400; mass = rng.uniform(0.1, 1.3, n); Teff = np.interp(mass, *seq("Msun", "Teff")); rad = np.interp(Teff, *seq("Teff", "R_Rsun"))
    MG = np.interp(Teff, *seq("Teff", "M_G")); MJ = np.interp(Teff, *seq("Teff", "M_J")); JH = np.interp(Teff, *seq("Teff", "JH")); MK = np.interp(Teff, *seq("Teff", "M_Ks"))
    logg = np.log10(F.G*(mass*F.Msun)/(rad*F.Rsun)**2); dmod = rng.uniform(8, 13, n)
    pd.DataFrame({"Mact": mass, "logg": logg, "logTe": np.log10(Teff), "[M/H]": 0.0, "TESS": MG+dmod-0.4, "J": MJ+dmod, "H": MJ-JH+dmod, "Ks": MK+dmod}).to_csv("bg_demo.csv", index=False)
    t0 = _t.time(); T = run(ph, fl, err, P, stars, "bg_demo.csv", N=int(sys.argv[1]) if len(sys.argv) > 1 else 20000, ID=120686465)
    print(T.probs[["scenario", "M_s", "R_s", "b", "R_p", "prob"]].to_string()); print("FPP %.4g NFPP %.4g  (%.0f s)" % (T.FPP, T.NFPP, _t.time()-t0))
