#!/usr/bin/env python
"""Dilution/contamination check: does the bright neighbour (TPF host KIC 9655236) show the same
6.7035-d signal in its own KBONUS PSF light curve?  Same from-scratch loader/detrend/box code.
The neighbour FITS is a fresh download (untrusted data; read with astropy only)."""
import sys, os, json, numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rederive_kic009655224 as R
NB = sys.argv[1]; TAG = sys.argv[2] if len(sys.argv) > 2 else "neighbour"
OUT = R.OUT
res = json.load(open(os.path.join(OUT, "rederive_results.json")))
P = res["ephemeris"]["period_d"]; t0 = res["ephemeris"]["t0_bjd"] - 2454833.0; DUR = res["ephemeris"]["duration_h"] / 24
S = {}
with fits.open(NB) as f:
    h = f[0].header
    S["header"] = {k: h.get(k) for k in ("OBJECT", "KEPLERID", "GAIAID", "RA_OBJ", "DEC_OBJ", "KEPMAG", "GMAG")}
    S["quarters_in_file"] = [f[i].header["QUARTER"] for i in range(2, len(f), 2)]
    S["per_quarter"] = {int(f[i].header["QUARTER"]): dict(PSFFRAC=f[i].header.get("PSFFRAC"), CROWDSAP=f[i].header.get("CROWDSAP"), FLFRCSAP=f[i].header.get("FLFRCSAP"), NPIXSAP=f[i].header.get("NPIXSAP"),
                                                           median_FLUX=float(np.nanmedian(f[i].data["FLUX"])), median_SAP=float(np.nanmedian(f[i].data["SAP_FLUX"]))) for i in range(2, len(f), 2)}
print("neighbour header:", S["header"]); print("quarters:", S["quarters_in_file"])
D = R.load(NB)
ta = D["t"]
# separation from target
ra_t, dec_t = 294.2826791587251, 46.34150883424167
sep = np.hypot((S["header"]["RA_OBJ"] - ra_t) * np.cos(np.radians(dec_t)), S["header"]["DEC_OBJ"] - dec_t) * 3600
S["separation_arcsec_from_target"] = float(sep); print(f"separation {sep:.2f} arcsec")
out = {}
for key, lab in (("fn", "PSF"), ("sapn", "SAP")):
    flat, keep, mad = R.detrend(ta, D[key], window=1.0)
    t, y = ta[keep], flat[keep]
    d, nin = R.box_depth(t, y, P, t0, DUR)
    rp = R.random_phase_depths(t, y, P, t0, DUR, n=1000)
    e = float(np.std(rp))
    med = float(np.median(D["f" if key == "fn" else "sap"]))
    out[lab] = dict(depth_ppm=float(d), err_ppm=e, sigma=float(d / e), n_in=int(nin), mad_ppm=float(mad * 1e6), n_good=int(len(t)),
                    median_flux_e_s=med, deficit_e_s=float(d * 1e-6 * med), deficit_err_e_s=float(e * 1e-6 * med))
    print(f"neighbour {lab}: depth at target ephemeris {d:.0f} +/- {e:.0f} ppm ({d/e:.1f} sigma), N_in={nin}, flux {med:.0f} e/s -> deficit {d*1e-6*med:.2f} +/- {e*1e-6*med:.2f} e/s")
    if key == "fn":
        # restrict to the same quarters as the target for a like-for-like comparison
        s = np.isin(D["q"][keep], res["data"]["quarters"])
        d2, n2 = R.box_depth(t[s], y[s], P, t0, DUR); e2 = float(np.std(R.random_phase_depths(t[s], y[s], P, t0, DUR, n=600)))
        out[lab]["same_quarters_as_target"] = dict(depth_ppm=float(d2), err_ppm=e2, sigma=float(d2 / e2), n_in=int(n2))
        print(f"  (same quarters as target: {d2:.0f} +/- {e2:.0f} ppm)")
        # narrow BLS around the period: is there any peak at 6.7035 in the neighbour's PSF flux?
        dy = np.full_like(y, mad)
        pg = np.linspace(P * 0.98, P * 1.02, 6001)
        r = BoxLeastSquares(t, y, dy).power(pg, [DUR], objective="likelihood", oversample=10)
        i = int(np.argmax(r.power))
        out[lab]["narrow_bls"] = dict(best_P_d=float(pg[i]), dlogL=float(r.power[i]), depth_ppm=float(r.depth[i] * 1e6), dlogL_at_target_P=float(np.interp(P, pg, r.power)), mean_dlogL=float(r.power.mean()), std_dlogL=float(r.power.std()))
        print("  narrow BLS +/-2%:", out[lab]["narrow_bls"])
        # per-quarter at ephemeris
        pq = {}
        for qq in np.unique(D["q"][keep]):
            s = D["q"][keep] == qq
            dq, nq = R.box_depth(t[s], y[s], P, t0, DUR); eq = float(np.std(R.random_phase_depths(t[s], y[s], P, t0, DUR, n=200)))
            pq[int(qq)] = dict(depth_ppm=float(dq), err_ppm=eq, n_in=int(nq))
        out[lab]["per_quarter_depth"] = pq
        print("  per-quarter:", {k: (round(v["depth_ppm"]), round(v["err_ppm"])) for k, v in pq.items()})
S["at_target_ephemeris"] = out
# expected depth on the neighbour if the eclipse were on the neighbour and leaked into the target:
# target absolute deficit (stitched scale) from main results
tgt = res["depth"]["psf_depth_ppm"] * 1e-6 * res["data"]["median_psf_flux_e_s"]
S["target_deficit_e_s"] = float(tgt)
S["note"] = ("If the signal originated on this neighbour, its own PSF flux should show a deficit at least as large in e/s as the "
             "target's (the PSF fit assigns a source's flux change mostly to that source). The target's deficit is ~1.07 e/s.")
json.dump(S, open(os.path.join(OUT, f"neighbour_check_{TAG}.json"), "w"), indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
print(f"wrote neighbour_check_{TAG}.json")
