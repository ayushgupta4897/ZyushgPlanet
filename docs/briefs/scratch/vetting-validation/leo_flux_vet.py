"""leo_flux_vet.py -- LEO-Vetter flux-level vetting of one TCE (Kepler/K2/TESS/KBONUS light curves).

Verified 2026-10-17 in an isolated venv: leo-vetter 1.2.0 (PyPI 2026-07-29), Python 3.13,
numpy 2.5, batman-package 2.5.1, lmfit 1.3.4. 121 metrics in ~2 s per TCE on one CPU.
Pixel-level vetting (leo_vetter.pixel) is TESS-only and needs:
    pip install git+https://github.com/stevepur/transit-diffImage  spiceypy
Reference: Kunimoto et al. 2025, AJ 170, 280 (arXiv:2509.10619); thresholds.py in the package.
Inputs: time (days), raw flux (normalised, undetrended), detrended flux (normalised to 1),
flux_err (same units), period (d), epoch (same time system as `time`), duration (d).
`star` dict keys actually read by the code: mass, e_mass, rad, e_rad, Teff, e_Teff, u1, u2
(rho/tic/ra/dec/sector/cam/ccd are used by the pixel module).
"""
import numpy as np, json, sys
from leo_vetter.main import TCELightCurve
from leo_vetter.thresholds import check_thresholds

def vet(tic, time, raw, flux, flux_err, period, epoch, dur_days, star, label="TCE"):
    tlc = TCELightCurve(int(tic), time, raw, flux, flux_err, period, epoch, dur_days)
    tlc.compute_flux_metrics(star, verbose=False)   # keyword defaults: cap_b=True, max_RpRs=1, frac=0.7, gap=0.3, chases=0.01, rubble=0.75, redchi2=5, A=0.3
    m = tlc.metrics
    fa = check_thresholds(m, "FA"); fp = check_thresholds(m, "FP")
    # replicate the published default sub-tests so the log says WHICH test fired
    sub = {
        "weak(MES<6.2)": m["MES"] < 6.2,
        "invalid_transits(new_MES<6.2|new_N<3)": (m["new_MES"] < 6.2) or (m["new_N_transit"] < 3),
        "bad_shape(SHP>0.6)": m["SHP"] > 0.6,
        "non_unique(MS1<0.2|MS2<0.8|MS3<0.8)": (m["sig_pri"]/m["Fred"] - m["FA1"] < 0.2) or (m["sig_pri"]-m["sig_ter"]-m["FA2"] < 0.8) or (m["sig_pri"]-m["sig_pos"]-m["FA2"] < 0.8),
        "chases(N<=5 & mean<0.78)": (m["N_transit"] <= 5) and (m["mean_chases"] < 0.78),
        "dmm(>1.5)": m["DMM"] > 1.5,
        "single_event(maxSES/MES>0.8 & N<=10)": (m["max_SES"]/m["MES"] > 0.8) and (m["N_transit"] <= 10),
        "sinusoidal(sine_sig>15 & P<10)": (m["sine_sig"] > 15) and (m["per"] < 10),
        "unphysical_duration": (m["transit_aRs"] < 1.5) or (m["aRs"] < 2) or (m["q"]/m["trap_qtran"] < 0.6) or (not np.isfinite(m["sig_sec"])) or (m["trap_qtran"] > 0.5),
        "chi(CHI<7.8)": m["CHI"] < 7.8,
        "data_gapped(>=50%)": m["N_gap_2.0"]/m["N_transit"] >= 0.5,
        "FP:odd_even": (m["sig_dep"] > 3 and (m["trap_sig_dep"] > 3 or m["transit_sig_dep"] > 3)) or m["transit_sig_epo"] > 10 or m["trap_sig_epo"] > 10,
        "FP:vshaped(b+RpRs>1.5)": m["transit_b"] + m["transit_RpRs"] > 1.5,
        "FP:large(Rp>22Re)": m["Rp"] > 22,
    }
    fired = [k for k, v in sub.items() if bool(v) and np.isfinite(float(v))]
    print(f"[{label}] MES={m['MES']:.2f} N={int(m['N_transit'])} dep={m['dep']*1e6:.0f}ppm Rp={m['Rp']:.2f}Re  FA={fa} FP={fp}  fired={fired}")
    return m, fa, fp, fired

if __name__ == "__main__":
    # self-test on a synthetic 13.77-d, 2100-ppm transit with 1500-ppm white noise, 30-min cadence, 360 d
    import batman
    from scipy.ndimage import median_filter
    rng = np.random.default_rng(1); t = np.arange(0, 360, 29.4/1440.)
    P, T0 = 13.769448, 5.3
    pm = batman.TransitParams(); pm.t0=T0; pm.per=P; pm.rp=0.046; pm.a=39.7; pm.inc=89.3; pm.ecc=0; pm.w=90; pm.u=[0.4,0.25]; pm.limb_dark="quadratic"
    raw = batman.TransitModel(pm, t).light_curve(pm)*(1+0.002*np.sin(2*np.pi*t/27.3)) + rng.normal(0, 1500e-6, t.size)
    flux = raw/median_filter(raw, size=97, mode="nearest")     # edge-safe demo detrend; use wotan biweight (window ~3x duration) in production
    err = np.full(t.size, 1500e-6)
    star = dict(mass=0.505, e_mass=0.02, rad=0.507, e_rad=0.015, Teff=3633., e_Teff=70., u1=0.4, u2=0.25)
    vet(120686465, t, raw, flux, err, P, T0, 2.3/24, star, "REAL")
    vet(120686465, t, 2-raw, 2-flux, err, P, T0, 2.3/24, star, "INVERTED")
