"""Transit-model fit (batman + emcee) of a stage-2 signal, with a stellar-density consistency check against Gaia.

    python scripts/transit_fit.py poolA_main kic-006699887 [--steps 4000] [--walkers 40]

Data: the same light curve stage 2 uses (KBONUS PSF FLUX, quality 0, per-quarter normalised, biweight re-detrend
with the transit masked), cut to +-3 transit durations around every epoch. Model: circular orbit, quadratic limb
darkening with Kipping (2013) q1, q2 free in [0, 1], batman supersampled 15x over the 29.4-min long cadence.
Free parameters: t0, P, rp/rs, b, log10(rho_star [g/cc]), q1, q2, log(jitter). Priors uniform; b in [0, 1 + rp/rs).

Output (runs/<run>/transit_fit/<fname>/): fit.json (medians, 16/84 percentiles, derived a/R*, T14, Rp with the Gaia
FLAME radius, and the density comparison), chain.npz, fit.png.

The density test: a transit on the target must give rho_star consistent with the host's Gaia DR3 FLAME density.
A transit on an unresolved companion or background star of different size generally does not. A match is supporting
evidence, not proof (eccentric orbits shift the transit density by ((1 + e sin w) / sqrt(1 - e^2))^3).
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.batch_search import fetch_kbonus, read_kbonus, _phase_mask  # noqa: E402
from pipeline.vet_batch import _normalise, _flatten  # noqa: E402

G = 6.674e-8
CAD = 29.4244 / 1440.0
RSUN_REARTH = 109.076

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("fname")
ap.add_argument("--steps", type=int, default=4000)
ap.add_argument("--walkers", type=int, default=40)
ap.add_argument("--burn", type=float, default=0.5)
ap.add_argument("--host-params", default="runs/poolA_main/host_params/gaia_dr3_astrophysical_parameters.csv")
ap.add_argument("--config", default="pipeline/prereg_poolA.json")
ap.add_argument("--vet-json", default=None)
ap.add_argument("--label", default="")
a = ap.parse_args()

import batman  # noqa: E402
import emcee  # noqa: E402

cfg = json.load(open(a.config))
vet = json.load(open(a.vet_json or f"runs/{a.run}/vet/{a.fname}/vet.json"))
P0, D, T00 = float(vet["P"]), float(vet["dur_h"]) / 24.0, float(vet["t0_bkjd"])
depth0 = float(vet["depth_pass2_ppm"]) * 1e-6
out = f"runs/{a.run}/transit_fit/{a.fname}{('_' + a.label) if a.label else ''}"
os.makedirs(out, exist_ok=True)

path = fetch_kbonus(a.fname, cache_dir="data/kbonus_batch", keep=True)
t, f, fe, q, quarter, meta = read_kbonus(path, flux_col=cfg.get("flux_col", "FLUX"))
t, y, e, quarter = _normalise(t, f, fe, q, quarter)
mask = _phase_mask(t, P0, 2.0 * D, T00)
y2, _ = _flatten(t, y, cfg.get("window_length_d", 1.0), mask=mask)
g = np.isfinite(y2)
t, y2 = t[g], y2[g]
# choose a reference epoch near the middle of the baseline to decorrelate t0 and P
nmid = np.round((np.median(t) - T00) / P0)
T0ref = T00 + nmid * P0
ph = ((t - T0ref + 0.5 * P0) % P0) - 0.5 * P0
w = np.abs(ph) < 3 * D
t, y2 = t[w], y2[w]
oot = np.abs((((t - T0ref + 0.5 * P0) % P0) - 0.5 * P0)) > 0.75 * D
sig0 = float(np.std(y2[oot]))
print(f"{a.fname}: {len(t)} cadences within 3 durations; OOT std {sig0*1e6:.0f} ppm; t0 ref {T0ref:.5f} BKJD (epoch {int(nmid)})")

params = batman.TransitParams()
params.t0, params.per, params.rp, params.a, params.inc, params.ecc, params.w = T0ref, P0, np.sqrt(depth0), 20.0, 89.0, 0.0, 90.0
params.u, params.limb_dark = [0.4, 0.25], "quadratic"
m = batman.TransitModel(params, t, supersample_factor=15, exp_time=CAD)


def unpack(th):
    t0, P, k, b, lrho, q1, q2, lj = th
    rho = 10 ** lrho
    aR = ((G * rho * (P * 86400.0) ** 2) / (3 * np.pi)) ** (1 / 3)
    return t0, P, k, b, rho, aR, q1, q2, lj


def model(th):
    t0, P, k, b, rho, aR, q1, q2, lj = unpack(th)
    cosi = b / aR
    if cosi >= 1:
        return None
    params.t0, params.per, params.rp, params.a, params.inc = t0, P, k, aR, np.degrees(np.arccos(cosi))
    sq = np.sqrt(q1)
    params.u = [2 * sq * q2, sq * (1 - 2 * q2)]
    return m.light_curve(params)


def lnprob(th):
    t0, P, k, b, lrho, q1, q2, lj = th
    if not (0.001 < k < 0.3 and 0 <= b < 1 + k and -2 < lrho < 2 and 0 < q1 < 1 and 0 < q2 < 1 and -12 < lj < -4
            and abs(t0 - T0ref) < 0.1 and abs(P - P0) < 0.01):
        return -np.inf
    mu = model(th)
    if mu is None:
        return -np.inf
    s2 = sig0 ** 2 + np.exp(2 * lj)
    r = y2 - mu
    return -0.5 * np.sum(r * r / s2 + np.log(2 * np.pi * s2))


# starting point: a short least-squares-ish scan over b and rho at fixed depth
best, bl = None, -np.inf
for b0 in (0.1, 0.4, 0.7, 0.85):
    for lr in np.linspace(-0.8, 1.0, 19):
        th = np.array([T0ref, P0, np.sqrt(depth0) * 1.05, b0, lr, 0.4, 0.35, np.log(sig0 * 0.1)])
        lp = lnprob(th)
        if lp > bl:
            bl, best = lp, th
ndim = len(best)
rng = np.random.default_rng(42)
p0 = best + rng.normal(0, 1, (a.walkers, ndim)) * np.array([2e-4, 2e-6, 2e-3, 0.03, 0.03, 0.05, 0.05, 0.2])
p0[:, 3] = np.abs(p0[:, 3]); p0[:, 5:7] = np.clip(p0[:, 5:7], 0.01, 0.99)
sampler = emcee.EnsembleSampler(a.walkers, ndim, lnprob)
sampler.run_mcmc(p0, a.steps, progress=False)
chain = sampler.get_chain(discard=int(a.burn * a.steps), flat=True)
lnp = sampler.get_log_prob(discard=int(a.burn * a.steps), flat=True)
try:
    tau = sampler.get_autocorr_time(tol=0).tolist()
except Exception:  # noqa: BLE001
    tau = None
names = ["t0_bkjd", "P_d", "rp_rs", "b", "log10_rho_star_gcc", "q1", "q2", "ln_jitter"]


def pct(x):
    lo, me, hi = np.percentile(x, [16, 50, 84])
    return {"median": float(me), "minus": float(me - lo), "plus": float(hi - me)}


res = {n: pct(chain[:, i]) for i, n in enumerate(names)}
rho = 10 ** chain[:, 4]
aR = ((G * rho * (chain[:, 1] * 86400.0) ** 2) / (3 * np.pi)) ** (1 / 3)
k, b = chain[:, 2], chain[:, 3]
inc = np.arccos(np.clip(b / aR, -1, 1))
T14 = chain[:, 1] / np.pi * np.arcsin(np.clip(np.sqrt((1 + k) ** 2 - b ** 2) / (aR * np.sin(inc)), 0, 1)) * 24
res["rho_star_gcc"] = pct(rho)
res["a_over_rstar"] = pct(aR)
res["T14_h"] = pct(T14)
res["depth_ppm_k2"] = pct(k ** 2 * 1e6)
res["t0_bjd"] = {kk: (v + 2454833.0 if kk == "median" else v) for kk, v in res["t0_bkjd"].items()}
hp = pd.read_csv(a.host_params) if os.path.exists(a.host_params) else None
dens = None
if hp is not None and (hp.fname == a.fname).any():
    r = hp[hp.fname == a.fname].iloc[0]
    R, Rl, Ru, M, Ml, Mu = r.radius_flame, r.radius_flame_lower, r.radius_flame_upper, r.mass_flame, r.mass_flame_lower, r.mass_flame_upper
    if np.isfinite(R) and np.isfinite(M):
        Rs = rng.normal(R, (Ru - Rl) / 2, len(chain)); Ms = rng.normal(M, (Mu - Ml) / 2, len(chain))
        rho_g = Ms * 1.989e33 / (4 / 3 * np.pi * (Rs * 6.957e10) ** 3)
        lr_t, lr_g = np.log10(rho), np.log10(rho_g[np.isfinite(rho_g) & (rho_g > 0)])
        diff = np.median(lr_t) - np.median(lr_g)
        sd = np.hypot(np.std(lr_t), np.std(lr_g))
        Rp = k * rng.normal(R, (Ru - Rl) / 2, len(chain)) * RSUN_REARTH
        dens = {"rho_gaia_flame_gcc": pct(rho_g), "log10_ratio_transit_over_gaia": float(diff), "tension_sigma": float(diff / sd),
                "fraction_of_transit_posterior_within_gaia_68pct": float(np.mean((rho > np.percentile(rho_g, 16)) & (rho < np.percentile(rho_g, 84)))),
                "Rp_earth": pct(Rp), "host_R_flame": float(R), "host_M_flame": float(M)}
res["density_check"] = dens
res.update({"fname": a.fname, "n_points": int(len(t)), "oot_std_ppm": sig0 * 1e6, "walkers": a.walkers, "steps": a.steps, "burn_frac": a.burn,
            "autocorr_time": tau, "acceptance": float(np.mean(sampler.acceptance_fraction)), "max_lnprob": float(np.max(lnp)),
            "model": "batman quadratic LD (Kipping q1,q2), circular orbit, 15x supersampled 29.4-min cadence, white noise + jitter",
            "data": f"KBONUS FLUX quality 0, per-quarter normalised, biweight {cfg.get('window_length_d', 1.0)} d re-detrend with transit masked; +-3 durations",
            "ephemeris_source": a.vet_json or f"runs/{a.run}/vet/{a.fname}/vet.json"})
json.dump(res, open(os.path.join(out, "fit.json"), "w"), indent=1)
np.savez_compressed(os.path.join(out, "chain.npz"), chain=chain.astype(np.float32), names=np.array(names))
print(json.dumps({k: res[k] for k in ("rp_rs", "b", "rho_star_gcc", "T14_h", "depth_ppm_k2", "a_over_rstar")}, indent=0))
print("density check:", json.dumps(dens, indent=0) if dens else None)
print("acceptance", round(res["acceptance"], 3), "autocorr", None if tau is None else [round(x) for x in tau])
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    med = np.median(chain, axis=0)
    ph = ((t - med[0] + 0.5 * med[1]) % med[1]) - 0.5 * med[1]
    o = np.argsort(ph)
    mu = model(med)
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    bins = np.linspace(-3 * D, 3 * D, 61)
    idx = np.digitize(ph, bins)
    bx = [0.5 * (bins[i - 1] + bins[i]) for i in range(1, len(bins)) if (idx == i).sum() > 3]
    by = [np.mean(y2[idx == i]) for i in range(1, len(bins)) if (idx == i).sum() > 3]
    be = [np.std(y2[idx == i]) / np.sqrt((idx == i).sum()) for i in range(1, len(bins)) if (idx == i).sum() > 3]
    ax[0].plot(ph * 24, (y2 - 1) * 1e6, ",", color="0.7", alpha=0.4)
    ax[0].errorbar(np.array(bx) * 24, (np.array(by) - 1) * 1e6, np.array(be) * 1e6, fmt="o", ms=3, color="k")
    ax[0].plot(ph[o] * 24, (mu[o] - 1) * 1e6, color="C3", lw=2)
    ax[0].set_ylim(-4 * depth0 * 1e6, 2 * depth0 * 1e6); ax[0].set_ylabel("ppm")
    ax[0].set_title(f"{a.fname}: P {med[1]:.5f} d, rp/rs {med[2]:.4f}, b {med[3]:.2f}, rho* {10**med[4]:.2f} g/cc")
    rb = [np.mean((y2 - mu)[idx == i]) for i in range(1, len(bins)) if (idx == i).sum() > 3]
    ax[1].errorbar(np.array(bx) * 24, np.array(rb) * 1e6, np.array(be) * 1e6, fmt="o", ms=3, color="k"); ax[1].axhline(0, color="C3")
    ax[1].set_xlabel("hours from mid-transit"); ax[1].set_ylabel("resid ppm")
    fig.tight_layout(); fig.savefig(os.path.join(out, "fit.png"), dpi=120)
except Exception as ex:  # noqa: BLE001
    print("plot failed:", ex)
