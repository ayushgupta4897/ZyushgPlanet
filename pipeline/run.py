"""End-to-end run for one KBONUS-BKG source.

Example (the KIC 5689121 reproduction):
  . .venv/bin/activate
  python -m pipeline.run --gaia 2103792461826460288 --ra 285.907827 --dec 40.941988 \
      --neighbour 2103792427466721792 --out runs/kic5689121 \
      --compare-period 13.769448251 --compare-t0-bjd 2455630.2454306 --compare-duration-h 2.3

Flow:  fetch -> describe FITS -> stitch/normalise -> pass-1 detrend (blind) -> noise metrics
       -> BLS 1-60 d -> TLS 1-30 d -> adopt ephemeris -> pass-2 detrend with transits masked
       -> BLS refine + vetting (box depth, odd/even, secondary, per-transit, empirical noise)
       -> inverted-flux null (BLS + TLS on the masked-detrended LC, flux -> 2 - flux)
       -> neighbour: same blind pipeline + depth at the target ephemeris -> JSON + figures.
Outputs in --out: results.json, fits_structure*.json, lightcurve_*.csv, 01..09 *.png, run.log
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from . import fetch, prep, search, vet, report
from .fetch import BKJD_OFFSET


class Tee:
    def __init__(self, path):
        self.f = open(path, "a"); self.o = sys.stdout
    def write(self, s):
        self.o.write(s); self.f.write(s); self.f.flush()
    def flush(self):
        self.o.flush(); self.f.flush()


def _save_lc(lc, path):
    np.savetxt(path, np.c_[lc.time, lc.time + BKJD_OFFSET, lc.flux, lc.flat, lc.flux_err, lc.quarter, lc.sap_flux],
               delimiter=",", header="time_bkjd,time_bjd_tdb,flux_norm,flux_detrended,flux_err,quarter,sap_flux_norm", comments="")


def vetting(lc, P, T0, D):
    V = {"box_depth": vet.box_depth(lc.time, lc.flat, P, T0, D),
         "odd_even": vet.odd_even(lc.time, lc.flat, P, T0, D),
         "secondary": vet.secondary_eclipse(lc.time, lc.flat, P, T0, D),
         "per_transit": vet.per_transit_depths(lc.time, lc.flat, lc.quarter, P, T0, D),
         "empirical_noise": vet.empirical_noise(lc.time, lc.flat, P, T0, D)}
    print(f"[vet] box depth {V['box_depth']['depth']*1e6:.0f} +/- {V['box_depth']['depth_err']*1e6:.0f} ppm; "
          f"odd {V['odd_even']['odd']['depth']*1e6:.0f}+/-{V['odd_even']['odd']['depth_err']*1e6:.0f}, "
          f"even {V['odd_even']['even']['depth']*1e6:.0f}+/-{V['odd_even']['even']['depth_err']*1e6:.0f} ppm "
          f"({V['odd_even']['difference_sigma']:.2f} sigma); secondary @0.5 {V['secondary']['phase_0p5']['depth']*1e6:.0f}+/-"
          f"{V['secondary']['phase_0p5']['depth_err']*1e6:.0f} ppm; max other dip {V['secondary']['scan_max_dip']['sigma']:.1f} sigma at phase "
          f"{V['secondary']['scan_max_dip']['phase']:.3f}; events {V['per_transit']['summary'].get('n_events')} "
          f"({V['per_transit']['summary'].get('n_positive')} positive); empirical SNR {V['empirical_noise']['empirical_snr']:.1f} "
          f"(random-phase std {V['empirical_noise']['random_depth_std_ppm']:.0f} ppm)")
    return V


def analyse_source(gaia, ra, dec, args, label, out: Path, tls=True, is_neighbour=False):
    """fetch -> describe -> prep (pass 1) -> noise -> BLS -> TLS. Returns (results, lc1, quarters, bls, tls)."""
    sfx = "_neighbour" if is_neighbour else ""
    R = {"gaia_source_id": str(gaia), "label": label}
    path = fetch.download(gaia, ra, dec)
    R["fits_file"] = str(path)
    struct = fetch.describe_fits(path, verbose=not is_neighbour)
    report.write_json(struct, out / f"fits_structure{sfx}.json")
    R["primary_header"] = struct["primary"]
    quarters, _ = fetch.load_quarters(path, min_psffrac=args.min_psffrac)
    lc = prep.prepare(quarters, quality_mode=args.quality, window_length=args.window, break_tolerance=0.5,
                      n_sigma_up=args.clip_sigma)
    R["prep"] = {"quality_mode": args.quality, "window_days": args.window or prep.default_window_days(),
                 "upward_clip_sigma": args.clip_sigma, "min_psffrac": args.min_psffrac, "break_tolerance_days": 0.5}
    R["noise"] = prep.noise_metrics(lc)
    _save_lc(lc, out / f"lightcurve_pass1{sfx}.csv")
    b = search.bls_search(lc.time, lc.flat, lc.flux_err, pmin=args.bls_pmin, pmax=args.bls_pmax, oversample=args.bls_oversample)
    R["bls"] = b.summary()
    T = None
    if tls:
        T = search.tls_search(lc.time, lc.flat, lc.flux_err, pmin=args.tls_pmin, pmax=args.tls_pmax,
                              oversampling_factor=args.tls_osf, duration_grid_step=args.tls_dgs, use_threads=args.threads)
        R["tls"] = T.summary()
    return R, lc, quarters, b, T


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gaia", required=True, help="Gaia DR3 source id of the target")
    ap.add_argument("--ra", type=float); ap.add_argument("--dec", type=float)
    ap.add_argument("--neighbour", help="Gaia DR3 id of a neighbour to run the same pipeline on (control)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--quality", default="hard", choices=["hard", "default", "none"])
    ap.add_argument("--min-psffrac", type=float, default=0.0)
    ap.add_argument("--window", type=float, default=None, help="wotan window (d); default prep.default_window_days() = 0.496")
    ap.add_argument("--clip-sigma", type=float, default=4.0)
    ap.add_argument("--bls-pmin", type=float, default=1.0); ap.add_argument("--bls-pmax", type=float, default=60.0)
    ap.add_argument("--bls-oversample", type=float, default=3.0)
    ap.add_argument("--tls-pmin", type=float, default=1.0); ap.add_argument("--tls-pmax", type=float, default=30.0)
    ap.add_argument("--tls-osf", type=int, default=3); ap.add_argument("--tls-dgs", type=float, default=1.1)
    ap.add_argument("--threads", type=int, default=3)
    ap.add_argument("--no-tls", action="store_true"); ap.add_argument("--no-tls-null", action="store_true")
    ap.add_argument("--tls-neighbour", action="store_true", help="also run TLS on the neighbour (default: BLS only)")
    ap.add_argument("--compare-period", type=float); ap.add_argument("--compare-t0-bjd", type=float)
    ap.add_argument("--compare-duration-h", type=float, default=None)
    args = ap.parse_args(argv)

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    sys.stdout = Tee(out / "run.log")
    t_start = time.time()
    print(f"[run] {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} args={vars(args)}")

    # ---------------- pass 1: blind search
    R, lc1, quarters, b, T = analyse_source(args.gaia, args.ra, args.dec, args, "target", out, tls=not args.no_tls)
    P1 = T.period if T else b.period; T01 = T.t0 if T else b.t0; D1 = b.duration
    R["pass1_ephemeris"] = {"period_d": P1, "t0_bkjd": T01, "duration_h": D1 * 24,
                            "source": "TLS period/T0 + BLS box duration" if T else "BLS"}
    R["pass1_box_depth"] = vet.box_depth(lc1.time, lc1.flat, P1, T01, D1)
    print(f"[pass1] box depth on blind-detrended LC: {R['pass1_box_depth']['depth']*1e6:.0f} +/- {R['pass1_box_depth']['depth_err']*1e6:.0f} ppm (biased low by in-transit trend pull)")

    # ---------------- pass 2: masked re-detrend, refine, vet
    lc2 = prep.prepare(quarters, quality_mode=args.quality, window_length=args.window, break_tolerance=0.5,
                       n_sigma_up=args.clip_sigma, mask_ephemeris=(P1, T01, D1))
    _save_lc(lc2, out / "lightcurve_pass2_masked.csv")
    R["noise_pass2"] = prep.noise_metrics(lc2, verbose=False)
    ref = search.bls_refine(lc2.time, lc2.flat, lc2.flux_err, P1)
    R["bls_refine_pass2"] = ref
    P, T0, D = ref["period_d"], ref["t0_bkjd"], ref["duration_h"] / 24
    R["adopted_ephemeris"] = {"period_d": P, "t0_bkjd": T0, "t0_bjd_tdb": T0 + BKJD_OFFSET, "duration_h": D * 24,
                              "source": "BLS refine on transit-masked detrended light curve"}
    R["vetting"] = vetting(lc2, P, T0, D)

    # ---------------- inverted-flux null on the masked-detrended LC (flux -> 2 - flux)
    inv = prep.invert(lc2)
    bi = search.bls_search(inv.time, inv.flat, inv.flux_err, pmin=args.bls_pmin, pmax=args.bls_pmax,
                           oversample=args.bls_oversample, refine=False)
    R["inverted_null"] = {"bls": {k: v for k, v in bi.summary().items() if k != "stats"},
                          "bls_excluding_candidate_band": search.null_sde_excluding(bi.periods, bi.powers, P)}
    ti = None
    if not args.no_tls and not args.no_tls_null:
        ti = search.tls_search(inv.time, inv.flat, inv.flux_err, pmin=args.tls_pmin, pmax=args.tls_pmax,
                               oversampling_factor=args.tls_osf, duration_grid_step=args.tls_dgs, use_threads=args.threads)
        R["inverted_null"]["tls"] = {k: v for k, v in ti.summary().items()
                                     if k not in ("transit_times_bkjd", "transit_depths", "per_transit_count", "snr_per_transit")}
        R["inverted_null"]["tls_excluding_candidate_band"] = search.null_sde_excluding(ti.periods, ti.powers, P)
    print(f"[null] inverted BLS best P={bi.period:.4f} d power={bi.power:.2f} SDE={bi.sde:.2f}" +
          (f"; inverted TLS best P={ti.period:.4f} SDE={ti.sde:.2f} SNR={ti.snr:.2f}" if ti is not None else ""))

    # ---------------- comparison with an external (published) ephemeris
    if args.compare_period and args.compare_t0_bjd:
        mid = 0.5 * (lc2.time.min() + lc2.time.max())
        Dc = (args.compare_duration_h / 24) if args.compare_duration_h else D
        R["comparison"] = vet.compare_ephemeris(P, T0, args.compare_period, args.compare_t0_bjd, mid)
        R["comparison"]["depth_at_reference_ephemeris_pass2"] = vet.ephemeris_depth(lc2.time, lc2.flat, args.compare_period,
                                                                                    args.compare_t0_bjd - BKJD_OFFSET, Dc)
        R["comparison"]["reference_duration_h_used"] = Dc * 24
        print(f"[compare] dP = {R['comparison']['period_diff_d']*86400:.2f} s ({R['comparison']['period_diff_ppm']:.1f} ppm); "
              f"T0 offset at epoch {R['comparison']['epoch_mid_baseline']} (mid-baseline): {R['comparison']['t0_offset_min_at_mid_baseline']:.1f} min; "
              f"depth at reference ephemeris ({Dc*24:.1f} h box): {R['comparison']['depth_at_reference_ephemeris_pass2']['depth']:.0f} +/- "
              f"{R['comparison']['depth_at_reference_ephemeris_pass2']['depth_err']:.0f} ppm")

    # ---------------- neighbour control
    N = None
    if args.neighbour:
        N, nlc1, nq, nb, nT = analyse_source(args.neighbour, None, None, args, "neighbour", out,
                                             tls=(not args.no_tls) and args.tls_neighbour, is_neighbour=True)
        nlc2 = prep.prepare(nq, quality_mode=args.quality, window_length=args.window, break_tolerance=0.5,
                            n_sigma_up=args.clip_sigma, mask_ephemeris=(P, T0, D), verbose=False)
        _save_lc(nlc2, out / "lightcurve_pass2_masked_neighbour.csv")
        N["depth_at_target_ephemeris_pass2"] = vet.ephemeris_depth(nlc2.time, nlc2.flat, P, T0, D)
        N["depth_at_target_ephemeris_pass1"] = vet.ephemeris_depth(nlc1.time, nlc1.flat, P, T0, D)
        N["odd_even_at_target_ephemeris"] = vet.odd_even(nlc2.time, nlc2.flat, P, T0, D)
        N["empirical_noise_at_target_ephemeris"] = vet.empirical_noise(nlc2.time, nlc2.flat, P, T0, D)
        N["bls_excluding_candidate_band"] = search.null_sde_excluding(nb.periods, nb.powers, P)
        gmag_t, gmag_n = R["primary_header"].get("GMAG"), N["primary_header"].get("GMAG")
        if gmag_t and gmag_n:
            ratio = 10 ** (-0.4 * (gmag_n - gmag_t))
            N["flux_ratio_neighbour_over_target_G"] = ratio
            N["depth_required_if_hosted_by_neighbour_ppm"] = R["vetting"]["box_depth"]["depth"] * 1e6 / ratio
        print(f"[neighbour] own best BLS P={nb.period:.4f} d power={nb.power:.2f} SDE={nb.sde:.2f}; depth at target ephemeris "
              f"{N['depth_at_target_ephemeris_pass2']['depth']:.0f} +/- {N['depth_at_target_ephemeris_pass2']['depth_err']:.0f} ppm "
              f"(would need ~{N.get('depth_required_if_hosted_by_neighbour_ppm', float('nan')):.0f} ppm to host the signal)")
        R["neighbour"] = N
        report.fig_fold(nlc2.time, nlc2.flat, P, T0, D, out / "08_neighbour_fold.png",
                        f"Neighbour Gaia DR3 {args.neighbour} (KIC {N['primary_header'].get('KEPLERID')}) at the TARGET ephemeris: "
                        f"{N['depth_at_target_ephemeris_pass2']['depth']:.0f} +/- {N['depth_at_target_ephemeris_pass2']['depth_err']:.0f} ppm",
                        depth_ppm=R["vetting"]["box_depth"]["depth"] * 1e6, color=report.C["neigh"], label="neighbour detrended flux")

    # ---------------- figures for the target
    name = f"Gaia DR3 {args.gaia} (KIC {R['primary_header'].get('KEPLERID')})"
    V = R["vetting"]
    report.fig_raw(lc2, out / "01_raw_stitched.png", f"{name}: KBONUS-BKG PSF flux, per-quarter normalised (trend = pass-2 masked biweight)")
    report.fig_detrended(lc2, out / "02_detrended.png", f"{name}: detrended (pass 2), P={P:.6f} d", P, T0, D)
    report.fig_periodogram(b.periods, b.powers, b.period, out / "03_bls_periodogram.png",
                           f"{name}: BLS 1-60 d (objective=snr), best P={b.period:.6f} d, SDE={b.sde:.1f}; purple = inverted-flux null",
                           "BLS SNR", null=(bi.periods, bi.powers))
    if T is not None:
        report.fig_periodogram(T.periods, T.powers, T.period, out / "04_tls_periodogram.png",
                               f"{name}: TLS {args.tls_pmin:g}-{args.tls_pmax:g} d, best P={T.period:.6f} d, SDE={T.sde:.1f}" +
                               ("; purple = inverted-flux null" if ti is not None else ""), "TLS SDE",
                               null=((ti.periods, ti.powers) if ti is not None else None))
    report.fig_fold(lc2.time, lc2.flat, P, T0, D, out / "05_phase_fold.png",
                    f"{name}: P={P:.6f} d, box depth {V['box_depth']['depth']*1e6:.0f}+/-{V['box_depth']['depth_err']*1e6:.0f} ppm, dur {D*24:.2f} h",
                    model=((T.model_phase - 0.5, T.model_y) if T is not None else None), depth_ppm=V["box_depth"]["depth"] * 1e6)
    report.fig_odd_even(lc2.time, lc2.flat, P, T0, D, V["odd_even"], out / "06_odd_even.png", name)
    s = V["secondary"]["phase_0p5"]
    report.fig_fold(lc2.time, lc2.flat, P, T0, D, out / "07_secondary.png",
                    f"{name}: phase 0.5, depth {s['depth']*1e6:.0f} +/- {s['depth_err']*1e6:.0f} ppm", center=0.5,
                    depth_ppm=V["box_depth"]["depth"] * 1e6)
    report.fig_per_transit(V["per_transit"], out / "09_per_transit.png", f"{name}: individual events (pass-2 LC)")

    R["runtime_total_s"] = time.time() - t_start
    R["software"] = _versions()
    report.write_json(R, out / "results.json")
    print(f"[run] done in {R['runtime_total_s']/60:.1f} min -> {out/'results.json'}")
    return R


def _versions():
    import importlib.metadata as m
    v = {"python": sys.version.split()[0]}
    for p in ("lightkurve", "astropy", "astroquery", "transitleastsquares", "wotan", "numpy", "scipy", "matplotlib", "batman-package"):
        try:
            v[p] = m.version(p)
        except Exception:
            v[p] = None
    return v


if __name__ == "__main__":
    main()
