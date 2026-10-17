"""
Triage a batch run: apply the pre-registered first screen, compare real vs inverted
pass rates, run ephemeris/known-object matching on every real survivor, and write
runs/<run-id>/survivors.csv + triage.md.

    python -m pipeline.triage --run-id poolA_pilot --targets targets/poolA_kbonus_bg_dwarfs_ranked.csv
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.batch_search import load_config, passes_first_screen  # noqa: E402
from pipeline.ephemeris_match import match_signal  # noqa: E402


def flat_signal(prefix: str, s: dict | None) -> dict:
    if not s:
        return {}
    keys = ("period_d", "t0_bjd", "depth_ppm", "depth_err_ppm", "duration_h", "sde", "snr_fit", "snr_emp",
            "n_transits", "max_single_event_frac", "odd_even_sigma", "secondary_sigma", "n_periods")
    return {f"{prefix}_{k}": s.get(k) for k in keys}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--targets", default="targets/poolA_kbonus_bg_dwarfs_ranked.csv")
    ap.add_argument("--config", default=None, help="defaults to runs/<run-id>/config_used.json")
    a = ap.parse_args()
    run_dir = os.path.join("runs", a.run_id)
    cfg = load_config(a.config or os.path.join(run_dir, "config_used.json"))
    rows = [json.loads(l) for l in open(os.path.join(run_dir, "results.jsonl"))]
    tg = pd.read_csv(a.targets).set_index("fname") if os.path.exists(a.targets) else None

    recs = []
    for r in rows:
        rec = {"fname": r["fname"], "ok": r.get("ok"), "error": r.get("error"), "skipped": r.get("skipped"),
               "runtime_s": r.get("runtime_s"), "n_pts": r.get("n_after_quality"), "n_quarters": r.get("n_quarters"),
               "p2p_ppm": r.get("p2p_ppm"), "ra": r.get("ra"), "dec": r.get("dec"), "tpforg": r.get("tpforg"),
               "kepmag": r.get("kepmag")}
        if r.get("ok") and "real" in r:
            rs1, is1 = r["real"]["s1"], r.get("inverted", {}).get("s1")
            pr, fr = passes_first_screen(rs1, cfg)
            pi, fi = passes_first_screen(is1, cfg) if is1 else (False, ["none"])
            rec.update({"pass_real": pr, "fail_real": ",".join(fr), "pass_inv": pi, "fail_inv": ",".join(fi)})
            rec.update(flat_signal("real", rs1))
            rec.update(flat_signal("real2", r["real"].get("s2")))
            rec.update(flat_signal("inv", is1))
            if r["real"].get("s2"):
                p2, f2 = passes_first_screen(r["real"]["s2"], cfg)
                rec.update({"pass_real2": p2, "fail_real2": ",".join(f2)})
        if tg is not None and r["fname"] in tg.index:
            t = tg.loc[r["fname"]]
            rec.update({"rank": int(t["rank"]), "gaia": t["gaia_designation"], "G": t["phot_g_mean_mag"],
                        "sptype": t["sptype"], "dist_pc": t["dist_pc"], "Rstar": t["Rstar"], "RUWE": t["RUWE"],
                        "kic": t["kic"]})
        recs.append(rec)
    df = pd.DataFrame(recs)
    df.to_csv(os.path.join(run_dir, "triage_all.csv"), index=False)

    ok = df[df.ok == True]  # noqa: E712
    n_ok = len(ok)
    n_real = int(ok.pass_real.fillna(False).sum())
    n_inv = int(ok.pass_inv.fillna(False).sum())
    lines = [f"# Triage: run `{a.run_id}`", "",
             f"- stars with results: {len(df)} (ok {n_ok}, errors {int(df.error.notna().sum())}, skipped {int(df.skipped.notna().sum())})",
             f"- first-screen passes: **real {n_real}**, **inverted {n_inv}** (pre-registered thresholds: {json.dumps(cfg['first_screen'])})",
             f"- mean runtime {df.runtime_s.mean():.0f} s/star; median point-to-point scatter {df.p2p_ppm.median():.0f} ppm", ""]
    if n_ok:
        lines.append("- What the inverted null does and does not calibrate: inverting the flux (2 - f) turns dips into bumps, so the inverted pass rate bounds the rate of SIGN-SYMMETRIC noise and variability passing the screen; it cannot mimic dip-shaped (sign-asymmetric) systematics such as post-gap thermal ramps. A real-minus-inverted excess is therefore NOT a count of astrophysical signals; compare per period bin:")
        bins = [0.5, 2.0, 10.0, 20.0, 40.1]
        for lo, hi in zip(bins[:-1], bins[1:]):
            nr = int(((ok.pass_real.fillna(False)) & (ok.real_period_d >= lo) & (ok.real_period_d < hi)).sum())
            ni = int(((ok.pass_inv.fillna(False)) & (ok.inv_period_d >= lo) & (ok.inv_period_d < hi)).sum())
            lines.append(f"  - P {lo:g}-{hi:g} d: real {nr}, inverted {ni}")
        lines.append("")
        # SDE / SNR distributions real vs inverted
        for col in ("sde", "snr_emp"):
            rv, iv = ok[f"real_{col}"].dropna(), ok[f"inv_{col}"].dropna()
            if len(rv) and len(iv):
                qs = [50, 90, 99, 100]
                lines.append(f"- {col}: real p50/p90/p99/max = {np.percentile(rv, qs).round(1).tolist()}; inverted = {np.percentile(iv, qs).round(1).tolist()}")
        lines.append("")

    surv = ok[ok.pass_real.fillna(False)].copy()
    surv_rows = []
    for _, s in surv.iterrows():
        hits = match_signal(s.ra, s.dec, s.real_period_d, tpforg_kic=int(s.tpforg) if pd.notna(s.tpforg) else None,
                            kic=int(s.kic) if pd.notna(s.get("kic")) else None, radius_arcmin=3.0,
                            t0_bkjd=(s.real_t0_bjd - 2454833.0) if pd.notna(s.real_t0_bjd) else None)
        note = "; ".join(f"{h.catalog}:{h.name} {h.ratio} sep={h.sep_arcmin:.2f}' {h.disposition}"
                         + (f" same_tpf" if h.same_tpf else "") + (f" same_kic" if h.same_kic else "")
                         + (f" dphase={h.epoch_offset_phase:.3f}" if pd.notna(h.epoch_offset_phase) else "")
                         for _, h in hits.iterrows()) if len(hits) else ""
        d = s.to_dict()
        # both-polarity veto (variability): same period or 1:2 / 2:1 harmonic found in the inverted run
        pr, pi = s.real_period_d, s.get("inv_period_d")
        d["both_polarity"] = bool(pd.notna(pi) and any(abs(pr - k * pi) / (k * pi) < 0.01 for k in (1.0, 2.0, 0.5)))
        d["known_matches"] = note
        d["known_match_flag"] = "SAME_STAR" if any(hits.same_kic) else ("SAME_TPF" if any(hits.same_tpf) else ("NEARBY" if len(hits) else "none"))
        surv_rows.append(d)
    sv = pd.DataFrame(surv_rows)
    if len(sv):
        sv = sv.sort_values("real_snr_emp", ascending=False)
        sv.to_csv(os.path.join(run_dir, "survivors.csv"), index=False)
        lines.append(f"## Real first-screen survivors ({len(sv)})\n")
        lines.append("| rank | fname | G | type | d (pc) | P (d) | depth (ppm) | dur (h) | SDE | SNR_emp | Ntr | o/e σ | sec σ | maxSE | inv SDE | both-pol | known |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for _, s in sv.iterrows():
            lines.append(f"| {s.get('rank','')} | {s.fname} | {s.get('G', float('nan')):.1f} | {s.get('sptype','')} | {s.get('dist_pc', float('nan')):.0f} | {s.real_period_d:.4f} | {s.real_depth_ppm:.0f}±{s.real_depth_err_ppm:.0f} | {s.real_duration_h:.1f} | {s.real_sde:.1f} | {s.real_snr_emp:.1f} | {s.real_n_transits} | {s.real_odd_even_sigma if s.real_odd_even_sigma is None else round(s.real_odd_even_sigma,1)} | {s.real_secondary_sigma if s.real_secondary_sigma is None else round(s.real_secondary_sigma,1)} | {s.real_max_single_event_frac:.2f} | {s.inv_sde:.1f} | {'YES' if s.both_polarity else ''} | {s.known_match_flag}: {s.known_matches[:90]} |")
    inv = ok[ok.pass_inv.fillna(False)]
    if len(inv):
        lines.append(f"\n## Inverted (null) passes ({len(inv)}) — the false-alarm floor\n")
        lines.append("| fname | P (d) | depth (ppm) | SDE | SNR_emp | Ntr |")
        lines.append("|---|---|---|---|---|---|")
        for _, s in inv.iterrows():
            lines.append(f"| {s.fname} | {s.inv_period_d:.4f} | {s.inv_depth_ppm:.0f} | {s.inv_sde:.1f} | {s.inv_snr_emp:.1f} | {s.inv_n_transits} |")
    with open(os.path.join(run_dir, "triage.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
