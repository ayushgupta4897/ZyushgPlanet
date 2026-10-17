"""Assemble the evidence pack for one signal under candidates/<fname>/ with a SHA-256 manifest.

    python scripts/build_evidence_pack.py poolA_main kic-006699887 [--deep runs/kic6699887_deep] [--second-signal P47]

Copies (never re-derives) the pipeline products: stage-1 row, stage-2 vet.json/events/plot, pixel verdict and
difference images, independent re-derivation and red-team verdict, TRICERATOPS summary and per-seed tables, the
deep single-target run and any post-hoc second-signal vetting, the prereg file and the known-catalogue manifest.
Writes README.md (the key numbers with the file each comes from), manifest.json and SHA256SUMS.

Language: the README calls the object a *candidate* only if pixel verdict is on-target, no ephemeris match,
red team "survives", re-derivation agrees, and a TRICERATOPS summary exists; otherwise "signal".
"""
import argparse
import datetime as dt
import glob
import hashlib
import json
import os
import shutil
import subprocess

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("fname")
ap.add_argument("--deep", default=None, help="directory of a pipeline.run deep run to include")
ap.add_argument("--second-signal", default=None, help="label of a runs/<run>/vet_posthoc/<fname>_<label> second-signal vetting")
ap.add_argument("--out", default=None)
a = ap.parse_args()

R = f"runs/{a.run}"
out = a.out or f"candidates/{a.fname}"
os.makedirs(out, exist_ok=True)


def cp(src, dst_sub):
    if src and os.path.exists(src):
        d = os.path.join(out, dst_sub)
        os.makedirs(os.path.dirname(d), exist_ok=True)
        if os.path.isdir(src):
            shutil.copytree(src, d, dirs_exist_ok=True)
        else:
            shutil.copy2(src, d)
        return dst_sub
    return None


def jload(p):
    return json.load(open(p)) if p and os.path.exists(p) else None


# ---- gather
vet = jload(f"{R}/vet/{a.fname}/vet.json")
pix = jload(f"{R}/pixel_v16/{a.fname}/pixel.json")
rt = jload(f"{R}/adversarial/{a.fname}/redteam_verdict.json")
rd = jload(f"{R}/adversarial/{a.fname}/rederive_result.json")
fpp = jload(f"{R}/fpp/{a.fname}/fpp_summary.json")
prereg = jload("pipeline/prereg_poolA.json")
known = jload("data/known/manifest.json")
s1 = None
with open(f"{R}/results.jsonl") as fh:
    for line in fh:
        if f'"fname": "{a.fname}"' in line:
            s1 = json.loads(line)
            break
if s1 is None and os.path.exists("runs/poolA_pilot/results.jsonl"):
    with open("runs/poolA_pilot/results.jsonl") as fh:
        for line in fh:
            if f'"fname": "{a.fname}"' in line:
                s1 = json.loads(line)
                break
tg = pd.read_csv("targets/poolA_kbonus_bg_dwarfs_ranked.csv")
star = tg[tg.fname == a.fname].iloc[0].to_dict()
second = jload(f"{R}/vet_posthoc/{a.fname}_{a.second_signal}/vet.json") if a.second_signal else None
second_pix = jload(f"{R}/pixel_v16_posthoc/{a.fname}_{a.second_signal}/pixel.json") if a.second_signal else None
deep = jload(os.path.join(a.deep, "results.json")) if a.deep else None

# ---- copy
files = {}
files["stage1_row.json"] = "stage1_row.json"
json.dump(s1, open(os.path.join(out, "stage1_row.json"), "w"), indent=1)
files["target_row.json"] = "target_row.json"
json.dump(star, open(os.path.join(out, "target_row.json"), "w"), indent=1, default=str)
for src, dst in [(f"{R}/vet/{a.fname}/vet.json", "stage2/vet.json"), (f"{R}/vet/{a.fname}/events.csv", "stage2/events.csv"),
                 (f"{R}/vet/{a.fname}/fold_events.png", "stage2/fold_events.png"),
                 (f"{R}/pixel_v16/{a.fname}/pixel.json", "pixels/pixel_v16.json"), (f"{R}/pixel_v16/{a.fname}/diffimages.png", "pixels/diffimages_v16.png"),
                 (f"{R}/pixel_v15/{a.fname}/pixel.json", "pixels/pixel_v15_descriptive.json"),
                 (f"{R}/pixel/{a.fname}/pixel.json", "pixels/pixel_firstpass_descriptive.json"),
                 (f"{R}/adversarial/{a.fname}", "independent"),
                 (f"{R}/fpp/{a.fname}", "fpp"),
                 (f"{R}/ephemeris_fieldwide/{a.fname}.json", "ephemeris_fieldwide.json"),
                 ("pipeline/prereg_poolA.json", "prereg_poolA.json"), ("data/known/manifest.json", "known_catalogues_manifest.json"),
                 (f"{R}/triage.md", "pool_triage.md"), (f"{R}/modal_status.json", "pool_modal_status.json"), (f"{R}/config_used.json", "pool_config_used.json")]:
    r = cp(src, dst)
    if r:
        files[dst] = r
if a.deep:
    cp(a.deep, "deep_run")
if a.second_signal:
    cp(f"{R}/vet_posthoc/{a.fname}_{a.second_signal}", f"second_signal_{a.second_signal}/stage2")
    cp(f"{R}/pixel_v16_posthoc/{a.fname}_{a.second_signal}", f"second_signal_{a.second_signal}/pixels")
# drop bulky scratch binaries from the copy
for pat in ["**/*.npz", "**/*.fits", "**/__pycache__", "deep_run/lightcurve_*.csv", "deep_run/nohup.log"]:
    for p in glob.glob(os.path.join(out, pat), recursive=True):
        (shutil.rmtree if os.path.isdir(p) else os.remove)(p)

# ---- status
pix_ok = bool(pix and pix.get("summary", {}).get("verdict") == "on-target")
eph_ok = bool(vet and not any(h.get("same_kic") or h.get("same_tpf") for h in vet.get("known_matches", [])))
rt_ok = bool(rt and rt.get("verdict") == "survives")
rd_ok = bool(rd and rd.get("agrees_with_pipeline"))
fpp_ok = bool(fpp and fpp.get("seeds_done", 0) >= 1)
s2_ok = bool(vet and vet.get("stage2_pass"))
label = "candidate" if (pix_ok and eph_ok and rt_ok and rd_ok and fpp_ok and s2_ok) else "signal"


def f(x, nd=1):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:  # noqa: BLE001
        return "n/a"


git = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
tls = (vet or {}).get("tls") or {}
real = ((s1 or {}).get("real") or {}).get("s1") or {}
inv = ((s1 or {}).get("inverted") or {}).get("s1") or {}
real2 = ((s1 or {}).get("real") or {}).get("s2") or {}
ps = pix["summary"] if pix else {}
L = []
L.append(f"# {a.fname} ({star.get('gaia_designation')}, KIC {star.get('kic')}) — evidence pack")
L.append("")
L.append(f"_Built {dt.datetime.utcnow().isoformat(timespec='seconds')}Z from run `{a.run}` at commit `{git[:12]}`. "
         f"Status under `CLAUDE.md` language: **{label}**. Nothing here is a planet; nothing is 'discovered'. "
         "Search, vetting, localisation, re-derivation, red team and FPP were run by Claude Code (Anthropic) agents under human direction; "
         "every number below names the file it comes from._")
L.append("")
L.append("## Key numbers")
L.append("")
L.append("| quantity | value | statistic / subset | file |")
L.append("|---|---|---|---|")
if vet:
    L.append(f"| Period | {f(vet['P'], 6)} d | stage-2 BLS refinement, all quarters | stage2/vet.json `P` |")
    L.append(f"| Mid-transit (BJD) | {f(vet['t0_bjd'], 5)} | stage-2 | stage2/vet.json `t0_bjd` |")
    L.append(f"| Box duration | {f(vet['dur_h'], 2)} h | stage-2 box | stage2/vet.json `dur_h` |")
    L.append(f"| Depth | {f(vet['depth_pass2_ppm'], 0)} ppm (formal ± {f(vet.get('depth_pass2_formal_err_ppm'), 0)}) | masked re-detrend, unweighted box mean | stage2/vet.json `depth_pass2_ppm` |")
    L.append(f"| Empirical SNR | {f(vet['snr_emp_pass2'])} | depth / scatter of 80 random off-transit boxes | stage2/vet.json `snr_emp_pass2` |")
    L.append(f"| Events | {vet.get('n_events')} ({f(100*(vet.get('frac_events_positive') or 0), 0)} % positive) | per-event local depths | stage2/events.csv |")
    L.append(f"| Per-event χ²/dof | {f(vet.get('event_chi2_const'))} / {vet.get('event_dof')} | formal errors understate red noise (controls 1.5–2.0) | stage2/vet.json |")
    L.append(f"| Odd/even | {f(vet.get('odd_even_sigma'), 2)} σ | stage-2 | stage2/vet.json `odd_even_sigma` |")
    L.append(f"| Secondary at phase 0.5 | {f(vet.get('secondary_depth_ppm'), 0)} ppm ({f(vet.get('secondary_sigma'), 2)} σ) | two-sided test | stage2/vet.json |")
    L.append(f"| TLS (±2 % window) | SDE {f(tls.get('SDE'))}, SNR {f(tls.get('snr'))}, {tls.get('distinct_transit_count')} transits | not a calibrated SDE (see prereg v1.5) | stage2/vet.json `tls` |")
    L.append(f"| R_p if on target | {f(vet.get('rp_earth_box'), 2)} R⊕ (box), {f(tls.get('rp_earth'), 2)} R⊕ (TLS) | R★ {f(vet.get('rstar_used'), 2)} R☉ colour-based | stage2/vet.json |")
    L.append(f"| Stage-2 verdict | {'PASS' if vet.get('stage2_pass') else 'FAIL ' + str(vet.get('stage2_fails'))} ({vet.get('stage2_rules')}); v1.4: {'PASS' if vet.get('stage2_pass_v14') else 'FAIL ' + str(vet.get('stage2_fails_v14'))} | post-hoc relative to stage 1 | stage2/vet.json |")
if real:
    L.append(f"| Stage-1 BLS | P {f(real.get('period_d'), 5)} d, SDE {f(real.get('sde'))}, SNR_emp {f(real.get('snr_emp'))}, {real.get('n_transits')} transits, max single-event fraction {f(real.get('max_single_event_frac'), 2)}; inverted (2−f) top peak SDE {f(inv.get('sde'))} at {f(inv.get('period_d'), 3)} d | blind screen, config {s1.get('config_sha')} | stage1_row.json |")
    if real2:
        L.append(f"| Stage-1 second signal (after masking the first) | P {f(real2.get('period_d'), 4)} d, depth {f(real2.get('depth_ppm'), 0)} ppm, SDE {f(real2.get('sde'))}, SNR_emp {f(real2.get('snr_emp'))} | 0.5–40 d window | stage1_row.json `real.s2` |")
if ps:
    L.append(f"| Pixel localisation (v1.6 χ²) | **{ps.get('verdict')}**: Σ(χ²null−χ²target) {f(ps.get('sum_dchi2_null_minus_target'), 0)}, Σ(χ²best-alt−χ²target) {f(ps.get('sum_dchi2_bestalt_minus_target'), 0)} over {ps.get('n_quarters_chi2')} quarters; strongest alternative {ps.get('strongest_alternative')} | post-hoc rule, passed control KIC 5689121 | pixels/pixel_v16.json |")
if rd:
    L.append(f"| Independent re-derivation | {'agrees' if rd.get('agrees_with_pipeline') else 'DISAGREES'}: P {f(rd.get('period_d'), 6)} d, depth {f(rd.get('depth_ppm'), 0)} ± {f(rd.get('depth_err_ppm'), 0)} ppm, SNR {f(rd.get('empirical_snr'))}, {rd.get('n_transits')} events | no shared code | independent/rederive_result.json |")
if rt:
    L.append(f"| Red team | **{rt.get('verdict')}**; kill reasons: {rt.get('kill_reasons') or 'none'} | | independent/redteam_verdict.json, independent/redteam.md |")
if fpp:
    L.append(f"| TRICERATOPS FPP | {f(fpp.get('FPP_mean'), 4)} ± {f(fpp.get('FPP_std'), 4)} ({fpp.get('seeds_done')} seeds × {fpp.get('N_draws')} draws) | Kepler mode, TIC stellar parameters, no contrast curve | fpp/fpp_summary.json |")
    L.append(f"| TRICERATOPS NFPP | {f(fpp.get('NFPP_mean'), 5)} ± {f(fpp.get('NFPP_std'), 5)} | nearby stars excluded only by their own light curves (fpp/fpp_summary.json `neighbour_tests`) | fpp/fpp_summary.json |")
L.append(f"| Host (Gaia DR3) | G {f(star.get('phot_g_mean_mag'), 2)}, BP−RP {f(star.get('BP-RP'), 2)}, Teff {f(star.get('Teff'), 0)} K, R★ {f(star.get('Rstar'), 2)} R☉, {f(star.get('dist_pc'), 0)} pc, RUWE {f(star.get('RUWE'), 2)} | colour-based radius; TIC values in fpp/stars_tic.csv | target_row.json |")
if vet:
    km = vet.get("known_matches") or []
    L.append(f"| Ephemeris matches (KOI/TCE/EB/confirmed within 3′, ratios 1:1,2:1,1:2,3:1,1:3, 0.2 %) | {len(km)} ({'none on the same star or TPF' if eph_ok else 'MATCH ON SAME STAR/TPF'}) | catalogues dated {known.get('fetched_utc') if known else 'n/a'} | stage2/vet.json `known_matches`, known_catalogues_manifest.json |")
fw = jload(f"{R}/ephemeris_fieldwide/{a.fname}.json")
if fw:
    L.append(f"| Field-wide ephemeris check (whole Kepler field, any channel; KOI, DR25 TCE, Kepler EB; harmonics 1:1–4:1) | {fw['n_period_coincidences']} period coincidences within 0.1 %, {fw['n_epoch_matches']} with a matching epoch, **{fw['n_viable_contaminants']} viable** (1:1 twins, or harmonics not excluded by this signal's odd/even and secondary tests) | a 4-year baseline separates period twins whose epochs differ; harmonic twins need an odd/even or phase-0.5 signature this signal lacks | ephemeris_fieldwide.json |")
L.append("")
if second:
    L.append(f"## Second signal (post-hoc, label {a.second_signal})")
    L.append("")
    L.append(f"Found by the red team outside the registered 0.5–40 d window. Stage 2 under current rules: "
             f"P {f(second['P'], 5)} d, depth {f(second['depth_pass2_ppm'], 0)} ppm, SNR {f(second['snr_emp_pass2'])}, {second.get('n_events')} events "
             f"({f(100*(second.get('frac_events_positive') or 0), 0)} % positive), odd/even {f(second.get('odd_even_sigma'), 2)} σ, secondary {f(second.get('secondary_sigma'), 2)} σ, "
             f"verdict {'PASS' if second.get('stage2_pass') else 'FAIL ' + str(second.get('stage2_fails'))}; known matches {len(second.get('known_matches') or [])}. "
             f"Pixel v1.6: {(second_pix or {}).get('summary', {}).get('verdict', 'not run')} ({(second_pix or {}).get('summary', {}).get('verdict_reason', '')}). Files: second_signal_{a.second_signal}/.")
    L.append("")
if deep:
    L.append("## Deep single-target run (`pipeline/run.py`, `deep_run/results.json`)")
    L.append("")
    ae = deep.get("adopted_ephemeris", {}); dv = deep.get("vetting", {}); cmpz = deep.get("comparison", {}); nul = deep.get("inverted_null", {})
    bd = dv.get("box_depth", {}); oe = dv.get("odd_even", {}); sc = dv.get("secondary", {}); en = dv.get("empirical_noise", {})
    nb = deep.get("neighbour", {}) or {}
    L.append("| quantity | value | note |")
    L.append("|---|---|---|")
    L.append(f"| Blind recovery | BLS 1–60 d best P {f((deep.get('bls') or {}).get('period'), 6)} d, SDE {f((deep.get('bls') or {}).get('sde'))}; TLS 1–50 d P {f((deep.get('tls') or {}).get('period'), 6)} d, SDE {f((deep.get('tls') or {}).get('SDE'))}, SNR {f((deep.get('tls') or {}).get('snr'))} | independent detrend window and grids from the batch pipeline |")
    L.append(f"| Adopted ephemeris | P {f(ae.get('period_d'), 6)} d, t0 BJD {f(ae.get('t0_bjd_tdb'), 5)}, duration {f(ae.get('duration_h'), 2)} h | {ae.get('source')} |")
    L.append(f"| Agreement with stage 2 | ΔP {f(cmpz.get('period_diff_ppm'), 1)} ppm, mid-baseline t0 offset {f(cmpz.get('t0_offset_min_at_mid_baseline'), 1)} min; depth at the stage-2 ephemeris {f((cmpz.get('depth_at_reference_ephemeris_pass2') or {}).get('depth'), 0)} ± {f((cmpz.get('depth_at_reference_ephemeris_pass2') or {}).get('depth_err'), 0)} ppm | |")
    L.append(f"| Box depth (pass 2) | {f(1e6*(bd.get('depth') or 0), 0)} ± {f(1e6*(bd.get('depth_err') or 0), 0)} ppm (formal), {bd.get('n_in')} in-transit cadences | empirical SNR {f(en.get('empirical_snr'))} (random-phase boxes: {en.get('n_random')}, std {f(en.get('random_depth_std_ppm'), 0)} ppm) |")
    L.append(f"| Odd / even | {f(1e6*((oe.get('odd') or {}).get('depth') or 0), 0)} / {f(1e6*((oe.get('even') or {}).get('depth') or 0), 0)} ppm, difference {f(oe.get('difference_sigma'), 2)} σ | |")
    L.append(f"| Secondary | phase 0.5: {f(1e6*((sc.get('phase_0p5') or {}).get('depth') or 0), 0)} ppm ({f((sc.get('phase_0p5') or {}).get('sigma'), 2)} σ); deepest other dip in a {sc.get('scan_n_phases')}-phase scan {f((sc.get('scan_max_dip') or {}).get('sigma'), 1)} σ | scan σ scatter {f(sc.get('scan_sigma_std'), 2)} |")
    L.append(f"| Inverted-flux null (masked LC) | BLS best SDE {f((nul.get('bls') or {}).get('sde'))}; TLS best SDE {f((nul.get('tls_excluding_candidate_band') or {}).get('best_sde_overall'))} at {f((nul.get('tls_excluding_candidate_band') or {}).get('period_overall_d'), 3)} d | full-range TLS SDE, not the stage-1 statistic |")
    nbp = nb.get('depth_at_target_ephemeris_pass2') or {}
    L.append(f"| TPF host ({(nb.get('primary_header') or {}).get('OBJECT')}) at the ephemeris | {f(nbp.get('depth'), 0)} ± {f(nbp.get('depth_err'), 0)} ppm ({f(nbp.get('sigma'), 1)} σ) | would need {f(nb.get('depth_required_if_hosted_by_neighbour_ppm'), 0)} ppm to host the signal |")
    L.append("")
L.append("## What this pack does not show")
L.append("")
L.append("- No high-resolution imaging, no radial velocities, no contrast curve: the unresolved-companion scenarios in the FPP carry their prior weight only.")
L.append("- The pixel model is a Gaussian point source, not the Kepler PRF; many background targets sit at the stamp edge (quarters flagged).")
L.append("- Stage-2 rules (v1.3–v1.7) and the pixel rule (v1.6) were registered after results had been seen; every such result is post-hoc relative to stage 1 (prereg_poolA.json `revision_log`).")
L.append("- 'Validated' requires FPP < 0.015 and NFPP < 0.001 in a refereed paper; this pack is a pipeline product.")
L.append("- Novelty holds for the catalogues listed in known_catalogues_manifest.json on their download date; the Kepler Bonus team's own candidate list is unpublished.")
L.append("")
L.append("## Files")
L.append("")
L.append("`SHA256SUMS` lists every file; `manifest.json` adds sizes and the git commit.")
open(os.path.join(out, "README.md"), "w").write("\n".join(L) + "\n")

# ---- manifest
entries = []
for root, _, fs in os.walk(out):
    for fn in fs:
        if fn in ("SHA256SUMS", "manifest.json"):
            continue
        p = os.path.join(root, fn)
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        entries.append({"path": os.path.relpath(p, out), "bytes": os.path.getsize(p), "sha256": h})
entries.sort(key=lambda e: e["path"])
open(os.path.join(out, "SHA256SUMS"), "w").write("".join(f"{e['sha256']}  {e['path']}\n" for e in entries))
json.dump({"fname": a.fname, "run_id": a.run, "built_utc": dt.datetime.utcnow().isoformat() + "Z", "git_commit": git, "status_label": label,
           "checks": {"stage2_pass": s2_ok, "pixel_on_target": pix_ok, "no_ephemeris_match_same_star_or_tpf": eph_ok, "redteam_survives": rt_ok,
                      "rederivation_agrees": rd_ok, "fpp_run_exists": fpp_ok},
           "prereg_sha256": hashlib.sha256(open("pipeline/prereg_poolA.json", "rb").read()).hexdigest(),
           "files": entries}, open(os.path.join(out, "manifest.json"), "w"), indent=1)
print(f"{out}: {len(entries)} files, status {label}; checks:", json.dumps({"stage2": s2_ok, "pixel": pix_ok, "ephemeris": eph_ok, "redteam": rt_ok, "rederive": rd_ok, "fpp": fpp_ok}))
