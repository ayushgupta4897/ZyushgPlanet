"""Build the Zenodo data package the RNAAS note cites -> papers/rnaas/zenodo/ (zip + .zenodo.json + SHA256SUMS).

    python papers/rnaas/build_zenodo.py [--author "Name"] [--orcid 0000-...] [--affiliation "..."]

Contents: the three evidence packs, report/data.json, the note source and figure, the pipeline and scripts at this commit,
the frozen target list and pre-registration, pool-level run summaries (first-screen results with the inverted-flux
null, triage, stage-2 table), and the three Kepler Bonus light-curve FITS files (MAST HLSP, CC BY 4.0, DOI 10.17909/7jbr-w430).
The zip is not committed to git; upload it and .zenodo.json to Zenodo from your own account.
"""
import argparse
import hashlib
import json
import os
import subprocess
import zipfile

ap = argparse.ArgumentParser()
ap.add_argument("--author", default="[AUTHOR NAME]")
ap.add_argument("--orcid", default="")
ap.add_argument("--affiliation", default="Independent researcher")
a = ap.parse_args()
root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
out = os.path.join(root, "papers", "rnaas", "zenodo")
os.makedirs(out, exist_ok=True)
commit = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
K = ["kic-006699887", "kic-003831969", "kic-011138423"]
files = []
for d in ["candidates", "pipeline", "scripts", "papers/rnaas"]:
    for r, _, fs in os.walk(os.path.join(root, d)):
        if "/zenodo" in r or "__pycache__" in r or "/build" in r:
            continue
        files += [os.path.join(r, f) for f in fs if not f.endswith((".cls",))]
files += [os.path.join(root, p) for p in ["report/data.json", "pipeline/prereg_poolA.json", "targets/poolA_kbonus_bg_dwarfs_ranked.csv",
                                          "runs/poolA_main/results.jsonl", "runs/poolA_pilot/results.jsonl", "runs/poolA_main/triage.md",
                                          "runs/poolA_main/vet_summary.csv", "runs/poolA_main/survivors_stage2.csv", "runs/poolA_main/modal_status.json",
                                          "runs/poolA_main/config_used.json", "data/known/manifest.json", "CLAUDE.md"]]
files += [os.path.join(root, "data/kbonus_batch", f"hlsp_kbonus-bkg_kepler_kepler_{k}_kepler_v1.0_lc.fits") for k in K]
files = sorted({f for f in files if os.path.isfile(f)})
zpath = os.path.join(out, f"zyushgplanet-kbonus-candidates-{commit[:8]}.zip")
sums = []
with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
    for f in files:
        rel = os.path.relpath(f, root)
        z.write(f, rel)
        sums.append(f"{hashlib.sha256(open(f, 'rb').read()).hexdigest()}  {rel}")
    z.writestr("SHA256SUMS", "\n".join(sums) + "\n")
    z.writestr("README.txt", f"ZyushgPlanet Kepler Bonus candidate package, git commit {commit}.\nStart with candidates/<fname>/README.md and report/data.json.\n"
                             "Kepler Bonus light curves: Martinez-Palomera, Hedges & Dotson 2023, AJ 166, 265; MAST HLSP DOI 10.17909/7jbr-w430 (CC BY 4.0).\n"
                             "Analysis code, re-derivations, reviews and drafts were produced with Claude Code (Anthropic) under the author's direction.\n")
open(os.path.join(out, "SHA256SUMS"), "w").write("\n".join(sums) + "\n")
creator = {"name": a.author, "affiliation": a.affiliation}
if a.orcid:
    creator["orcid"] = a.orcid
meta = {"upload_type": "dataset", "title": "Data, code and evidence for three sub-Neptune transit candidates around Kepler Bonus background stars",
        "creators": [creator], "license": "cc-by-4.0", "version": commit[:8], "keywords": ["exoplanets", "transits", "Kepler", "Kepler Bonus", "planet candidates"],
        "description": "Evidence packs, run outputs, code and light curves for KIC 6699887, KIC 3831969 and KIC 11138423, three transit candidates from a pre-registered search of 6,465 Kepler Bonus background dwarfs. Candidates, not validated planets. Analysis performed with Claude Code (Anthropic) under the author's direction.",
        "related_identifiers": [{"identifier": "10.17909/7jbr-w430", "relation": "isDerivedFrom", "scheme": "doi"}]}
json.dump(meta, open(os.path.join(out, ".zenodo.json"), "w"), indent=1)
print(f"{zpath}: {len(files)} files, {os.path.getsize(zpath) / 1e6:.1f} MB")
