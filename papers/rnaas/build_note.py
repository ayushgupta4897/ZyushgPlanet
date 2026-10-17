"""Fill papers/rnaas/note.tex.tpl from the run outputs, compile it, and count words the way RNAAS does.

    python papers/rnaas/build_note.py [--author "Name"] [--orcid 0000-...] [--affiliation "..."] [--doi 10.5281/zenodo.N]

Every number comes from report/data.json (exported from runs/ by scripts/report_data.py) or runs/poolA_main/transit_fit/*/fit.json;
the filled values and their sources are written to papers/rnaas/numbers.json.
"""
import argparse
import json
import os
import re
import subprocess

ap = argparse.ArgumentParser()
ap.add_argument("--author", default="[AUTHOR NAME]")
ap.add_argument("--orcid", default="0000-0000-0000-0000")
ap.add_argument("--affiliation", default="Independent researcher")
ap.add_argument("--email", default="[author email]")
ap.add_argument("--doi", default="[Zenodo DOI to be added]")
ap.add_argument("--tectonic", default="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/tex/tectonic")
ap.add_argument("--aastex", default="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/tex/aastex/aastex701.cls")
a = ap.parse_args()
here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", ".."))
D = json.load(open(os.path.join(root, "report", "data.json")))
K = ["kic-006699887", "kic-003831969", "kic-011138423"]
C = [D["candidates"][k] for k in K]
F = [json.load(open(os.path.join(root, f"runs/poolA_main/transit_fit/{k}/fit.json"))) for k in K]
f4 = json.load(open(os.path.join(root, "runs/poolA_main/fpp/kic-003963519/fpp_summary.json")))


def pm(d, nd, scale=1.0):
    m, p, n = d["median"] * scale, d["plus"] * scale, d["minus"] * scale
    return f"{m:.{nd}f}^{{+{p:.{nd}f}}}_{{-{n:.{nd}f}}}"


V = {"AUTHOR": a.author, "EMAIL": a.email, "ORCID": a.orcid, "AFFILIATION": a.affiliation, "DATA_DOI": a.doi,
     "n_targets": f"{D['funnel']['targets']:,}", "n_s1_real": D["funnel"]["first_screen_real"], "n_s1_inv": D["funnel"]["first_screen_inverted"],
     "w_pix": ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven"][D["funnel"]["path_pixel_on_target"]],
     "w_pix_l": ["zero", "one", "two", "three", "four", "five", "six", "seven"][D["funnel"]["path_pixel_on_target"]],
     "w_rt": ["zero", "one", "two", "three", "four", "five", "six", "seven"][D["funnel"]["path_redteam_survives"]],
     "n_stage2": D["funnel"]["path_stage2"], "n_pix": D["funnel"]["path_pixel_on_target"], "n_rt": D["funnel"]["path_redteam_survives"],
     "host_excl_min": f"{min(c['fpp']['host_excluded_sigma'] for c in C):.0f}",
     "ruwe_min": f"{min(c['host']['ruwe'] for c in C):.2f}", "ruwe_max": f"{max(c['host']['ruwe'] for c in C):.2f}",
     "prior_def_min": f"{100 * min(c['fpp']['companion_prior_default'] for c in C):.0f}", "prior_def_max": f"{100 * max(c['fpp']['companion_prior_default'] for c in C):.0f}",
     "prior_mol_min": f"{100 * min(c['fpp']['companion_prior_molusc'] for c in C):.0f}", "prior_mol_max": f"{100 * max(c['fpp']['companion_prior_molusc'] for c in C):.0f}",
     "fppmin": f"{100 * min(c['fpp']['molusc']['mean'] for c in C):.1f}", "fppmax": f"{100 * max(c['fpp']['molusc']['mean'] for c in C):.1f}",
     "fpp4": f"{100 * f4['FPP_mean']:.0f}", "dens_max": f"{max(abs(c['density']['tension_sigma']) for c in C):.1f}",
     "P47": f"{D['second_signal']['period_d']:.4f}", "P47s": f"{D['second_signal']['period_d']:.2f}", "dep47": f"{D['second_signal']['depth_ppm']:,.0f}", "n47": D["second_signal"]["n_events"]}
for i, (c, f) in enumerate(zip(C, F), 1):
    V.update({f"P{i}": pm(f["P_d"], 5), f"P{i}s": f"{c['period_d']:.3f}", f"T0{i}": pm(f["t0_bjd"], 4), f"dep{i}": f"{c['depth_ppm']:,.0f}",
              f"T14{i}": pm(f["T14_h"], 2), f"rp{i}": pm(c["rp_earth"], 2), f"rp{i}s": f"{c['rp_earth']['median']:.1f}",
              f"fppb{i}": f"{100 * c['fpp']['baseline']['mean']:.1f}", f"fppm{i}": f"{100 * c['fpp']['molusc']['mean']:.2f}", f"fpps{i}": f"{100 * c['fpp']['molusc']['std']:.2f}",
              f"gaia{i}": c["gaia_dr3"], f"tic{i}": c["tic"], f"teff{i}": f"{c['host']['teff_gspphot']:,.0f}", f"dist{i}": f"{c['host']['dist_pc']:,.0f}"})
tpl = open(os.path.join(here, "note.tex.tpl")).read()
missing = sorted(set(re.findall(r"<<(\w+)>>", tpl)) - set(V))
assert not missing, missing
tex = re.sub(r"<<(\w+)>>", lambda m: str(V[m.group(1)]), tpl)
open(os.path.join(here, "note.tex"), "w").write(tex)
json.dump({"values": {k: str(v) for k, v in V.items()}, "sources": ["report/data.json (scripts/report_data.py; sources listed inside)"] +
           [f"runs/poolA_main/transit_fit/{k}/fit.json" for k in K] + ["runs/poolA_main/fpp/kic-003963519/fpp_summary.json"],
           "data_json_commit": D["generated_from_commit"]}, open(os.path.join(here, "numbers.json"), "w"), indent=1)

# RNAAS-style word count: title, abstract, headers, body, caption and references; LaTeX commands and math removed
body = tex.split("\\begin{document}")[1].split("\\end{document}")[0]
abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", body, re.S).group(1)


def words(s):
    s = re.sub(r"%.*", "", s)
    s = re.sub(r"\\(author|affiliation|label|plotone|bibitem)\[[^\]]*\]\{[^}]*\}|\\(affiliation|label|plotone)\{[^}]*\}", " ", s)
    s = re.sub(r"\$[^$]*\$", " X ", s)
    s = re.sub(r"\\citep\{[^}]*\}", " (Ref. 2000) ", s)
    s = re.sub(r"\\citeyear\{[^}]*\}", " 2000 ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", s)
    s = re.sub(r"[{}~\\]", " ", s)
    return len([w for w in s.split() if re.search(r"[A-Za-z0-9]", w)])


n_total, n_abs = words(body), words(abstract)
print(f"words: total {n_total} (limit 1,500 incl. title, headers, caption, references), abstract {n_abs} (limit 150)")
os.makedirs(os.path.join(here, "build"), exist_ok=True)
import shutil
shutil.copy(a.aastex, os.path.join(here, "aastex701.cls"))
r = subprocess.run([a.tectonic, "--outdir", os.path.join(here, "build"), os.path.join(here, "note.tex")], capture_output=True, text=True, cwd=here)
print(r.stdout[-400:], r.stderr[-1500:])
json.dump({"total_words": n_total, "abstract_words": n_abs}, open(os.path.join(here, "word_count.json"), "w"))
