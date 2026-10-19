"""Target-selection cut chain with counts, Kepler Bonus source catalogue -> the frozen Pool A list.

    python scripts/target_selection_counts.py      # -> docs/data/target_selection_counts.json

Inputs: data/kbonus-catalogs/source_catalog_annotated.parquet (KBONUS-BKG v1.0 source catalogue with is_bg and the
Mehta-searched flag), data/known/kepler_target_kepids.csv (pipeline.fetch_known) and the frozen target list. The
last cut (Gaia DR3 dwarf, prereg target_list_rule) is checked as a subset relation against the frozen list.
"""
import json

import pandas as pd

a = pd.read_parquet("data/kbonus-catalogs/source_catalog_annotated.parquet")
kt = set(pd.read_csv("data/known/kepler_target_kepids.csv").iloc[:, 0].astype(int))
tg = pd.read_csv("targets/poolA_kbonus_bg_dwarfs_ranked.csv")
cuts = [("background source (kic != tpforg)", lambda d: d.is_bg), ("PSF fraction > 0.5", lambda d: d.psffrac > 0.5),
        ("crowding metric > 0.7", lambda d: d.crowdsap > 0.7), (">= 8 PSF quarters", lambda d: d.nq_psf >= 8),
        ("G < 17.5", lambda d: d.phot_g_mean_mag < 17.5), ("never a Kepler target", lambda d: ~d.kic.fillna(-1).astype(int).isin(kt)),
        ("not among the 1,800 sources searched by Mehta (2026)", lambda d: ~d.kuber_searched)]
m, chain = a, [["Kepler Bonus sources", len(a)]]
for name, f in cuts:
    m = m[f(m)]
    chain.append([name, len(m)])
assert set(tg.gaia_designation) <= set(m.gaia_designation)
chain.append(["Gaia DR3 dwarf (frozen list)", len(tg)])
json.dump({"chain": chain}, open("docs/data/target_selection_counts.json", "w"), indent=1)
for n, c in chain:
    print(f"{c:>8,}  {n}")
