"""Run the stage-2 vetting battery on a hand-supplied ephemeris (POST-HOC; not part of the registered screen).

    python scripts/vet_posthoc_signal.py poolA_main kic-006699887 47.6906 4.5 135.7185 --label P47

Used for signals that the registered 0.5-40 d search could not reach (e.g. the 47.69-d second signal on
kic-006699887 found by the red team at its P/2 alias). Output goes to runs/<run>/vet_posthoc/<fname>_<label>/
and is labelled post-hoc in every write-up; the stage-2 verdict is reported under the current prereg rules.
"""
import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.batch_search import load_config  # noqa: E402
from pipeline.vet_batch import vet_one  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("fname")
ap.add_argument("P", type=float)
ap.add_argument("dur_h", type=float)
ap.add_argument("t0_bkjd", type=float)
ap.add_argument("--label", default="posthoc")
ap.add_argument("--config", default="pipeline/prereg_poolA.json")
ap.add_argument("--targets", default="targets/poolA_kbonus_bg_dwarfs_ranked.csv")
a = ap.parse_args()

cfg = load_config(a.config)
tg = pd.read_csv(a.targets)
star = tg[tg.fname == a.fname].iloc[0].to_dict()
out = f"runs/{a.run}/vet_posthoc/{a.fname}_{a.label}"
res = vet_one(a.fname, a.P, a.dur_h, a.t0_bkjd, star, cfg, out, make_plots=True)
res["post_hoc"] = True
res["post_hoc_note"] = "Ephemeris supplied by hand (outside the registered 0.5-40 d search); stage-2 rules applied as registered at run time."
with open(os.path.join(out, "vet.json"), "w") as fh:
    json.dump(res, fh, indent=1, default=str)
keys = ["P", "dur_h", "t0_bjd", "depth_pass2_ppm", "snr_emp_pass2", "n_events", "frac_events_positive", "event_chi2_const", "event_dof",
        "odd_even_sigma", "secondary_sigma", "core_to_box_depth_ratio", "gap_frac_1d", "duration_ratio", "rp_earth_box", "stage2_pass", "stage2_fails", "stage2_rules"]
print(json.dumps({k: res.get(k) for k in keys}, indent=1, default=str))
tls = res.get("tls") or {}
print("TLS:", {k: tls.get(k) for k in ("period", "SDE", "snr", "depth", "duration", "distinct_transit_count")})
print("known_matches:", res.get("known_matches"))
