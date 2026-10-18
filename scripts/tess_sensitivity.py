"""Can TESS test the three candidates? Upper-bound transit S/N in all TESS sectors to date.

    python scripts/tess_sensitivity.py      # -> runs/poolA_main/tess_sensitivity.json

Tmag from the TIC (MAST), sectors from the TESScut sector service, one-hour noise from the Stassun et al. (2018)
TIC noise model (ticgen 1.0.4), which assumes an uncrowded star and an ideal aperture. S/N = depth / sigma_1h *
sqrt(T14 [h] * N), with N = (sectors x 27.4 d) / P, i.e. full coverage and no dilution by neighbours inside TESS's
21" pixels: every simplification raises the S/N, so the result is an upper bound.
"""
import datetime as dt
import json
import math
import urllib.request

from astroquery.mast import Catalogs
from ticgen import Star

D = json.load(open("report/data.json"))
out = {"computed_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "noise_model": "Stassun et al. 2018 (ticgen 1.0.4)",
       "assumptions": "full sector coverage (27.4 d per sector), no dilution, ideal aperture: S/N is an upper bound", "candidates": {}}
for k in ["kic-006699887", "kic-003831969", "kic-011138423"]:
    c = D["candidates"][k]
    fit = json.load(open(f"runs/poolA_main/transit_fit/{k}/fit.json"))
    tmag = float(Catalogs.query_criteria(catalog="Tic", ID=int(c["tic"]))["Tmag"][0])
    secs = json.load(urllib.request.urlopen(f"https://mast.stsci.edu/tesscut/api/v0.1/sector?ra={c['ra']}&dec={c['dec']}&radius=0", timeout=60))
    sectors = sorted(int(s["sector"]) for s in secs["results"])
    sig1h = float(Star(Tmag=tmag).oneSigmaNoise)  # ppm in one hour (integration default 60 min)
    t14 = fit["T14_h"]["median"]
    n_tr = len(sectors) * 27.4 / c["period_d"]
    snr = c["depth_ppm"] / sig1h * math.sqrt(t14 * n_tr)
    out["candidates"][k] = {"tic": c["tic"], "Tmag": round(tmag, 3), "sectors": sectors, "n_sectors": len(sectors), "sigma_1h_ppm": round(sig1h),
                            "depth_ppm": round(c["depth_ppm"]), "T14_h": round(t14, 2), "n_transits_max": round(n_tr, 1),
                            "snr_upper_bound": round(snr, 1)}
    print(k, out["candidates"][k])
json.dump(out, open("runs/poolA_main/tess_sensitivity.json", "w"), indent=1)
