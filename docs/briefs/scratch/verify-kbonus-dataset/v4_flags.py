import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
src=pd.read_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/source_catalog.parquet")
kep=set(pd.read_csv("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/keplerstellar_kepids.csv").kepid.astype(int))
A=np.array([list(s) for s in src.psf_avail]); F=np.array([list(s) for s in src.psf_fraction_flag]); S=np.array([list(s) for s in src.sap_avail])
print("psf_avail x psf_fraction_flag joint counts (per quarter-slot):")
print(pd.crosstab(pd.Series(A.ravel(),name='psf_avail'), pd.Series(F.ravel(),name='frac_flag')))
print("sap_avail x psf_avail:"); print(pd.crosstab(pd.Series(S.ravel(),name='sap_avail'), pd.Series(A.ravel(),name='psf_avail')))
t=src[src.gaia_designation.str.contains('2103792461826460288')].iloc[0]
print("TARGET psf_avail", t.psf_avail, "sap_avail", t.sap_avail, "flag", t.psf_fraction_flag, "psffrac", t.psffrac, "crowdsap", t.crowdsap, "psf_mean_flux", t.psf_mean_flux, "kic", t.kic, "tpforg", t.tpforg, "kic_sep", t.kic_sep, "kep_mag", t.kep_mag)
print("QDETECT from FITS was 011110111011101110")
n=src[src.gaia_designation.str.contains('2103792427466721792')].iloc[0]
print("NEIGHBOUR psf_avail", n.psf_avail, "flag", n.psf_fraction_flag, "psffrac", n.psffrac, "crowdsap", n.crowdsap, "G", n.phot_g_mean_mag, "kic", n.kic, "tpforg", n.tpforg)
# Kepler target definitions
tp=set(src.tpforg.unique()); print("distinct tpforg values (TPFs processed):", len(tp))
src['is_ks']=src.kic.isin(kep); src['is_tp']=src.kic.isin(tp)
print("kic in keplerstellar:", src.is_ks.sum(), "| kic in tpforg-set:", src.is_tp.sum(), "| union:", (src.is_ks|src.is_tp).sum(), "| paper says 200,352 targets / 406,548 background; background by union:", (~(src.is_ks|src.is_tp)).sum())
print("kic in tpforg-set but not keplerstellar:", (src.is_tp&~src.is_ks).sum(), "; in keplerstellar but not tpforg-set:", (src.is_ks&~src.is_tp).sum())
print("kic>0 but neither:", ((src.kic>0)&~src.is_ks&~src.is_tp).sum(), " (KIC stars that were never Kepler targets)")
# psf_mean_flux NaN
print("psf_mean_flux NaN:", src.psf_mean_flux.isna().sum(), "sap_mean_flux NaN:", src.sap_mean_flux.isna().sum())
# M-stars with union definition
ms=pd.read_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/mstars_joined.parquet")
ms['nq']=ms.psf_avail.str.count('1')
e=ms[(ms.psffrac>0.5)&(ms.crowdsap>0.7)&(ms.nq>=8)&(ms.kic!=ms.tpforg)]
print("M-star chain after kic!=tpforg:", len(e), "; minus keplerstellar:", (~e.kic.isin(kep)).sum(), "; minus union:", (~(e.kic.isin(kep)|e.kic.isin(tp))).sum())
print("M-star entries in tpforg-set but not keplerstellar (and kic!=tpforg):", (e.kic.isin(tp)&~e.kic.isin(kep)).sum())
# Background quality pool with union
nq=src.psf_avail.str.count('1'); ng=src.psf_fraction_flag.str.count('1')
bkg=~(src.is_ks|src.is_tp)
for name,q in [("nq_any>=8",nq>=8),("nq_good>=8",ng>=8)]:
    print(f"background(union) & psffrac>0.5 & crowdsap>0.7 & {name}:", (bkg&(src.psffrac>0.5)&(src.crowdsap>0.7)&q).sum(), " G<17:", (bkg&(src.psffrac>0.5)&(src.crowdsap>0.7)&q&(src.phot_g_mean_mag<17)).sum())
