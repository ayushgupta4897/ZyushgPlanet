import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from astropy.io import fits
from astropy.table import Table
C="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/catalogs/"
src=Table.read(C+"hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits").to_pandas()
for c in src.columns:
    if src[c].dtype==object: src[c]=src[c].str.decode('utf-8') if isinstance(src[c].iloc[0],bytes) else src[c]
print("SOURCE CATALOG rows", len(src)); print(src.dtypes); print(src.head(3).to_string())
print("\nexample psf_avail / sap_avail / psf_fraction_flag:", src.psf_avail.iloc[0], src.sap_avail.iloc[0], src.psf_fraction_flag.iloc[0], "len", len(src.psf_avail.iloc[0]))
bkg = src.kic!=src.tpforg
print("\nbackground (kic != tpforg):", bkg.sum(), " kepler-target (kic==tpforg):", (~bkg).sum())
print("kic==0 or NaN count:", (src.kic<=0).sum(), src.kic.isna().sum())
print("\nG mag describe all:\n", src.phot_g_mean_mag.describe())
print("\nG hist (all / background):")
bins=np.arange(9,19.5,1)
h1,_=np.histogram(src.phot_g_mean_mag,bins); h2,_=np.histogram(src.phot_g_mean_mag[bkg],bins)
for b,a,c in zip(bins[:-1],h1,h2): print(f"  G {b:4.0f}-{b+1:4.0f}: all {a:7d}  bkg {c:7d}")
nq = src.psf_avail.str.count('1')
print("\nN PSF quarters available distribution (all):"); print(nq.value_counts().sort_index().to_string())
print("\nN PSF quarters (background only):"); print(nq[bkg].value_counts().sort_index().to_string())
print("\npsffrac describe:", src.psffrac.describe().to_string()); print("crowdsap describe:", src.crowdsap.describe().to_string()); print("flfrcsap describe:", src.flfrcsap.describe().to_string())
print("\nfraction psffrac>0.5:", (src.psffrac>0.5).mean(), " crowdsap>0.7:", (src.crowdsap>0.7).mean(), " both:", ((src.psffrac>0.5)&(src.crowdsap>0.7)).mean())
print("background & psffrac>0.5 & crowdsap>0.7 & nq>=8:", (bkg&(src.psffrac>0.5)&(src.crowdsap>0.7)&(nq>=8)).sum())
print("background & psffrac>0.5 & crowdsap>0.7 & nq>=8 & G<17:", (bkg&(src.psffrac>0.5)&(src.crowdsap>0.7)&(nq>=8)&(src.phot_g_mean_mag<17)).sum())
# per quarter availability counts
qa=np.array([[int(ch) for ch in s] for s in src.psf_avail])
print("\nper-quarter PSF availability counts Q0..Q17:", qa.sum(0).tolist())
qb=np.array([[int(ch) for ch in s] for s in src.sap_avail]); print("per-quarter SAP availability Q0..Q17:", qb.sum(0).tolist())
# colour
src['bp_rp']=src.phot_bp_mean_mag-src.phot_rp_mean_mag
print("\nbp_rp describe:", src.bp_rp.describe().to_string())
# mstars
ms=Table.read(C+"hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits").to_pandas()
for c in ms.columns:
    if ms[c].dtype==object and isinstance(ms[c].iloc[0],bytes): ms[c]=ms[c].str.decode('utf-8')
print("\nMSTARS rows", len(ms), ms.columns.tolist()); print(ms.head(3).to_string())
m=src.merge(ms[['gaia_designation']],on='gaia_designation')
print("mstars matched in source cat:", len(m))
mb = m.kic!=m.tpforg
print("mstars background:", mb.sum(), " kepler-target mstars:", (~mb).sum())
print("mstars G describe:\n", m.phot_g_mean_mag.describe().to_string())
print("mstars bp_rp describe:\n", m.bp_rp.describe().to_string())
print("mstars psffrac>0.5:", (m.psffrac>0.5).sum(), " +crowdsap>0.7:", ((m.psffrac>0.5)&(m.crowdsap>0.7)).sum(), " +nq>=8:", ((m.psffrac>0.5)&(m.crowdsap>0.7)&(m.psf_avail.str.count('1')>=8)).sum(), " +bkg:", ((m.psffrac>0.5)&(m.crowdsap>0.7)&(m.psf_avail.str.count('1')>=8)&mb).sum())
kn=Table.read(C+"hlsp_kbonus-bkg_kepler_kepler_kois-nns-catalog_kepler_v1.0_cat.fits").to_pandas(); print("\nKOIS-NNS rows", len(kn), kn.columns.tolist())
wd=Table.read(C+"hlsp_kbonus-bkg_kepler_kepler_wd-catalog_kepler_v1.0_cat.fits").to_pandas(); print("WD rows", len(wd), wd.columns.tolist())
# target row
t=src[src.gaia_designation.str.contains('2103792461826460288')]
print("\nTARGET ROW:\n", t.T.to_string())
n=src[src.gaia_designation.str.contains('2103792427466721792')]
print("\nNEIGHBOUR ROW:\n", n.T.to_string())
# file names pattern check
print("\nfile_name examples:", src.fname.iloc[0], src.fname[bkg].iloc[0], src[src.kic<=0].fname.iloc[0] if (src.kic<=0).any() else None)
# save parquet for later
src.to_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/source_catalog.parquet")
m.to_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/mstars_joined.parquet")
