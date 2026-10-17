import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from astropy.table import Table
from astropy.io import fits
C="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/catalogs/"
for f in ["source-catalog","mstars-catalog","kois-nns-catalog","wd-catalog"]:
    p=C+f"hlsp_kbonus-bkg_kepler_kepler_{f}_kepler_v1.0_cat.fits"
    with fits.open(p) as h:
        print(f, "rows", h[1].header['NAXIS2'], "cols", h[1].columns.names)
src=pd.read_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/source_catalog.parquet")
print(src.dtypes)
print("G min/max", src.phot_g_mean_mag.min(), src.phot_g_mean_mag.max())
print("kic>0 count", (src.kic>0).sum(), "kic==0", (src.kic==0).sum())
nq=src.psf_avail.str.count('1'); ng=src.psf_fraction_flag.str.count('1')
print("psf_avail quarter count dist (1,17,18):", (nq==1).sum(), (nq==17).sum(), (nq==18).sum())
print("psf_fraction_flag chars:", sorted(set(''.join(src.psf_fraction_flag.sample(5000, random_state=0)))))
print("psffrac>0.5 frac", (src.psffrac>0.5).mean(), "crowdsap>0.7", (src.crowdsap>0.7).mean(), "both", ((src.psffrac>0.5)&(src.crowdsap>0.7)).mean())
qa=np.array([[int(ch) for ch in s] for s in src.psf_avail]); print("per-quarter PSF availability:", qa.sum(0).tolist())
# psf_avail vs psf_fraction_flag consistency
print("rows where psf_avail count != (flag count of 1 + flag count of 2):", (nq != (ng + src.psf_fraction_flag.str.count('2'))).sum())
print("rows where nq_any>=8 but nq_good<8:", ((nq>=8)&(ng<8)).sum())
# M stars chain
ms=pd.read_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/mstars_joined.parquet")
kep=set(pd.read_csv("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/keplerstellar_kepids.csv").kepid.astype(int))
print("keplerstellar distinct kepids:", len(kep))
ms['nq']=ms.psf_avail.str.count('1'); ms['ng']=ms.psf_fraction_flag.str.count('1')
a=ms[ms.psffrac>0.5]; b=a[a.crowdsap>0.7]; c=b[b.nq>=8]; d=c[c.kic!=c.tpforg]; e=d[~d.kic.isin(kep)]
print("chain:", len(ms), len(a), len(b), len(c), len(d), len(e))
c2=b[b.ng>=8]; print("if >=8 counted with psf_fraction_flag=='1' instead:", len(c2), len(c2[(c2.kic!=c2.tpforg)&~c2.kic.isin(kep)]))
e=e.sort_values(['phot_g_mean_mag','gaia_designation'])
print("G of rank1,64,576,1800:", e.phot_g_mean_mag.iloc[[0,63,575,len(e)-1]].tolist())
print("first 5:", e[['gaia_designation','kic','tpforg','phot_g_mean_mag','psffrac','crowdsap','nq']].head(5).to_string())
print("target in eligible? rank:", int(np.where(e.gaia_designation.str.contains('2103792461826460288'))[0][0])+1 if e.gaia_designation.str.contains('2103792461826460288').any() else None)
# M-star G stats
print("mstars G median", ms.phot_g_mean_mag.median(), "quantiles", ms.phot_g_mean_mag.quantile([.25,.5,.75]).tolist())
print("mstars kic==tpforg:", (ms.kic==ms.tpforg).sum(), "in keplerstellar:", ms.kic.isin(kep).sum())
# whole catalogue background accounting
src['is_kt']=src.kic.isin(kep)
print("Kepler targets (keplerstellar):", src.is_kt.sum(), "non", (~src.is_kt).sum(), "; kic==tpforg", (src.kic==src.tpforg).sum(), "; kic!=tpforg & is_kt", ((src.kic!=src.tpforg)&src.is_kt).sum(), "; kic==tpforg & ~is_kt", ((src.kic==src.tpforg)&~src.is_kt).sum())
pool=(~src.is_kt)&(src.psffrac>0.5)&(src.crowdsap>0.7)&(nq>=8)
print("non-target quality pool (nq_any>=8):", pool.sum(), "; with ng>=8:", ((~src.is_kt)&(src.psffrac>0.5)&(src.crowdsap>0.7)&(ng>=8)).sum())
print("pool G<17:", (pool&(src.phot_g_mean_mag<17)).sum())
