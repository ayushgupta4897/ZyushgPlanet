import numpy as np, pandas as pd
src=pd.read_parquet("source_catalog.parquet"); g=pd.read_parquet("mstars_gaia.parquet")
kep=set(pd.read_csv("keplerstellar_kepids.csv").kepid.astype(int))
src['is_kepler_target']=src.kic.isin(kep)
src['nq']=src.psf_avail.str.count('1')
print("Kepler-target (kic in keplerstellar):", src.is_kepler_target.sum(), " not Kepler target:", (~src.is_kepler_target).sum(), " (paper: 200,352 / 406,548)")
print("kic==tpforg:", (src.kic==src.tpforg).sum(), "; kic!=tpforg but kic in keplerstellar:", ((src.kic!=src.tpforg)&src.is_kepler_target).sum(), "; kic==tpforg but NOT in keplerstellar:", ((src.kic==src.tpforg)&~src.is_kepler_target).sum())
bkg=~src.is_kepler_target
q=(src.psffrac>0.5)&(src.crowdsap>0.7)&(src.nq>=8)
print("\nBackground & psffrac>0.5 & crowdsap>0.7 & nq>=8:", (bkg&q).sum())
src['bp_rp']=src.phot_bp_mean_mag-src.phot_rp_mean_mag
src.loc[src.bp_rp.abs()>50,'bp_rp']=np.nan
print("  of which with valid bp_rp:", (bkg&q&src.bp_rp.notna()).sum())
print("\nBackground sources passing quality cuts, by G bin and colour class (bp_rp: <0.8 'blue/FG', 0.8-1.4 'GK', 1.4-1.84 'late K', >=1.84 'M' per Gaia DR3 approx):")
bins=[10,12,13,14,15,16,17,18,19]
cc=pd.cut(src.bp_rp,[-5,0.8,1.4,1.84,10],labels=['<0.8','0.8-1.4','1.4-1.84','>=1.84'])
tab=pd.crosstab(pd.cut(src.phot_g_mean_mag[bkg&q],bins), cc[bkg&q], margins=True)
print(tab.to_string())
print("\nSame, ALL quality (background only, any psffrac/crowdsap/nq):")
print(pd.crosstab(pd.cut(src.phot_g_mean_mag[bkg],bins), cc[bkg], margins=True).to_string())
# M-star catalogue: Kuber's eligible set; G range
ms=pd.read_parquet("mstars_joined.parquet"); ms['nq']=ms.psf_avail.str.count('1'); ms['sid']=ms.gaia_designation.str.split().str[-1].astype(np.int64)
ms=ms.merge(g.rename(columns={'phot_g_mean_mag':'G_gaia'}),left_on='sid',right_on='source_id',how='left')
ms['is_kepler_target']=ms.kic.isin(kep)
elig=ms[(ms.psffrac>0.5)&(ms.crowdsap>0.7)&(ms.nq>=8)&(ms.kic!=ms.tpforg)&(~ms.is_kepler_target)].sort_values(['phot_g_mean_mag','gaia_designation'])
print("\nM-star eligible (Kuber cuts, primary excluded via keplerstellar):", len(elig), " (Kuber: 1800)")
print("  G range of eligible:", elig.phot_g_mean_mag.min(), elig.phot_g_mean_mag.max(), "; G at rank 64:", elig.phot_g_mean_mag.iloc[63], "; rank 576:", elig.phot_g_mean_mag.iloc[575], "; rank 1736:", elig.phot_g_mean_mag.iloc[min(1735,len(elig)-1)])
print("  eligible by G bin:\n", pd.cut(elig.phot_g_mean_mag,[11,14,15,16,17,17.5,18,18.5,19]).value_counts().sort_index().to_string())
print("  eligible distance (plx/err>5): <100 pc:", ((elig.parallax>10)&(elig.parallax/elig.parallax_error>5)).sum(), " <200 pc:", ((elig.parallax>5)&(elig.parallax/elig.parallax_error>5)).sum(), " <300 pc:", ((elig.parallax>3.333)&(elig.parallax/elig.parallax_error>5)).sum())
# M-star remainder accounting
tot=len(ms); print("\nM-star catalogue remainder accounting (29,800):")
print("  psffrac<=0.5:", (ms.psffrac<=0.5).sum())
print("  psffrac>0.5 & crowdsap<=0.7:", ((ms.psffrac>0.5)&(ms.crowdsap<=0.7)).sum())
print("  psffrac>0.5 & crowdsap>0.7 & nq<8:", ((ms.psffrac>0.5)&(ms.crowdsap>0.7)&(ms.nq<8)).sum())
print("  ... & nq>=8 & Kepler target (kic==tpforg or in keplerstellar):", ((ms.psffrac>0.5)&(ms.crowdsap>0.7)&(ms.nq>=8)&((ms.kic==ms.tpforg)|ms.is_kepler_target)).sum())
print("  M stars within 100 pc (plx>10, plx/err>5) regardless of quality:", ((ms.parallax>10)&(ms.parallax/ms.parallax_error>5)).sum(), "; of those psffrac>0.5:", ((ms.parallax>10)&(ms.parallax/ms.parallax_error>5)&(ms.psffrac>0.5)).sum(), "; and nq>=8:", ((ms.parallax>10)&(ms.parallax/ms.parallax_error>5)&(ms.psffrac>0.5)&(ms.nq>=8)).sum(), "; not Kepler target:", ((ms.parallax>10)&(ms.parallax/ms.parallax_error>5)&(ms.psffrac>0.5)&(ms.nq>=8)&~ms.is_kepler_target).sum())
elig.to_csv("mstar_eligible_ranked.csv",index=False)
src[bkg&q][['gaia_designation','ra','dec','kic','tpforg','phot_g_mean_mag','bp_rp','psffrac','crowdsap','flfrcsap','nq','psf_avail','fname']].to_csv("background_quality_pool.csv",index=False)
print("\nwrote mstar_eligible_ranked.csv and background_quality_pool.csv")
