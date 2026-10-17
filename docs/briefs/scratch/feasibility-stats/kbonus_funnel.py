import pandas as pd, numpy as np
src=pd.read_parquet('dl/kbonus_source_catalog.parquet')
ms=pd.read_csv('dl/kbonus_mstars.csv'); ms['gaia_designation']=ms['gaia_designation'].astype(str).str.replace(r"^b'",'',regex=True).str.replace(r"'$",'',regex=True).str.strip()
print("source catalog", len(src), "mstars", len(ms))
m=ms.merge(src, on='gaia_designation', suffixes=('_m',''))
print("joined M stars:", len(m))
m['n_psf_q']=m['psf_avail'].astype(str).str.count('1')
m['background']=(m['kic']!=m['tpforg'])
steps=[("M-star catalogue", np.ones(len(m),bool))]
c1=m['psffrac']>0.5; steps.append(("PSFFRAC>0.5", c1))
c2=c1&(m['crowdsap']>0.7); steps.append(("CROWDSAP>0.7", c2))
c3=c2&(m['n_psf_q']>=8); steps.append((">=8 PSF quarters", c3))
c4=c3&m['background']; steps.append(("kic != tpforg (background)", c4))
for lab,c in steps: print(f"{lab:32s} {c.sum():6d}")
elig=m[c4].copy()
print("eligible G mag quantiles:", elig['phot_g_mean_mag'].quantile([0,0.05,0.25,0.5,0.75,0.95,1]).round(2).to_dict())
bins=[10,13,14,15,16,16.5,17,18,19.5]
print("\nCounts by G bin: all M-star cat | background (kic!=tpforg) | eligible (full funnel)")
tab=pd.DataFrame({'all':pd.cut(m['phot_g_mean_mag'],bins).value_counts().sort_index(),
                  'background':pd.cut(m.loc[m.background,'phot_g_mean_mag'],bins).value_counts().sort_index(),
                  'eligible':pd.cut(elig['phot_g_mean_mag'],bins).value_counts().sort_index()})
print(tab)
print("\nKepler-target (kic==tpforg) M stars in catalogue:", (~m.background).sum(), "; with psffrac>0.5&crowdsap>0.7&>=8q:", (c3&~m.background).sum())
# rank by G and save eligible list with ranks (approximate reproduction of Kuber's ranking; he also excluded 36 primary targets)
elig=elig.sort_values(['phot_g_mean_mag','gaia_designation']).reset_index(drop=True); elig['rank']=elig.index+1
elig.to_csv('kbonus_eligible_background_mdwarfs.csv',index=False)
print("rank 1-5:", elig[['rank','gaia_designation','kic','phot_g_mean_mag','psffrac','crowdsap','n_psf_q']].head(5).to_string(index=False))
# where does KIC 5689121 rank?
k=elig[elig.kic==5689121]; print("KIC 5689121 in eligible list:", k[['rank','phot_g_mean_mag','psffrac','crowdsap','n_psf_q','kep_mag','tpforg']].to_string(index=False))
# noise-measurement sample: 2 eligible stars nearest each target G
targets=[13.5,14.5,15.5,16.5,17.5,18.5]
sel=[]
for g in targets:
    cand=elig.iloc[(elig.phot_g_mean_mag-g).abs().argsort()[:2]]
    sel.append(cand)
sel=pd.concat(sel); sel.to_csv('kbonus_noise_sample.csv',index=False)
print(sel[['rank','gaia_designation','kic','phot_g_mean_mag','psffrac','crowdsap','n_psf_q','fname']].to_string(index=False))
