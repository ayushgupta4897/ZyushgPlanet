import pandas as pd, numpy as np
src=pd.read_parquet("source_catalog.parquet"); src['nq']=src.psf_avail.str.count('1')
kep=set(pd.read_csv("keplerstellar_kepids.csv").kepid.astype(int))
pool=src[(~src.kic.isin(kep))&(src.psffrac>0.7)&(src.crowdsap>0.8)&(src.nq>=12)].copy()
rng=np.random.default_rng(42); rows=[]
for lo in [12,13,14,15,16,17,18]:
    b=pool[(pool.phot_g_mean_mag>=lo)&(pool.phot_g_mean_mag<lo+1)]
    rows.append(b.sample(min(4,len(b)),random_state=1))
s=pd.concat(rows)
def url(r):
    f=r.fname
    if f.startswith('kic-'):
        k=f[4:]; return f"https://archive.stsci.edu/hlsps/kbonus-bkg/lcs/{k[:4]}/{k}/hlsp_kbonus-bkg_kepler_kepler_{f}_kepler_v1.0_lc.fits"
    else:
        gid=f.replace('gaia-dr3-',''); return f"https://archive.stsci.edu/hlsps/kbonus-bkg/lcs/{gid[:4]}/{gid}/hlsp_kbonus-bkg_kepler_kepler_{f}_kepler_v1.0_lc.fits"
s['url']=s.apply(url,axis=1)
s[['gaia_designation','fname','phot_g_mean_mag','psffrac','crowdsap','nq','url']].to_csv("sample/sample_list.csv",index=False)
print(s[['fname','phot_g_mean_mag','psffrac','crowdsap','nq']].to_string(index=False))
