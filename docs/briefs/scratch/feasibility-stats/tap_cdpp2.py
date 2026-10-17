import requests, io, pandas as pd, sys, time
TAP="https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
def q(sql):
    r=requests.get(TAP, params=dict(query=sql, format='csv'), timeout=600)
    if r.status_code!=200: print("ERR", sql[:80], r.text[:300]); return None
    return pd.read_csv(io.StringIO(r.text))
# Try one aggregate query first
t=q("select count(*) as n, avg(rrmscdpp06p0) as cdpp6 from Q1_Q17_DR25_KS where teff<4000 and logg>4 and kepmag>=15 and kepmag<16 and rrmscdpp06p0>0")
print(t)
rows=[]
for lo in [10,11,12,13,14,15,16,17]:
    for teffcut,lab in [("teff<4000","cool"),("teff>=4000 and teff<7000","FGK")]:
        d=q(f"select count(*) as n, avg(rrmscdpp03p0) as cdpp3, avg(rrmscdpp06p0) as cdpp6, avg(rrmscdpp12p0) as cdpp12 from Q1_Q17_DR25_KS where {teffcut} and logg>4 and kepmag>={lo} and kepmag<{lo+1} and rrmscdpp06p0>0")
        if d is not None: rows.append(dict(kp_bin=f"{lo}-{lo+1}", sample=lab, **d.iloc[0].to_dict()))
df=pd.DataFrame(rows); pd.set_option('display.width',200); print(df.round(1).to_string(index=False)); df.to_csv('dr25_cdpp_by_kp.csv',index=False)
# Cool dwarf (M-dwarf proxy) target counts by Kp bin, all DR25 stellar
d=q("select count(*) as n from Q1_Q17_DR25_KS where teff<4000 and logg>4"); print("DR25 cool dwarfs total:", d.n[0])
d=q("select count(*) as n from Q1_Q17_DR25_KS where teff<4000 and logg>4 and kepmag<16.5"); print("DR25 cool dwarfs Kp<16.5:", d.n[0])
d=q("select count(*) as n from Q1_Q17_DR25_KS where teff<4000 and logg>4 and st_quarters is not null"); print("DR25 cool dwarfs with quarters info:", d.n[0])
# medians are not available in ADQL; download kepmag+cdpp for cool dwarfs to compute medians locally
d=q("select kepid,kepmag,teff,logg,radius,rrmscdpp03p0,rrmscdpp06p0,rrmscdpp12p0,st_quarters from Q1_Q17_DR25_KS where teff<4000 and logg>4")
d.to_csv('dl/dr25_cool_dwarfs.csv',index=False); print("downloaded cool dwarfs rows", len(d))
d['kpbin']=pd.cut(d.kepmag,[9,10,11,12,13,14,15,16,17,20])
g=d[d.rrmscdpp06p0>0].groupby('kpbin',observed=True)
print(pd.DataFrame({'n':g.size(),'median_cdpp3':g.rrmscdpp03p0.median().round(0),'median_cdpp6':g.rrmscdpp06p0.median().round(0),'median_cdpp12':g.rrmscdpp12p0.median().round(0),'p10_cdpp6':g.rrmscdpp06p0.quantile(.1).round(0),'median_radius':g.radius.median().round(2)}))
# K2 targets
c=q("select column_name from TAP_SCHEMA.columns where table_name='K2TARGETS'"); print("K2TARGETS cols:", c.column_name.tolist()[:60])
