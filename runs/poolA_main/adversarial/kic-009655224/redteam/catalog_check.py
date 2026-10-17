import pandas as pd, numpy as np, json, sys
K='/home/user/ZyushgPlanet/data/known/'
ra0,dec0=294.2826791587251,46.34150883424167
P=6.703464; t0=265.18316
out={}
def sep(ra,dec): return np.hypot((ra-ra0)*np.cos(np.radians(dec0)),dec-dec0)*3600
koi=pd.read_csv(K+'koi_cumulative.csv'); tce=pd.read_csv(K+'dr25_tce.csv'); ps=pd.read_csv(K+'ps_keplerfield.csv')
eb=pd.read_csv(K+'J_AJ_151_68_catalog.csv'); ks=pd.read_csv(K+'keplerstellar.csv')
print('EB cols', list(eb.columns)[:20])
# 1. direct KIC matches
for name,df,col in [('koi',koi,'kepid'),('tce',tce,'kepid'),('ebs',eb,[c for c in eb.columns if 'KIC' in c.upper()][0]),('keplerstellar',ks,'kepid')]:
    m=df[df[col].isin([9655224,9655236,9655226])]
    print(name,'direct KIC rows:',len(m)); print(m.to_string()[:1500]) if len(m) else None
# 2. within 3 arcmin
for name,df in [('koi',koi),('tce',tce),('ps',ps)]:
    rc=[c for c in df.columns if c.lower() in ('ra','ra_obj')][0]; dc=[c for c in df.columns if c.lower() in ('dec','dec_obj')][0]
    s=sep(df[rc].values,df[dc].values); m=df[s<180]; print(name,'within 3 arcmin:',len(m)); print(m.assign(sep=s[s<180]).to_string()[:2000]) if len(m) else None
# EB catalogue coordinates?
rc=[c for c in eb.columns if c.lower().startswith('ra')]; dc=[c for c in eb.columns if c.lower().startswith('de')]
print('eb coord cols',rc,dc)
if rc and dc:
    s=sep(eb[rc[0]].values,eb[dc[0]].values); m=eb[s<180]; print('ebs within 3 arcmin:',len(m)); print(m.to_string()[:1500]) if len(m) else None
# 3. period matches anywhere in the field (harmonics), with epoch
def epoch_match(Pc,T0c,n_h):
    # phase of our t0 in their ephemeris, in hours
    ph=((t0-T0c)/Pc+0.5)%1-0.5; return ph*Pc*24
print('\n=== period matches (|dP/P|<2e-3) at P*{1/3,1/2,1,2,3} ===')
for name,df,pc,tc,dcol,ic in [('koi',koi,'koi_period','koi_time0bk','koi_depth','kepoi_name'),('tce',tce,'tce_period','tce_time0bk','tce_depth','kepid')]:
    rc=[c for c in df.columns if c.lower()=='ra'][0]; dcc=[c for c in df.columns if c.lower()=='dec'][0]
    for h in [1/3,0.5,1,2,3]:
        Pc=P*h; m=df[np.abs(df[pc]/Pc-1)<2e-3].copy()
        if len(m):
            m['sep_arcmin']=sep(m[rc].values,m[dcc].values)/60; m['dt_hours']=[epoch_match(p,t,h) for p,t in zip(m[pc],m[tc])]
            print(f'{name} harmonic {h}: {len(m)} rows'); print(m[[ic,pc,tc,dcol,'sep_arcmin','dt_hours']].sort_values('sep_arcmin').head(15).to_string())
pcol=[c for c in eb.columns if c.lower() in ('period','per','p')]; print('eb period col',pcol)
if pcol:
    for h in [1/3,0.5,1,2,3]:
        m=eb[np.abs(eb[pcol[0]]/(P*h)-1)<2e-3]
        if len(m): print('EB harmonic',h,len(m)); print(m.to_string()[:1500])
