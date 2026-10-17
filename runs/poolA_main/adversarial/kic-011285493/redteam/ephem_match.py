import pandas as pd, numpy as np
K='/home/user/ZyushgPlanet/data/known/'
ra0,dec0=285.1259461685635,49.01660406817707
P=10.790158; t0=843.73888  # BKJD
def sep(ra,dec): return 3600*np.degrees(np.arccos(np.clip(np.sin(np.radians(dec0))*np.sin(np.radians(dec))+np.cos(np.radians(dec0))*np.cos(np.radians(dec))*np.cos(np.radians(ra-ra0)),-1,1)))
koi=pd.read_csv(K+'koi_cumulative.csv'); tce=pd.read_csv(K+'dr25_tce.csv'); eb=pd.read_csv(K+'J_AJ_151_68_catalog.csv'); ps=pd.read_csv(K+'ps_keplerfield.csv')
print('EB cols',list(eb.columns)[:15])
for name,df,pc,tc,idc in [('KOI',koi,'koi_period','koi_time0bk','kepoi_name'),('TCE',tce,'tce_period','tce_time0bk','kepid'),('PS',ps,'pl_orbper','pl_tranmid','pl_name')]:
    df=df.copy(); df['sep']=sep(df['ra'].values,df['dec'].values)
    near=df[df['sep']<180]
    print(f'--- {name}: {len(near)} within 3 arcmin'); print(near[[idc,pc,'sep']].to_string() if len(near) else '')
    for kic in [11285493,11285496]:
        if 'kepid' in df: print(name,'rows for KIC',kic,':',(df['kepid']==kic).sum())
    hits=[]
    for mult in [1,2,3,4,0.5,1/3,0.25,1.5,2/3]:
        Pt=P*mult
        s=np.abs(df[pc]-Pt)/Pt<2e-3
        for _,row in df[s].iterrows():
            tt=row[tc]; 
            if name=='PS': tt=tt-2454833
            Pm=row[pc]
            # phase offset of our t0 relative to their ephemeris in units of their period
            dphi=((t0-tt)/Pm+0.5)%1-0.5
            hits.append((row[idc],round(Pm,5),round((Pm-Pt)/Pt*1e3,3),mult,round(dphi,3),round(abs(dphi)*Pm*24,1),round(row['sep']/60,1)))
    hits.sort(key=lambda x:abs(x[5]))
    print(f'{name} period matches |dP/P|<2e-3 (id, P, dP/P 1e-3, mult, dphase, |dt| h, sep arcmin):')
    for h in hits[:25]: print('  ',h)
    print('  total',len(hits))
# EB catalogue
pcol=[c for c in eb.columns if c.lower() in ('period','per','p')][0]
print('EB period col',pcol)
rc=[c for c in eb.columns if c.lower().startswith('ra')][0]; dc=[c for c in eb.columns if c.lower().startswith('de')][0]
eb['sep']=sep(eb[rc].values,eb[dc].values)
print('EBs within 3 arcmin:',eb[eb['sep']<180].to_string())
for mult in [1,2,3,0.5,1/3,1.5]:
    s=np.abs(eb[pcol]-P*mult)/(P*mult)<2e-3
    for _,row in eb[s].iterrows(): print('EB match',mult,row[[c for c in eb.columns if 'KIC' in c.upper()][0]],row[pcol],'sep arcmin',round(row['sep']/60,1), {c:row[c] for c in eb.columns if 'bjd' in c.lower() or 't0' in c.lower()})
