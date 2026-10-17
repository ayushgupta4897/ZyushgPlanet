import warnings; warnings.filterwarnings('ignore')
import lightkurve as lk, numpy as np, pandas as pd, collections, random, time
# Kuber's star: everything
sr=lk.search_lightcurve('TIC 120686465'); tb=sr.table
d=collections.defaultdict(set)
for a,s in zip(tb['author'],tb['sequence_number']):
    if not np.ma.is_masked(s): d[a].add(int(s))
print('TIC 120686465 all products:',{k:sorted(v) for k,v in d.items()})
# QLP latest sector availability
from astroquery.mast import Observations
for s in [104,105]:
    try:
        o=Observations.query_criteria(obs_collection='HLSP',provenance_name='QLP',sequence_number=s,target_name='TIC 150428135')
        print('QLP S',s,'TIC150428135 n',len(o))
    except Exception as e: print('QLP',s,'ERR',e)
# random sample of faint nearby M dwarfs: which authors have products?
df=pd.read_csv('/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/other-datasets/tic_mdwarfs_100pc_faint.csv')
print('100pc file rows',len(df),'cols',list(df.columns)[:12])
tcol=[c for c in df.columns if c.lower()=='tmag'][0]; dcol=[c for c in df.columns if c.lower()=='d'][0]; idcol=[c for c in df.columns if c.lower() in ('id','ticid','tic_id')][0]
print('T bins',np.histogram(df[tcol],bins=[13.5,14,14.5,15,15.5,16,16.5])[0],'median T',df[tcol].median(),'d<50',(df[dcol]<50).sum(),'d<50&T<16',((df[dcol]<50)&(df[tcol]<16)).sum(),'d<30&T<16',((df[dcol]<30)&(df[tcol]<16)).sum())
random.seed(1); ids=random.sample(list(df[idcol]),60)
cnt=collections.Counter(); nsec=[]; none=0; spoc2=0
t0=time.time()
for i in ids:
    try:
        sr=lk.search_lightcurve(f'TIC {i}'); tb=sr.table
        auths=set(tb['author']); 
        for a in auths: cnt[a]+=1
        if 'SPOC' in auths: spoc2+=1
        if len(tb)==0: none+=1
    except Exception as e: print('err',i,e)
print('sample of',len(ids),'faint nearby M dwarfs: stars with any product by author',dict(cnt),'no product at all',none,'with SPOC 2-min',spoc2,'elapsed',round(time.time()-t0))
