import warnings; warnings.filterwarnings('ignore')
import lightkurve as lk, astropy.units as u, numpy as np
from astropy.coordinates import SkyCoord
print('lightkurve',lk.__version__)
def summ(sr,label):
    try:
        tb=sr.table
        secs=[int(x) for x in tb['sequence_number'] if not np.ma.is_masked(x) and x is not None]
        print(label,'n',len(sr),'sectors',sorted(set(secs)), 'authors',sorted(set(tb['author'])))
    except Exception as e: print(label,'ERR',e)
r=lk.search_lightcurve(SkyCoord(285.907827,40.941988,unit='deg'),author='KBONUS-BKG',radius=20*u.arcsec)
print('KBONUS',len(r)); print(r.table['target_name','author','mission','distance'][:5] if len(r) else 'none')
for tic,auths in [('TIC 120686465',['TARS','QLP','TESS-SPOC','SPOC','TGLC','GSFC-ELEANOR-LITE']),('TIC 142960957',None),('TIC 130644400',None),('TIC 150428135',['QLP','TESS-SPOC'])]:
    if auths is None:
        sr=lk.search_lightcurve(tic)
        tb=sr.table
        import collections
        d=collections.defaultdict(set)
        for a,s in zip(tb['author'],tb['sequence_number']):
            if not np.ma.is_masked(s): d[a].add(int(s))
        print(tic,{k:sorted(v) for k,v in d.items()})
    else:
        for a in auths:
            summ(lk.search_lightcurve(tic,author=a),f'{tic} {a}')
