import requests, io, pandas as pd, sys, time
TAP="https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
def dl(sql, fn):
    t0=time.time()
    r=requests.get(TAP, params=dict(query=sql, format='csv'), timeout=900)
    if r.status_code!=200:
        print(fn, "HTTP", r.status_code, r.text[:800]); return None
    open(fn,'w').write(r.text)
    df=pd.read_csv(fn)
    print(f"{fn}: {len(df)} rows, {time.time()-t0:.0f}s, cols={list(df.columns)[:12]}"); sys.stdout.flush()
    return df
dl("select kepid,koi_disposition,koi_pdisposition,koi_score,koi_period,koi_prad,koi_depth,koi_duration,koi_kepmag,koi_steff,koi_srad,koi_smass,koi_model_snr,koi_max_mult_ev,koi_num_transits from cumulative", "dl/koi_cumulative.csv")
dl("select tid,toi,tfopwg_disp,pl_orbper,pl_rade,pl_trandep,pl_trandurh,st_tmag,st_teff,st_rad,st_dist from toi", "dl/toi.csv")
dl("select epic_hostname,pl_name,disposition,pl_orbper,pl_rade,st_teff,st_rad,sy_kepmag,sy_dist from k2pandc", "dl/k2pandc.csv")
dl("select kepid,kepmag,teff,logg,radius,mass,dist,rrmscdpp03p0,rrmscdpp06p0,rrmscdpp12p0,dutycycle,dataspan from q1_q17_dr25_stellar", "dl/dr25_stellar.csv")
