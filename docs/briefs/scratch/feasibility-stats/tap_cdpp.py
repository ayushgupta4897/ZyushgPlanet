import requests, io, pandas as pd, sys, time
TAP="https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
def q(sql, label):
    t0=time.time()
    r=requests.get(TAP, params=dict(query=sql, format='csv'), timeout=600)
    if r.status_code!=200: print(label, r.status_code, r.text[:500]); return None
    df=pd.read_csv(io.StringIO(r.text))
    print(f"--- {label} ({time.time()-t0:.1f}s)"); print(df.to_string(index=False)); sys.stdout.flush()
    return df
# DR25 stellar table: median 6-hr CDPP by Kp bin for cool dwarfs and for all dwarfs
for teffcut,label in [("teff<4000","cool"),("teff between 4000 and 7000","FGK")]:
    for lo in range(10,18):
        q(f"""select {lo} as kp_lo, count(*) as n, avg(rrmscdpp06p0) as mean_cdpp6, avg(rrmscdpp03p0) as mean_cdpp3, avg(rrmscdpp12p0) as mean_cdpp12
        from q1_q17_dr25_stellar where {teffcut} and logg>4 and kepmag>={lo} and kepmag<{lo+1} and rrmscdpp06p0 is not null and rrmscdpp06p0>0""", f"DR25 {label} dwarfs Kp {lo}-{lo+1}")
