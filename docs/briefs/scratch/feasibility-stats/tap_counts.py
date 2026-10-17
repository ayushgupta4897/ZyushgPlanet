import requests, io, pandas as pd, sys, time
TAP="https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
def q(sql, label):
    t0=time.time()
    r=requests.get(TAP, params=dict(query=sql, format='csv'), timeout=300)
    r.raise_for_status()
    df=pd.read_csv(io.StringIO(r.text))
    print(f"--- {label} ({time.time()-t0:.1f}s)"); print(df.to_string(index=False)); sys.stdout.flush()
    return df
# KOI cumulative: candidates+confirmed around cool hosts (Teff<4000), by radius & period bins, by Kp bin
q("""select koi_disposition, count(*) as n from cumulative where koi_steff<4000 and koi_disposition in ('CANDIDATE','CONFIRMED') group by koi_disposition""","KOI cumulative Teff<4000 by disposition")
q("""select case when koi_prad<1.5 then 'a_1.0-1.5' when koi_prad<2 then 'b_1.5-2' when koi_prad<4 then 'c_2-4' else 'd_>4' end as rbin,
 case when koi_period<10 then 'P<10' when koi_period<50 then 'P10-50' when koi_period<200 then 'P50-200' else 'P>200' end as pbin, count(*) as n
 from cumulative where koi_steff<4000 and koi_disposition in ('CANDIDATE','CONFIRMED') and koi_prad>=1 group by rbin,pbin order by rbin,pbin""","KOI cool hosts by Rp and P bin")
q("""select case when koi_kepmag<14 then 'Kp<14' when koi_kepmag<15 then 'Kp14-15' when koi_kepmag<16 then 'Kp15-16' when koi_kepmag<17 then 'Kp16-17' else 'Kp>17' end as kbin, count(*) as n
 from cumulative where koi_steff<4000 and koi_disposition in ('CANDIDATE','CONFIRMED') group by kbin order by kbin""","KOI cool hosts by Kp bin")
q("""select count(*) as n_kepler_targets_teff_lt_4000 from keplerstellar where teff<4000""","keplerstellar Teff<4000 count")
q("""select case when kepmag<14 then 'Kp<14' when kepmag<15 then 'Kp14-15' when kepmag<16 then 'Kp15-16' when kepmag<17 then 'Kp16-17' else 'Kp>17' end as kbin, count(*) as n
 from keplerstellar where teff<4000 group by kbin order by kbin""","keplerstellar Teff<4000 by Kp")
# TOI table: cool hosts
q("""select tfopwg_disp, count(*) as n from toi where st_teff<4000 group by tfopwg_disp""","TOI st_teff<4000 by disposition")
q("""select case when pl_rade<1.5 then 'a_1.0-1.5' when pl_rade<2 then 'b_1.5-2' when pl_rade<4 then 'c_2-4' else 'd_>4' end as rbin,
 case when pl_orbper<10 then 'P<10' when pl_orbper<50 then 'P10-50' else 'P>50' end as pbin, count(*) as n
 from toi where st_teff<4000 and pl_rade>=1 and tfopwg_disp not in ('FP','FA') group by rbin,pbin order by rbin,pbin""","TOI cool hosts by Rp,P (excluding FP/FA)")
q("""select case when st_tmag<10 then 'T<10' when st_tmag<12 then 'T10-12' when st_tmag<13 then 'T12-13' when st_tmag<14 then 'T13-14' when st_tmag<15 then 'T14-15' else 'T>15' end as tbin, count(*) as n
 from toi where st_teff<4000 and tfopwg_disp not in ('FP','FA') group by tbin order by tbin""","TOI cool hosts by Tmag")
q("""select case when st_dist<30 then 'd<30' when st_dist<60 then 'd30-60' when st_dist<100 then 'd60-100' else 'd>100' end as dbin, count(*) as n
 from toi where st_teff<4000 and tfopwg_disp not in ('FP','FA') group by dbin order by dbin""","TOI cool hosts by distance")
# K2 candidates cool hosts
q("""select disposition, count(*) as n from k2pandc where st_teff<4000 group by disposition""","K2 pandc Teff<4000 by disposition")
# confirmed planets around M dwarfs discovered by transit, by facility
q("""select disc_facility, count(*) as n from pscomppars where st_teff<4000 and discoverymethod='Transit' group by disc_facility order by n desc""","pscomppars M-dwarf transiting planets by facility")
