"""Back-of-envelope yield model for a solo archival transit search around M dwarfs.
All inputs are listed at the top with their sources; outputs are printed as markdown tables."""
import numpy as np, pandas as pd, json, os
pd.set_option('display.width',250)
G_MSUN_AU_YR = 1.0  # Kepler's 3rd law in AU, yr, Msun
RSUN_AU = 0.00465047
REARTH_RSUN = 1/109.2
# ---------------------------------------------------------------- 1. occurrence (DC15 Table 5, Huber-radii block; planets per 100 stars, cumulative in period)
# rows: radius bin; cols: P<10, P<50, P<100, P<150, P<200 d
DC15_T5 = {
 (1.0,1.5):[21.18,56.25,65.09,82.17,88.74],
 (1.5,2.0):[10.76,45.95,57.35,66.58,69.72],
 (2.0,2.5):[ 8.68,35.62,49.81,54.59,55.72],
 (2.5,3.0):[ 3.89,15.46,21.97,23.19,23.45],
 (3.0,3.5):[ 1.73, 4.65, 7.15, 7.38, 7.41],
 (3.5,4.0):[ 0.72, 1.63, 2.11, 2.27, 2.27]}
PEDGES=[0.5,10,50,100,150,200]
def occ_density(rlo,rhi):
    """piecewise-constant occurrence density per ln P (planets per star per ln P) from the cumulative table"""
    c=np.array(DC15_T5[(rlo,rhi)])/100.
    d=np.diff(np.concatenate([[0],c]))               # planets per star in each period interval
    width=np.diff(np.log(PEDGES))
    return d/width                                   # per ln P
# fine grid in ln P
lnP=np.linspace(np.log(0.5),np.log(200),400); P=np.exp(lnP); dlnP=lnP[1]-lnP[0]
def dens_on_grid(rlo,rhi):
    dens=occ_density(rlo,rhi); idx=np.clip(np.searchsorted(PEDGES,P,side='right')-1,0,len(dens)-1)
    return dens[idx]
# ---------------------------------------------------------------- 2. star and geometry
MSTAR, RSTAR = 0.5, 0.5
def a_AU(P_d): return (MSTAR*(P_d/365.25)**2)**(1/3)
def p_transit(P_d): return 1.08*RSTAR*RSUN_AU/a_AU(P_d)   # 1.08 = eccentricity boost (Kipping 2013) as in DC15
def t14_hours(P_d):  # mean over impact parameter (pi/4 of central) of the full duration for Rp<<R*
    return (P_d*24/np.pi)*(RSTAR*RSUN_AU/a_AU(P_d))*(np.pi/4)
# ---------------------------------------------------------------- 3. noise models (ppm on the transit-duration timescale)
# Kepler official: mean 3/6/12-h rrmsCDPP of DR25 cool dwarfs (Teff<4000, logg>4) by Kp bin, Exoplanet Archive Q1_Q17_DR25_KS (queried 2026-10-17)
_c=pd.read_csv('dl/dr25_cool_dwarfs.csv'); _c=_c[_c.rrmscdpp06p0>0]; _c['bin']=pd.cut(_c.kepmag,[10,12,13,14,15,16,17,20])
_g=_c.groupby('bin',observed=True); cdpp6_kep=_g.rrmscdpp06p0.median().values; cdpp3_kep=_g.rrmscdpp03p0.median().values; cdpp12_kep=_g.rrmscdpp12p0.median().values
kp_centres=np.array([11.5,12.5,13.5,14.5,15.5,16.5,17.5])
print('DR25 cool-dwarf median CDPP by Kp bin (3h/6h/12h ppm, n):', list(zip(kp_centres, cdpp3_kep.round(0), cdpp6_kep.round(0), cdpp12_kep.round(0), _g.size().values)))
alpha_kep=np.log(cdpp3_kep/cdpp6_kep).mean()/np.log(2)   # CDPP(T) ~ T^-alpha
def cdpp_kepler(kp, T_h): return np.interp(kp, kp_centres, cdpp6_kep)*(6.0/T_h)**alpha_kep
# KBONUS PSF photometry: measured in this sandbox (kbonus_noise_vs_G.csv) on eligible background M dwarfs; fallback = Kepler x 1.5
kb=pd.read_csv('kbonus_noise_vs_G.csv')
_c2=pd.read_csv('dl/dr25_cool_dwarfs.csv'); _c2=_c2[_c2.rrmscdpp06p0>0]
_ratios=[]
for _,r in kb.iterrows():
    mm=_c2[(_c2.kepmag>r.G-0.25)&(_c2.kepmag<r.G+0.25)]
    if len(mm)>5: _ratios.append(r['flux_cdpp6h']/mm.rrmscdpp06p0.median())
KB_RATIO=float(np.median(_ratios))
KB_RATIO_PESSIMISTIC=1.5
def cdpp_kbonus(G, T_h): return KB_RATIO*cdpp_kepler(G, T_h)
KB_SOURCE=f"KBONUS PSF 6-h noise / Kepler DR25 cool-dwarf median CDPP at same magnitude: median ratio {KB_RATIO:.2f} over {len(_ratios)} measured background M dwarfs (range {min(_ratios):.2f}-{max(_ratios):.2f}); model = {KB_RATIO:.2f} x Kepler CDPP(G); pessimistic alt = 1.5x"
# TESS: Stassun et al. 2018 (arXiv:1706.00495) 1-hour noise polynomial as implemented in ticgen 1.0.4; scaled as white noise to duration
def tess_1h(T):
    lnA=3.29685004771; B=0.850021465753; C=-0.285041632435; D=0.0395908321369; E=-0.0022308015894; F=4.73508403525e-5
    return np.exp(lnA+B*T+C*T**2+D*T**3+E*T**4+F*T**5)
def cdpp_tess(T, T_h): return tess_1h(T)/np.sqrt(T_h)
# K2: EVEREST 2.0 recovers Kepler precision to Kp~14-15 (Luger+2018); assume Kepler CDPP x1.0 to Kp 14.5, degrading linearly to x2.0 at Kp 17 (ASSUMPTION)
def cdpp_k2(kp, T_h): return cdpp_kepler(kp,T_h)*np.interp(kp,[14.5,17.0],[1.0,2.0])
# ---------------------------------------------------------------- 4. detectability
SNR_THRESH=7.1; MIN_TRANSITS=3
import sys
REALISTIC = '--realistic' in sys.argv   # logistic completeness vs step
def pdet(snr):
    if REALISTIC: return 1/(1+np.exp(-(snr-9.0)/1.2))   # ~0.17 at 7.1, 0.5 at 9, 0.92 at 12 (shape mimics Christiansen+2020 DR25 gamma-CDF)
    return (snr>=SNR_THRESH).astype(float)
def detect_frac(noise_fn, mag, baseline_d, duty, rlo, rhi, pmask, dilution=1.0, nR=15):
    """occurrence-weighted fraction of transiting planets in [rlo,rhi] x pmask with SNR>=7.1 and >=3 transits (radius uniform in ln R within each DC15 sub-bin)"""
    Pm=P[pmask]; T_h=t14_hours(Pm); ntr=duty*baseline_d/Pm; noise=noise_fn(mag, T_h)
    tot=0; det=0
    for (a,b) in DC15_T5:
        if a<rlo or b>rhi: continue
        dens=dens_on_grid(a,b)[pmask]
        for r in np.exp(np.linspace(np.log(a),np.log(b),nR)):
            depth=(r*REARTH_RSUN/RSTAR)**2*1e6*dilution
            snr=depth/noise*np.sqrt(ntr)
            ok=pdet(snr)*(ntr>=MIN_TRANSITS)
            w=dens*p_transit(Pm)*dlnP
            tot+=w.sum(); det+=(w*ok).sum()
    return det/tot if tot>0 else 0.0
def expected_transiting(rlo,rhi,pmask):
    dens=dens_on_grid(rlo,rhi)[pmask]; Pm=P[pmask]
    return (dens*p_transit(Pm)*dlnP).sum(), (dens*dlnP).sum()
RBINS=[(1.0,1.5),(1.5,2.0),(2.0,4.0)]
PBINS={'<10 d':(P<10),'10-50 d':(P>=10)&(P<50),'50-200 d':(P>=50)}
def occ_bin(rlo,rhi,pmask):
    tot=0
    for (a,b) in DC15_T5:
        if a>=rlo and b<=rhi: tot+=(dens_on_grid(a,b)[pmask]*dlnP).sum()
    return tot
def ptr_bin(rlo,rhi,pmask):
    tot=0
    for (a,b) in DC15_T5:
        if a>=rlo and b<=rhi: tot+=(dens_on_grid(a,b)[pmask]*p_transit(P[pmask])*dlnP).sum()
    return tot
print("## 1. Occurrence (DC15 Table 5, Huber-radii block), transit probability, transiting planets per 1000 M dwarfs (M*=0.5, R*=0.5)\n")
rows=[]
for (rlo,rhi) in RBINS:
    for pl,pm in PBINS.items():
        f=occ_bin(rlo,rhi,pm); ft=ptr_bin(rlo,rhi,pm)
        rows.append(dict(Rp=f"{rlo}-{rhi} Re", P=pl, planets_per_star=round(f,3), transiting_per_1000_stars=round(1000*ft,2), mean_ptransit=round(ft/f,4) if f>0 else np.nan))
tab1=pd.DataFrame(rows); print(tab1.to_markdown(index=False))
f_all=sum(occ_bin(a,b,P>0) for a,b in RBINS); ft_all=sum(ptr_bin(a,b,P>0) for a,b in RBINS)
print(f"\nTotal 1-4 Re, 0.5-200 d: {f_all:.2f} planets/star; transiting: {1000*ft_all:.1f} per 1000 stars")
m50=P<50; f50=sum(occ_bin(a,b,m50) for a,b in [(2.0,4.0)]); ft50=ptr_bin(2.0,4.0,m50)
print(f"P<50 d, Rp 2-4 Re: {f50:.3f} planets/star; transiting {1000*ft50:.1f} per 1000 M dwarfs")
f50b=sum(occ_bin(a,b,m50) for a,b in RBINS); ft50b=sum(ptr_bin(a,b,m50) for a,b in RBINS)
print(f"P<50 d, Rp 1-4 Re: {f50b:.2f} planets/star; transiting {1000*ft50b:.1f} per 1000 M dwarfs")
print("\nGeometric transit probability R*/a (x1.08 ecc. boost) for M*=0.5 Msun, R*=0.5 Rsun:")
print(pd.DataFrame({'P_d':[1,3,5,10,13.77,30,50,100,200],'a_AU':[round(a_AU(p),4) for p in [1,3,5,10,13.77,30,50,100,200]],'p_transit':[round(p_transit(p),4) for p in [1,3,5,10,13.77,30,50,100,200]],'T14_mean_h':[round(t14_hours(p),2) for p in [1,3,5,10,13.77,30,50,100,200]]}).to_markdown(index=False))
# ---------------------------------------------------------------- 5. noise tables and minimum detectable radius
print("\n## 2. Noise vs magnitude (ppm) and minimum detectable radius (SNR>=7.1) for a 0.5 Rsun / 0.5 Msun M dwarf\n")
print(f"Kepler CDPP duration scaling exponent alpha (from DR25 cool-dwarf 3h/6h means): {alpha_kep:.2f}; KBONUS noise model: {KB_SOURCE}")
rows=[]
for mag in [12,13,14,15,16,17,18]:
    r=dict(mag=mag)
    r['Kepler_CDPP6h']=round(cdpp_kepler(mag,6)); r['KBONUS_CDPP6h']=round(cdpp_kbonus(mag,6)); r['K2_CDPP6h']=round(cdpp_k2(mag,6)); r['TESS_1h(Tmag=mag)']=round(tess_1h(mag))
    rows.append(r)
print(pd.DataFrame(rows).to_markdown(index=False))
def rmin(noise_fn, mag, baseline, duty, P_d, dilution=1.0):
    T_h=t14_hours(P_d); ntr=duty*baseline/P_d
    if ntr<MIN_TRANSITS: return np.nan
    depth_min=SNR_THRESH*noise_fn(mag,T_h)/np.sqrt(ntr)
    return RSTAR/REARTH_RSUN*np.sqrt(depth_min/1e6/dilution)
print("\nMinimum detectable planet radius (Re) at P = 5 / 13.77 / 30 d:")
rows=[]
for mag in [12,13,14,15,16,17,18]:
    r=dict(mag=mag)
    for lab,fn,base,duty in [('Kepler 4yr',cdpp_kepler,1460,0.92),('KBONUS 13q',cdpp_kbonus,1170,0.9),('K2 80d',cdpp_k2,80,0.95),('TESS 2 sec',cdpp_tess,54.8,0.9),('TESS 4 sec',cdpp_tess,109.6,0.9),('TESS 13 sec',cdpp_tess,356,0.9)]:
        for pd_ in [5,13.77,30]:
            v=rmin(fn,mag,base,duty,pd_); r[f'{lab} P{pd_}']=round(v,2) if np.isfinite(v) else '-'
    rows.append(r)
tabR=pd.DataFrame(rows); print(tabR.to_markdown(index=False))
tabR.to_csv('min_radius_table.csv',index=False)
# ---------------------------------------------------------------- 6. stars per dataset and magnitude bin
print("\n## 3. Yield\n")
# Kepler primary cool dwarfs (DR25, Teff<4000, logg>4) per Kp bin
cool=pd.read_csv('dl/dr25_cool_dwarfs.csv'); cool['bin']=pd.cut(cool.kepmag,[9,12,13,14,15,16,17,20],labels=['<12','12-13','13-14','14-15','15-16','16-17','>17'])
kep_counts=cool.groupby('bin',observed=True).size()
# KBONUS background M dwarfs (eligible funnel) per G bin
elig=pd.read_csv('kbonus_eligible_background_mdwarfs.csv'); elig['bin']=pd.cut(elig.phot_g_mean_mag,[9,12,13,14,15,16,17,20],labels=['<12','12-13','13-14','14-15','15-16','16-17','>17'])
kb_counts=elig.groupby('bin',observed=True).size()
# relaxed KBONUS: PSFFRAC>0.5 & CROWDSAP>0.7 & background, any number of quarters
src=pd.read_parquet('dl/kbonus_source_catalog.parquet'); ms=pd.read_csv('dl/kbonus_mstars.csv'); ms['gaia_designation']=ms.gaia_designation.astype(str).str.replace(r"^b'",'',regex=True).str.replace(r"'$",'',regex=True)
m=ms[['gaia_designation']].merge(src,on='gaia_designation'); m['n_psf_q']=m.psf_avail.astype(str).str.count('1'); relaxed=m[(m.psffrac>0.5)&(m.crowdsap>0.7)&(m.kic!=m.tpforg)&(m.n_psf_q>=4)]
relaxed_counts=relaxed.assign(bin=pd.cut(relaxed.phot_g_mean_mag,[9,12,13,14,15,16,17,20],labels=['<12','12-13','13-14','14-15','15-16','16-17','>17'])).groupby('bin',observed=True).size()
# GCNS M dwarfs within 100 pc by approx Tmag bin (T ~ G - 0.7 for BP-RP~2.3): counts queried 2026-10-17 (gaia_gcns3.py)
gcns_cum={'<12':11610,'<13':36700,'<14':74331,'<15':118970}
tess_counts={'<12':11610,'12-13':36700-11610,'13-14':74331-36700,'14-15':118970-74331}
known_koi={'<12':0,'12-13':0,'13-14':8,'14-15':24,'15-16':104,'16-17':19,'>17':2}  # cool-host KOIs (cand+conf) by Kp, all P/R (known_bins.py)
toi=pd.read_csv('dl/toi.csv'); t=toi[(toi.st_teff<4000)&~toi.tfopwg_disp.isin(['FP','FA'])&(toi.st_dist<100)&(toi.pl_orbper<50)&(toi.pl_rade>=1)&(toi.pl_rade<=4)]
t['bin']=pd.cut(t.st_tmag,[0,12,13,14,15,20],labels=['<12','12-13','13-14','14-15','>15']); toi_counts=t.groupby('bin',observed=True).size()
koi=pd.read_csv('dl/koi_cumulative.csv'); k=koi[(koi.koi_steff<4000)&koi.koi_disposition.isin(['CANDIDATE','CONFIRMED'])&(koi.koi_period<50)&(koi.koi_prad>=1)&(koi.koi_prad<=4)]
k['bin']=pd.cut(k.koi_kepmag,[9,12,13,14,15,16,17,20],labels=['<12','12-13','13-14','14-15','15-16','16-17','>17']); koi_counts=k.groupby('bin',observed=True).size()
m50=P<50
def yield_table(label, counts, noise_fn, baseline, duty, magmap, known=None, dilution=1.0):
    rows=[]; tot_exp=0; tot_known=0
    for b,n in counts.items():
        mag=magmap[b]
        exp=0; det_small=0
        for (rlo,rhi) in RBINS:
            ft=ptr_bin(rlo,rhi,m50)  # transiting per star
            C=detect_frac(noise_fn,mag,baseline,duty,rlo,rhi,m50,dilution)
            exp+=n*ft*C
        exp2=0
        for (rlo,rhi) in [(2.0,4.0)]:
            exp2+=n*ptr_bin(rlo,rhi,m50)*detect_frac(noise_fn,mag,baseline,duty,rlo,rhi,m50,dilution)
        kn=int(known.get(b,0)) if known is not None else 0
        rows.append(dict(mag_bin=b,n_stars=int(n),transiting_1_4Re_P50=round(n*sum(ptr_bin(a,c,m50) for a,c in RBINS),1),detectable_1_4Re=round(exp,1),detectable_2_4Re=round(exp2,1),known_1_4Re=kn,remaining_est=round(max(exp-kn,0),1)))
        tot_exp+=exp; tot_known+=kn
    df=pd.DataFrame(rows); print(f"\n### {label}\n"); print(df.to_markdown(index=False)); print(f"Totals: detectable {tot_exp:.1f}, known {tot_known}, remaining ~{max(tot_exp-tot_known,0):.1f}")
    return df
kepmap={'<12':11.5,'12-13':12.5,'13-14':13.5,'14-15':14.5,'15-16':15.5,'16-17':16.5,'>17':17.5}
yield_table("Kepler primary cool dwarfs (DR25 Teff<4000, logg>4; official PDC photometry; 4-yr; known = cool-host KOIs cand+conf with P<50 d, 1-4 Re)", kep_counts, cdpp_kepler, 1460, 0.92, kepmap, koi_counts.to_dict())
yield_table("KBONUS eligible background M dwarfs (Kuber funnel, 1,836 stars; PSF photometry; ~13 quarters; known = 0 catalogued, 1 Kuber candidate)", kb_counts, cdpp_kbonus, 1170, 0.9, kepmap, None)
yield_table("KBONUS relaxed background M dwarfs (PSFFRAC>0.5, CROWDSAP>0.7, background, >=4 PSF quarters; ~9 quarters assumed)", relaxed_counts, cdpp_kbonus, 810, 0.9, kepmap, None)
tessmap={'<12':11.5,'12-13':12.5,'13-14':13.5,'14-15':14.5}
for nsec in [2,4,8]:
    yield_table(f"TESS FFI: GCNS M dwarfs within 100 pc by approx Tmag (T~G-0.7), {nsec} sectors assumed, dilution 0.85; known = cool-host TOIs d<100 pc, P<50 d, 1-4 Re (excl FP/FA)", pd.Series(tess_counts), cdpp_tess, 27.4*nsec, 0.9, tessmap, toi_counts.to_dict(), dilution=0.85)
# K2: approximate M-dwarf count via colour cuts is uncertain; use the literature figure of ~2,500-4,000 M dwarfs per Hardegree-Ullman 2020 style samples -> run with 3,000 at Kp 14.5 and 15.5 split
k2_counts=pd.Series({'<14':2316,'14-15':3843,'15-16':4805,'16-16.5':2927})
k2known={'<14':40,'14-15':40,'15-16':45,'16-16.5':7}
yield_table("K2: M-dwarf proxy from K2TARGETS (J-K 0.78-1.05, Kp-J>2.3, PM>20 mas/yr, long cadence; 13,891 stars Kp<16.5; 80-d campaigns; EVEREST-like noise; known = 132 K2 cool-host planets P<50 d 1-4 Re spread by Kp)", k2_counts, cdpp_k2, 80, 0.95, {'<14':13.5,'14-15':14.5,'15-16':15.5,'16-16.5':16.25}, k2known)
print("\nMODE:", "REALISTIC logistic completeness (50% at SNR 9)" if REALISTIC else "STEP completeness at SNR 7.1 (upper envelope)")
json.dump(dict(alpha_kep=alpha_kep, kb_source=KB_SOURCE), open('yield_model_meta.json','w'))
