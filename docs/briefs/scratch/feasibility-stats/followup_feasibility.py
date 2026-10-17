import numpy as np, pandas as pd
pd.set_option('display.width',250)
# --- Case A: Kuber-like Kepler-field candidate. Case B: Rabtsevich-like nearby TESS candidate.
cases={'A: Kp16 M1.5V @300pc (KIC 5689121-like)':dict(G=16.28, BP_RP=2.19, dist=298, Rstar=0.507, Mstar=0.505, Rp=2.7, P=13.77, Teff=3633),
       'B: Tmag12 mid-M @30pc (TIC 4206066-like)':dict(G=12.8, BP_RP=2.6, dist=30, Rstar=0.33, Mstar=0.32, Rp=1.4, P=3.18, Teff=3300)}
def mass_from_radius(rp):   # Chen & Kipping 2017 broken power law (approx): Mp ~ R^3.58 (R<1.23), ~ R^1.70 (1.23<R<14.3) with continuity
    if rp<1.23: return rp**3.58
    return 1.23**3.58*(rp/1.23)**1.70
rows=[]
for name,c in cases.items():
    mp=mass_from_radius(c['Rp'])                        # Earth masses
    a=(c['Mstar']*(c['P']/365.25)**2)**(1/3)
    K=28.4329*(mp/317.83)*(c['P']/365.25)**(-1/3)*c['Mstar']**(-2/3)  # m/s (Lovis & Fischer 2010 eq.)
    depth=(c['Rp']/109.2/c['Rstar'])**2*1e6
    T14=(c['P']*24/np.pi)*(c['Rstar']*0.00465047/a)*np.sqrt(1-0.25)  # b=0.5
    # approximate V, Ic, J from Gaia G and BP-RP (Riello+2021 / Pecaut & Mamajek colour tables; rough)
    V=c['G']+ (0.3 if c['BP_RP']<2.3 else 0.6)+0.3*(c['BP_RP']-2.0)
    Ic=c['G']-1.1-0.3*(c['BP_RP']-2.2)
    J=c['G']-2.6-0.5*(c['BP_RP']-2.2)
    # ground-based photon noise per transit for a 1-m and a 0.4-m telescope in Ic: zero point ~470 ph/s/cm2/A at Ic=0, 1500 A band, throughput 0.4
    def photon_noise_ppm(D_m, mag=Ic, hours=T14, sky_mag_arcsec2=19.5, ap_radius_arcsec=3.0, throughput=0.4):
        area=np.pi*(D_m*100/2)**2
        rate0=470*1500*area*throughput
        star=rate0*10**(-0.4*mag)*hours*3600
        sky=rate0*10**(-0.4*(sky_mag_arcsec2-2.5*np.log10(np.pi*ap_radius_arcsec**2)))*hours*3600
        return 1e6*np.sqrt(star+sky)/star*np.sqrt(2)  # x sqrt2 for in- vs out-of-transit comparison
    rows.append(dict(case=name, V=round(V,1), Ic=round(Ic,1), J=round(J,1), depth_ppm=round(depth), T14_h=round(T14,2), Mp_Me=round(mp,1), K_ms=round(K,2),
                     photon_noise_1m_per_transit_ppm=round(photon_noise_ppm(1.0)), photon_noise_0p4m_per_transit_ppm=round(photon_noise_ppm(0.4)),
                     sigma_per_transit_1m_incl_500ppm_floor=round(np.sqrt(photon_noise_ppm(1.0)**2+500**2)),
                     transit_SNR_1m_single=round(depth/np.sqrt(photon_noise_ppm(1.0)**2+500**2),1),
                     transit_SNR_0p4m_single=round(depth/np.sqrt(photon_noise_ppm(0.4)**2+800**2),1),
                     sep_0p1arcsec_au=round(0.1*c['dist'],1), sep_1arcsec_au=round(1.0*c['dist'],0)))
df=pd.DataFrame(rows); print(df.T.to_markdown())
print("\nNotes: V/Ic/J are rough transformations from Gaia G and BP-RP (+-0.3 mag). Photon noise: Ic band, 1500 A, 40% throughput, 3-arcsec aperture, dark sky 19.5 mag/arcsec^2; real systematics floors of 500 (1-m) and 800 (0.4-m) ppm per transit-length bin are ASSUMED. RV K uses Chen & Kipping 2017 M-R relation.")
