import numpy as np, warnings; warnings.filterwarnings("ignore")
from astropy.io import fits
import lightkurve as lk
D="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/"
T=D+"mastDownload/HLSP/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits"
N=D+"mastDownload/HLSP/hlsp_kbonus-bkg_kepler_kepler_kic-005689116_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689116_kepler_v1.0_lc.fits"
def cdpp6(t,f,nbin=13):
    """6-hr CDPP a la lightkurve: savgol flatten (window ~2 d), sigma clip, bin 13 cadences, std*1e6."""
    lc=lk.LightCurve(time=t,flux=f).remove_nans().normalize()
    flat=lc.flatten(window_length=101, polyorder=2).remove_outliers(sigma=5)
    b=flat.bin(bins=len(flat)//nbin) if False else flat
    # manual binning in groups of 13 consecutive cadences
    fl=np.asarray(b.flux.value); n=len(fl)//nbin*nbin
    binned=fl[:n].reshape(-1,nbin).mean(1)
    return 1e6*np.std(binned), 1e6*np.std(np.asarray(flat.flux.value)), float(flat.estimate_cdpp(transit_duration=nbin))
def per_quarter(path,label):
    print(f"\n=== {label} ===")
    with fits.open(path) as h:
        d=h[1].data; g=(d['SAP_QUALITY']==0)
        for col in ['FLUX','SAP_FLUX','PSF_FLAT_FLUX']:
            f=d[col][g]; t=d['TIME'][g]
            c6,pp,lkc=cdpp6(t,f)
            print(f"  stitched {col:14s}: 6-h CDPP (manual bin) = {c6:6.0f} ppm ; per-cadence rms after flatten = {pp:6.0f} ppm ; lightkurve estimate_cdpp(13) = {lkc:6.0f} ppm")
        print("  per quarter 6-h CDPP [PSF FLUX / SAP_FLUX]:")
        out=[]
        for hdu in h:
            if hdu.name.startswith('LIGHTCURVE_Q'):
                q=hdu.data; gg=(q['SAP_QUALITY']==0)&np.isfinite(q['FLUX'])&np.isfinite(q['SAP_FLUX'])
                if gg.sum()<500: continue
                a=cdpp6(q['TIME'][gg],q['FLUX'][gg])[0]; b=cdpp6(q['TIME'][gg],q['SAP_FLUX'][gg])[0]
                out.append((hdu.header['QUARTER'],a,b))
        print("   "+"  ".join(f"Q{q}:{a:.0f}/{b:.0f}" for q,a,b in out))
        return out
per_quarter(T,"KBONUS target KIC 5689121 (G=16.28)")
per_quarter(N,"KBONUS neighbour KIC 5689116 (G=15.48), a real Kepler target")
# official Kepler pipeline for the neighbour
print("\n=== Official Kepler pipeline (author='Kepler') for KIC 5689116 ===")
sr=lk.search_lightcurve("KIC 5689116", author="Kepler", cadence="long")
print(sr)
lcs=sr.download_all(download_dir=D)
res=[]
for lc in lcs:
    q=lc.meta['QUARTER']
    for col in ['pdcsap_flux','sap_flux']:
        l=lc.copy(); l.flux=l[col]; l=l[l.quality==0].remove_nans()
        if len(l)<500: continue
        c=cdpp6(l.time.value,l.flux.value)[0]
        res.append((q,col,c))
import collections
byq=collections.defaultdict(dict)
for q,col,c in res: byq[q][col]=c
print("   per quarter 6-h CDPP [PDCSAP / SAP]:")
print("   "+"  ".join(f"Q{q}:{v.get('pdcsap_flux',np.nan):.0f}/{v.get('sap_flux',np.nan):.0f}" for q,v in sorted(byq.items())))
print("   header CDPP6_0 keyword (pipeline's own 6-h CDPP) per quarter:", {lc.meta['QUARTER']:lc.meta.get('CDPP6_0') for lc in lcs})
# stitched pdcsap
st=lcs.stitch(corrector_func=lambda x: x[x.quality==0].remove_nans().normalize())
c6,pp,lkc=cdpp6(st.time.value, st.flux.value)
print(f"   stitched PDCSAP 6-h CDPP (manual) = {c6:.0f} ppm; lightkurve estimate = {lkc:.0f} ppm")
