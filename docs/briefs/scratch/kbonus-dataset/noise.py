import numpy as np, warnings; warnings.filterwarnings("ignore")
from astropy.io import fits
import lightkurve as lk
D="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/mastDownload/HLSP/"
T=D+"hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits"
N=D+"hlsp_kbonus-bkg_kepler_kepler_kic-005689116_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689116_kepler_v1.0_lc.fits"
def cdpp(t,f,hours=6.0):
    # lightkurve estimate_cdpp on a normalized LC (savgol flatten + binning)
    lc=lk.LightCurve(time=t,flux=f).remove_nans().normalize()
    return float(lc.estimate_cdpp(transit_duration=int(round(hours*2))))  # 2 cadences/hour at 29.4 min... 6h ~ 13 cadences
for label,path in [("TARGET KIC 5689121",T),("NEIGHBOUR KIC 5689116",N)]:
    print("\n======", label, "======")
    with fits.open(path) as h:
        hd=h[0].header; print("GMAG",hd['GMAG'],"KEPMAG",hd['KEPMAG'],"TPFORG",hd['TPFORG'],"QDETECT",hd['QDETECT'],"PARALLAX",hd['PARALLAX'])
        d=h[1].data
        print("stitched cols:", d.columns.names)
        print("SAP_QUALITY unique (first 20):", np.unique(d['SAP_QUALITY'])[:20], "n nonzero", (d['SAP_QUALITY']!=0).sum(), "of", len(d))
        print("FLATTEN_MASK unique:", np.unique(d['FLATTEN_MASK']), "n nonzero", (d['FLATTEN_MASK']!=0).sum())
        good=(d['SAP_QUALITY']==0)&np.isfinite(d['FLUX'])&np.isfinite(d['TIME'])
        print("finite FLUX", np.isfinite(d['FLUX']).sum(), "finite SAP", np.isfinite(d['SAP_FLUX']).sum(), "finite PSF_FLAT", np.isfinite(d['PSF_FLAT_FLUX']).sum())
        t=d['TIME'][good]
        print("time first/last BKJD:", t.min(), t.max(), " -> BJD", t.min()+2454833, t.max()+2454833, "; cadence median diff [min]", np.median(np.diff(t))*1440)
        for col in ['FLUX','SAP_FLUX','PSF_FLAT_FLUX']:
            f=d[col][good]; ok=np.isfinite(f)
            print(f"  {col}: median {np.nanmedian(f):.1f} e-/s; 6-h CDPP (lightkurve) = {cdpp(t[ok],f[ok]):.0f} ppm; point-to-point rms = {1e6*np.nanstd(np.diff(f[ok]))/np.sqrt(2)/np.nanmedian(f):.0f} ppm; median FLUX_ERR/FLUX = {1e6*np.nanmedian(d[col+'_ERR'][good][ok]/f[ok]):.0f} ppm")
        # per quarter medians (offsets)
        print("  per-quarter: Q, module/channel, season, N, median PSF flux, median SAP flux, PSFFRAC, CROWDSAP, FLFRCSAP, NPIXSAP, PERTRATI, PERTSTD, SAP mode, 6h CDPP PSF, 6h CDPP SAP, median RED_CHI2")
        for hdu in h:
            if hdu.name.startswith('LIGHTCURVE_Q'):
                q=hdu.data; hh=hdu.header; g=(q['SAP_QUALITY']==0)&np.isfinite(q['FLUX'])
                try: c1=cdpp(q['TIME'][g],q['FLUX'][g]); c2=cdpp(q['TIME'][g],q['SAP_FLUX'][g])
                except Exception as e: c1=c2=np.nan
                print(f"   Q{hh['QUARTER']:2d} mod{hh['MODULE']}.{hh['OUTPUT']} ch{hh['CHANNEL']:2d} s{hh['SEASON']} N={g.sum():4d} PSF={np.nanmedian(q['FLUX'][g]):8.1f} SAP={np.nanmedian(q['SAP_FLUX'][g]):8.1f} PSFFRAC={hh['PSFFRAC']:.3f} CROWD={hh['CROWDSAP']:.3f} FLFRC={hh['FLFRCSAP']:.3f} NPIX={hh['NPIXSAP']} PERTRATI={hh['PERTRATI']:.4f} PERTSTD={hh['PERTSTD']:.2e} SAP={hh['SAP']} CDPP_PSF={c1:5.0f} CDPP_SAP={c2:5.0f} redchi2={np.nanmedian(q['RED_CHI2'][g]):.2f} cen=({np.nanmedian(q['CENTROID_COLUMN'][g]):.3f},{np.nanmedian(q['CENTROID_ROW'][g]):.3f}) cen_std=({np.nanstd(q['CENTROID_COLUMN'][g]):.2e},{np.nanstd(q['CENTROID_ROW'][g]):.2e})")
