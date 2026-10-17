import numpy as np, glob, warnings, pandas as pd; warnings.filterwarnings("ignore")
from astropy.io import fits
import lightkurve as lk
def cdpp6(t,f,nbin=13):
    lc=lk.LightCurve(time=t,flux=f).remove_nans().normalize()
    flat=lc.flatten(window_length=101, polyorder=2).remove_outliers(sigma=5)
    fl=np.asarray(flat.flux.value); n=len(fl)//nbin*nbin
    return 1e6*np.std(fl[:n].reshape(-1,nbin).mean(1)), 1e6*np.std(fl)
rows=[]
for f in sorted(glob.glob("sample/*.fits")):
    with fits.open(f) as h:
        hd=h[0].header; d=h[1].data; g=(d['SAP_QUALITY']==0)&np.isfinite(d['FLUX'])&np.isfinite(d['SAP_FLUX'])
        # per-quarter CDPP then median across quarters (avoids inter-quarter stitching artefacts)
        qc=[];qs=[]
        for hdu in h:
            if hdu.name.startswith('LIGHTCURVE_Q'):
                q=hdu.data; gg=(q['SAP_QUALITY']==0)&np.isfinite(q['FLUX'])&np.isfinite(q['SAP_FLUX'])
                if gg.sum()<800: continue
                qc.append(cdpp6(q['TIME'][gg],q['FLUX'][gg])[0]); qs.append(cdpp6(q['TIME'][gg],q['SAP_FLUX'][gg])[0])
        c_psf,pp=cdpp6(d['TIME'][g],d['FLUX'][g]); c_sap,_=cdpp6(d['TIME'][g],d['SAP_FLUX'][g]); c_flat,_=cdpp6(d['TIME'][g],d['PSF_FLAT_FLUX'][g])
        rows.append(dict(obj=hd['OBJECT'],G=hd['GMAG'],Kp=hd.get('KEPMAG'),bp_rp=hd['BPMAG']-hd['RPMAG'],nq=len(qc),med_err_ppm=1e6*np.nanmedian(d['FLUX_ERR'][g]/d['FLUX'][g]),
                         cdpp6_psf_stitched=c_psf,cdpp6_sap_stitched=c_sap,cdpp6_psfflat_stitched=c_flat,cdpp6_psf_quarter_median=np.median(qc),cdpp6_sap_quarter_median=np.median(qs),rms_percadence=pp))
df=pd.DataFrame(rows).sort_values('G')
pd.set_option('display.width',250)
print(df.round(0).to_string(index=False))
df.to_csv("sample/sample_cdpp.csv",index=False)
print("\nBinned by G (median of per-quarter-median 6-h CDPP, PSF flux):")
df['Gbin']=pd.cut(df.G,[12,13,14,15,16,17,18,19])
print(df.groupby('Gbin',observed=True)[['cdpp6_psf_quarter_median','cdpp6_sap_quarter_median','med_err_ppm','rms_percadence']].median().round(0).to_string())
