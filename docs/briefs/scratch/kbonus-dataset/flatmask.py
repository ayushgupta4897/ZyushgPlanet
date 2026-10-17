import numpy as np, warnings; warnings.filterwarnings("ignore")
from astropy.io import fits
T="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/mastDownload/HLSP/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0/hlsp_kbonus-bkg_kepler_kepler_kic-005689121_kepler_v1.0_lc.fits"
P=13.769448; T0_bjd=2455630.245431; T0=T0_bjd-2454833  # BKJD
with fits.open(T) as h:
    d=h[1].data
    t=d['TIME']; fm=d['FLATTEN_MASK']==1; q=d['SAP_QUALITY']
    ph=((t-T0+P/2)%P)-P/2
    intr=np.abs(ph)<(2.3/24)/2
    print("cadences in transit (2.3h box):", intr.sum(), " flagged by FLATTEN_MASK:", fm.sum(), " overlap:", (intr&fm).sum())
    print("FLATTEN_MASK flagged cadences: fraction with SAP_QUALITY!=0:", (q[fm]!=0).mean())
    f=d['FLUX']; pf=d['PSF_FLAT_FLUX']
    # are flagged cadences outliers?
    med=np.nanmedian(pf); mad=1.4826*np.nanmedian(np.abs(pf-med))
    z=(pf-med)/mad
    print("flagged: median |z| of PSF_FLAT_FLUX =", np.nanmedian(np.abs(z[fm])), "; unflagged:", np.nanmedian(np.abs(z[~fm])))
    print("flagged z sign: frac positive", (z[fm]>0).mean())
    print("time clustering of flagged: first 30 times", np.round(t[fm][:30],2))
    # Compare depth at the ephemeris in FLUX (per-quarter normalized) vs PSF_FLAT_FLUX
    good=(q==0)&np.isfinite(f)
    qs=d['QUARTER']
    fn=f.copy()
    for Q in np.unique(qs):
        m=qs==Q; fn[m]=f[m]/np.nanmedian(f[m&good])
    pfn=pf/np.nanmedian(pf[good])
    oot=(np.abs(ph)>0.1)&(np.abs(ph)<0.5)&good
    for name,arr in [("FLUX per-quarter-normalized",fn),("PSF_FLAT_FLUX",pfn)]:
        dep=1-np.nanmean(arr[intr&good]); err=np.nanstd(arr[oot])/np.sqrt((intr&good).sum())
        print(f"  box depth at Kuber ephemeris in {name}: {dep*1e6:.0f} +/- {err*1e6:.0f} ppm (n_in={ (intr&good).sum()})")
    # what if we also exclude FLATTEN_MASK cadences
    g2=good&~fm
    dep=1-np.nanmean(pfn[intr&g2]); print(f"  PSF_FLAT_FLUX depth excluding FLATTEN_MASK: {dep*1e6:.0f} ppm (n_in={(intr&g2).sum()})")
    # SAP depth
    sap=d['SAP_FLUX']; sn=sap.copy()
    for Q in np.unique(qs):
        m=qs==Q; sn[m]=sap[m]/np.nanmedian(sap[m&good])
    dep=1-np.nanmean(sn[intr&good]); print(f"  SAP_FLUX per-quarter-normalized depth: {dep*1e6:.0f} ppm")
    # quarter-to-quarter medians of FLUX relative to overall
    print("  per-quarter median FLUX / global median:", {int(Q): round(float(np.nanmedian(f[(qs==Q)&good])/np.nanmedian(f[good])),3) for Q in np.unique(qs)})
