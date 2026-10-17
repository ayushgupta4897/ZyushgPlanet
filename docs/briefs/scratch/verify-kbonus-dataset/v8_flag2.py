import numpy as np, pandas as pd, glob, warnings, subprocess, os; warnings.filterwarnings("ignore")
from astropy.io import fits
src=pd.read_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/source_catalog.parquet")
S="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/sample/"
local={fits.getheader(f)['GAIAID']:f for f in glob.glob(S+"*.fits")}
rows=src[src.gaia_designation.isin(local)]
print("local sample stars with a '2' in flag:", rows.psf_fraction_flag.str.contains('2').sum())
# pick a catalogue star with mixed flag (both 1 and 2) and avail=0&flag=1 slots, small file: few quarters
cand=src[(src.psf_fraction_flag.str.contains('2'))&(src.psf_fraction_flag.str.contains('1'))&(src.psf_avail.str.count('1')<=6)&(src.phot_g_mean_mag<15)]
print("candidates:", len(cand)); r=cand.iloc[0]; print(r[['gaia_designation','kic','tpforg','fname','psf_avail','sap_avail','psf_fraction_flag','psffrac','phot_g_mean_mag']].to_string())
fn=r.fname; sid=fn.split('-')[-1] if fn.startswith('gaia') else fn.split('-')[-1]
d4=sid[:4]; url=f"https://archive.stsci.edu/hlsps/kbonus-bkg/lcs/{d4}/{sid}/hlsp_kbonus-bkg_kepler_kepler_{fn}_kepler_v1.0_lc.fits"
out="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/verify-kbonus-dataset/dl/"+os.path.basename(url); os.makedirs(os.path.dirname(out),exist_ok=True)
print(url); subprocess.run(["curl","-sS","-L","--max-time","120","-o",out,url],check=True); print("size", os.path.getsize(out))
with fits.open(out) as h:
    print("QDETECT", h[0].header['QDETECT'])
    hd=[x for x in h if x.name.startswith('LIGHTCURVE_Q')]
    print("HDUs:", [(x.header['QUARTER'], round(x.header['PSFFRAC'],3), int(np.isfinite(x.data['FLUX']).sum()), int(np.isfinite(x.data['SAP_FLUX']).sum()), x.header.get('ITERNEG')) for x in hd])
    print("stitched quarters:", sorted(set(h[1].data['QUARTER'].tolist())), "stitched finite FLUX per quarter:", {int(q): int(np.isfinite(h[1].data['FLUX'][h[1].data['QUARTER']==q]).sum()) for q in np.unique(h[1].data['QUARTER'])})
# Also: a star with avail=0 & flag=1: check FLUX NaN in that quarter
g="Gaia DR3 2102612715913245696"; f=local[g]; row=src[src.gaia_designation==g].iloc[0]
print("\n",g,"avail",row.psf_avail,"flag",row.psf_fraction_flag)
with fits.open(f) as h:
    for x in h:
        if x.name.startswith('LIGHTCURVE_Q'):
            q=x.header['QUARTER']; print(f"  Q{q} avail={row.psf_avail[q]} flag={row.psf_fraction_flag[q]} finiteFLUX={int(np.isfinite(x.data['FLUX']).sum())}/{len(x.data)} finiteSAP={int(np.isfinite(x.data['SAP_FLUX']).sum())} PSFFRAC={x.header['PSFFRAC']:.3f} ITERNEG={x.header.get('ITERNEG')} finite PSF_FLUX_NOVA={int(np.isfinite(x.data['PSF_FLUX_NOVA']).sum())}")
    st=h[1].data
    print("  stitched finite FLUX by quarter:", {int(q): int(np.isfinite(st['FLUX'][st['QUARTER']==q]).sum()) for q in np.unique(st['QUARTER'])})
# Count usable quarters correctly across the catalogue
A=np.array([list(s) for s in src.psf_avail])=='1'; F=np.array([list(s) for s in src.psf_fraction_flag])=='1'
src['nq_usable']=(A&F).sum(1); src['nq_any']=A.sum(1); src['nq_flag1']=F.sum(1)
print("\nusable (avail==1 & flag==1) quarter count distribution head:", src.nq_usable.value_counts().sort_index().to_dict())
kep=set(pd.read_csv("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/keplerstellar_kepids.csv").kepid.astype(int)); tp=set(src.tpforg.unique())
bkg=~(src.kic.isin(kep)|src.kic.isin(tp))
q=(src.psffrac>0.5)&(src.crowdsap>0.7)
for name,col in [("nq_any",src.nq_any),("nq_flag1",src.nq_flag1),("nq_usable",src.nq_usable)]:
    print(f"background & psffrac>0.5 & crowdsap>0.7 & {name}>=8: {(bkg&q&(col>=8)).sum()}  (G<17: {(bkg&q&(col>=8)&(src.phot_g_mean_mag<17)).sum()})")
ms=pd.read_parquet("/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/mstars_joined.parquet")
m=ms.merge(src[['gaia_designation','nq_usable','nq_any','nq_flag1']],on='gaia_designation')
e=m[(m.psffrac>0.5)&(m.crowdsap>0.7)&(m.kic!=m.tpforg)&~m.kic.isin(kep)]
print("M-star eligible by nq_any>=8:", (e.nq_any>=8).sum(), " by nq_usable>=8:", (e.nq_usable>=8).sum(), " nq_any>=8 but nq_usable<8:", ((e.nq_any>=8)&(e.nq_usable<8)).sum())
print("note: psffrac (catalogue) is the MINIMUM over detected quarters, so psffrac>0.5 implies every detected quarter has PSFFRAC>0.5 -> flag has no '2'? check:", (src[src.psffrac>0.5].psf_fraction_flag.str.contains('2')).sum(), "rows with psffrac>0.5 yet a '2' in flag")
src[['gaia_designation','nq_usable','nq_any','nq_flag1']].to_parquet("nq_counts.parquet")
