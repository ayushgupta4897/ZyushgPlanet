# SAP vs PSF depth ratio using MEASURED flux levels instead of header CROWDSAP.
# If the dip is on the target with true depth d, SAP depth = d * (F_target_in_SAP / F_SAP) = d * FLFRCSAP * F_psf / F_sap.
# A dip on a resolved neighbour (outside the PSF core) would give SAP/PSF depth ratio > F_psf/F_sap... i.e. > ~0.8 and typically >> 1.
import numpy as np, json
from astropy.io import fits
f='/home/user/ZyushgPlanet/data/kbonus_batch/hlsp_kbonus-bkg_kepler_kepler_kic-003831969_kepler_v1.0_lc.fits'
h=fits.open(f)
lc=json.load(open('/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-003831969/redteam_work/lc_tests.json'))
seasons={'S_Q1,5,9,13,17':[1,5,9,13,17],'S_Q2,6,10,14':[2,6,10,14],'S_Q3,7,11,15':[3,7,11,15],'S_Q4,8,12,16':[4,8,12,16]}
out={}
rows=[]
for sname,qs in seasons.items():
    fp=[];fs=[];fl=[];cr=[]
    for q in qs:
        e=h[f'LIGHTCURVE_Q{q}']; d=e.data; good=(d['SAP_QUALITY']==0)&np.isfinite(d['FLUX'])&np.isfinite(d['SAP_FLUX'])
        fp.append(np.nanmedian(d['FLUX'][good])); fs.append(np.nanmedian(d['SAP_FLUX'][good])); fl.append(e.header['FLFRCSAP']); cr.append(e.header['CROWDSAP'])
    Fp=np.mean(fp); Fs=np.mean(fs); FL=np.mean(fl); CR=np.mean(cr)
    dp,ep,_=lc['target']['FLUX']['per_season'][sname]; ds,es,_=lc['target']['SAP_FLUX']['per_season'][sname]
    r=ds/dp; er=r*np.hypot(ep/dp,es/ds)
    exp_fluxlevel=FL*Fp/Fs      # expected SAP/PSF depth ratio if dip on target, from measured flux levels
    exp_header=CR               # expected if header CROWDSAP correct and PSF depth = true depth
    rows.append(dict(season=sname,F_psf=round(Fp,1),F_sap=round(Fs,1),FLFRCSAP=round(FL,3),CROWDSAP_hdr=round(CR,3),
        depth_psf=dp,depth_sap=ds,ratio=round(r,3),ratio_err=round(er,3),expect_fluxlevel=round(exp_fluxlevel,3),expect_header=round(exp_header,3),
        sigma_vs_fluxlevel=round((exp_fluxlevel-r)/er,2),sigma_vs_header=round((exp_header-r)/er,2)))
# global
dp=lc['target']['FLUX']['depth_ppm']; ep=lc['target']['FLUX']['err_randbox_ppm']; ds=lc['target']['SAP_FLUX']['depth_ppm']; es=lc['target']['SAP_FLUX']['err_randbox_ppm']
r=ds/dp; er=r*np.hypot(ep/dp,es/ds)
Fp=lc['target']['FLUX']['median_flux_e_s']; Fs=lc['target']['SAP_FLUX']['median_flux_e_s']
FLg=np.mean([h[f'LIGHTCURVE_Q{q}'].header['FLFRCSAP'] for q in range(1,18)])
glob=dict(ratio=round(r,3),ratio_err=round(er,3),F_psf=round(Fp,1),F_sap=round(Fs,1),FLFRCSAP_mean=round(FLg,3),
    expect_fluxlevel=round(FLg*Fp/Fs,3),expect_fluxlevel_FL1=round(Fp/Fs,3),
    sigma_vs_fluxlevel=round((FLg*Fp/Fs-r)/er,2),sigma_vs_fluxlevel_FL1=round((Fp/Fs-r)/er,2),
    note='expected SAP/PSF depth ratio for an on-target dip = FLFRCSAP*F_psf/F_sap; ratio>F_psf/F_sap would require the eclipsing source to contribute more to SAP than to the PSF-fit flux, i.e. a resolved neighbour')
out=dict(per_season=rows,global_=glob)
json.dump(out,open('/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-003831969/redteam_work/sap_psf_fluxlevel.json','w'),indent=1)
for rr in rows: print(rr)
print(glob)
