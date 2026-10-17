import sys
from astropy.io import fits
import numpy as np
d='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/other-datasets/kbonus/'
for name in ['mstars','source']:
    t=fits.open(d+f'hlsp_kbonus-bkg_kepler_kepler_{name}-catalog_kepler_v1.0_cat.fits')[1].data
    print(name,'rows',len(t)); 
    if name=='mstars': print([c for c in t.columns.names])
    cols=t.columns.names
    psf=t['psffrac']; crowd=t['crowdsap']
    pa=t['psf_avail']
    print('psf_avail dtype',pa.dtype,'example',pa[0])
    nq=np.array([str(x).count('1') for x in pa])
    kic=t['kic']; tpf=t['tpforg']
    s1=psf>0.5; s2=s1&(crowd>0.7); s3=s2&(nq>=8); s4=s3&(kic!=tpf)
    print(name,'funnel',s1.sum(),s2.sum(),s3.sum(),s4.sum())
    if name=='source':
        print('kic!=tpforg total',(kic!=tpf).sum(),'kic==tpforg',(kic==tpf).sum())
        G=t['phot_g_mean_mag']
        print('non-M? need to exclude mstars ids')
        m=fits.open(d+'hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits')[1].data
        # key: gaia designation column?
        print('cols source',cols[:40])
        key=[c for c in cols if 'gaia' in c.lower() or 'designation' in c.lower()]
        print('key cols',key)
        k=key[0]
        mset=set(m[k])
        isM=np.array([x in mset for x in t[k]])
        print('M in source',isM.sum())
        nm=~isM
        f1=nm&s1; f2=nm&s2; f3=nm&s3; f4=nm&s4
        print('nonM funnel',nm.sum(),f1.sum(),f2.sum(),f3.sum(),f4.sum())
        g=G[f4]; print('nonM G<16',(g<16).sum(),'G<17',(g<17).sum(),'G<15',(g<15).sum(),'G<14',(g<14).sum())
        print('hist',np.histogram(g,bins=[0,13,14,15,16,17,18,19,25])[0])
        gm=G[isM&s4]; print('M bkg G hist',np.histogram(gm,bins=[0,13,14,15,16,17,18,19,25])[0], 'G<16',(gm<16).sum())
        print('primary (kic==tpforg) passing cuts',(s3&(kic==tpf)).sum())
        print('median kic_sep bkg survivors',np.nanmedian(t['kic_sep'][s4]))
        bprp=t['phot_bp_mean_mag'][f4]-t['phot_rp_mean_mag'][f4]; print('median bp-rp nonM surv',np.nanmedian(bprp))
