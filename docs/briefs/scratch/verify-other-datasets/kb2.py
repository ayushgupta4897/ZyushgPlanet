from astropy.io import fits
import numpy as np
d='/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/other-datasets/kbonus/'
m=fits.open(d+'hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits')[1].data
s=fits.open(d+'hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits')[1].data
print('source cols',s.columns.names)
print('source rows',len(s),'mstars rows',len(m))
idx={g.strip():i for i,g in enumerate(s['gaia_designation'])}
mi=np.array([idx.get(g.strip(),-1) for g in m['gaia_designation']])
print('mstars matched to source',(mi>=0).sum())
isM=np.zeros(len(s),bool); isM[mi[mi>=0]]=True
psf=s['psffrac']; crowd=s['crowdsap']; pa=s['psf_avail']
print('psf_avail dtype',pa.dtype,'ex',repr(pa[0]),'sap_avail ex',repr(s['sap_avail'][0]))
nq=np.array([str(x).count('1') for x in pa])
kic=s['kic']; tpf=s['tpforg']
print('kic dtype',kic.dtype,'tpforg dtype',tpf.dtype, 'kic sample',kic[:3],tpf[:3])
c1=psf>0.5; c2=c1&(crowd>0.7); c3=c2&(nq>=8); bkg=(kic!=tpf); c4=c3&bkg
for lab,mask in [('M',isM),('nonM',~isM),('all',np.ones(len(s),bool))]:
    print(lab,mask.sum(),(mask&c1).sum(),(mask&c2).sum(),(mask&c3).sum(),(mask&c4).sum())
print('bkg total',bkg.sum(),'primary',(~bkg).sum())
G=s['phot_g_mean_mag']
for lab,mask in [('M bkg surv',isM&c4),('nonM bkg surv',(~isM)&c4)]:
    g=G[mask]; print(lab,'G hist [<13,13-14,14-15,15-16,16-17,17-18,18-19,>19]',np.histogram(g,bins=[0,13,14,15,16,17,18,19,30])[0],'G<16',(g<16).sum(),'G<17',(g<17).sum())
print('primary passing c3',(c3&~bkg).sum())
print('median kic_sep bkg surv',np.nanmedian(s['kic_sep'][c4]))
print('median bp-rp nonM surv',np.nanmedian((s['phot_bp_mean_mag']-s['phot_rp_mean_mag'])[(~isM)&c4]))
# Kuber's star
k=idx.get('Gaia DR3 2103792461826460288')
if k is not None:
    for c in ['kic','tpforg','psffrac','crowdsap','flfrcsap','npixsap','psf_avail','sap_avail','kic_sep','kep_mag','phot_g_mean_mag','psf_fraction_flag','pertrati','pertstd']:
        if c in s.columns.names: print(' ',c,s[c][k])
# M stars catalogue ambiguity: count gaia designations in mstars but also how many are KIC
print('mstars with kic>0',(m['kic']>0).sum() if m['kic'].dtype.kind in 'iu' else 'n/a')
