import numpy as np, json, glob
from astropy.io import fits
from astropy.wcs import WCS
np.seterr(all='ignore')
C='/root/.cache/lightkurve/mastDownload/Kepler/kplr006515734_lc_Q111100100011001000/'
D='/home/user/ZyushgPlanet/data/kbonus_batch/'
OUT={}
BITS={1:'AttitudeTweak',2:'SafeMode',4:'CoarsePoint',8:'EarthPoint',16:'ZeroCrossing',32:'Desat',64:'Argabrightening_ch',128:'CosmicRayAp',256:'ManualExclude',1024:'SuddenSensDropout',2048:'ImpulsiveOutlier',4096:'Argabrightening_field',8192:'CosmicRayCollat',16384:'DetectorAnomaly',32768:'NoFinePoint',65536:'NoData',131072:'RollingBand_ap',262144:'RollingBand_mask'}
bits=lambda v:[n for b,n in BITS.items() if v&b]
tgt=(289.74983644414,41.95718176903)
# 1. other PSF-dip cadences feeding the 36.2075 d ephemeris: Q3 272.75, Q11 1033.06/1033.14, 1069.3
tpfs={}
for fn in glob.glob(C+'*.fits.gz'):
    h=fits.open(fn); tpfs[int(h[0].header['QUARTER'])]=h
def tpfcheck(qq,tt):
    h=tpfs[qq]; d=h[1].data; t=d['TIME']; fl=d['FLUX']; qual=d['QUALITY']
    j=np.nanargmin(abs(t-tt))
    w=np.arange(max(0,j-24),min(len(t),j+25)); w=w[(w!=j)&(qual[w]==0)&np.isfinite(t[w])]
    ref=np.nanmedian(fl[w],axis=0); sig=1.4826*np.nanmedian(abs(fl[w]-ref),axis=0); z=(fl[j]-ref)/sig
    px=WCS(h[2].header).all_world2pix([list(tgt)],0)[0]
    r,c=int(round(px[1])),int(round(px[0]))
    pc=[float((d['POS_CORR%d'%k][j]-np.nanmedian(d['POS_CORR%d'%k][w]))/(1.4826*np.nanmedian(abs(d['POS_CORR%d'%k][w]-np.nanmedian(d['POS_CORR%d'%k][w]))))) for k in (1,2)]
    return dict(time=float(t[j]),cadenceno=int(d['CADENCENO'][j]),quality=int(qual[j]),bits=bits(int(qual[j])),quality_pm2=[int(x) for x in qual[j-2:j+3]],channel=int(h[0].header['CHANNEL']),
        whole_stamp_rel=float(np.nansum(fl[j])/np.nansum(ref)-1),target_3x3_rel=float(np.nansum((fl[j]-ref)[max(0,r-1):r+2,max(0,c-1):c+2])/np.nansum(ref[max(0,r-1):r+2,max(0,c-1):c+2])),
        n_px_below_m3=int(np.nansum(z<-3)),n_px_above_3=int(np.nansum(z>3)),min_z=float(np.nanmin(z)),max_z=float(np.nanmax(z)),pos_corr_z=pc,
        cadences_in_quarter=int(len(t)), is_first_or_last_good=bool(j<=w.min() or j>=w.max()))
OUT['other_event_cadences']={}
for qq,tt in ((3,272.75),(3,272.73),(11,1033.0632),(11,1033.1449),(11,1069.29),(11,1069.32),(11,1001.229),(3,349.495)):
    try: OUT['other_event_cadences']['Q%d_%.3f'%(qq,tt)]=tpfcheck(qq,tt)
    except Exception as e: OUT['other_event_cadences']['Q%d_%.3f'%(qq,tt)]=str(e)
# 2. control baseline for the cross-source test: random Q2 cadences
rng=np.random.default_rng(2)
files=sorted(glob.glob(D+'*.fits'))
data=[]
for fn in files[:200]:
    try:
        g=fits.open(fn); d=g['LIGHTCURVE_STITCHED'].data
        m=(d['QUARTER']==2)&(d['SAP_QUALITY']==0)&np.isfinite(d['FLUX'])
        if m.sum()<500: continue
        data.append((d['TIME'][m],d['FLUX'][m],d['CADENCENO'][m])); g.close()
    except: pass
cands=[4473]+list(rng.choice(np.arange(3000,6500),40,replace=False))
res={}
for cn in cands:
    zs=[]
    for tt,yy,cc in data:
        j=np.where(cc==cn)[0]
        if len(j)==0: continue
        j=j[0]; w=(tt>tt[j]-0.5)&(tt<tt[j]+0.5); w[j]=False
        if w.sum()<10: continue
        base=np.median(yy[w]); sc=1.4826*np.median(abs(yy[w]-base)); zs.append((yy[j]-base)/sc)
    zs=np.array(zs); res[int(cn)]=dict(n=int(len(zs)),frac_below_m3=float(np.mean(zs<-3)),n_below_m10=int(np.sum(zs<-10)),median_z=float(np.median(zs)))
ctrl=[v for k,v in res.items() if k!=4473]
OUT['cross_source_control_Q2']=dict(cadence_4473=res[4473],control_n_cadences=len(ctrl),control_frac_below_m3_median=float(np.median([v['frac_below_m3'] for v in ctrl])),control_frac_below_m3_max=float(np.max([v['frac_below_m3'] for v in ctrl])),control_n_below_m10_max=int(np.max([v['n_below_m10'] for v in ctrl])),control_n_below_m10_mean=float(np.mean([v['n_below_m10'] for v in ctrl])))
# 3. per-quarter heavy tails for this source
f=fits.open(D+'hlsp_kbonus-bkg_kepler_kepler_kic-006515744_kepler_v1.0_lc.fits'); s=f['LIGHTCURVE_STITCHED'].data
ok=(s['SAP_QUALITY']==0)&np.isfinite(s['FLUX']); t=s['TIME'][ok]; fl=s['FLUX'][ok]; q=s['QUARTER'][ok]
tails={}
for qq in np.unique(q):
    m=q==qq; y=fl[m]/np.median(fl[m])-1; mad=1.4826*np.median(abs(y-np.median(y)))
    tails[int(qq)]=dict(n=int(m.sum()),mad_ppm=float(mad*1e6),n_below_m10000ppm=int(np.sum(y<-0.01)),n_above_10000ppm=int(np.sum(y>0.01)),n_below_m30000ppm=int(np.sum(y<-0.03)))
OUT['per_quarter_tails']=tails
json.dump(OUT,open('/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006515744/redteam_extra.json','w'),indent=1)
print(json.dumps(OUT,indent=1))
