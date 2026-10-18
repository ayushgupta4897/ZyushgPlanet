import numpy as np, json, glob
from astropy.io import fits
np.seterr(all='ignore')
C='/root/.cache/lightkurve/mastDownload/Kepler/kplr006515734_lc_Q111100100011001000/'
OUT={}
# Kepler QUALITY bits
BITS={1:'AttitudeTweak',2:'SafeMode',4:'CoarsePoint',8:'EarthPoint',16:'ZeroCrossing',32:'Desat',64:'Argabrightening_ch',128:'CosmicRayAp',256:'ManualExclude',1024:'SuddenSensDropout',2048:'ImpulsiveOutlier',4096:'Argabrightening_field',8192:'CosmicRayCollat',16384:'Detector anomaly',32768:'NoFinePoint',65536:'NoData',131072:'RollingBand_ap',262144:'RollingBand_mask',524288:'PossibleThrusterFiring',1048576:'ThrusterFiring'}
def bits(v):
    return [n for b,n in BITS.items() if v & b]
tgt_ra,tgt_dec=289.74983644414,41.95718176903
for fn,cadno,label in ((C+'kplr006515734-2009259160929_lpd-targ.fits.gz',4473,'Q2_cad4473'),(C+'kplr006515734-2012277125453_lpd-targ.fits.gz',61720,'Q14_cad61720')):
    h=fits.open(fn); d=h[1].data; hd=h[1].header; ap=h[2].data
    j=np.where(d['CADENCENO']==cadno)[0][0]
    t=d['TIME']; fl=d['FLUX']; qual=d['QUALITY']
    w=np.arange(max(0,j-24),min(len(t),j+25)); w=w[(w!=j)&(qual[w]==0)&np.isfinite(t[w])]
    ref=np.nanmedian(fl[w],axis=0); sig=1.4826*np.nanmedian(abs(fl[w]-ref),axis=0)
    diff=fl[j]-ref; z=diff/sig
    # target pixel position via WCS
    from astropy.wcs import WCS
    wcs=WCS(h[2].header)
    px=wcs.all_world2pix([[tgt_ra,tgt_dec]],0)[0]
    tot_ref=np.nansum(ref); tot_j=np.nansum(fl[j])
    OUT[label]=dict(file=fn.split('/')[-1],time=float(t[j]),quality=int(qual[j]),quality_bits=bits(int(qual[j])),quality_pm3=[int(x) for x in qual[j-3:j+4]],
        module=int(hd.get('MODULE',h[0].header.get('MODULE'))),output=int(hd.get('OUTPUT',h[0].header.get('OUTPUT'))),channel=int(hd.get('CHANNEL',h[0].header.get('CHANNEL'))),
        stamp_shape=list(fl[j].shape),target_px_colrow=[float(px[0]),float(px[1])],
        whole_stamp_rel_change=float(tot_j/tot_ref-1),
        n_pix_z_below_m3=int(np.nansum(z<-3)),n_pix_z_above_3=int(np.nansum(z>3)),
        min_z=float(np.nanmin(z)),min_z_px_rowcol=[int(x) for x in np.unravel_index(np.nanargmin(z),z.shape)],
        max_z=float(np.nanmax(z)),max_z_px_rowcol=[int(x) for x in np.unravel_index(np.nanargmax(z),z.shape)],
        diff_e_s_sum=float(np.nansum(diff)),
        diff_e_s_at_target_3x3=float(np.nansum(diff[max(0,int(round(px[1]))-1):int(round(px[1]))+2,max(0,int(round(px[0]))-1):int(round(px[0]))+2])),
        ref_e_s_at_target_3x3=float(np.nansum(ref[max(0,int(round(px[1]))-1):int(round(px[1]))+2,max(0,int(round(px[0]))-1):int(round(px[0]))+2])),
        zmap=np.round(np.nan_to_num(z),1).tolist(),
        # also: how many cadences in this quarter have the same QUALITY bits as this cadence (context)
        n_cad_with_same_quality=int(np.sum(qual==qual[j])),
        # time-series of whole-stamp sum around the cadence
        stamp_sum_pm5=[float(np.nansum(fl[k])/tot_ref-1) for k in range(j-5,j+6)])
    # original Kepler SAP for the host at this cadence
    pos=(d['POS_CORR1'][j],d['POS_CORR2'][j]) if 'POS_CORR1' in d.names else None
    OUT[label]['pos_corr']=[float(x) for x in pos] if pos is not None else None
    OUT[label]['pos_corr_z']=[float((d['POS_CORR%d'%k][j]-np.nanmedian(d['POS_CORR%d'%k][w]))/(1.4826*np.nanmedian(abs(d['POS_CORR%d'%k][w]-np.nanmedian(d['POS_CORR%d'%k][w]))))) for k in (1,2)] if pos is not None else None
json.dump(OUT,open('/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-006515744/redteam_tpf_check.json','w'),indent=1)
for k,v in OUT.items():
    print('==',k); 
    for a,b in v.items():
        if a!='zmap': print(' ',a,b)
    print('  zmap (rows x cols):'); 
    for r in v['zmap']: print('   ',r)
