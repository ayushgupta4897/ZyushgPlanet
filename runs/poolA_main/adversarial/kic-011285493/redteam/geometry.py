import warnings; warnings.filterwarnings('ignore')
import glob, numpy as np, lightkurve as lk
from astropy.io import fits
from astropy.wcs import WCS
ra_t,dec_t=285.1259332116789,49.016594587983434   # target Gaia
ra_n,dec_n=285.12944145127454,49.01770673400924    # Kp13.6 neighbour
ra_e,dec_e=285.21619,49.08947                      # KIC 11285625 EB
print('--- our TPF (KIC 11285496): channel, ref col/row, target CCD col/row')
for f in sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr011285496_lc_*/*targ.fits.gz')):
    h0=fits.getheader(f,0); h2=fits.getheader(f,2); w=WCS(h2)
    xt,yt=w.all_world2pix(ra_t,dec_t,0); xn,yn=w.all_world2pix(ra_n,dec_n,0); xe,ye=w.all_world2pix(ra_e,dec_e,0)
    c0=h2['CRVAL1P']; r0=h2['CRVAL2P']
    print(f"Q{h0['QUARTER']:2d} s{h0['SEASON']} ch{h0['CHANNEL']:3d} mod{h0['MODULE']}.{h0['OUTPUT']} stamp col0={c0} row0={r0}  target col={c0+xt:.1f} row={r0+yt:.1f} | nbr col={c0+xn:.1f} row={r0+yn:.1f} | EB(ext) col={c0+xe:.1f} row={r0+ye:.1f}")
print('--- EB TPF (KIC 11285625)')
sr=lk.search_targetpixelfile('KIC 11285625',mission='Kepler',cadence='long'); print(len(sr),'TPFs')
want={}
for i in range(len(sr)):
    q=int(str(sr.mission[i]).split('Q')[-1]) if 'Q' in str(sr.mission[i]) else i
    want.setdefault(q,i)
for q in sorted(want):
    try:
        tpf=sr[want[q]].download(quality_bitmask='default')
    except Exception as e: print('Q',q,'dl err',e); continue
    h0=tpf.hdu[0].header; h2=tpf.hdu[2].header; w=WCS(h2)
    xe,ye=w.all_world2pix(ra_e,dec_e,0); c0=h2['CRVAL1P']; r0=h2['CRVAL2P']
    ap=tpf.pipeline_mask; cols=np.where(ap.any(axis=0))[0]+c0; rows=np.where(ap.any(axis=1))[0]+r0
    print(f"Q{h0['QUARTER']:2d} s{h0['SEASON']} ch{h0['CHANNEL']:3d} mod{h0['MODULE']}.{h0['OUTPUT']} Kp{h0['KEPMAG']} stamp {tpf.shape[1]}x{tpf.shape[2]} col0={c0} row0={r0} EB col={c0+xe:.1f} row={r0+ye:.1f} aperture cols {cols.min()}-{cols.max()} rows {rows.min()}-{rows.max()}  medflux={np.nanmedian(np.nansum(tpf.flux.value[:,ap],axis=1)):.0f} e/s")
    if q in (2,3,4,5): continue
