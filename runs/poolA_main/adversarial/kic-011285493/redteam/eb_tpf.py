import warnings; warnings.filterwarnings('ignore')
import numpy as np, lightkurve as lk
from astropy.wcs import WCS
ra_e,dec_e=285.21619,49.08947
sr=lk.search_targetpixelfile('KIC 11285625',mission='Kepler',cadence='long')
for want in ['Quarter 01','Quarter 02','Quarter 03','Quarter 04']:
    i=[k for k in range(len(sr)) if want in str(sr.mission[k])][0]
    tpf=sr[i].download(quality_bitmask='default')
    h0=tpf.hdu[0].header; h2=tpf.hdu[2].header; w=WCS(h2)
    xe,ye=w.all_world2pix(ra_e,dec_e,0); c0=h2['CRVAL1P']; r0=h2['CRVAL2P']
    ap=tpf.pipeline_mask; cols=np.where(ap.any(axis=0))[0]+c0; rows=np.where(ap.any(axis=1))[0]+r0
    print(f"EB Q{h0['QUARTER']} s{h0['SEASON']} ch{h0['CHANNEL']} mod{h0['MODULE']}.{h0['OUTPUT']} Kp{h0['KEPMAG']} stamp {tpf.shape[1]}x{tpf.shape[2]} col0={c0} row0={r0} EB col={c0+xe:.1f} row={r0+ye:.1f} aper cols {cols.min()}-{cols.max()} rows {rows.min()}-{rows.max()} medflux={np.nanmedian(np.nansum(tpf.flux.value[:,ap],axis=1)):.0f} e/s")
