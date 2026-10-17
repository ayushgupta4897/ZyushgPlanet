import numpy as np, lightkurve as lk, warnings, time
from wotan import flatten
warnings.filterwarnings('ignore')
# Measure 1-hr CDPP-like noise on real TESS-SPOC FFI light curves: KIC 5689121 (TIC 120686465, T~15), Kepler-169 (KIC 5689351 ~Kp14.5), and a brighter Kepler star
targets=[("TIC 120686465","KIC 5689121 (Kuber target)"),("TIC 120686496","Kepler-169 host? (check)"),("Kepler-169","Kepler-169"),("Kepler-138","Kepler-138 (Kp 12.9 M dwarf)"),("Kepler-42","Kepler-42 (Kp 15.9 M dwarf)")]
for name,desc in targets:
    t0=time.time()
    try:
        sr=lk.search_lightcurve(name, author='TESS-SPOC')
        if len(sr)==0: print(name, "no TESS-SPOC products"); continue
        print(name, desc, "sectors:", sorted(set(int(s) for s in sr.table['sequence_number'])), "Tmag?", end=' ')
        lc=sr[-1].download(download_dir='lkcache')
        tm=lc.meta.get('TESSMAG'); exptime=sr.table['exptime'][-1]
        tt=np.asarray(lc.time.value,float); fl=np.asarray(lc['pdcsap_flux'].value,float); ok=np.isfinite(fl)&np.isfinite(tt)&(lc['quality'].value==0)
        fn_=fl[ok]/np.nanmedian(fl[ok]); tt=tt[ok]
        flat=flatten(tt,fn_,method='biweight',window_length=0.75); flat=flat[np.isfinite(flat)]
        cad=np.nanmedian(np.diff(tt))*1440
        n=max(1,int(round(60/cad))); cs=np.convolve(flat-1,np.ones(n)/n,mode='valid')
        cdpp1h=1.4826*np.nanmedian(np.abs(cs-np.nanmedian(cs)))*1e6
        n2=max(1,int(round(138/cad))); cs2=np.convolve(flat-1,np.ones(n2)/n2,mode='valid')
        cdpp23=1.4826*np.nanmedian(np.abs(cs2-np.nanmedian(cs2)))*1e6
        p2p=1.4826*np.nanmedian(np.abs(np.diff(fn_)-np.nanmedian(np.diff(fn_))))/np.sqrt(2)*1e6
        print(f"Tmag={tm}, sector {lc.meta.get('SECTOR')}, cadence {cad:.1f} min, N={ok.sum()}, p2p {p2p:.0f} ppm, 1-h CDPP-like {cdpp1h:.0f} ppm, 2.3-h {cdpp23:.0f} ppm, CROWDSAP {lc.meta.get('CROWDSAP')}, {time.time()-t0:.0f}s", flush=True)
    except Exception as e:
        print(name, "err", repr(e)[:200], flush=True)
