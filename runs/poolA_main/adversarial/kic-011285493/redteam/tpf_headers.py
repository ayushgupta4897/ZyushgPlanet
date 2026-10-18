import glob
from astropy.io import fits
fs=sorted(glob.glob('/root/.cache/lightkurve/mastDownload/Kepler/kplr011285496_lc_*/*targ.fits.gz'))
print(len(fs),'TPFs cached')
for f in fs:
    h=fits.getheader(f,0); h1=fits.getheader(f,1)
    print(f.split('/')[-1], 'Q',h.get('QUARTER'),'season',h.get('SEASON'),'mod',h.get('MODULE'),'out',h.get('OUTPUT'),'ch',h.get('CHANNEL'),'Kp',h.get('KEPMAG'),'Teff',h.get('TEFF'),'logg',h.get('LOGG'),'R',h.get('RADIUS'),'crowd',h1.get('CROWDSAP'),'flfrc',h1.get('FLFRCSAP'),'ncols',h1.get('NAXIS2'))
