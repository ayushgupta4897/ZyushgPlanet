import lightkurve as lk, sys, warnings; warnings.filterwarnings('ignore')
out=sys.argv[1]
for qq in (5,13,1):
    try:
        sr=lk.search_targetpixelfile('KIC 5950394', author='Kepler', cadence='long', quarter=qq)
        print(qq, sr)
        tpf=sr.download(download_dir=out)
        print('ok', qq, tpf.path, tpf.shape)
    except Exception as e:
        print('fail', qq, repr(e)[:300])
