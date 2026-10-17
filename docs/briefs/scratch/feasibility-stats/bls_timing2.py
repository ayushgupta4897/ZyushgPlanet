import numpy as np, time, sys, os
from astropy.timeseries import BoxLeastSquares
rng = np.random.default_rng(1)
def ofir_grid(S, pmin, pmax, dur_min, OS=3):
    fmin, fmax = 1/pmax, 1/pmin
    n = int(np.ceil(OS*S/dur_min*np.log(fmax/fmin)))
    f = fmin*np.exp(np.arange(n)*dur_min/(OS*S))
    return 1/f[::-1]
durations = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 6.0])/24.
print("durations (h):", durations*24, flush=True)
def run(label, t, sigma, pmin, pmax, OS):
    f = 1 + rng.normal(0, sigma, t.size)
    bls = BoxLeastSquares(t, f, dy=np.full(t.size, sigma))
    periods = ofir_grid(t.max()-t.min(), pmin, pmax, durations.min(), OS)
    t0 = time.time(); res = bls.power(periods, durations, objective='snr'); dt = time.time()-t0
    pw = res.power; i = np.argmax(pw); sde = (pw.max()-pw.mean())/pw.std()
    print(f"{label}: N={t.size}, span={t.max()-t.min():.0f} d, grid {pmin}-{pmax} d OS={OS}: {len(periods)} periods x {len(durations)} durations -> {dt:.1f} s; max SNR {pw.max():.2f} SDE {sde:.2f} at P={res.period[i]:.4f}", flush=True)
    return dt
t = np.arange(0, 1460, 29.4/1440.); t = t[rng.random(t.size) < (46000/t.size)]
t1 = t[t<90]
t2 = np.arange(0, 54, 10/1440.); t2 = t2[rng.random(t2.size)<0.9]
starts = np.array([0, 330, 710, 740, 1070, 1100, 1650, 1680, 1830, 1860])
t3 = np.concatenate([s+np.arange(0,27,10/1440.) for s in starts]); t3=t3[rng.random(t3.size)<0.9]
run("Kepler 1 quarter OS3", t1, 1500e-6, 0.5, 25.0, 3)
run("TESS 2 sectors 10-min OS3", t2, 3000e-6, 0.5, 25.0, 3)
run("TESS 10 sectors over 5 yr OS1", t3, 3000e-6, 0.5, 25.0, 1)
run("TESS 10 sectors over 5 yr OS3", t3, 3000e-6, 0.5, 25.0, 3)
run("Kepler 4yr OS1", t, 1500e-6, 1.0, 30.0, 1)
run("Kepler 4yr OS2", t, 1500e-6, 1.0, 30.0, 2)
run("Kepler 4yr 0.5-25d OS2", t, 1500e-6, 0.5, 25.0, 2)
print("--- white-noise null MC (max BLS 'snr' objective over grid) ---", flush=True)
for label, tt, sig, OS, nrep in [("TESS 2 sectors", t2, 3000e-6, 3, 8), ("Kepler 1 quarter", t1, 1500e-6, 3, 8), ("Kepler 4yr OS1", t, 1500e-6, 1, 3)]:
    vals=[]
    for k in range(nrep):
        f = 1 + rng.normal(0, sig, tt.size)
        bls = BoxLeastSquares(tt, f, dy=np.full(tt.size, sig))
        periods = ofir_grid(tt.max()-tt.min(), 0.5, 25.0, durations.min(), OS)
        res = bls.power(periods, durations, objective='snr'); pw=res.power
        vals.append((pw.max(), (pw.max()-pw.mean())/pw.std()))
    v=np.array(vals)
    print(f"{label}: {nrep} realisations, {len(periods)} periods: max-SNR mean {v[:,0].mean():.2f} (range {v[:,0].min():.2f}-{v[:,0].max():.2f}); SDE mean {v[:,1].mean():.2f} (range {v[:,1].min():.2f}-{v[:,1].max():.2f})", flush=True)
