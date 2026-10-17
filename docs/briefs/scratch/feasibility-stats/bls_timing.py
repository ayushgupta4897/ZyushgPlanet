import numpy as np, time, sys
from astropy.timeseries import BoxLeastSquares
# synthetic Kepler-like LC: 4 yr, 30-min cadence, ~46k points (with gaps)
rng = np.random.default_rng(1)
t = np.arange(0, 1460, 29.4/1440.)
# emulate ~70% duty cycle -> 46k points
mask = rng.random(t.size) < (46000/t.size)
t = t[mask]
f = 1 + rng.normal(0, 1500e-6, t.size)
print("N points:", t.size)
# Period grid: Ofir 2014-style frequency spacing; use astropy autoperiod with coarse oversampling
bls = BoxLeastSquares(t, f, dy=np.full(t.size, 1500e-6))
durations = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 6.0])/24.
for (pmin, pmax, freq_factor) in [(1.0, 30.0, 1.0), (0.5, 25.0, 1.0), (1.0, 30.0, 0.5)]:
    periods = bls.autoperiod(durations, minimum_period=pmin, maximum_period=pmax, frequency_factor=freq_factor)
    t0 = time.time()
    res = bls.power(periods, durations, objective='snr')
    dt = time.time()-t0
    print(f"P {pmin}-{pmax} d, freq_factor {freq_factor}: {len(periods)} periods x {len(durations)} durations -> {dt:.1f} s; peak power {res.power.max():.2f} at P={res.period[np.argmax(res.power)]:.4f}")
    sys.stdout.flush()
# TESS-like: 2 sectors 10-min cadence ~ 7000 points, and 10 sectors ~ 36k
for npts, baseline, label in [(7000, 54, 'TESS 2 sectors 10-min'), (36000, 2200, 'TESS 10 sectors (sparse, spanning 6 yr)')]:
    t2 = np.sort(rng.choice(np.arange(0, baseline, 10/1440.), size=npts, replace=False)) if baseline*144 > npts else np.arange(0, baseline, 10/1440.)
    f2 = 1 + rng.normal(0, 3000e-6, t2.size)
    bls2 = BoxLeastSquares(t2, f2, dy=np.full(t2.size, 3000e-6))
    periods = bls2.autoperiod(durations, minimum_period=0.5, maximum_period=25.0, frequency_factor=1.0)
    t0 = time.time()
    res = bls2.power(periods, durations, objective='snr')
    print(f"{label}: N={t2.size}, {len(periods)} periods -> {time.time()-t0:.1f} s")
    sys.stdout.flush()
