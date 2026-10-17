import numpy as np
try:
    import ticgen
    print("ticgen version", getattr(ticgen,'__version__','?'))
    print("Tmag  1hr_noise_ppm  (ticgen Star(Tmag=...).oneHourNoiseLnsigma)")
    for T in [8,9,10,11,12,13,14,15,16]:
        s = ticgen.Star(Tmag=T)
        # ticgen returns ln(sigma) per 1 hour in ppm? inspect
        print(T, getattr(s,'oneHourNoiseLnsigma',None), getattr(s,'oneSigmaNoise',None))
    import inspect
    print(inspect.getsource(ticgen.Star))
except Exception as e:
    print("ticgen failed:", e)
