import numpy as np, glob, warnings; warnings.filterwarnings("ignore")
import lightkurve as lk
D="/tmp/claude-0/-home-user-ZyushgPlanet/0a31cad6-114a-5cca-8219-924182316e90/scratchpad/kbonus-dataset/downloads/"
P=13.769448; T0=2455630.245431-2454833; dur=2.3/24
files=sorted(glob.glob(D+"mastDownload/Kepler/kplr005689116_lc_Q*/*llc.fits")); print(len(files),"official LC files")
for col in ['pdcsap_flux','sap_flux']:
    deps=[];
    for f in files:
        lc=lk.read(f, quality_bitmask='hardest'); lc=lc.select_flux(col).remove_nans().normalize()
        fl=lc.flatten(window_length=301, polyorder=2)  # ~6 d window; transits only 2.3 h
        t=fl.time.value; ph=((t-T0+P/2)%P)-P/2
        intr=np.abs(ph)<0.4*dur; oot=(np.abs(ph)>dur/2+0.05)&(np.abs(ph)<0.5)
        if intr.sum()<5: continue
        deps.append((1-np.mean(fl.flux.value[intr]))*1e6); 
    # pooled estimate over all quarters
    allf=[];allph=[]
    for f in files:
        lc=lk.read(f, quality_bitmask='hardest').select_flux(col).remove_nans().normalize().flatten(window_length=301)
        t=lc.time.value; allph.append(((t-T0+P/2)%P)-P/2); allf.append(lc.flux.value)
    ph=np.concatenate(allph); fl=np.concatenate(allf); intr=np.abs(ph)<0.4*dur; oot=(np.abs(ph)>dur/2+0.05)&(np.abs(ph)<0.5)
    dep=(1-np.mean(fl[intr]))*1e6; err=np.std(fl[oot])/np.sqrt(intr.sum())*1e6
    print(f"KIC 5689116 official {col}: box depth at Kuber ephemeris = {dep:.0f} +/- {err:.0f} ppm (n_in={intr.sum()}); per-quarter depths median {np.median(deps):.0f}")
