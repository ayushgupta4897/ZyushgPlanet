import numpy as np, pandas as pd
from scipy import stats
pd.set_option('display.width',250)
print("## Analytic trial counts and Gaussian-tail expectations (white noise)")
def ofir_n(S, pmin, pmax, dur_min, OS): return int(np.ceil(OS*S/dur_min*np.log(pmax/pmin)))
rows=[]
for lab,S,pmin,pmax in [("Kepler 4 yr",1460,0.5,25),("Kepler 4 yr 1-30 d",1460,1,30),("TESS 2 sectors",54,0.5,25),("TESS 10 sectors/5 yr",1887,0.5,25),("K2 80 d",80,0.5,25)]:
    for OS in [1,3]:
        nper=ofir_n(S,pmin,pmax,1/24,OS)
        # independent phase trials per period ~ P/duration; use geometric mean over log-uniform period grid
        Pg=np.exp(np.linspace(np.log(pmin),np.log(pmax),200)); nphase=np.mean(Pg/(2/24))
        neff=nper*6*nphase/OS   # OS-oversampled periods are correlated -> divide back; 6 durations partly correlated (over-estimate)
        exp_max=np.sqrt(2*np.log(neff))  # expected max of neff standard normals (leading order)
        thr_1_in_10k=stats.norm.isf(1e-4/neff)
        rows.append(dict(dataset=lab, OS=OS, n_periods=nper, n_eff_trials=f"{neff:.1e}", expected_max_SNR_white=round(exp_max,2), SNR_for_lt1_FA_per_10k_stars=round(thr_1_in_10k,2), P_white_SNR_ge9=f"{neff*stats.norm.sf(9):.1e}"))
print(pd.DataFrame(rows).to_markdown(index=False))
print("\nTLS/BLS single-star FAP vs SDE in white noise (TLS docs, >10,000 noise-only runs): SDE 5.7 (FAP 10%), 6.1 (5%), 7.0 (1%), 8.3 (0.1%), 9.1 (0.01%).")
print("Expected number of SDE>=9.1 white-noise false alarms: per 1,000 stars = 0.1; per 10,000 stars = 1; per 100,000 = 10. For <1 FA per 10k stars one needs FAP<1e-4 per star, i.e. SDE > ~9.1 (white) and empirically higher under red noise.")
print("\n## Kuber's matched inverted null (editorial audit 2026-10-13)")
def pois_ci(k, cl=0.9):
    lo=0.5*stats.chi2.ppf((1-cl)/2, 2*k) if k>0 else 0.0; hi=0.5*stats.chi2.ppf(1-(1-cl)/2, 2*k+2); return lo,hi
N=1736
for lab,k in [("first-screen inverted signals (BLS SNR>=9 etc.)",32),("inverted stars passing",28),("inverted stars at candidate strength (P>1.5 d, BLS SNR>=26.7, holdout SNR>=15.5)",2)]:
    lo,hi=pois_ci(k); print(f"{lab}: {k}/{N} -> per-star rate {k/N:.2e} (90% CI {lo/N:.1e}-{hi/N:.1e}); expected in a 10k-star search: {1e4*k/N:.0f} (CI {1e4*lo/N:.0f}-{1e4*hi/N:.0f})")
real_first=31; inv_first=32
print(f"Real vs inverted first-screen passes: {real_first} vs {inv_first}; Poisson two-sample test p = {stats.binomtest(real_first, real_first+inv_first, 0.5).pvalue:.2f} -> no excess of real passes.")
k=2; print(f"Probability that >=1 of {N} *inverted* stars reaches the candidate's scores, given the measured rate: 1 - exp(-{k}) = {1-np.exp(-k):.2f}. The two observed inverted passes are themselves the evidence.")
print(f"The candidate's detection statistic alone therefore carries a per-search false-alarm expectation of ~{k} (90% CI {pois_ci(2)[0]:.1f}-{pois_ci(2)[1]:.1f}) equivalent events; its credibility must come from localization/dilution/pixel tests, not from BLS SNR.")
print(f"Resolving power of a 1,736-star inverted null: the smallest per-star FAR it can bound (0 passes, 95%%) is {3/N:.1e}; a 1e-4 per-star FAR (1 FA/10k stars) needs >= {int(3/1e-4):,} inverted stars with zero passes.")
print("\n## Chopra (docs/GUIDE.md)")
n=526
print(f"743 real vs 752 flipped signals passing search cuts (SDE>=7, SNR>=6, >=3 transits): binomial p = {stats.binomtest(743,743+752,0.5).pvalue:.2f} -> search passes are noise-dominated.")
print(f"Flag-free: 10 real vs 17 flipped -> real excess? p(one-sided) = {stats.binomtest(10,27,0.5,alternative='greater').pvalue:.2f} (none).")
print(f"0 flag-free flipped signals with SDE>=9 in {n} stars -> 95% upper limit on per-star FAR {3/n:.1e}; expected flag-free SDE>=9 fakes in the survey < 3 (as GUIDE.md states).")
print(f"Max flipped flag-free SDE 8.4 vs real 11.3 and 10.3: with white-noise FAP(SDE) table, SDE 11.3 ~ 1e-6-ish per star, but the flipped maximum of 8.4 over 526 stars implies the red-noise tail is already ~1.4x the white-noise expectation (white: {stats.norm.isf(1/n):.1f} sigma-equivalent).")
print(f"Scrambles: best fake SDE 8.5 and 9.3 in 38 scrambles each -> p < 1/38 = {1/38:.3f} per candidate (not corrected for 526 stars: Sidak-corrected p ~ {1-(1-1/38)**(1/526):.1e} per-star equivalent needed).")
print("\n## eumemic (README): inverted-star calibration on 798,769 stars, scaled to 4,039,203 searched stars")
rows=[dict(cut="SNR>=9, clean",expected_FA=96,observed_real=9938-2462+19, note="19 inverted found; 96 after scaling"),dict(cut="SNR>=9, total incl. neighbour/centroid classes",expected_FA=2462,observed_real=9938,note="FA fraction ~25% of real signals at SNR>=9"),dict(cut="SNR 7.5-9 total",expected_FA=2756,observed_real=3482,note="FA fraction ~79%")]
print(pd.DataFrame(rows).to_markdown(index=False))
print("Scaling: 798,769 inverted stars resolve per-star FARs down to ~4e-6 (3/798769), i.e. <1 FA per 10k stars is testable only at that scale.")
