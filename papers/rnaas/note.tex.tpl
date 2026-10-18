\documentclass[RNAAS]{aastex701}
\usepackage{amsmath}

\begin{document}

\title{Three Sub-Neptune Transit Candidates Around Kepler Bonus Background Stars}

<<AUTHORLINE>>
\affiliation{<<AFFILIATION>>}
\email{<<EMAIL>>}

\begin{abstract}
I report three transiting sub-Neptune candidates from a pre-registered search of <<n_targets>> dwarf stars that fell on Kepler target-pixel files but were never Kepler targets, using Kepler Bonus photometry. KIC 6699887 ($P = <<P1s>>$\,d, <<rp1s>>\,$R_\oplus$), KIC 3831969 (<<P2s>>\,d, <<rp2s>>\,$R_\oplus$) and KIC 11138423 (<<P3s>>\,d, <<rp3s>>\,$R_\oplus$) show transits that localize on the target in difference images, no odd--even or secondary signal, no ephemeris match to any known Kepler object, and transit-derived stellar densities consistent with Gaia DR3. Using TRICERATOPS with companion constraints from Gaia astrometry and imaging, their false-positive probabilities are <<fppmin>>--<<fppmax>>\% and their nearby-star false-positive probabilities are zero. KIC 6699887 also shows a second on-target signal at <<P47s>>\,d. High-resolution imaging is needed before these can be considered validated.
\end{abstract}

\section{Sample and search}
Kepler Bonus \citep{MartinezPalomera2023} provides point-spread-function photometry for 606,900 sources on Kepler target-pixel files, most of them background stars. Bienias \& Szab\'o (\citeyear{Bienias2025}) reported hot-Jupiter candidates among them, and the Kepler Bonus team is running its own search. I froze a list of <<n_targets>> sources before searching: PSF fraction $>0.5$, crowding metric $>0.7$, $\geq8$ quarters, $G<17.5$, dwarfs by Gaia DR3 \citep{Gaia2023}, never observed as Kepler targets. I normalized each quarter, detrended with a 1-day biweight \citep{Hippke2019}, and ran box least squares \citep{Kovacs2002} over 0.5--40\,d with thresholds registered in advance. The same search ran on every light curve inverted ($2-f$): <<n_s1_real>> stars passed on real flux and <<n_s1_inv>> on inverted flux. Two positive controls (KIC 5689121 and Kepler-138) passed before any pool result was read.

\section{Vetting}
Of these, <<n_stage2>> passed a second stage (masked re-detrending, per-transit depths, gap proximity, variability, duration). Some second-stage criteria and both pixel tests were defined after an early look at the results; every revision is dated in the public log, and I treat everything after the first screen as post hoc. For localization I differenced in- and out-of-transit images from the originating target-pixel files and fitted point sources at the Gaia positions of the target and of every neighbor within two pixels, plus a scaled stamp image, per quarter and per Kepler season (Figure~\ref{fig}). The test places KIC 5689121 on target and two known off-target signals on other stars. <<w_pix>> of the <<n_stage2>> localize on target. Applied to the <<n_s1_inv>> inverted signals, the second stage passes <<inv_s2>> and the pixel test places <<inv_pix>> on target. Each signal then received an independent re-derivation with separate code and an adversarial review, and <<w_rt>> of the <<w_pix_l>> survived. Across all reviewed signals, rejections included post-downlink thermal ramps, single-cadence artifacts, a saturated neighbor, and an eclipsing binary $5.6'$ away on the same detector channel. No candidate matches the period (within 0.1\%, harmonics 1:4--4:1) and epoch of any KOI, DR25 TCE \citep{Thompson2018} or Kepler eclipsing binary \citep{Kirk2016} in the Kepler field, and the brighter star on each stamp shows no dip ($\geq<<host_excl_min>>\sigma$).

\section{False-positive probability}
I ran TRICERATOPS \citep{Giacalone2021} in Kepler mode, eight runs of $10^6$ draws each, removing a neighbor only when its own light curve excludes the eclipse it would need. With default priors the false-positive probabilities (FPP) are <<fppb1>>, <<fppb2>> and <<fppb3>>\%, almost entirely a planet orbiting an unresolved bound companion. The hosts have Gaia RUWE <<ruwe_min>>--<<ruwe_max>>. MOLUSC \citep{Wood2021}, using RUWE and Gaia imaging and a multiplicity fraction of 0.46, reduces the companion prior from <<prior_def_min>>--<<prior_def_max>>\% to <<prior_mol_min>>--<<prior_mol_max>>\%, giving FPP $= <<fppm1>> \pm <<fpps1>>$, $<<fppm2>> \pm <<fpps2>>$ and $<<fppm3>> \pm <<fpps3>>$\% (mean and scatter over runs) and zero nearby-star FPP. A fourth survivor, KIC 3963519, has FPP $\approx$ <<fpp4>>\% because of an unmodeled $G=18$ neighbor and is not included. Transit fits \citep{Kreidberg2015,ForemanMackey2013} give stellar densities within $<<dens_max>>\sigma$ of Gaia FLAME values \citep{Creevey2023}.

\section{Candidates}
Ephemerides are BJD$_{\rm TDB}$; radii use Gaia FLAME stellar radii; S/N is the depth over the scatter of same-duration boxes at random phases, and SDE is the BLS peak significance (in parentheses, the strongest peak on inverted flux).
\begin{itemize}
\item KIC 6699887 (Gaia DR3 <<gaia1>>, TIC <<tic1>>; $T_{\rm eff}=<<teff1>>$\,K, <<dist1>>\,pc): $P=<<P1>>$\,d, $T_0=<<T01>>$, depth <<dep1>>\,ppm, $T_{14}=<<T141>>$\,h, $R_p=<<rp1>>\,R_\oplus$, S/N <<snr1>> over <<nev1>> transits, SDE <<sde1>> (<<sdei1>>). A second signal ($P=<<P47>>$\,d, <<dep47>>\,ppm, <<n47>> transits) passes the same vetting and the seasonal pixel test.
\item KIC 3831969 (Gaia DR3 <<gaia2>>, TIC <<tic2>>; <<teff2>>\,K, <<dist2>>\,pc): $P=<<P2>>$\,d, $T_0=<<T02>>$, depth <<dep2>>\,ppm, $T_{14}=<<T142>>$\,h, $R_p=<<rp2>>\,R_\oplus$, S/N <<snr2>> over <<nev2>> transits, SDE <<sde2>> (<<sdei2>>).
\item KIC 11138423 (Gaia DR3 <<gaia3>>, TIC <<tic3>>; <<teff3>>\,K, <<dist3>>\,pc): $P=<<P3>>$\,d, $T_0=<<T03>>$, depth <<dep3>>\,ppm, $T_{14}=<<T143>>$\,h, $R_p=<<rp3>>\,R_\oplus$, S/N <<snr3>> over <<nev3>> transits, SDE <<sde3>> (<<sdei3>>).
\end{itemize}
None appears in public candidate lists as of 2026 October 17. Statistical limits on companions within $0.5''$ are the main remaining uncertainty; speckle or adaptive-optics imaging of these $G=16$--17 hosts would test them. TESS observed each host in <<tess_nsec>> sectors, but at $T=<<tess_tmin>>$--<<tess_tmax>> even an undiluted, ideal-aperture estimate \citep{Stassun2018} gives a combined transit S/N of at most <<tess_snr_max>>, so TESS can neither confirm nor exclude these signals. I used Kepler Bonus light curves (\dataset[doi:10.17909/7jbr-w430]{https://doi.org/10.17909/7jbr-w430}) and Kepler long-cadence data, Q0--Q17 (\dataset[doi:10.17909/T9488N]{https://doi.org/10.17909/T9488N}), from MAST, and the KOI cumulative, DR25 TCE, Kepler stellar and Planetary Systems tables of the NASA Exoplanet Archive (\dataset[doi:10.26133/NEA4]{https://doi.org/10.26133/NEA4}, \dataset[doi:10.26133/NEA15]{https://doi.org/10.26133/NEA15}, \dataset[doi:10.26133/NEA6]{https://doi.org/10.26133/NEA6}, \dataset[doi:10.26133/NEA12]{https://doi.org/10.26133/NEA12}), accessed <<known_date>>. Code, light curves, run outputs and evidence packs are archived at <<DATA_DOI>>.

\begin{figure}[ht!]
\plotone{figure1.pdf}
\caption{Top: phase-folded Kepler Bonus photometry (gray), binned (black) with the transit model (red). Bottom: missing flux during transit (S/N) for one Kepler season; red circle, target Gaia position; blue crosses, other Gaia sources. $\Delta\chi^2$ values compare a source at the target with no signal and with the best alternative, summed over four seasons; the season shown is labeled.\label{fig}}
\end{figure}

\begin{acknowledgments}
This paper includes data collected by the Kepler mission and obtained from the MAST data archive at the Space Telescope Science Institute (STScI); funding for the Kepler mission is provided by the NASA Science Mission Directorate, and STScI is operated by the Association of Universities for Research in Astronomy, Inc., under NASA contract NAS 5--26555. It uses data from the European Space Agency mission Gaia, processed by the Gaia Data Processing and Analysis Consortium, and has made use of the NASA Exoplanet Archive, operated by the California Institute of Technology under contract with NASA under the Exoplanet Exploration Program. I used Claude (Anthropic) to help write analysis code and draft the text; I reviewed all results and the text and take responsibility for them.
\end{acknowledgments}

\facilities{Kepler, Gaia, Exoplanet Archive}

\begin{thebibliography}{}
\bibitem[Bienias \& Szab\'o(2025)]{Bienias2025} Bienias, D., \& Szab\'o, R. 2025, A\&A, 694, A209
\bibitem[Creevey et al.(2023)]{Creevey2023} Creevey, O.~L., et al. 2023, A\&A, 674, A26
\bibitem[Foreman-Mackey et al.(2013)]{ForemanMackey2013} Foreman-Mackey, D., et al. 2013, PASP, 125, 306
\bibitem[Gaia Collaboration et al.(2023)]{Gaia2023} Gaia Collaboration, Vallenari, A., et al. 2023, A\&A, 674, A1
\bibitem[Giacalone et al.(2021)]{Giacalone2021} Giacalone, S., et al. 2021, AJ, 161, 24
\bibitem[Hippke et al.(2019)]{Hippke2019} Hippke, M., et al. 2019, AJ, 158, 143
\bibitem[Kirk et al.(2016)]{Kirk2016} Kirk, B., et al. 2016, AJ, 151, 68
\bibitem[Kov\'acs et al.(2002)]{Kovacs2002} Kov\'acs, G., Zucker, S., \& Mazeh, T. 2002, A\&A, 391, 369
\bibitem[Kreidberg(2015)]{Kreidberg2015} Kreidberg, L. 2015, PASP, 127, 1161
\bibitem[Mart\'inez-Palomera et al.(2023)]{MartinezPalomera2023} Mart\'inez-Palomera, J., Hedges, C., \& Dotson, J. 2023, AJ, 166, 265
\bibitem[Stassun et al.(2018)]{Stassun2018} Stassun, K.~G., et al. 2018, AJ, 156, 102
\bibitem[Thompson et al.(2018)]{Thompson2018} Thompson, S.~E., et al. 2018, ApJS, 235, 38
\bibitem[Wood et al.(2021)]{Wood2021} Wood, M.~L., Mann, A.~W., \& Kraus, A.~L. 2021, AJ, 162, 128
\end{thebibliography}

\end{document}
