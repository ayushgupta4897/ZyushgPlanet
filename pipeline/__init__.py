"""ZyushgPlanet transit-search pipeline for Kepler Bonus (KBONUS-BKG) light curves.

Modules
-------
fetch   : locate + download all KBONUS-BKG quarters for a Gaia DR3 id / sky position, describe the FITS
prep    : stitch quarters, quality mask, per-quarter normalisation, upward-only clipping, wotan detrending
search  : astropy BoxLeastSquares (1-60 d) and transitleastsquares (restricted grid) searches
vet     : fold, odd/even, secondary, per-transit depths, empirical noise, inverted null, neighbour check
report  : JSON results + PNG figures
run     : command-line entry point (python -m pipeline.run ...)
"""
__version__ = "0.1.0"
