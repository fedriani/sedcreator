*********
SedFluxer
*********

Introduction
------------

SedFluxer is a class with a number of convinience functions to measure fluxes from images
in fits format input by the user. get_flux() reads the header and does units transformation, provided
that the header is correct, else get_raw_flux() does not do any units transformation and give the
results in raw units so the user can perform its own units transformation.

Starting the SedFluxer object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First of all we need to initialise an SedFluxer object with a loaded fits file, i.e., with a 2D array that is the image
and a header. There are various ways to load in a fits file, either from your local disk using, e.g. astropy.io.fits or
from archives loaded (and downloaded if desired) into a python variable, e.g., using astroquery, in particular,
ESASky or SkyView.