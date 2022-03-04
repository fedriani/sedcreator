*********
SedFluxer
*********

Introduction
------------

SedFluxer is a class with a number of convenience functions to measure fluxes from images
in fits format input by the user. get_flux() reads the header and does units transformation, provided
that the header is correct, else get_raw_flux() does not do any units transformation and gives the
results in raw units so the user can perform their own units transformation.

Starting the SedFluxer object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First of all we need to initialise an SedFluxer object with a loaded fits file, i.e., with a 2D array that is the image
and a header. There are various ways to load in a fits file, either from your local disk using, e.g. astropy.io.fits or
from archives loaded (and downloaded if desired) into a python variable, e.g., using astroquery, in particular,
ESASky or SkyView::

    >>> from sedcreator import SedFluxer
    >>> fluxer = SedFluxer(hdu)

Once the SedFluxer is initialised, we can now measure fluxes performing aperture photometry by either using get_flux() for automatic unit transformation based on the header or get_raw_fux() to make our own unit transformation.
These functions need central coordinates, which are astropy.SkyCoord type, the aperture radius, inner annulus radius,
and outer annulus radius given in arcsec::

    >>> flux = fluxer.get_flux(central_coords,aper_rad,inner_annu,outer_annu)
    >>> #or
    >>> flux = fluxer.get_raw_flux(central_coords,aper_rad,inner_annu,outer_annu)

With this flux object, one can retrieve the values for the background subtracted flux and the total flux.
One can also print useful information from the header and parameters used and plot the image used together
with the aperture and annulus defined::

    >>> flux_bkg_sub,flux = flux.get_value()
    >>> flux.plot()
    >>> flux.info

Below is given a short example::

    >>> from sedcreator import SedFluxer
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> from astropy.io import fits as pyfits
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u

    >>> filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    >>> hdu = pyfits.open(filename)[0]

    >>> GC = SkyCoord(l='00:00:00.00', b='00:00:00.00', unit=(u.deg, u.deg),frame='galactic')
    >>> GC_fluxer = SedFluxer(hdu)

    >>> aperture_radius = 300.0#arcsec
    >>> GC_flux = GC_fluxer.get_raw_flux(central_coords=GC,
    ...                                  aper_rad=aperture_radius,inner_annu=1.0*aperture_radius,
    ...                                  outer_annu=2.0*aperture_radius)

The plot should look like::

    >>> GC_flux.plot(cmap='rainbow')

.. plot::
   :context: close-figs
   :format: doctest
   :include-source: False

    >>> from sedcreator import SedFluxer
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> from astropy.io import fits as pyfits
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u

    >>> filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    >>> hdu = pyfits.open(filename)[0]

    >>> GC = SkyCoord(l='00:00:00.00', b='00:00:00.00', unit=(u.deg, u.deg),frame='galactic')
    >>> GC_fluxer = SedFluxer(hdu)

    >>> aperture_radius = 300.0#arcsec
    >>> GC_flux = GC_fluxer.get_raw_flux(central_coords=GC,
    ...                                  aper_rad=aperture_radius,inner_annu=1.0*aperture_radius,
    ...                                  outer_annu=2.0*aperture_radius)
    >>> GC_flux.plot(cmap='rainbow')

And printing the info::

    >>> GC_flux.info  # doctest: +FLOAT_CMP
    The aperture used is 300.0 arcsec
    pixel scale is 24.0 arcsec/pixel
    ~ 12.5 pixels are used for the aperture radius
    units in the image are: W/m^2-sr
    Regarding observing time:
    You are probably using HERSCHEL, look at the first extension of the header
    Regarding wavelength:
    You are probably using HERSCHEL or ALMA, look at the first extension of the header
    ############################
    Flux bkg sub 0.04420493760862438 unitless
    Flux         0.05151128117223151 unitless
    Background   0.007306343563607133 unitless
    ############################

And finally retrieve the fluxes (background subtracted, and without background subtraction) for our own calculations
(and use with the SedFitter class)::

    >>> flux_bkg_sub,flux = GC_flux.get_value()
    >>> print(flux_bkg_sub,flux)  # doctest: +FLOAT_CMP
    0.04420493760862438 0.05151128117223151


	