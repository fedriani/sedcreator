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

The data used in this example can be accessed at this `link <https://github.com/fedriani/sedcreator/tree/main/examples/AFGL2591_data>`

Below is given a short example::

    >>> from sedcreator import SedFluxer
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> from astropy.io import fits as pyfits
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u

    >>> filename = '../examples/AFGL2591_data/AFGL2591_Herschel_70.fits'

    >>> hdu = pyfits.open(filename)[1]

    >>> AFGL2591_coords = SkyCoord(ra='20h29m24.8916s', dec='+40d11m19.388s', frame='fk5')
    >>> AFGL2591_fluxer = SedFluxer(hdu)

    >>> aperture_radius = 25.0#arcsec
    >>> AFGL2591_H70_gt = AFGL2591_fluxer.get_flux(central_coords=GC,
    ...                                  	       aper_rad=aperture_radius,inner_annu=1.0*aperture_radius,
    ...                                  	       outer_annu=2.0*aperture_radius)

The plot should look like::

    >>> AFGL2591_H70_gt.plot(cmap='rainbow')

.. plot::
   :context: close-figs
   :format: doctest
   :include-source: False

    >>> from sedcreator import SedFluxer
    >>> from astropy.io import fits as pyfits
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u

    >>> filename = '../examples/AFGL2591_data/AFGL2591_Herschel_70.fits'

    >>> hdu = pyfits.open(filename)[1]

    >>> AFGL2591_coords = SkyCoord(ra='20h29m24.8916s', dec='+40d11m19.388s', frame='fk5')
    >>> AFGL2591_fluxer = SedFluxer(hdu)

    >>> aperture_radius = 25.0#arcsec
    >>> AFGL2591_H70_gt = AFGL2591_fluxer.get_flux(central_coords=AFGL2591_coords,
    ...                                  	       aper_rad=aperture_radius,inner_annu=1.0*aperture_radius,
    ...                                  	       outer_annu=2.0*aperture_radius)

    >>> AFGL2591_H70_gt.plot(cmap='rainbow')

And printing the info::

    >>> AFGL2591_H70_gt.info  # doctest: +FLOAT_CMP
    The aperture used is 25.0 arcsec
    pixel scale is 3.2 arcsec/pixel
    ~ 7.812 pixels are used for the aperture radius
    units in the image are: Jy/pixel
    Regarding observing time:
    You are probably using HERSCHEL, look at the first extension of the header
    Regarding wavelength:
    You are probably using HERSCHEL or ALMA, look at the first extension of the header
    ############################
    Flux bkg sub 4920.804203524983 Jy
    Flux         5060.075065156685 Jy
    Background   139.27086163170225 Jy
    ############################

And finally retrieve the fluxes (background subtracted, and without background subtraction) for our own calculations
(and use with the SedFitter class)::

    >>> flux_bkg_sub,flux = AFGL2591_H70_gt.get_value()
    >>> print(flux_bkg_sub,flux)  # doctest: +FLOAT_CMP
    4920.804203524983 5060.075065156685


	
