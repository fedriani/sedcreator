************
SedFitter
************

Introduction
------------

SedFitter is a class with a number of convinience functions to fit the measured fluxes input by the user.


Starting the SedFitter object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First of all we need to initialise an SedFitter object with an extinction law (default 'kmh', `Kim et al 1994 <https://ui.adsabs.harvard.edu/abs/1994ApJ...422..164K/abstract>`__), and the arrays to be used in the fit, i.e. wavelength, flux, error flux (in percentage), identify upper limits, and the filter names for each wavelength used::

    >>> from sedcreator import SedFitter
    >>> source_sed = SedFitter(extinction_law,source_wave,source_flux,source_flux_err,source_uplim,source_filter)

With this class initialised we can now fit the observed SED by setting the distance to the source (dist), the maximum visual extinction one wants to explore, and the fitting method

    >>> source_sed_results = source_sed.sed_fit(dist=source_dist,AV_max=AV_max,method=method)

In this way, the results from the SED fit will be stored in source_sed_results.
From this object we can retrieve the tables with physical parameters, together with its chisq, for each model.

    >>> source_models_3p = source_sed_results.get_model_info(keys=['mcore','sigma','mstar'],
    ... tablename='table_432models.txt')    
    >>> source_models_4p = source_sed_results.get_model_info(keys=['mcore','sigma','mstar','theta_view'],
    ... tablename='table_8640models.txt')

Below one can find an example: