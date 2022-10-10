****************************
Filters: adding and removing
****************************

Introduction
------------

It may well happen that the filter response of an instrument of your interest is not in the default database. sedcreator, within the class ``SedFitter``, provides a number of functions to add, plot, and remove custom filter responses. These filters responses are used to convolve the model flux in order to compare it to observations.

The are four main functions to deal with filters in sedcreator

* add_filter()

* add_square_filter()

* remove_filter()

* plot_filter()

Let's start by checking the default filters added to the database::

    >>> SedFitter().print_default_filters

    filter wavelength   instrument  
    ------ ---------- --------------
        2J        1.2          2MASS
        2H        1.6          2MASS
        2K        2.2          2MASS
        I1        3.6   Spitzer_IRAC
        I2        4.5   Spitzer_IRAC
        I3        5.6   Spitzer_IRAC
        I4        8.0   Spitzer_IRAC
        M1       24.0   Spitzer_MIPS
        M2       70.0   Spitzer_MIPS
        M3      160.0   Spitzer_MIPS
        F1        5.4  SOFIA_FORCAST
        F2        6.4  SOFIA_FORCAST
        F3        6.6  SOFIA_FORCAST
        F4        7.7  SOFIA_FORCAST
        F5        8.6  SOFIA_FORCAST
        F6       11.1  SOFIA_FORCAST
        F7       11.3  SOFIA_FORCAST
        F8       19.2  SOFIA_FORCAST
        F9       24.2  SOFIA_FORCAST
        L1       31.5  SOFIA_FORCAST
        L2       33.6  SOFIA_FORCAST
        L3       34.8  SOFIA_FORCAST
        L4       37.1  SOFIA_FORCAST
        P1       70.0  Herschel_PACS
        P2      100.0  Herschel_PACS
        P3      160.0  Herschel_PACS
        P4      250.0 Herschel_SPIRE
        P5      350.0 Herschel_SPIRE
        P6      500.0 Herschel_SPIRE
        R1       12.0           IRAS
        R2       25.0           IRAS
        R3       60.0           IRAS
        R4      100.0           IRAS
        W1        3.4           WISE
        W2        4.6           WISE
        W3       12.0           WISE
        W4       22.0           WISE
        S1      450.0          Scuba
        S2      850.0          Scuba

Add filter
----------

If you have a filter response for your instrument, you can add it to the sedcreator database and use it afterwards for the fitting by using ``add_filter()``. The inputs needed for this function are:

* filter_name
* instrument
* lambda_array
* response_array

Take into account that the array for the wavelength (lambda_array) needs to be in microns. The units of the filter response are not important since it normalised by the peak. It is recommended that the filter name is informative but at the same time not too long, e.g., use one letter and one number or two letters and two numbers. This does not really matter because then you will need to call the filter name in the same way as you defined it for the ``SedFitter``. The wavelength shown in the database will be the median of the lambda_array, but the entire filter response will be used when convolving the model flux.

We recommend using the following database to download the filter response for any given instrument: http://svo2.cab.inta-csic.es/theory/fps/

Add square filter
-----------------

On the other hand, if your instrument has no filter response, but it can be well approximated by a top hat profile, you can use the function ``add_square_filter()``, that assumes that the full filter response is flat to 1 throught the entire width of the filter. The inputs for this function are:

* filter_name
* instrument
* filter_lambda
* filter_width

Remember that the central lambda and the width of the filter have to be given in microns.

Plot filter
-----------

If you are interested, you can also plot the filter response to make sure everything is in order.

Remove filter
-------------

Finally, you can also remove filters in case you made a mistake or are not interested in the default filters. To do that the function ``remove_filter()`` will help you.

See example notebook here
