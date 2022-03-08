************
SedFitter
************

Introduction
------------

SedFitter is a class with a number of convenience functions to fit the measured fluxes input by the user.


Starting the SedFitter object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First of all we need to initialise an ``SedFitter`` object with an extinction law (default 'kmh', `Kim et al 1994 <https://ui.adsabs.harvard.edu/abs/1994ApJ...422..164K/abstract>`__), and the arrays to be used in the fit, i.e. wavelength, flux, error flux (in percentage), identify upper limits, and the filter names for each wavelength used::

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

We first need to define the arrays to be input in the ``SedFitter`` class. That is 

* wavelength array in micrometers

* flux array in Jy

* error_flux array in percentage with respect to flux

* upp_limit as a boolean array, 1 (or True) for upper limit and 0 (or False) for non upper limit.

* filter_name array of string with the name of the filter (see below to retrieve the filter names)

 We are going to use the data for the SOMA source AFGL437::


    >>> from sedcreator import SedFitter
    >>> import numpy as np

    >>> wavelength = np.array([3.6,4.5,5.8,8.0,
    ... 7.7,19.2,31.5,37.1,
    ... 70.0,160.0,250.0,350.0,500.0]) #micron
    >>> flux = np.array([1.8835,2.0301,9.8080,14.3600,
    ... 28.4727,266.4529,718.3385,865.5905,
    ... 1126.7534,603.8391,197.4299,68.2472,18.4542]) #Jy
    >>> error_flux = np.array([0.1,0.1,0.1,0.1052,
    ... 0.1,0.1,0.1,0.1,
    ... 0.1,0.2117,0.3970,0.6091,0.9186]) #percent
    >>> upp_limit = np.array([1,1,1,1,
    ... 1,0,0,0,
    ... 0,0,0,0,0],dtype=bool)
    >>> filter_name = np.array(['I1','I2','I3','I4',
    ... 'F4','F8','L1','L4',
    ... 'P1','P3','P4','P5','P6'])

    >>> source_sed = SedFitter(extc_law='kmh',lambda_array=wavelength,flux_array=flux,err_flux_array=error_flux,
    ... upper_limit_array=upp_limit,filter_array=filter_name)

We know need to call the ``sed_fit()`` function to fit out observations, but first we need to define also the distance to the source in pc, the maximum visual extinction (AV) in mag we are going to allow the fit, and the method (the default and recomended is 'minimize'::

    >>> distance = 2000.0 #pc
    >>> AV_max = 1000.0 #mag

    >>> source_sed_results = source_sed.sed_fit(dist=distance,AV_max=AV_max,method='minimize')

Once the fit is done, we can retrieve the best models ordered by chisq using the function ``get_model_info()``. One can save the table setting tablename (e.g. tablename='se_results.txt'. We can either retrieve the 432 physical models or the full 8640 model grid::

    >>> source_models_3p = source_sed_results.get_model_info(keys=['mcore','sigma','mstar'],
    ... tablename=None)
    >>> source_models_4p = source_sed_results.get_model_info(keys=['mcore','sigma','mstar','theta_view'],
    ... tablename=None)
    >>> print(source_models_3p)

     SED_number   chisq    chisq_nonlim mcore ...   lbol   lbol_iso lbol_av  t_now  
    ----------- ---------- ------------ ----- ... -------- -------- ------- --------
    10_01_07_12     1.0734      1.55047 160.0 ...  33070.0  15218.5 15218.5 351681.0
    10_01_08_19    1.98664      2.86959 160.0 ...  78272.7  18702.0 15570.6 445867.0
    11_01_07_10    2.09091      3.02021 200.0 ...  34235.3  19413.9 18612.9 329404.0
    11_01_06_05    2.12947      3.46039 200.0 ...  19867.8  14905.2 14905.2 282819.0
    04_04_06_05    2.16712      3.52156  40.0 ...  49108.3  31242.6 17867.0  33382.6
    12_01_06_04    2.23936      2.64651 240.0 ...  20164.6  16053.3 16053.3 270376.0
    05_04_05_03     2.3662      3.07606  50.0 ...  16605.8  16419.7 15170.6  24991.8
    05_03_06_05    2.76284      3.99077  50.0 ...  46086.1  24153.8 17875.8  73957.4
    07_02_06_04    2.82566      3.33941  80.0 ...  33793.1  21093.9 17340.2 155958.0
    07_03_08_10    2.97906      4.30308  80.0 ... 119307.0  15346.0 15346.0  96353.8
            ...        ...          ...   ... ...      ...      ...     ...      ...
    11_02_02_01  943.50503   1533.19568 200.0 ...  522.287   654.49  654.49  33909.6
    08_01_01_01 1005.12119   1633.32193 100.0 ...  91.4825  136.731 136.731  67546.6
    06_02_01_01 1013.60566   1647.10919  60.0 ...  183.603  302.278 302.278  32471.2
    01_01_03_01 1020.75287   1658.72341  10.0 ...  127.848  386.266 386.266 265950.0
    07_02_01_01 1029.43729    1672.8356  80.0 ...  256.351  366.457 366.457  30140.7
    05_03_01_01  1070.3813   1739.36962  50.0 ...  397.096  626.247 626.247  14314.1
    09_01_01_01 1082.49042   1759.04693 120.0 ...  88.0842  119.239 119.239  64535.4
    11_01_01_01 1090.83179   1772.60166 200.0 ...  131.924  154.212 154.212  56715.7
    10_01_01_01 1132.13915   1839.72613 160.0 ...   98.071  121.715 121.715  59976.8
    08_02_01_01 1189.68553   1933.23899 100.0 ...  242.939  314.297 314.297  28469.4
    06_03_01_01 1211.29374   1968.35232  60.0 ...   398.87  569.117 569.117  13683.8
    Length = 432 rows

Now, we can generate very interesting plots to show our data and the best models. To do that we need first to initilise the ``ModelPlotter`` class with the object from the sed_fit::

    >>> from sedcreator import ModelPlotter
    >>> md = ModelPlotter(source_sed_results)

It is very simple then to plot, for example the best 5 SEDs from the 432 physical models::

        >>> md.plot_multiple_seds(source_models_3p[0:5],xlim=[1e0,1e3],ylim=[1e-12,1e-6],title='Best 5 SEDs models',marker='rs',cmap='gray',colorbar=False,figname=None)


.. plot::
   :context: close-figs
   :format: doctest
   :include-source: False

    >>> from sedcreator import SedFitter,ModelPlotter
    >>> import numpy as np

    >>> wavelength = np.array([3.6,4.5,5.8,8.0,
    ... 7.7,19.2,31.5,37.1,
    ... 70.0,160.0,250.0,350.0,500.0]) #micron
    >>> flux = np.array([1.8835,2.0301,9.8080,14.3600,
    ... 28.4727,266.4529,718.3385,865.5905,
    ... 1126.7534,603.8391,197.4299,68.2472,18.4542]) #Jy
    >>> error_flux = np.array([0.1,0.1,0.1,0.1052,
    ... 0.1,0.1,0.1,0.1,
    ... 0.1,0.2117,0.3970,0.6091,0.9186]) #percent
    >>> upp_limit = np.array([1,1,1,1,
    ... 1,0,0,0,
    ... 0,0,0,0,0],dtype=bool)
    >>> filter_name = np.array(['I1','I2','I3','I4',
    ... 'F4','F8','L1','L4',
    ... 'P1','P3','P4','P5','P6'])

    >>> source_sed = SedFitter(extc_law='kmh',lambda_array=wavelength,flux_array=flux,err_flux_array=error_flux,
    ... upper_limit_array=upp_limit,filter_array=filter_name)

    >>> distance = 2000.0 #pc
    >>> AV_max = 1000.0 #mag

    >>> source_sed_results = source_sed.sed_fit(dist=distance,AV_max=AV_max,method='minimize')

    >>> source_models_3p = source_sed_results.get_model_info(keys=['mcore','sigma','mstar'],
    ... tablename=None)

    >>> md = ModelPlotter(source_sed_results)

    >>> md.plot_multiple_seds(source_models_3p[0:5],xlim=[1e0,1e3],ylim=[1e-12,1e-6],title='Best 5 SEDs models',marker='rs',cmap='gray',colorbar=False,figname=None)

Let's also do a more colorful plot by plotting all SED with a chisq<50, considering this time the 8640 models::

        >>> md.plot_multiple_seds(source_models_4p[source_models_4p['chisq']<50.0],xlim=[1e0,1e3],ylim=[1e-12,1e-6],title=r'SEDs with $\chi^2<50$',marker='ks',cmap='rainbow_r',colorbar=True,figname=None)

.. plot::
   :context: close-figs
   :format: doctest
   :include-source: False

    >>> from sedcreator import SedFitter,ModelPlotter
    >>> import numpy as np

    >>> wavelength = np.array([3.6,4.5,5.8,8.0,
    ... 7.7,19.2,31.5,37.1,
    ... 70.0,160.0,250.0,350.0,500.0]) #micron
    >>> flux = np.array([1.8835,2.0301,9.8080,14.3600,
    ... 28.4727,266.4529,718.3385,865.5905,
    ... 1126.7534,603.8391,197.4299,68.2472,18.4542]) #Jy
    >>> error_flux = np.array([0.1,0.1,0.1,0.1052,
    ... 0.1,0.1,0.1,0.1,
    ... 0.1,0.2117,0.3970,0.6091,0.9186]) #percent
    >>> upp_limit = np.array([1,1,1,1,
    ... 1,0,0,0,
    ... 0,0,0,0,0],dtype=bool)
    >>> filter_name = np.array(['I1','I2','I3','I4',
    ... 'F4','F8','L1','L4',
    ... 'P1','P3','P4','P5','P6'])

    >>> source_sed = SedFitter(extc_law='kmh',lambda_array=wavelength,flux_array=flux,err_flux_array=error_flux,
    ... upper_limit_array=upp_limit,filter_array=filter_name)

    >>> distance = 2000.0 #pc
    >>> AV_max = 1000.0 #mag

    >>> source_sed_results = source_sed.sed_fit(dist=distance,AV_max=AV_max,method='minimize')

    >>> md = ModelPlotter(source_sed_results)

    >>> source_models_4p = source_sed_results.get_model_info(keys=['mcore','sigma','mstar','theta_view'],
    ... tablename=None)

        >>> md.plot_multiple_seds(source_models_4p[source_models_4p['chisq']<50.0],xlim=[1e0,1e3],ylim=[1e-12,1e-6],title=r'SEDs with $\chi^2<50$',marker='ks',cmap='rainbow_r',colorbar=True,figname=None)

It is also interesting to plot the 2D distribution of the 3 main parameters of the model, i.e., m*, sigma_cl, and M_c::

    >>> md.plot2d(source_models_4p[source_models_4p['chisq']<=50.0],title=None,figname=None)

.. plot::
   :context: close-figs
   :format: doctest
   :include-source: False

    >>> from sedcreator import SedFitter,ModelPlotter
    >>> import numpy as np

    >>> wavelength = np.array([3.6,4.5,5.8,8.0,
    ... 7.7,19.2,31.5,37.1,
    ... 70.0,160.0,250.0,350.0,500.0]) #micron
    >>> flux = np.array([1.8835,2.0301,9.8080,14.3600,
    ... 28.4727,266.4529,718.3385,865.5905,
    ... 1126.7534,603.8391,197.4299,68.2472,18.4542]) #Jy
    >>> error_flux = np.array([0.1,0.1,0.1,0.1052,
    ... 0.1,0.1,0.1,0.1,
    ... 0.1,0.2117,0.3970,0.6091,0.9186]) #percent
    >>> upp_limit = np.array([1,1,1,1,
    ... 1,0,0,0,
    ... 0,0,0,0,0],dtype=bool)
    >>> filter_name = np.array(['I1','I2','I3','I4',
    ... 'F4','F8','L1','L4',
    ... 'P1','P3','P4','P5','P6'])

    >>> source_sed = SedFitter(extc_law='kmh',lambda_array=wavelength,flux_array=flux,err_flux_array=error_flux,
    ... upper_limit_array=upp_limit,filter_array=filter_name)

    >>> distance = 2000.0 #pc
    >>> AV_max = 1000.0 #mag

    >>> source_sed_results = source_sed.sed_fit(dist=distance,AV_max=AV_max,method='minimize')

    >>> source_models_4p = source_sed_results.get_model_info(keys=['mcore','sigma','mstar','theta_view'],
    ... tablename=None)

    >>> md = ModelPlotter(source_sed_results)

    >>> md.plot2d(source_models_4p[source_models_4p['chisq']<=50.0],title=None,figname=None)

To check the name of the filter::

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
