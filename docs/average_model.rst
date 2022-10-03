**********************
Average Model
**********************

Introduction
------------

``get_average_model()`` is a function to average the models considered in the input table.

To obtain the average model of the input table, the geometric mean, i.e, in log space, is considered in all parameters but the visual extinction (av), inclination viewing angle (theta_view) and jet opening angle (theta_w_esc), where the arithmetic mean, i.e. in linear space, is used. In a similar way, the dispersion is calculated via geometric standard deviation for all the parameters but av, theta_view, and theta_w_esc that the standard deviation is used. To perform this geometric mean and dispersion, ``get_average_model()`` uses the functions gmean and gstd ``from scipy.stats``. More information.

The usage is very simple and following the example from ``SedFitter`` we first need to fit our measured fluxes and retrieve the models:

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
    >>> print(source_models_4p)

     SED_number   chisq   chisq_nonlim mcore ...   lbol    lbol_iso lbol_av  t_now  
    ----------- --------- ------------ ----- ... -------- --------- ------- --------
    10_01_07_13   1.03751      1.49863 160.0 ...  33070.0   14990.9 14990.9 351681.0
    10_01_07_12   1.07123      1.54733 160.0 ...  33070.0   15218.5 15218.5 351681.0
    10_01_07_14   1.10807      1.80062 160.0 ...  33070.0   14797.9 14797.9 351681.0
    10_01_07_15   1.19599      1.94349 160.0 ...  33070.0   14629.7 14629.7 351681.0
    10_01_07_11    1.2377      1.78779 160.0 ...  33070.0   15490.7 15045.6 351681.0
    10_01_07_16   1.26143      2.04983 160.0 ...  33070.0   14491.3 14491.3 351681.0
    10_01_07_17   1.33213      2.16472 160.0 ...  33070.0   14379.6 14379.6 351681.0
    03_04_06_08   1.33344      2.16684  30.0 ...  49398.0   14558.3 14558.3  37037.1
    10_01_07_18   1.37847      2.24001 160.0 ...  33070.0   14293.7 14293.7 351681.0
    10_01_07_19   1.40774      2.28758 160.0 ...  33070.0   14238.9 14238.9 351681.0
            ...       ...          ...   ... ...      ...       ...     ...      ...
    15_04_08_10 419.90125    682.33954 480.0 ... 293984.0  248772.0 34094.9  23942.2
    15_04_08_12 420.09126    682.64829 480.0 ... 293984.0  246621.0 34039.5  23942.2
    15_04_08_09 421.55693    685.03001 480.0 ... 293984.0  250185.0 34139.1  23942.2
    15_04_08_08 421.66362    685.20339 480.0 ... 293984.0  251789.0 34187.2  23942.2
    15_04_08_07 423.00606    687.38485 480.0 ... 293984.0  253969.0 34263.6  23942.2
    15_04_08_06  424.0662    689.10757 480.0 ... 293984.0  256449.0 34336.3  23942.2
    15_04_08_05 425.34389    691.18382 480.0 ... 293984.0  259815.0 34427.2  23942.2
    15_04_08_04 427.06071    693.97366 480.0 ... 293984.0  264673.0 34578.3  23942.2
    15_04_08_03 430.73338    699.94175 480.0 ... 293984.0  272536.0 34776.6  23942.2
    15_04_08_02 434.95745    706.80586 480.0 ... 293984.0  290514.0 35083.1  23942.2
    15_04_08_01 450.09227    731.39994 480.0 ... 293984.0 1085010.0 35999.7  23942.2
    Length = 8640 rows

Now with the model table we can calculate the average model based on some desired criteria, for example, let's take all models below chisq 10 and with a core radius below 30 arcsec. The function will convert the input in arcsec to pc considering the input distance in order to compare it with the core radius that is given in pc in the model grid.

   >>> average_table = SedFitter().get_average_model(models=source_models_4p,number_of_models=5,
   ...chisq_cut=5,core_radius_cut=30)

   >>> print(average_table)
