# sedcreator
sedcreator is a package that has two main classes, SedFluxer and SedFitter. SedFluxer performs aperture photometry on a given image, coordinates and aperture size. It has a number of functions to print useful information and to plot the image together with apertures. SedFitter fits observations to a grid of models following the Zhang&Tan radiative transfer models.

![alt text](https://github.com/fedriani/sedcreator/raw/main/sedcreator_logo.png)


## Installation

type the following in the terminal to install:

.. code:: bash
    $ pip install sedcreator

## usage

.. code:: python

    from sedcreator import SedFluxer,SedFitter

    YSO_image = SedFluxer(image)
    results = YSO_image.get_flux(central_coords,aper_rad,inner_annu,outer_annu)

    results.value #get the values from the aperture photometry
    results.info #get useful information
    results.plot() #plots the image togethet with apertures

    YSO_fit = SedFitter(extinction_law,lambda_array,flux_array,err_flux_array,upper_limit_array,filter_array)

    YSO_fit.sed_fit(dist,AV_max)

## Example


An example can be found in the notebook 'working_example_sed_creator.ipynb'. Here it is copied two simple example for the use of SedFluxer and SedFitter.

.. code:: python

    #few packages needed for the example
    import numpy as np
    from astroquery.skyview import SkyView
    from astropy.coordinates import SkyCoord

    #main sed_creator package
    from sedcreator import SedFluxer,SedFitter

    #START main block for SedFluxer example
    #defining some coordinates
    G33_coord = SkyCoord(ra='18h52m50.273s', dec='00d55m29.564s', frame='icrs')

    #Downloading WISE 3.4 data around these coordinates
    G33_image = SkyView.get_images(position=G33_coord, survey=['WISE 3.4'],pixels=500)

    #Initialasing the class
    G33 = SedFluxer(G33_image)

    #Calling SedFluxer get_flux() to perform aperture photometry on the image given the coords with an aperture radius of 10 arcsec and annulus 1 and 2 times the aperture radius to subtract the background

    results_G33 = G33.get_flux(G33_coord,aper_rad=10.0,inner_annu=1.0 * 10.0,outer_annu=2.0 * 10.0)

    #get the values of the flux background subtracted and flux with background
    G33_flux_bkg,G33_flux = results_G33.value

    #to print useful info do
    results_G33.info

    #to plot the image and the apertures use do
    results_G33.plot()
    #END main block for SedFluxer example

    #START main block for SedFitter example

    #defining some data to run the function

    G35N_filt_wav_arr1 = np.array([3.5999999,4.5,5.8000002,8.0,19.7000008,31.5,
                              37.0999985,70.0,160.0,350.0,500.0]) #micron

    G35N_flux_arr1 = np.array([4.970000e-01,1.240000e+00,1.840000e+00,3.220000e+00,6.818000e+01,5.528000e+02,
                         1.192500e+03,2.627600e+03,2.385600e+03,4.286000e+02,1.274000e+02]) #Jy

    G35N_errup_arr1 = np.array([0.2160000,0.1000000,0.4340000,0.5220000,0.1800000,0.1000000,
                             0.1000000,0.1000000,0.1760000,0.3870000,0.5360000]) #error percentage

    G35N_filter_arr1 = np.array(['I1','I2','I3','I4','F8','L1','L4','P1','P3','P5','P6']) #filter names

    G35N_dist = 2200.0 #distance to the source in pc

    G35N_upper_limit_idx = np.array([1,1,0,0,0,0,0,0,0,1,1],dtype=bool) #setting the indexes for upper limits. IMPORTANT TO SET IT TO BOOL

    #initialasing the SedFitter class

    YSO_fit =SedFitter('khm',G35N_filt_wav_arr1,G35N_flux_arr1,G35N_errup_arr1,G35N_upper_limit_idx,G35N_filter_arr1)

    #performing the SED grid fit
    YSO_fit.sed_fit(dist=G35N_dist,AV_max=100)

    #END main block for SedFitter example
