**********************
Error estimation
**********************

Introduction
------------

Here we give some further explanation on how ``SedFluxer get_flux`` function estimates the error.

Background Method
-----------------

Here we simply consider the median value within the annulus and multiply it by the area of the main aperture.

.. figure:: _static/bkg_method.png

Fluctuation Method
------------------

.. figure:: _static/fluctuation_method.png

Aliase over the full range

.. figure:: _static/aliase_flu_method.gif