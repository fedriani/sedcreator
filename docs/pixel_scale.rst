**********************
Pixel Scale Estimation
**********************

Introduction to pixel scale
---------------------------

In order to measure fluxes using the SedFluxer class in sedcreator, we need to understand the pixel scale to properly deal with apertures for photometry as well as units transformation. The pixel scale (or scale plate) in fits files is the translation between celestial coordinates and pixels, which in turn is what python uses.

In this section, we outline some history on how the pixel scale has been treated in FITS files. We have used as a reference the papers `Representation of world coordinates in FITS from Greisen and Calabreta (2002) <https://ui.adsabs.harvard.edu/abs/2002A%26A...395.1061G/abstract>`_ and `FITS: A FLEXIBLE IMAGE TRANSPORT SYSTEM from Wells, Greisen, and Harten (1981) <https://ui.adsabs.harvard.edu/abs/1981A%26AS...44..363W/abstract>`_.

In 1981 the main keywords for the World Coordinate System (WCS) in FITS headers were define as::

     CRVALn #coordinate value at reference point
     CRPIXn #array location of the reference point in pixels
     CDELTn #coordinate increment at reference point
     CTYPEn #axis type (8 characters)
     CROTAn #rotation from stated coordinate type.

In summary, the ``CDELT`` keyword gives the pixel scale and the CROTAn the Position angle of the celestial axis with respect to the xy plane. Therefore, if we want to know the equivalence between angular distance and pixels, we have retrieve the value of ``CDELT`` (which is usually given in degrees).

In 2002, the Linear transformation matrix was presented (although a general linear transformation matrix dates from `Hanisch and Wells 1988 <http://www.cv.nrao.edu/fits/wcs/wcs88.ps.Z>`_) that aimed to augment the information from ``CDELTi`` information. This how the ``CDi_j`` and ``PCi_j`` matrices were included in the FITS formalism. The main difference between the ``CD`` and ``PC`` matrices is that in the case of ``CD`` both rotation and scaling information is given whereas in the case of ``PC`` only rotation information is stored.

Therefore since we are only interested in the pixel scale (or plate scale) values to draw our aperture and to make units transformation, we are going to focus on ``CDELT`` and ``CD`` values.

Most, if not all, FITS header include a combination of ``CDELT1, CDELT1, CD1_1, CD1_2, CD2_1,`` and ``CD2_2`` keywords. The relation between them is (`link <https://danmoser.github.io/notes/gai_fits-imgs.html>`_)::


     CD1_1 = CDELT1 * cos (CROTA2)
     CD1_2 = -CDELT2 * sin (CROTA2)
     CD2_1 = CDELT1 * sin (CROTA2)
     CD2_2 = CDELT2 * cos (CROTA2)

In some headers we can also find the keyword ``PIXSCALEi``, but this is not always the case and we will not consider it in our routine.

Retrieving pixel scale in SedFluxer
-----------------------------------

The first thing we have to warn the user of sedcreator, in particular the class SedFluxer, is that the header of the FITS file is assumed to be right. The user should be responsible of inputting a checked and correct header. We also note that the pixels are assumed to be squares, which is usally the case. Nonetheless, there may be cases where the pixel is a rectangle and sedcreator does not deal with those.

As explained above, the pixel scale information is stored either in the ``CD`` matrix and/or in the ``CDELT`` value and are normally given in degrees (sedcreator assumes this in all cases). There are then three cases we need to consider:

* 1) Only ``CD`` information is available in header

* 2) Only ``CDELT`` information is available in header

* 3) Both ``CD`` and ``CDELT`` information is available in header

In case 1, the pixel scale information is stored in the ``CDi_j`` keywords. In this case, there are two subcases, when the position angle is zero (1.1) or different than zero (1.2). In the case 1.1 the pixel information is simply retrieved from ``CD1_1`` as ``pixel_scale=(header['CD1_1']**2)**0.5``. In the case 1.2 we need to consider the fact that the position angle is different than zero and therefore, using the relation stated in the previous section, the pixel scale would then be ``pixel_scale=(header['CD1_1']**2+header['CD2_1']**2)**0.5``.

The case 2 is very simple as ``CDELT`` keyword stores the pixel scale, and in turn it is ``pixel_scale = (header['CDELT1']**2)**0.5``.

The case 3 is trickier as in some headers either the ``CD`` or ``CDELT`` values are wrongly given. Of course there are cases where both keywords are correct and it does not really matters which one one gets. In this case, ``CD`` keywords take precedence over ``CDELT`` (see examples section below).

In summary, with the assumptions above, we calculate the pixel scale as follows::

    >>> if 'CD1_1' in header:
    >>>     if 'CD2_1' in header:
    >>>         pixel_scale_deg = (header['CD1_1']**2+header['CD2_1']**2)**0.5 #deg
    >>>     elif 'CD2_1' not in header:
    >>>         pixel_scale_deg = (header['CD1_1']**2)**0.5 #deg
    >>>     else:
    >>>         raise Exception('Problem with CD values, check image header')
    >>> elif 'CDELT1' in header:
    >>>     pixel_scale_deg = (header['CDELT1']**2)**0.5 #deg
    >>> else:
    >>>     raise Exception('No CD nor CDELT were found in the header, so no pixscale could be retrieved')
    >>> pixel_scale_arcsec = pixel_scale_deg*3600.0 #arcsec


Better solution?
----------------

If you have a better solution to deal in a robust manner with pixel scales in FITS files, please open an `Issue in the github repository <https://github.com/fedriani/sedcreator/issues>`_ and share it!