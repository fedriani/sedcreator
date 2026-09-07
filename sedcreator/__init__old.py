#packages needed for the sedcreator to work
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os
import pkg_resources
from tqdm import tqdm
from scipy.optimize import curve_fit, minimize, Bounds
from scipy.stats import gstd,gmean

from matplotlib import rcParams
rcParams['font.sans-serif'] = ['Times']
rcParams['pdf.use14corefonts'] = True
rcParams['font.family'] = 'serif'
#rcParams['text.usetex'] = True
#rcParams['ps.useafm'] = True
rcParams['font.size'] = 12
rcParams['mathtext.fontset'] = 'cm'

from astropy.io import fits as pyfits
from astropy.io import ascii
import astropy.units as u
from astropy.table import Table, unique
from astropy.wcs import WCS
from astropy.time import Time
from astropy.visualization import simple_norm
from astropy.stats import sigma_clipped_stats
from astropy.constants import c, m_e, m_n, m_p

from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus

#constants for SedFitter
pc2cm = u.pc.to(u.cm)
Lsun2erg_s = u.L_sun.to(u.erg*u.s**-1)
c_micron_s = c.to(u.micron*u.s**-1).value
Jy2erg_s_cm2 = u.Jy.to(u.erg*u.s**-1*u.cm**-2*u.Hz**-1)

class FluxerContainer():
    '''
    A class to store the results from the SedFluxer class
    '''
    def __init__(self,data=None,flux_bkgsub=None,flux=None,
                 central_coords=None,aper_rad=None,inner_annu=None,outer_annu=None,
                 x_source=None,y_source=None,aper_rad_pixel=None,wcs_header= None,
                 aperture=None,annulus_aperture=None,mask=None,flux_method=None):
        self.data = data
        self.flux_bkgsub = flux_bkgsub
        self.flux = flux
        self.central_coords = central_coords
        self.aper_rad = aper_rad
        self.inner_annu = inner_annu
        self.outer_annu = outer_annu
        self.x_source = x_source
        self.y_source = y_source
        self.aper_rad_pixel = aper_rad_pixel
        self.wcs_header = wcs_header
        self.aperture = aperture
        self.annulus_aperture = annulus_aperture
        self.mask = mask
        self.flux_method = flux_method
        self.__value = None
        self.__info = None

    
    @property
    def value(self):
        if self.__value is None:
            self.__value = self.get_value()
        return self.__value
    
    def get_value(self):
        '''
        Outputs the background-subtracted flux and the flux including the background
        
        Returns
        -------
        flux_bkgsub,flux: numpy array
        
        '''
        
        return(self.flux_bkgsub,self.flux)
    
    @property
    def info(self):
        if self.__info is None:
            self.__info = self.get_info()
        return self.__info
    
    def get_info(self):
        '''
        Prints the results and useful information of the image based on the header and aperture photometry
        '''
        
        #reading the data
        data,header = self.data
        flux_method = self.flux_method

        #retrieves the pixel scale or size from the header
        if 'CD1_1' in header:
            pixel_scale = np.absolute(header['CD1_1'])*3600.0
        elif 'CDELT1' in header:
            pixel_scale = np.absolute(header['CDELT1'])*3600.0
        else:
            raise Exception('Neither CD1_1 nor CDELT1 were found in the header')

        print('The aperture used is', round(self.aper_rad_pixel*pixel_scale,3), 'arcsec')
        print('pixel scale is', round(pixel_scale,3), 'arcsec/pixel')
        print('~',round(self.aper_rad_pixel,3),'pixels are used for the aperture radius')
        if 'BUNIT' in header:
            print('units in the image are:', header['BUNIT'])
        elif 'FUNITS' in header:
                print('units in the image are:', header['FUNITS'])
        elif 'COMMENT' in header and len(header['COMMENT'])>=17: #work around to print date in WISE data
            print('units in the image are:',header['COMMENT'][17])
        else:
            print('No BUNIT nor FUNITS key in the header')

        #to print the observing time
        #I add this to know the date in case we are interested in variable sources
        if 'MJDSTART' in header:
            time_MJD = Time(header['MJDSTART'],format='mjd')
            print('Observing date:',time_MJD.to_value(format='iso'))
        elif 'DATE-OBS' in header:
            print('Observing date:',header['DATE-OBS'])
        elif 'COMMENT' in header and len(header['COMMENT'])>=21: #work around to print date in WISE data
            print('Observing date:',header['COMMENT'][21])
        elif 'MIDOBS' in header: #work around to print date in WISE data
            print('Observing date:',header['MIDOBS'])
        else:
            print('Regarding observing time:\nYou are probably using HERSCHEL, look at the first extension of the header')
        #to print the wavelength
        if 'WAVELEN' in header:
            print('Wavelength:',header['WAVELEN'])
        elif 'WAVELNTH' in header:
            print('Wavelength:',header['WAVELNTH'])
        elif 'FILTER' in header:
            print('Wavelength:',header['FILTER'])
        elif 'SURVEY' in header:
            print('Wavelength:',header['SURVEY'])
        #TODO: make sure the below works for other images
        #workaround to also print wavelengths in Spitzer IRAC cases
        elif 'CHNLNUM' in header:
            if header['CHNLNUM']==1:
                print('Wavelength:',3.6)
            if header['CHNLNUM']==2:
                print('Wavelength:',4.5)
            if header['CHNLNUM']==3:
                print('Wavelength:',5.8)
            if header['CHNLNUM']==4:
                print('Wavelength:',8.0)
        else:
            print('Regarding wavelength:\nYou are probably using HERSCHEL or ALMA, look at the first extension of the header')

        if flux_method=='get_flux':
            print('############################')
            print('Flux bkg sub',self.flux_bkgsub, 'Jy')
            print('Flux        ',self.flux, 'Jy')
            print('Background  ',self.flux-self.flux_bkgsub, 'Jy')
            print('############################')

        elif flux_method=='get_raw_flux':
            print('############################')
            print('Flux bkg sub',self.flux_bkgsub, 'unitless')
            print('Flux        ',self.flux, 'unitless')
            print('Background  ',self.flux-self.flux_bkgsub, 'unitless')
            print('############################')
            print('Please, perform your own units transformation')
        else:
            print('')

        
    def plot(self,cmap='gray',percent=100.0,stretch='log',colorbar=True,aperture_color='black',annulus_color='red',plot_mask=False,title=None,figname=None):
        '''
        Plots the used image together with the aperture used for the flux
        and the annulus for the background subtraction.
                 
        Parameters
        ----------
        cmap: str
            Color map from matplolib to be used in the image. Default is 'gray'
        stretch: {'linear', 'sqrt', 'power', 'log', 'asinh'}, optional
            The stretch function to apply to the image. The default is 'log'.
        percent: float, optional
            The percentage of the image values used to determine the pixel values
            of the minimum and maximum cut levels. The lower cut level will set at the
            (100 - percent) / 2 percentile, while the upper cut level will be set at the
            (100 + percent) / 2 percentile. The default is 100.0.
        aperture_color: str
                    Circular aperture color to show in the image. Default is 'black'
                    
        annulus_color: str
                    Annulus aperture color to show in the image. Default is 'red'

        plot_mask: bool, optional
                    Plots the mask used in when getting the flux. Default is False
                    
        title: str, optional
                    Set title for the image. Default is None

        figname: str, optional
                    A path, or a Python file-like object. Note that fname is used verbatim,
                    and there is no attempt to make the extension. Default is None.
                    Note that one can choose the format of the figure by changing the extension,
                    e.g., figure.pdf would generate the figure in PDF format.
        '''
        
        #reading the data
        data,header = self.data
        mask = self.mask
        
        plt.figure()
        plt.subplot(projection=self.wcs_header)
        #This is to get around no good stretch for SOFIA images
        if 'OBSERVAT' in header:
            if 'SOFIA' in header['OBSERVAT']:
                data_for_norm = data
            else:
                data_for_norm = data[int(self.y_source-5.0*self.aper_rad_pixel):int(self.y_source+5.0*self.aper_rad_pixel),
                                     int(self.x_source-5.0*self.aper_rad_pixel):int(self.x_source+5.0*self.aper_rad_pixel)]
        elif 'TELESCOP' in header:
            if 'SOFIA 2.5m' in header['TELESCOP']:
                data_for_norm = data
            else:
                data_for_norm = data[int(self.y_source-5.0*self.aper_rad_pixel):int(self.y_source+5.0*self.aper_rad_pixel),
                                     int(self.x_source-5.0*self.aper_rad_pixel):int(self.x_source+5.0*self.aper_rad_pixel)]
        else:
            data_for_norm = data[int(self.y_source-5.0*self.aper_rad_pixel):int(self.y_source+5.0*self.aper_rad_pixel),
                                 int(self.x_source-5.0*self.aper_rad_pixel):int(self.x_source+5.0*self.aper_rad_pixel)]

        norm = simple_norm(data_for_norm[data_for_norm>0], stretch=stretch, percent=percent)
        plt.imshow(data, cmap=cmap, origin='lower', norm=norm)
        plt.plot(self.x_source,self.y_source,'rx')
        plt.xlim(self.x_source-5.0*self.aper_rad_pixel,self.x_source+5.0*self.aper_rad_pixel)
        plt.ylim(self.y_source-5.0*self.aper_rad_pixel,self.y_source+5.0*self.aper_rad_pixel)
        self.aperture.plot(color=aperture_color)
        self.annulus_aperture.plot(color=annulus_color,ls='--')
        #in case we are dealing with galactic coordinates
        if 'GLON' in self.wcs_header.axis_type_names or 'glon' in self.wcs_header.axis_type_names:
            plt.xlabel('GLON')
            plt.ylabel('GLAT')
        else:
            plt.xlabel('RA (J2000)')
            plt.ylabel('Dec (J2000)')
        
        if colorbar:
            cbar_ticks = np.around(np.linspace(norm.vmin,norm.vmax,num=5))
            if 'BUNIT' in header:
                cbar = plt.colorbar(label='PixelUnits: {0}'.format(header['BUNIT']), pad=0.01)
            elif 'FUNITS' in header:
                cbar = plt.colorbar(label='PixelUnits: {0}'.format(header['FUNITS']), pad=0.01)
            elif 'COMMENT' in header and len(header['COMMENT'])>=17: #work around to print date in WISE data
                cbar = plt.colorbar(label='{0}'.format(header['COMMENT'][17]), pad=0.01)
            else:
                cbar = plt.colorbar(label='PixelUnits: check header', pad=0.01)
            cbar.set_ticks(cbar_ticks)
            cbar.set_ticklabels(cbar_ticks)
        if title is not None:
            plt.title(title)
            
        if plot_mask:
            mask_to_plot = np.array(mask,dtype=int)
            mask_cmap = plt.cm.Reds_r
            mask_cmap.set_bad(color='white',alpha=0)
            plt.imshow(mask,vmin=0,vmax=1,cmap=mask_cmap,alpha=0.1)
            
        if figname is not None:
            plt.savefig(figname,dpi=300,bbox_inches='tight')
            print('Image saved in ',figname)
        plt.show()

class SedFluxer:
    '''
    A class used to measure flux
    
    Attributes
    ----------
    get_data: function
        Get the data to measure flux
    '''
    def __init__(self,image):
        self.image = image
        self.data = self.get_data()

    def get_data(self):
        try:
            data = self.image.data
            header = self.image.header
            #this is to deal with ALMA and VLA data that has multiple extensions
            if len(np.shape(self.image.data))!=2:
                data = self.image.data[0][0]
        except:
            data = self.image[0][0].data
            header = self.image[0][0].header
            #this is to deal with ALMA and VLA data that has multiple extensions
            if len(np.shape(self.image[0][0].data))!=2:
                data = self.image[0][0].data[0][0]
        return(data,header)
        
        
    def get_flux(self,central_coords,aper_rad,inner_annu,outer_annu,mask=None):
        '''
        Performs circular aperture photometry for a given image, specifying the central coordinates and
        the aperture radius. It returns the background subtracted flux and the flux including the background in units of Jy.
        Considering the header of the image, the function makes the units transformations.
        
        Current units supported are (BUNIT):
        MJy/sr; MJY/SR; Jy/beam; mJy/beam; Jy/pix
        
        If you are not sure about the units in your header, use get_raw_flux() function and make separately the transformation.
        
        Parameters
        ----------

        central_coords: `astropy.coordinates.SkyCoord`
                    Central coordinates where the apertures will be centred on the data image.
                    This must be a SkyCoord statement.

        aper_rad: float
                Aperture radius given in arcsec. The script will take care to convert it into pixels.
                
        inner_annu: float
                Inner annulus radius given in arcsec. The script will take care to convert it into pixels.
                Usually this is set to 1*aper_rad.

        outer_annu: float
                Outer annulus radius given in arcsec. The script will take care to convert it into pixels.
                Usually this is set to 2*aper_rad.
                
        mask: 2D array bool or int
                mask with the same shape of the image where the flux wants to be calculated.
                True values (or 1) will be ignore in the aperture_photometry.

        Returns
        -------

        array: float,float
            with flux with background subtracted and flux including background in units of Jy
        '''

        data,header = self.data

                
        wcs_header = WCS(header).celestial
        x_source,y_source = wcs_header.world_to_pixel(central_coords)
        
        #retrieves the pixel scale or size from the header
        if 'CD1_1' in header:
            pixel_scale = np.absolute(header['CD1_1'])*3600.0
        elif 'CDELT1' in header:
            pixel_scale = np.absolute(header['CDELT1'])*3600.0
        else:
            raise Exception('Neither CD1_1 nor CDELT1 were found in the header')

        #defines the aperture size in pixels
        aper_rad_pixel = aper_rad/pixel_scale
        inner_annu_pixel = inner_annu/pixel_scale
        outer_annu_pixel = outer_annu/pixel_scale
        aperture = CircularAperture([[x_source,y_source]], r=aper_rad_pixel)
        annulus_aperture = CircularAnnulus([[x_source,y_source]],
                                           r_in=inner_annu_pixel, r_out=outer_annu_pixel)

        #here is the main code for the aperture photometry
        #--> START of the aperture photometry block
        annulus_masks = annulus_aperture.to_mask(method='center')

        error = 0.1*data
        bkg_median = []
        for annu_mask in annulus_masks:
            annulus_data = annu_mask.multiply(data)
            annulus_data_1d = annulus_data[annu_mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        ap_phot = aperture_photometry(data, aperture,error=error,mask=mask) #circular aperture sum
        ap_phot['annulus_median'] = bkg_median
        ap_phot['aper_bkg'] = bkg_median * aperture.area #consider a median value for the background subtraction
        ap_phot['aper_sum_bkgsub'] = ap_phot['aperture_sum'] - ap_phot['aper_bkg'] #subtracting the background
        #--> END of the aperture photometry block

        #This block is to deal with WISE DN units
        #The first if else is to deal with both WISE images the one downloaded from Skyview and from WISE archive
        if ('TELESCOP' in header) & ('BAND' in header):
            if header['TELESCOP']=='WISE':
                if ((header['TELESCOP']=='WISE') & (header['BAND']==1)):
                    M_0_inst = 20.752 #mag for W1
                    AC = 0.222 #mag for W1
                    F_nu_0 = 309.540 #Jy
                elif ((header['TELESCOP']=='WISE') & (header['BAND']==2)):
                    M_0_inst = 19.596 #mag for W2
                    AC = 0.280 #mag for W2
                    F_nu_0 = 171.787 #Jy
                elif ((header['TELESCOP']=='WISE') & (header['BAND']==3)):
                    M_0_inst = 17.800 #mag for W3
                    AC = 0.665 #mag for W3
                    F_nu_0 = 31.674 #Jy
                elif ((header['TELESCOP']=='WISE') & (header['BAND']==4)):
                    M_0_inst = 12.945 #mag for W4
                    AC = 0.616 #mag for W4
                    F_nu_0 = 8.363 #Jy
                    
                #magnitude transformation for WISE depending on the above constants (band dependent)
                Mcal_bkg = M_0_inst-2.5*np.log10(ap_phot['aper_sum_bkgsub'])-AC
                Mcal = M_0_inst-2.5*np.log10(ap_phot['aperture_sum'])-AC

                #flux conversion for WISE
                flux_bkgsub = F_nu_0*10.0**(-Mcal_bkg.data[0]/2.5) #Jy
                #DISCLAMER:Flux without bkgsub does not really makes sense, only here for completeness
                flux = F_nu_0*10.0**(-Mcal.data[0]/2.5) #Jy

            elif 'SOFIA 2.5m' in header['TELESCOP']:#To get around HAWC+ data
                #TODO: careful here, I am not sure it will work 100%.
                #TEMPORARY FIX to give some value
                if 'BUNIT' in header:
                    if 'Jy/pix' in header['BUNIT']:
                        if 'mJy/pix' in header['BUNIT']:
                            unit_factor_Jy = 0.001 #from mJy to Jy
                        else:
                            unit_factor_Jy = 1.0 #leave it in Jy
                        flux_bkgsub = unit_factor_Jy*ap_phot['aper_sum_bkgsub'].data[0] #Jy
                        flux = unit_factor_Jy*ap_phot['aperture_sum'].data[0] #Jy
                elif 'FUNITS' in header:
                    if 'Jy/pix' in header['FUNITS']:#This is mainly for SOFIA data that does not have BUNIT, it is FUNIT
                        if 'mJy/pix' in header['FUNITS']:
                            unit_factor_Jy = 0.001 #from mJy to Jy
                        else:
                            unit_factor_Jy = 1.0 #leave it in Jy
                        flux_bkgsub = unit_factor_Jy*ap_phot['aper_sum_bkgsub'].data[0] #Jy
                        flux = unit_factor_Jy*ap_phot['aperture_sum'].data[0] #Jy
                else:
                    raise Exception('Neither BUNIT nor FUNITS found in the header, use get_raw_flux() function and perform own units transformation')
                    
            else:
                raise Exception('No valid information found in the header, use get_raw_flux() function and perform own units transformation')



        elif 'SURVEY' in header:#This is when WISE is considered
            if header['SURVEY']=='WISE 3.4 micron allWISE release' or\
            header['SURVEY']=='WISE 4.6 micron allWISE release' or\
            header['SURVEY']=='WISE 12 micron allWISE release' or\
                header['SURVEY']=='WISE 22 micron allWISE release':
                if header['SURVEY']=='WISE 3.4 micron allWISE release':
                    M_0_inst = 20.752 #mag for W1
                    AC = 0.222 #mag for W1
                    F_nu_0 = 309.540 #Jy
                elif header['SURVEY']=='WISE 4.6 micron allWISE release':
                    M_0_inst = 19.596 #mag for W2
                    AC = 0.280 #mag for W2
                    F_nu_0 = 171.787 #Jy
                elif header['SURVEY']=='WISE 12 micron allWISE release':
                    M_0_inst = 17.800 #mag for W3
                    AC = 0.665 #mag for W3
                    F_nu_0 = 31.674 #Jy
                elif header['SURVEY']=='WISE 22 micron allWISE release':
                    M_0_inst = 12.945 #mag for W4
                    AC = 0.616 #mag for W4
                    F_nu_0 = 8.363 #Jy
                #magnitude transformation for WISE depending on the above constants (band dependent)
                Mcal_bkg = M_0_inst-2.5*np.log10(ap_phot['aper_sum_bkgsub'])-AC
                Mcal = M_0_inst-2.5*np.log10(ap_phot['aperture_sum'])-AC

                #flux conversion for WISE
                flux_bkgsub = F_nu_0*10.0**(-Mcal_bkg.data[0]/2.5) #Jy
                #DISCLAMER:Flux without bkgsub does not really makes sense, only here for completeness
                flux = F_nu_0*10.0**(-Mcal.data[0]/2.5) #Jy
                    
            else:
                raise Exception('No valid information found in the header, use get_raw_flux() function and perform own units transformation')

        #TODO: deal with Jy and mJy units
        else:
            # Mengyao's conversion from MJy/sr to Jy/pixel
            if 'BUNIT' in header:
                if 'MJy/sr' in header['BUNIT'] or 'MJY/SR' in header['BUNIT']: #mainly for Herschel and some Spitzer data (and IRAS)
                    if 'CD1_1' in header:
                        flux_bkgsub = ap_phot['aper_sum_bkgsub'].data[0]*304.6*(np.absolute(header['CD1_1']))**2 #Jy
                        flux = ap_phot['aperture_sum'].data[0]*304.6*(np.absolute(header['CD1_1']))**2 #Jy
                    elif 'CDELT1' in header:
                        flux_bkgsub = ap_phot['aper_sum_bkgsub'].data[0]*304.6*(np.absolute(header['CDELT1']))**2 #Jy
                        flux = ap_phot['aperture_sum'].data[0]*304.6*(np.absolute(header['CDELT1']))**2 #Jy
                    else:
                        raise Exception('Neither CD1_1 nor CDELT1 were found in the header')
                elif 'Jy/beam' in header['BUNIT']:#mainly for ALMA data or some other radio data
                    beam = np.pi/(4.0*np.log(2.0))*header['BMAJ']*header['BMIN']
                    if 'mJy/beam' in header['BUNIT']:
                        unit_factor_Jy = 0.001 #from mJy to Jy
                    else:
                        unit_factor_Jy = 1.0 #leave it in Jy
                    if 'CD1_1' in header:
                        flux_bkgsub = unit_factor_Jy*ap_phot['aper_sum_bkgsub'].data[0]/(beam/(np.absolute(header['CD1_1']))**2) #Jy
                        flux = unit_factor_Jy*ap_phot['aperture_sum'].data[0]/(beam/(np.absolute(header['CD1_1']))**2) #Jy
                    elif 'CDELT1' in header:
                        flux_bkgsub = unit_factor_Jy*ap_phot['aper_sum_bkgsub'].data[0]/(beam/(np.absolute(header['CDELT1']))**2) #Jy
                        flux = unit_factor_Jy*ap_phot['aperture_sum'].data[0]/(beam/(np.absolute(header['CDELT1']))**2) #Jy
                    else:
                        raise Exception('Neither CD1_1 nor CDELT1 were found in the header')
                elif 'Jy/pix' in header['BUNIT']:
                    if 'mJy/pix' in header['BUNIT']:
                        unit_factor_Jy = 0.001 #from mJy to Jy
                    else:
                        unit_factor_Jy = 1.0 #leave it in Jy
                    flux_bkgsub = unit_factor_Jy*ap_phot['aper_sum_bkgsub'].data[0] #Jy
                    flux = unit_factor_Jy*ap_phot['aperture_sum'].data[0] #Jy
                else:
                    raise Exception('BUNIT (',header['BUNIT'],') found in the header but units not yet supported, use get_raw_flux() function and perform own units transformation')
            #TODO: add mJy/pix here!
            elif 'FUNITS' in header:
                if 'Jy/pix' in header['FUNITS']:#This is mainly for SOFIA data that does not have BUNIT, it is FUNIT
                    if 'mJy/pix' in header['FUNITS']:
                        unit_factor_Jy = 0.001 #from mJy to Jy
                    else:
                        unit_factor_Jy = 1.0 #leave it in Jy
                    flux_bkgsub = unit_factor_Jy*ap_phot['aper_sum_bkgsub'].data[0] #Jy
                    flux = unit_factor_Jy*ap_phot['aperture_sum'].data[0] #Jy
                elif 'Jy/sq-arc' in header['FUNITS']:#This is mainly for SOFIA data
                    if 'mJy/sq-arc' in header['FUNITS']:
                        unit_factor_Jy = 0.001 #from mJy to Jy
                    else:
                        unit_factor_Jy = 1.0 #leave it in Jy
                    flux_bkgsub = unit_factor_Jy*ap_phot['aper_sum_bkgsub'].data[0]*pixel_scale**2 #Jy
                    flux = unit_factor_Jy*ap_phot['aperture_sum'].data[0]*pixel_scale**2 #Jy
                else:
                    raise Exception('FUNITS (',header['FUNITS'],') found in the header but units not yet supported, use get_raw_flux() function and perform own units transformation')
            else:
                raise Exception('Neither BUNIT nor FUNITS found in the header, use get_raw_flux() function and perform own units transformation')
                

        return FluxerContainer(data=self.data,flux_bkgsub=flux_bkgsub,flux=flux,
                 central_coords=central_coords,aper_rad=aper_rad,inner_annu=inner_annu,outer_annu=outer_annu,
                 x_source=x_source,y_source=y_source,aper_rad_pixel=aper_rad_pixel,wcs_header=wcs_header,
                 aperture=aperture,annulus_aperture=annulus_aperture,mask=mask,flux_method='get_flux')

    def get_raw_flux(self,central_coords,aper_rad,inner_annu,outer_annu,mask=None):
        '''
        Performs circular aperture photometry for a given image, specifying the central coordinates and
        the aperture radius. It returns the background subtracted flux and the flux including the background without
        any transformations of units. This functions is intended in case the user wants to apply a custom tranformation.
        
        
        Parameters
        ----------

        central_coords: `astropy.coordinates.SkyCoord`
                    Central coordinates where the apertures will be centred on the data image.
                    This must be a SkyCoord statement.

        aper_rad: float
                Aperture radius given in arcsec. The script will take care to convert it into pixels.
                
        inner_annu: float
                Inner annulus radius given in arcsec. The script will take care to convert it into pixels.
                Usually this is set to 1*aper_rad.

        outer_annu: float
                Outer annulus radius given in arcsec. The script will take care to convert it into pixels.
                Usually this is set to 2*aper_rad.

        mask: 2D array bool or int
                mask with the same shape of the image where the flux wants to be calculated.
                True values (or 1) will be ignore in the aperture_photometry.

        Returns
        -------

        array: float,float
            with flux with background subtracted and flux including background without any units transformations
        '''

        data,header = self.data

                
        wcs_header = WCS(header).celestial
        x_source,y_source = wcs_header.world_to_pixel(central_coords)
        
        #retrieves the pixel scale or size from the header
        if 'CD1_1' in header:
            pixel_scale = np.absolute(header['CD1_1'])*3600.0
        elif 'CDELT1' in header:
            pixel_scale = np.absolute(header['CDELT1'])*3600.0
        else:
            raise Exception('Neither CD1_1 nor CDELT1 were found in the header')

        #defines the aperture size in pixels
        aper_rad_pixel = aper_rad/pixel_scale
        inner_annu_pixel = inner_annu/pixel_scale
        outer_annu_pixel = outer_annu/pixel_scale
        aperture = CircularAperture([[x_source,y_source]], r=aper_rad_pixel)
        annulus_aperture = CircularAnnulus([[x_source,y_source]],
                                           r_in=inner_annu_pixel, r_out=outer_annu_pixel)

        #here is the main code for the aperture photometry
        #--> START of the aperture photometry block
        annulus_masks = annulus_aperture.to_mask(method='center')

        error = 0.1*data
        bkg_median = []
        for annu_mask in annulus_masks:
            annulus_data = annu_mask.multiply(data)
            annulus_data_1d = annulus_data[annu_mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        ap_phot = aperture_photometry(data, aperture,error=error,mask=mask) #circular aperture sum
        ap_phot['annulus_median'] = bkg_median
        ap_phot['aper_bkg'] = bkg_median * aperture.area #consider a median value for the background subtraction
        ap_phot['aper_sum_bkgsub'] = ap_phot['aperture_sum'] - ap_phot['aper_bkg'] #subtracting the background
        #--> END of the aperture photometry block
        
        flux_bkgsub = ap_phot['aper_sum_bkgsub'].data[0]
        flux = ap_phot['aperture_sum'].data[0]
        

        return FluxerContainer(data=self.data,flux_bkgsub=flux_bkgsub,flux=flux,
                 central_coords=central_coords,aper_rad=aper_rad,inner_annu=inner_annu,outer_annu=outer_annu,
                 x_source=x_source,y_source=y_source,aper_rad_pixel=aper_rad_pixel,wcs_header=wcs_header,
                 aperture=aperture,annulus_aperture=annulus_aperture,mask=mask,flux_method='get_raw_flux')
    
    def get_optimal_aperture(self,central_coords,ap_inner=5.0,ap_outer=60.0,step_size=0.25,aper_increase=1.3,threshold=1.1,profile_plot=False):
        '''
        Finds the optimal aperture based on the following method:
        EXPLAIN METHOD.
        In short: if one increases by 30% the aperture size, the flux does not increase by
        more than threshold %


        Parameters
        ----------
        image: fits file
            Image for which we would like to calculate the "optimal" aperture radius
            
        central_coords: `astropy.coordinates.SkyCoord`
                    Central coordinates where the apertures will be centred on the data image.
                    This must be a SkyCoord statement.
                    
        ap_inner: float
                Inner aperture radius given in arcsec. The script will take care to convert it into pixels.
                It the starting point to find the optimal aperture. Default is 5.0 arcsec.
                
        ap_outer: float
                Outer aperture radius given in arcsec. The script will take care to convert it into pixels.
                It the ending point to find the optimal aperture. Default is 60.0 arcsec.
                
        step_size: float
                Step size for sampling aperture radii. Default is 0.25arcsec
                
        aper_increase: float
            It is the increase in aperture to meet the condition
            aper_increase*optimal radius <= flux at opt.rad. * threshold.
            Default is 1.30, i.e. 30% increase in aperture.
            
        threshold: float
            It is the condition to be met such as, at the optimal aperture radius,
            1.3*optimal radius <= flux at opt.rad. * threshold. Default is 1.08, i.e. 10% increase in flux.
            
        profile_plot: bool
            Plots the flux profile versus the aperture radius size. Default is False
            
        Returns
        -------
        opt_rad: float
            Optimal aperture radius defined by this method given in arcsec

        '''
        
        FLUX_BKG = []
        FLUX = []


        #Step size for sampling aperture radii
        #Sample aperture radii
        APER_RAD = np.arange(ap_inner, ap_outer, step_size)#arcsec


        for radius in APER_RAD:
            try:
                flux_bkg, flux = self.get_flux(central_coords=central_coords,
                                               aper_rad=radius,
                                               inner_annu=1.0*radius,
                                               outer_annu=2.0*radius).value
                unit = 1 # to consider Jy in label
            except:
                flux_bkg, flux = self.get_raw_flux(central_coords=central_coords,
                                               aper_rad=radius,
                                               inner_annu=1.0*radius,
                                               outer_annu=2.0*radius).value
                unit = 0 # to consider Jy in label
                
            FLUX_BKG.append(flux_bkg)
            FLUX.append(flux)


        x = np.copy(APER_RAD)
        y = np.copy(FLUX_BKG)

        for i in range(len(x)-1):
            #ideal radius is the radius aper_increase (default 30%) past the current radius
            ideal_radius = aper_increase*x[i]
            #closest_radius_ind is the index of the radius we have sampled closest to this ideal
            closest_radius_ind = np.argmin(np.abs(x - ideal_radius))
            #flux at this closest radius
            flux = y[closest_radius_ind]

            if flux <= y[i]*threshold:
                opt_rad = x[i]
                break

        if profile_plot:

            plt.figure()
            plt.plot(x, y,'ro', markersize=2.5, label="Bkg-sub flux")
#             plt.plot(x, FLUX,'bo', markersize=2.5, label="Non-bkg-sub flux")

            plt.axvline(x=opt_rad, color="black", label="Optimal Aperture")
            plt.xlabel("Aperture radius (arcsec)",  fontsize=14)
            if unit:
                plt.ylabel('Flux (Jy)', fontsize=14)
            else:
                plt.ylabel('Flux (unitless)', fontsize=14)
            plt.legend()
            plt.show()

        return opt_rad
    
    def plot_aps(self, coord_list, aperture_list):
        
        '''
        Plots a circular aperture for each source ontop of the image.
        
        Parameters
        ----------
        image: fits file
            Image to plot
        coord_list: list of `astropy.coordinates.SkyCoord`
            List of coordinates on which each aperture is centered
        aperture_list: list of floats
            List of aperture radii to plot
        
        Returns
        -------
        Plots apertures over the image.
        '''
        
        data,header = self.data 
        
        wcs_header = WCS(self.image.header).celestial
        plt.figure(figsize=(6,6))
        plt.subplot(projection=wcs_header)
        
        if 'CD1_1' in header:
            pixel_scale = np.absolute(header['CD1_1'])*3600.0
        elif 'CDELT1' in header:
            pixel_scale = np.absolute(header['CDELT1'])*3600.0
        else:
            raise Exception('Neither CD1_1 nor CDELT1 were found in the header')

        
        x_source_cent, y_source_cent = wcs_header.world_to_pixel(coord_list[0])
        
        #plot image
        primary_ap_rad_pixel = np.max(aperture_list)/pixel_scale 
        data_for_norm = data[int(y_source_cent-5.0*primary_ap_rad_pixel):int(y_source_cent+5.0*primary_ap_rad_pixel),
                                 int(x_source_cent-5.0*primary_ap_rad_pixel):int(x_source_cent+5.0*primary_ap_rad_pixel)]
        norm = simple_norm(data_for_norm[data_for_norm>0], stretch='log', percent=99.5)
        #tr = scipy.ndimage.rotate(data, 45)
        plt.imshow(data, cmap='inferno', origin='lower', norm=norm)
        
        #set x and y limits and labels
        plt.xlim(x_source_cent-5.0*primary_ap_rad_pixel, x_source_cent+5.0*primary_ap_rad_pixel)
        plt.ylim(y_source_cent-5.0*primary_ap_rad_pixel, y_source_cent+5.0*primary_ap_rad_pixel)
        plt.xlabel('RA (J2000)')
        plt.ylabel('Dec (J2000)')
        
        #plot apertures
        for i in range(len(aperture_list)):
            x_source, y_source = wcs_header.world_to_pixel(coord_list[i])
            plt.plot(x_source, y_source,'kx')
            #defines the aperture size in pixels
            aper_rad_pixel = aperture_list[i]/pixel_scale
            aperture = CircularAperture([[x_source,y_source]], r=aper_rad_pixel)
            aperture.plot(color='white', linewidth=1.5)
            
        colorbar=True
        if colorbar:
            cbar_ticks = np.around(np.linspace(norm.vmin,norm.vmax,num=5))
            if 'BUNIT' in header:
                cbar = plt.colorbar(label='PixelUnits: {0}'.format(header['BUNIT']), pad=0.01)
            elif 'FUNITS' in header:
                cbar = plt.colorbar(label='PixelUnits: {0}'.format(header['FUNITS']), pad=0.01)
            elif 'COMMENT' in header and len(header['COMMENT'])>=17: #work around to print date in WISE data
                cbar = plt.colorbar(label='{0}'.format(header['COMMENT'][17]), pad=0.01)
            else:
                cbar = plt.colorbar(label='PixelUnits: check header', pad=0.01)
            cbar.set_ticks(cbar_ticks)
            cbar.set_ticklabels(cbar_ticks)

        
class FitterContainer():
    '''
    A class to store the results from the SedFluxer class
    '''
    
    def __init__(self,models,master_dir="./",extinction_law="kmh",
                 lambda_array=None,flux_array=None,err_flux_array=None,upper_limit_array=None,dist=None):
        self.models_array = models
        self.master_dir = master_dir
        self.extc_law = extinction_law
        self.lambda_array = lambda_array
        self.flux_array = flux_array
        self.err_flux_array = err_flux_array
        self.upper_limit_array = upper_limit_array
        self.dist = dist
        self.__best_model = None

    #TODO: DELETE all these functions as now everything goes through get_model_info()
    #I do not think this is necessary anymore
    def get_available_mc(self):
        return set(self.models_array[:,0])

    def get_available_sigma(self):
        return set(self.models_array[:,1])

    def get_available_ms(self):
        return set(self.models_array[:,2])

    def get_index_from_model(self,mc,sigma,mstar,theta_view):
        return mc,sigma,mstar,theta_view

    #TODO: Change mu variable by theta or theta_view
    def get_model_info(self,keys=['mcore','sigma','mstar','theta_view'],tablename=None):
        '''
        Gets the model information from the results of the SED fit. Automatically sorts the results by chisq.
        Columns units are given as `astropy.units`
                 
        Parameters
        ----------
        keys: str or array of str
            keys to be used to select the unique final table. It should be use either
            ['mcore','sigma','mstar'] for the 432 physical models considering the best (by chisq) viewing angle or
            ['mcore','sigma','mstar','theta_view'] for the 8640 models including the viewing angle
            Default is ['mcore','sigma','mstar','theta_view']

        tablename: str, optional
            A path, or a Python file-like object. Note that fname is used verbatim,
            and there is no attempt to make the extension. Default is None.
            Note that one can choose the format of the table by changing the extension.
                    
        Returns
        ----------
        table: `astropy.table` with all the model information based on the given keys
        '''
                        
        models = self.models_array
        master_dir = self.master_dir
        norm_etxc_law = self.extc_law
        dist = self.dist
        
        nmu=20
        mu_arr=np.arange(float(nmu))/float(nmu)+1.0/float(nmu)/2.0
        mu_arr=mu_arr[::-1] #reversing the array
        theta_arr=np.arccos(mu_arr)/np.pi*180.0
       
        #filtering the table before getting into the physical model information and luminosity calculations
        columns_names = ['mcore','sigma','mstar','theta_view','av','chisq','chisq_nonlim']
        
        models_table = Table(data=models,names = columns_names) #create astropy table
        models_table.sort('chisq') #sort by chisq
        models_table_filter = unique(models_table,keys=keys) #filter by keys
    
        model_info = []
        #TODO: get rid of the try except. Also the function should accept simply the idx of the models.
        #Temporary fix to accept just one line table or multiple lines table
        try:
            mc_idx,sigma_idx,mstar_idx,mu_idx,av,chisq,chisq_nonlim = models_table_filter['mcore','sigma','mstar','theta_view','av','chisq','chisq_nonlim']
            model_dat = np.loadtxt(master_dir+'/Model_SEDs/model_info/{0:0=2d}_{1:0=2d}_{2:0=2d}'.format(int(mc_idx),int(sigma_idx),
                                                                                                         int(mstar_idx))+'.dat',
                                   unpack=True,skiprows=1)
            
            sed = np.loadtxt(master_dir+'/Model_SEDs/sed/{0:0=2d}_{1:0=2d}_{2:0=2d}_{3:0=2d}'.format(int(mc_idx),int(sigma_idx),
                                                                                                     int(mstar_idx),int(mu_idx))+'.dat',unpack=True)
            lambda_model = sed[0] #micron
            nu_model=c_micron_s/lambda_model #s-1
            flux_model = sed[1] #Lsun
            lbol_iso = np.trapz(y=flux_model/nu_model,x=nu_model) #bolometric luminosity of the unextincted SED
            flux_model_extincted = flux_model*10.0**(-0.4*av*norm_etxc_law) #Lsun extincted
            lbol_av = np.trapz(y=flux_model_extincted/nu_model,x=nu_model) #bolometric luminosity of the extincted SED
            
            mcore,sigma,mstar,rcore,massenv,theta_w,rstar,lstar,tstar,mdisk,rdisk,mdotd,lbol,tnow = model_dat
            model_info.append(['{0:0=2d}_{1:0=2d}_{2:0=2d}'.format(int(mc_idx),int(sigma_idx),int(mstar_idx)),
                               chisq,chisq_nonlim,mcore,sigma,mstar,theta_arr[int(mu_idx)-1],
                               dist,av,rcore,massenv,theta_w,rstar,lstar,tstar,mdisk,rdisk,
                               mdotd,lbol,lbol_iso,lbol_av,tnow])
        except:
            for ind_mod in models_table_filter:
                mc_idx,sigma_idx,mstar_idx,mu_idx,av,chisq,chisq_nonlim = ind_mod['mcore','sigma','mstar','theta_view','av','chisq','chisq_nonlim']
                model_dat = np.loadtxt(master_dir+'/Model_SEDs/model_info/{0:0=2d}_{1:0=2d}_{2:0=2d}'.format(int(mc_idx),int(sigma_idx),
                                                                                                             int(mstar_idx))+'.dat',
                           unpack=True,skiprows=1)

                sed = np.loadtxt(master_dir+'/Model_SEDs/sed/{0:0=2d}_{1:0=2d}_{2:0=2d}_{3:0=2d}'.format(int(mc_idx),int(sigma_idx),
                                                                                                         int(mstar_idx),int(mu_idx))+'.dat',unpack=True)
                lambda_model = sed[0] #micron
                nu_model=c_micron_s/lambda_model #s-1
                flux_model = sed[1] #Lsun
                lbol_iso = np.trapz(y=flux_model/nu_model,x=nu_model) #bolometric luminosity of the unextincted SED
                flux_model_extincted = flux_model*10.0**(-0.4*av*norm_etxc_law) #Lsun extincted
                lbol_av = np.trapz(y=flux_model_extincted/nu_model,x=nu_model) #bolometric luminosity of the extincted SED
                
                mcore,sigma,mstar,rcore,massenv,theta_w,rstar,lstar,tstar,mdisk,rdisk,mdotd,lbol,tnow = model_dat
                model_info.append(['{0:0=2d}_{1:0=2d}_{2:0=2d}_{3:0=2d}'.format(int(mc_idx),int(sigma_idx),int(mstar_idx),int(mu_idx)),
                                   chisq,chisq_nonlim,mcore,sigma,mstar,theta_arr[int(mu_idx)-1],
                                   dist,av,rcore,massenv,theta_w,rstar,lstar,tstar,mdisk,rdisk,
                                   mdotd,lbol,lbol_iso,lbol_av,tnow])
                
        model_info = np.array(model_info,dtype=object)
        
        columns_names = ['SED_number','chisq','chisq_nonlim',
                         'mcore','sigma','mstar','theta_view',
                         'dist','av','rcore','massenv','theta_w_esc',
                         'rstar','lstar','tstar','mdisk','rdisk',
                         'mdotd','lbol','lbol_iso','lbol_av','t_now']

        units = [u.dimensionless_unscaled,u.dimensionless_unscaled,u.dimensionless_unscaled,
                 u.M_sun,u.g*u.cm**-2,u.M_sun,u.deg,
                 u.pc,u.mag,u.pc,u.M_sun,u.deg,
                 u.R_sun,u.L_sun,u.K,u.M_sun,u.au,
                 u.M_sun/u.yr,u.L_sun,u.L_sun,u.L_sun,u.yr]
        
        dtype = [str,float,float,
                 float,float,float,float,
                 float,float,float,float,float,
                 float,float,float,float,float,
                 float,float,float,float,float]
        
        table_model_info = Table(data=[model_info[:,0],model_info[:,1],model_info[:,2],
                                       model_info[:,3],model_info[:,4],model_info[:,5],
                                       model_info[:,6],model_info[:,7],model_info[:,8],
                                       model_info[:,9],model_info[:,10],model_info[:,11],
                                       model_info[:,12],model_info[:,13],model_info[:,14],
                                       model_info[:,15],model_info[:,16],model_info[:,17],
                                       model_info[:,18],model_info[:,19],model_info[:,20],
                                       model_info[:,21]],
                                 names = columns_names, units = units, dtype= dtype)
        
        #sort the table by chisq values and filter by unique values
        table_model_info.sort('chisq')
        table_model_unique = unique(table_model_info,keys=keys)
        table_model_unique.sort('chisq') #the previous step sort them by they keys
                
        #For format consistency
        table_model_unique['chisq'].info.format = '%.5f'
        table_model_unique['chisq_nonlim'].info.format = '%.5f'
        table_model_unique['theta_view'].info.format = '%.5f'
        table_model_unique['av'].info.format = '%.5f'
        table_model_unique['mdotd'].info.format = '%.5e'
        table_model_unique['lbol'].info.format = '%.5e'
        table_model_unique['lbol_iso'].info.format = '%.5e'
        table_model_unique['lbol_av'].info.format = '%.5e'
        table_model_unique['t_now'].info.format = '%.5e'
        
        if tablename is not None:
            ascii.write(table_model_unique,tablename)
            print('Table saved in ',tablename)
        
        return(table_model_unique)


    @property
    def best_model(self):
        if self.__best_model is None:
            self.__best_model = self.get_best_model()
        return self.__best_model
    
    def get_best_model(self):
        FULL_MODEL = self.models_array
        
        nmc=15
        mc_arr=np.array([10.0,20.0,30.0,40.0,50.0,60.0,80.0,100.0,120.0,160.0,200.0,240.0,320.0,400.0,480.0])

        nsigma=4
        sigma_arr=np.array([0.1,0.316,1.0,3.16])

        nms=14
        ms_arr=np.array([0.5,1.0,2.0,4.0,8.0,12.0,16.0,24.0,32.0,48.0,64.0,96.0,128.0,160.0])

        nmu=20
        mu_arr=np.arange(float(nmu))/float(nmu)+1.0/float(nmu)/2.0
        mu_arr=mu_arr[::-1] #reversing the array
        theta_arr=np.arccos(mu_arr)/np.pi*180.0
        
        self.best_mc_idx = int(FULL_MODEL[FULL_MODEL[:,5]==np.min(FULL_MODEL[:,5])][0][0])
        self.best_sigma_idx = int(FULL_MODEL[FULL_MODEL[:,5]==np.min(FULL_MODEL[:,5])][0][1])
        self.best_ms_idx = int(FULL_MODEL[FULL_MODEL[:,5]==np.min(FULL_MODEL[:,5])][0][2])
        self.best_theta_idx = int(FULL_MODEL[FULL_MODEL[:,5]==np.min(FULL_MODEL[:,5])][0][3])

        self.best_AV = FULL_MODEL[FULL_MODEL[:,5]==np.min(FULL_MODEL[:,5])][0][4]
        self.best_chisq = FULL_MODEL[FULL_MODEL[:,5]==np.min(FULL_MODEL[:,5])][0][5]

        print('best Mc',mc_arr[self.best_mc_idx-1])
        print('best sigma',sigma_arr[self.best_sigma_idx-1])
        print('best m*',ms_arr[self.best_ms_idx-1])
        print('best theta_view', theta_arr[self.best_theta_idx-1])

        print('best AV', self.best_AV)
        print('best chisq', self.best_chisq)
        
        return(self.best_mc_idx,self.best_sigma_idx,self.best_ms_idx,self.best_theta_idx,self.best_AV)


class SedFitter(object):
    '''
    A class used to fit the SED model grid
    '''
    def __init__(self, extc_law=None,
                 lambda_array=None,
                 flux_array=None,
                 err_flux_array=None,
                 upper_limit_array=None,
                 filter_array=None):


        self.lambda_array = lambda_array
        self.flux_array = flux_array
        self.err_flux_array = err_flux_array
        self.upper_limit_array = upper_limit_array
        self.filter_array = filter_array

        self.master_dir = self.get_master_dir()
        self.extc_law = self.get_norm_extinction_law(extc_law)
        self.model_data = self.get_model_data()
        self.default_filters = self.__get_default_filters()
        self.__print_default_filters = None
    
    def get_master_dir(self):
        master_dir = pkg_resources.resource_filename("sedcreator","/")
        return(master_dir)
        
    def get_model_data(self):
        master_dir = self.master_dir
        ALL_model_dat = sorted(os.listdir(master_dir+'/Model_SEDs/sed/'))
        
        
        ALL_model_idx = []
        for i in ALL_model_dat:
            ALL_model_idx.append([int(i[0:11][0:2]),int(i[0:11][3:5]),
                                  int(i[0:11][6:8]),int(i[0:11][9:11])])
        return(ALL_model_dat,ALL_model_idx)

    @property
    def print_default_filters(self):
        if self.__print_default_filters is None:
            self.__print_default_filters = self.__get_default_filters().pprint_all()
        return self.__print_default_filters
        
    def __get_default_filters(self):
        master_dir = self.master_dir
        default_filt = ascii.read(master_dir+'/Model_SEDs/parfiles/filter_default.dat')
        filter_name = default_filt['filter']
        filter_wavelength = default_filt['wavelength']
        instrument = default_filt['instrument']
        default_filt_table = Table(data=[filter_name,filter_wavelength,instrument])
        return(default_filt_table)
                
    def get_extinction_laws_available(self):
        '''
        print all the extinction laws available
        '''
        return 0
    
    def get_norm_extinction_law(self,extc_law='kmh'):
        '''
        Computes the normalised extinction law based on the extc_law file
        '''
        
        # available_ext = Utilities.get_available_ext_law()
        # if extc_law not in available_ext:
        #    raise ValueError("%s not available. \n  Use this %s"%(extc_law,available_ext) )
        if extc_law is None:
            return 0
        
        master_dir = self.master_dir
        
        #load a model to retrieve the lambda
        #this is to interpolate the extinction law in the same lambda as the model
        sed_test = np.loadtxt(master_dir+'/Model_SEDs/sed/01_01_01_01.dat',unpack=True)#load the first SED model
        lambda_model = sed_test[0] #micron

        #load the opacities and interpolate in the range of model wavelength
        extc_file = np.loadtxt(master_dir+'/Model_SEDs/parfiles/'+extc_law+'.par',unpack=True)

        klam = extc_file[0]
        kkap = extc_file[3]

        interp_kkap = np.interp(lambda_model, klam, kkap)#in the whole vector
        interp_kV = np.interp(0.55, klam, kkap)#in the visible
        
        norm_extc_law = interp_kkap/interp_kV
        return(norm_extc_law)
    
    def sed_convolution(self,filter_wave_resp,lambda_array_model,sed_lambda_model,flux_model_extincted):
        '''
        Performs the convolution of the extincted SED with the filter response
        '''
        
        flux_model_extincted_CONV = []
        for filt_conv,filt_wave in zip(filter_wave_resp,lambda_array_model):
            fwave = filt_conv[0]
            fresponse = filt_conv[1]
            dfnu=filt_conv[2]

            I_arr1=np.interp(fwave[::-1],sed_lambda_model[::-1],flux_model_extincted[::-1])
            I_arr1=I_arr1/c_micron_s*fwave[::-1]
            I_filt=np.sum(I_arr1*dfnu[::-1]*fresponse[::-1])
            I_filt=I_filt*c_micron_s/filt_wave
            flux_model_extincted_CONV.append(I_filt)
        flux_model_extincted_CONV = np.array(flux_model_extincted_CONV)
        
        return(flux_model_extincted_CONV)

    
    def sed_extinction(self,flux_model_log,av):
        '''
        Extincts the flux_model given av
        '''
        
        norm_extc_law = self.extc_law
        FILTER_wave_resp = self.FILTER_wave_resp
        lambda_array_model = self.lambda_array_model
        sed_lambda_model = self.sed_lambda_model

        flux_model_log_extc = flux_model_log-0.4*av*norm_extc_law
        flux_model_log_extc_conv = self.sed_convolution(FILTER_wave_resp,
                                                        lambda_array_model,
                                                        sed_lambda_model,
                                                        flux_model_log_extc)
        
        return(flux_model_log_extc_conv)
    

    def chisq(self,flux_model,av):
        '''
        chisq function as defined in Eq 4 of Z&T 2018.
        This equation is used to compute the chisq stored in the tables.
        '''

        flux_fit_log_arr=np.log10(self.flux_array,dtype=np.float64)
        errup_fit_log_arr=np.log10(1.+self.err_flux_array,dtype=np.float64)# this is absolute error in log space
        errlo_fit_log_arr=-np.log10(1.-self.err_flux_array,dtype=np.float64)# note 100% error means a infinite lower error in log

        #these two lines are to avoid singularities in the case of error=0
        errup_fit_log_arr[errup_fit_log_arr==0.0]=1.0e-33
        errlo_fit_log_arr[errlo_fit_log_arr==0.0]=1.0e-33
        
        #this line is to avoid nan in log10(negative_value), we set it to very high number
        errlo_fit_log_arr[np.isnan(errlo_fit_log_arr)]=1.0e33
        
        nfit = len(self.upper_limit_array) #total number of points
#         nfit_nonlimit = len(self.upper_limit_array[self.upper_limit_array==0]) #total number of points that are NON upper limits
        errup_fit_log_arr[self.upper_limit_array] = 1.0e33 #set upper limit errors to very high value
        errlo_fit_log_arr[self.upper_limit_array] = 1.0e33 #set lower limit errors to very high value


        flux_model_ext_conv = self.sed_extinction(flux_model,av) #extincted and then convoluted

        up_err_for_chisq = flux_model_ext_conv>=flux_fit_log_arr
        errup_fit_log_arr[up_err_for_chisq] = np.log10(1.0+self.err_flux_array[up_err_for_chisq],dtype=np.float64)
        chisq_up = np.sum((flux_model_ext_conv[up_err_for_chisq]-flux_fit_log_arr[up_err_for_chisq])**2.0/errup_fit_log_arr[up_err_for_chisq]**2.0)

        lo_err_for_chisq = flux_model_ext_conv<flux_fit_log_arr
        chisq_lo = np.sum((flux_model_ext_conv[lo_err_for_chisq]-flux_fit_log_arr[lo_err_for_chisq])**2.0/errlo_fit_log_arr[lo_err_for_chisq]**2.0)
        
        chisq = (chisq_up+chisq_lo)/float(nfit)

        nfit_nonlimit = len(errup_fit_log_arr[errup_fit_log_arr<1.0e30]) #updating the number of points considered to be upper limits
        chisq_nonlimit = chisq*float(nfit)/float(nfit_nonlimit)

        return([chisq,chisq_nonlimit])

    #chisq function to use with minimise (note change of inputs w.r.t. chisq)
    def chisq_to_minimize(self,av,flux_model):
        '''
        chisq function to be used with scipy.optmize.minimize function.
        Note the slight change in sysntax from chisq() function above.
        It short it is exactly the same, but the parameters are called in a different order.
        '''

        flux_fit_log_arr=np.log10(self.flux_array,dtype=np.float64)
        errup_fit_log_arr=np.log10(1.+self.err_flux_array,dtype=np.float64)# this is absolute error in log space
        errlo_fit_log_arr=-np.log10(1.-self.err_flux_array,dtype=np.float64)# note 100% error means a infinite lower error in log
        
        #these two lines are to avoid singularities in the case of error=0
        errup_fit_log_arr[errup_fit_log_arr==0.0]=1.0e-33
        errlo_fit_log_arr[errlo_fit_log_arr==0.0]=1.0e-33
        
        #this line is to avoid nan in log10(negative_value), we set it to very high number
        errlo_fit_log_arr[np.isnan(errlo_fit_log_arr)]=1.0e33

        nfit = len(self.upper_limit_array) #total number of points
#         nfit_nonlimit = len(self.upper_limit_array[self.upper_limit_array==0]) #total number of points that are NON upper limits
        errup_fit_log_arr[self.upper_limit_array] = 1.0e33 #set upper limit errors to very high value
        errlo_fit_log_arr[self.upper_limit_array] = 1.0e33 #set lower limit errors to very high value


        flux_model_ext_conv = self.sed_extinction(flux_model,av) #extincted and then convoluted

        up_err_for_chisq = flux_model_ext_conv>=flux_fit_log_arr
        errup_fit_log_arr[up_err_for_chisq] = np.log10(1.0+self.err_flux_array[up_err_for_chisq],dtype=np.float64)
        chisq_up = np.sum((flux_model_ext_conv[up_err_for_chisq]-flux_fit_log_arr[up_err_for_chisq])**2.0/errup_fit_log_arr[up_err_for_chisq]**2.0)

        lo_err_for_chisq = flux_model_ext_conv<flux_fit_log_arr
        chisq_lo = np.sum((flux_model_ext_conv[lo_err_for_chisq]-flux_fit_log_arr[lo_err_for_chisq])**2.0/errlo_fit_log_arr[lo_err_for_chisq]**2.0)

        chisq = (chisq_up+chisq_lo)/float(nfit)

        nfit_nonlimit = len(errup_fit_log_arr[errup_fit_log_arr<1.0e30]) #updating the number of points considered to be upper limits
        chisq_nonlimit = chisq*float(nfit)/float(nfit_nonlimit)

        return(chisq)
    
    def sed_fit(self,dist,AV_max=1000,method='minimize',avopt=0):
        #TODO: write proper function description.
        '''
        Fits the SED observations to the Z&T18 set of models
        
        Parameters
        ----------
        dist: float
            distance given in pc to the object.

        AV_max: float
            Maximum visual extunction value to consider in the fit. The fit is perform in the range [0,AV_max].
            Default is 1000.0

        method: {'minimize', 'grid_search', 'idl'}
            Method to perform the fit.
            'minimize' method uses the `scipy.optimize` minimize function over chisq_to_minimize()
            to find the best AV for each 8640 model.
            'grid_search' performs a for loop over theAV_array = np.arange(0.0,AV_max+1.0,1.0)
            and calculates the chi square as define in Z&T18 (See also De Buizer et al. 2017).
            'idl' is a translation of the IDL version that also performs a grid search,
            it keeps the compatibility with the previous version using the same constants and fits files.
                        
        avopt: int
            This is only relevant when choosing method = 'idl'. It sets the visual extunction option.
            0 is equally distributed AV point in the range of (0,AV_max)
            and 1 is define taking into account the mass surface density value.
            Default is 0 and should be the one used if wants to compare results with minimize or grid search.
            
        Returns
        ----------
        FitterContainer: class with needed information to use the functions. COMPLETE!
        '''
        
        if method not in ('grid_search', 'idl', 'minimize'):
            raise ValueError("'method' must be either 'minimize', 'grid_search' or 'idl'")
        
        self.dist = dist
        
        #preparing the fluxes and errors in log space
        flux_fit_log_arr=np.log10(self.flux_array,dtype=np.float64)
        errup_fit_log_arr=np.log10(1.+self.err_flux_array,dtype=np.float64)# this is absolute error in log space
        errlo_fit_log_arr=-np.log10(1.-self.err_flux_array,dtype=np.float64)# note 100% error means a infinite lower error in log

        #these two lines are to avoid singularities in the case of error=0
        errup_fit_log_arr[errup_fit_log_arr==0.0]=1.0e-33
        errlo_fit_log_arr[errlo_fit_log_arr==0.0]=1.0e-33
        
        #this line is to avoid nan in log10(negative_value), we set it to very high number
        errlo_fit_log_arr[np.isnan(errlo_fit_log_arr)]=1.0e33
        
        nfit = len(self.upper_limit_array) #total number of points
        nfit_nonlimit = len(self.upper_limit_array[self.upper_limit_array==0]) #total number of points that are NON upper limits
        errup_fit_log_arr[self.upper_limit_array] = 1.0e33 #set upper limit errors to very high value
        errlo_fit_log_arr[self.upper_limit_array] = 1.0e33 #set lower limit errors to very high value

        #Creating the AV array, simply from 0 to AV_max in steps of 1
        AV_array = np.arange(0.0,AV_max+1.0,1.0)

        #loading here SED model files, extinction law, default parameters
        norm_extc_law = self.extc_law
        MODEL_DATA, MODEL_IDX = self.model_data
        master_dir = self.master_dir
        default_filters_table = self.default_filters
        filter_name = default_filters_table['filter']
        filter_wavelength = default_filters_table['wavelength']

        filter_idx = []
        for filter_value in self.filter_array:
            filter_idx.append(np.where(filter_name==filter_value)[0][0])

        self.lambda_array_model = filter_wavelength[filter_idx] #this is use for the convolution
        filter_array_model = filter_name[filter_idx]

        #load the filters lambda and responses to make the convolution
        FILTER_wave_resp = []
        for filter_NAME in filter_array_model:
            filter_file = np.loadtxt(master_dir+'/Model_SEDs/parfiles/'+filter_NAME+'.txt',unpack=True)
            fwave = filter_file[0]
            fresponse = filter_file[1]

            fnu=c_micron_s/fwave
            nf=len(fnu)

            if fnu[0] > fnu[1]:
                fnu=fnu[::-1]
                fresponse=fresponse[::-1]

            dfnu=fnu[1:nf-1]-fnu[0:nf-2]
            fint=np.sum(0.5*(fresponse[0:nf-2]+fresponse[1:nf-1])*dfnu)
            fresponse=fresponse[1:nf-1]/fint
            fwave=fwave[1:nf-1]

            FILTER_wave_resp.append([fwave,fresponse,dfnu])
        self.FILTER_wave_resp = FILTER_wave_resp


        FULL_MODEL = []


        if method == 'minimize':
            for model_data,model_idx in tqdm(zip(MODEL_DATA,MODEL_IDX),total=len(MODEL_DATA)):
                sed_model = np.loadtxt(master_dir+'/Model_SEDs/sed/'+model_data,unpack=True)
                self.sed_lambda_model = sed_model[0] #micron
                sed_flux_model = sed_model[1]

                fnorm=1.0/(4.0*np.pi)/self.dist**2/pc2cm/pc2cm*Lsun2erg_s
                modelflux_arr1=sed_flux_model*fnorm # now in erg/s/cm^2. It was in Lsun
                flux_model_Jy=modelflux_arr1/c_micron_s*self.sed_lambda_model/Jy2erg_s_cm2# now in Jy
                flux_model_Jy_log=np.log10(flux_model_Jy,dtype=np.float64)

                #TEMPORARY FIX TO AVOID CURVE FIT TO DIE
                flux_model_Jy_log[flux_model_Jy_log==np.inf]=1.0e33
                flux_model_Jy_log[flux_model_Jy_log==-np.inf]=1.0e-33

                #fitting the best av for each model (8640)
                result = minimize(self.chisq_to_minimize,x0=np.array([AV_max/2.0]),args=(flux_model_Jy_log),
                                  bounds=Bounds(0.0,AV_max))

                chisq,chisq_nonlimit = self.chisq(flux_model_Jy_log,result.x[0])

                FULL_MODEL.append(model_idx+[result.x[0],chisq,chisq_nonlimit])
                            
                
        elif method == 'grid_search':
            for model_data,model_idx in tqdm(zip(MODEL_DATA,MODEL_IDX),total=len(MODEL_DATA)):
                sed_model = np.loadtxt(master_dir+'/Model_SEDs/sed/'+model_data,unpack=True)
                self.sed_lambda_model = sed_model[0] #micron
                sed_flux_model = sed_model[1]

                fnorm=1.0/(4.0*np.pi)/self.dist**2/pc2cm/pc2cm*Lsun2erg_s
                modelflux_arr1=sed_flux_model*fnorm # now in erg/s/cm^2. It was in Lsun
                flux_model_Jy=modelflux_arr1/c_micron_s*self.sed_lambda_model/Jy2erg_s_cm2# now in Jy
                flux_model_Jy_log=np.log10(flux_model_Jy,dtype=np.float64)
                
                #TEMPORARY FIX TO AVOID inf values in chisq
                flux_model_Jy_log[flux_model_Jy_log==np.inf]=1.0e33
                flux_model_Jy_log[flux_model_Jy_log==-np.inf]=1.0e-33

                
                for av in AV_array:
                    #extincting the model flux using the av value from AV_array
                    #and the normalised (in Vband) etxc_law
                    flux_model_Jy_log_extincted = flux_model_Jy_log-0.4*av*norm_extc_law

                    flux_model_Jy_extincted_log_CONV = self.sed_convolution(FILTER_wave_resp,
                                                                            self.lambda_array_model,self.sed_lambda_model,
                                                                            flux_model_Jy_log_extincted)

                    #preliminar consideration of low and up errors
                    up_err_for_chisq = flux_model_Jy_extincted_log_CONV>=flux_fit_log_arr
                    errup_fit_log_arr[up_err_for_chisq] = np.log10(1.0+self.err_flux_array[up_err_for_chisq],dtype=np.float64)
                    chisq_up = np.sum((flux_model_Jy_extincted_log_CONV[up_err_for_chisq]-flux_fit_log_arr[up_err_for_chisq])**2.0/errup_fit_log_arr[up_err_for_chisq]**2.0)

                    lo_err_for_chisq = flux_model_Jy_extincted_log_CONV<flux_fit_log_arr
                    chisq_lo = np.sum((flux_model_Jy_extincted_log_CONV[lo_err_for_chisq]-flux_fit_log_arr[lo_err_for_chisq])**2.0/errlo_fit_log_arr[lo_err_for_chisq]**2.0)

                    chisq = (chisq_up+chisq_lo)/float(nfit)
                    
                    nfit_nonlimit = len(errup_fit_log_arr[errup_fit_log_arr<1.0e30])
                    chisq_nonlimit = chisq*float(nfit)/float(nfit_nonlimit)


                    FULL_MODEL.append(model_idx+[av,chisq,chisq_nonlimit])
                    
                    
        elif method == 'idl':
            #load the opacities and interpolate in the range of model wavelength
            par_test = np.loadtxt(master_dir+'/Model_SEDs/parfiles/kmh.par',unpack=True)

            klam = par_test[0]
            kkap = par_test[3]

            fits_path = master_dir+'/Model_SEDs/flux_filt/'

            nmc=15
            mc_arr=np.array([10.0,20.0,30.0,40.0,50.0,60.0,80.0,100.0,120.0,160.0,200.0,240.0,320.0,400.0,480.0])

            nsigma=4
            sigma_arr=np.array([0.1,0.316,1.0,3.16])

            nms=14
            mstar_arr=np.array([0.5,1.0,2.0,4.0,8.0,12.0,16.0,24.0,32.0,48.0,64.0,96.0,128.0,160.0])

            nmu=20
            mu_arr=np.arange(float(nmu))/float(nmu)+1.0/float(nmu)/2.0
            mu_arr=mu_arr[::-1] #reversing the array
            theta_arr=np.arccos(mu_arr)/np.pi*180.0

            #these are the constants used in the the IDL code. In the new python version more accurate ones are used
            pc=3.0857e18
            lsun=3.845e33
            mH=1.6733e-24
            clight=2.9979e14

            for mc in tqdm(range(nmc)):
                for sigma in range(nsigma):
                    if avopt:
                        av_clump = sigma/1.6733e-24/1.8e21/2.0
                        av_min=0.0
                        av_max=av_clump*5.0
                        av_arr=av_min+np.arange(AV_max)/float(AV_max-1.0)*(av_max-av_min)
                    else:
                        av_min=0.0
                        av_max = AV_max
                        av_arr=av_min+np.arange(AV_max)/float(AV_max-1)*(av_max-av_min) #like idl
                    for ms in range(nms):
                        for mu in range(nmu):
                            filt_conv = []
                            for filter_NAME in filter_array_model:
                                test_fits = pyfits.open(fits_path+filter_NAME+'.fits')#consider putting out of the loop
                                filt_conv.append(test_fits[0].data[mu][ms][sigma][mc])
                            filt_conv = np.array(filt_conv)
                            if all(filt_conv==0.0):
                                continue
                            else:
                                filt_conv_idx = filt_conv > 0.0 #create boolean array to take only values >0
                                lambda_model = self.lambda_array[filt_conv_idx] #use basically the same lambda_Array
                                flux_fit_log_arr_gt0 = flux_fit_log_arr[filt_conv_idx]
                                err_flux_array_gt0 = self.err_flux_array[filt_conv_idx]
                                errup_fit_log_arr_gt0 = errup_fit_log_arr[filt_conv_idx]
                                errlo_fit_log_arr_gt0 = errlo_fit_log_arr[filt_conv_idx]
                                fnorm=1./(4.0*np.pi)/dist**2/pc/pc*lsun
                                modelflux_arr1=filt_conv[filt_conv_idx]*fnorm # now in erg/s/cm^2. it was in lsun
                                flux_model_Jy=modelflux_arr1/3.0e14*lambda_model*1.0e23 # now in Jy
                                modelflux_fit_log=np.log10(flux_model_Jy,dtype=np.float64)

                                interp_kkap = np.interp(lambda_model, klam, kkap)#in the whole vector
                                interp_kV = np.interp(0.55, klam, kkap)#in the visible

                                for av in av_arr:
                                    log_flux_model_Jy_extincted = modelflux_fit_log-0.4*av*interp_kkap/interp_kV

                                    #preliminar consideration of low and up errors
                                    up_err_for_chisq = log_flux_model_Jy_extincted>=flux_fit_log_arr_gt0
                                    #here considers some upper limits not to be upper limit
                                    errup_fit_log_arr_gt0[up_err_for_chisq] = np.log10(1.0+err_flux_array_gt0[up_err_for_chisq],dtype=np.float64)
                                    chisq_up = np.sum((log_flux_model_Jy_extincted[up_err_for_chisq]-flux_fit_log_arr_gt0[up_err_for_chisq])**2.0/errup_fit_log_arr_gt0[up_err_for_chisq]**2.0)

                                    lo_err_for_chisq = log_flux_model_Jy_extincted<flux_fit_log_arr_gt0
                                    chisq_lo = np.sum((log_flux_model_Jy_extincted[lo_err_for_chisq]-flux_fit_log_arr_gt0[lo_err_for_chisq])**2.0/errlo_fit_log_arr_gt0[lo_err_for_chisq]**2.0)

                                    nfit = len(lambda_model)
                                    chisq = (chisq_up+chisq_lo)/float(nfit)
                                        
                                    nfit_nonlimit = len(errup_fit_log_arr[errup_fit_log_arr<1.0e30])
                                    chisq_nonlimit = chisq*float(nfit)/float(nfit_nonlimit)

                                    FULL_MODEL.append([mc+1,sigma+1,ms+1,mu+1,av,chisq,chisq_nonlimit])

            
        else:
            raise TypeError('method must be either minimize, grid_search or idl')
            
            
        FULL_MODEL = np.array(FULL_MODEL)
        return FitterContainer(FULL_MODEL,master_dir = self.master_dir,
                              extinction_law = self.extc_law,
                              lambda_array = self.lambda_array,flux_array= self.flux_array,
                              err_flux_array = self.err_flux_array,upper_limit_array = self.upper_limit_array,
                              dist = self.dist)
    
    
              
    #TODO: Consider putting utilities in a class
    #UTILITIES
    
    def add_filter(self,filter_name,instrument,lambda_array,response_array):
        '''
        Adds an user defined filter to the database.
        It adds a txt file to the database with columns lambda_array and response_array.
        It also adds the fits file with the convolved fluxes for the use of the idl method.
        
        Parameters
        ----------
        filter_name: str
                string to define the name of the filter. Following the convention LLNN,
                where L is the letter of the Telescope/Instrument and N is the number refering to the wavelength

        instrument: str
                string to specify the name of the telescope_instrument

        lambda_array: array
                Wavelength array in microns where the response of the filter is defined
                
        response_array: array
                    Response array that corresponds to instrumental response of the filter.
                    The response array does not need to be normalised.
            
        Returns
        ----------
        None
        '''
        
        master_dir = self.master_dir

        existing_filters = os.listdir(master_dir+'/Model_SEDs/parfiles/')

        if filter_name+'.txt' in existing_filters:
            print('WARNING! The filter ' + filter_name + ' already exists in the database')
        else:
            header0='Filter name '+filter_name+'\n'
            header1='Instrument '+instrument+'\n'
            header2='Central wavelength '+str(np.around(np.median(lambda_array),decimals=1))+' micron'+'\nwavelength(micron) response(arbitrary_units)'
            HEADER = header0+header1+header2
            np.savetxt(master_dir+'/Model_SEDs/parfiles/'+filter_name+'.txt',
                       list(zip(lambda_array,response_array)),header=HEADER,fmt='%.5e',delimiter=' ',newline='\n')
            with open(master_dir+'/Model_SEDs/parfiles/filter_default.dat','ab') as filter_default:
                np.savetxt(filter_default,[np.array([filter_name,np.around(np.median(lambda_array),decimals=1),instrument],dtype=str)],fmt='%s', delimiter=' ')

            print(filter_name+'.txt succesfully saved in '+master_dir+'Model_SEDs/parfiles/')
            
        #For the IDL method we need the convolved fluxes for the new filter in a fits file
        #(this is a translation from the IDL script filtflux.pro)
        #Note that this is not needed for the new method as it convolves on the fly
        existing_FITS_filters = os.listdir(master_dir+'/Model_SEDs/flux_filt/')
        
        if filter_name+'.fits' in existing_FITS_filters:
            print('WARNING! The filter file ' + filter_name + '.fits already exists in the database')
        else:
            nu_array = 3.0e14/lambda_array #keep consistency in constants
            nf = len(nu_array)
            if nu_array[0] > nu_array[1]:
                nu_array = nu_array[::-1]
                response_array = response_array[::-1]
            dfnu= nu_array[1:nf-1]-nu_array[0:nf-2]
            fint=np.sum(0.5*(response_array[0:nf-2]+response_array[1:nf-1])*dfnu)
            response_array=response_array[1:nf-1]/fint
            lambda_array=lambda_array[1:nf-1]

            nmc=15
            nsigma=4
            nms=14
            nmu=20

            flux_model_conv = np.zeros([nmu,nms,nsigma,nmc])

            sed_list = os.listdir(master_dir+'/Model_SEDs/sed/')

            for imc in range(1,nmc+1,1):
                for isigma in range(1,nsigma+1,1):
                    for ims in range(1,nms+1,1):
                        for imu in range(1,nmu+1,1):
                            sed_file = '{0:0=2d}_{1:0=2d}_{2:0=2d}_{3:0=2d}'.format(int(imc),int(isigma),
                                                                                    int(ims),int(imu))+'.dat'
                            if sed_file in sed_list:
                                sed = np.loadtxt(master_dir+'/Model_SEDs/sed/'+sed_file,unpack=True)

                                lambda_model = sed[0] #micron
                                flux_model = sed[1] #Lsun
                                I_arr1=np.interp(lambda_array[::-1],lambda_model[::-1],flux_model[::-1])
                                I_arr1=I_arr1/3.0e14*lambda_array[::-1]
                                I_filt=np.sum(I_arr1*dfnu[::-1]*response_array[::-1])
                                #TODO: Double check this for given filter
                                I_filt=I_filt*3.0e14/np.median(lambda_array)
                                flux_model_conv[imu-1][ims-1][isigma-1][imc-1] = I_filt
                            else:
                                continue

            pyfits.writeto(master_dir+'Model_SEDs/flux_filt/'+filter_name+'.fits',flux_model_conv)

    
    def add_SQUARE_filter(self,filter_name,instrument,filter_lambda,filter_width):
        '''
        Adds square filter to the database given the central filter_lambda and the filter_width.
        It adds a txt file to the database with columns lambda_array and response_array.
        It also adds the fits file with the convolved fluxes for the use of the idl method.
        
        Parameters
        ----------
        filter_name: str
                string to define the name of the filter. Following the convention LLNN,
                where L is the letter of the Telescope/Instrument and N is the number refering to the wavelength

        instrument: str
                string to specify the name of the telescope_instrument

        filter_lambda: float
                    Wavelength array in microns where the response of the filter is defined
                
        filter_width: float
                Response array that corresponds to instrumental response of the filter.
                The response array does not need to be normalised.
            
        Returns
        ----------
        None
        '''
        
        master_dir = self.master_dir

        existing_filters = os.listdir(master_dir+'/Model_SEDs/parfiles/')

        dwave=filter_width/100.0
        lambda_array=(np.arange(110)-55.0)*dwave+filter_lambda
        lambda_array = lambda_array[::-1]
        response_array=np.zeros(110)
        a=np.absolute(lambda_array-filter_lambda) < filter_width/2.0
        response_array[a]=1.0

        if filter_name+'.txt' in existing_filters:
            print('WARNING! The filter ' + filter_name + ' already exists in the database')
        else:
            header0='Filter name '+filter_name+'\n'
            header1='Instrument '+instrument+'\n'
            header2='Central wavelength '+str(np.around(np.median(lambda_array),decimals=1))+' micron'+'\nwavelength(micron) response(arbitrary_units)'
            HEADER = header0+header1+header2
            np.savetxt(master_dir+'/Model_SEDs/parfiles/'+filter_name+'.txt',
                       list(zip(lambda_array,response_array)),header=HEADER,fmt='%.5e',delimiter=' ',newline='\n')
            with open(master_dir+'/Model_SEDs/parfiles/filter_default.dat','ab') as filter_default:
                np.savetxt(filter_default,[np.array([filter_name,filter_lambda,instrument],dtype=str)],fmt='%s', delimiter=' ')

                print(filter_name+'.txt succesfully saved in '+master_dir+'Model_SEDs/parfiles/')
        
        #For the IDL method we need the convolved fluxes for the new filter in a fits file
        #(this is a translation from the IDL script filtflux.pro)
        #Note that this is not needed for the new method as it convolves on the fly
                existing_FITS_filters = os.listdir(master_dir+'/Model_SEDs/flux_filt/')
        
        if filter_name+'.fits' in existing_FITS_filters:
            print('WARNING! The filter file ' + filter_name + '.fits already exists in the database')
        else:
            nu_array = 3.0e14/lambda_array #keep consistency in constants
            nf = len(nu_array)
            if nu_array[0] > nu_array[1]:
                nu_array = nu_array[::-1]
                response_array = response_array[::-1]
            dfnu= nu_array[1:nf-1]-nu_array[0:nf-2]
            fint=np.sum(0.5*(response_array[0:nf-2]+response_array[1:nf-1])*dfnu)
            response_array=response_array[1:nf-1]/fint
            lambda_array=lambda_array[1:nf-1]

            nmc=15
            nsigma=4
            nms=14
            nmu=20

            flux_model_conv = np.zeros([nmu,nms,nsigma,nmc])

            sed_list = os.listdir(master_dir+'/Model_SEDs/sed/')

            for imc in range(1,nmc+1,1):
                for isigma in range(1,nsigma+1,1):
                    for ims in range(1,nms+1,1):
                        for imu in range(1,nmu+1,1):
                            sed_file = '{0:0=2d}_{1:0=2d}_{2:0=2d}_{3:0=2d}'.format(int(imc),int(isigma),
                                                                                    int(ims),int(imu))+'.dat'
                            if sed_file in sed_list:
                                sed = np.loadtxt(master_dir+'/Model_SEDs/sed/'+sed_file,unpack=True)

                                lambda_model = sed[0] #micron
                                flux_model = sed[1] #Lsun
                                I_arr1=np.interp(lambda_array[::-1],lambda_model[::-1],flux_model[::-1])
                                I_arr1=I_arr1/3.0e14*lambda_array[::-1]
                                I_filt=np.sum(I_arr1*dfnu[::-1]*response_array[::-1])
                                I_filt=I_filt*3.0e14/filter_lambda
                                flux_model_conv[imu-1][ims-1][isigma-1][imc-1] = I_filt
                            else:
                                continue

            pyfits.writeto(master_dir+'Model_SEDs/flux_filt/'+filter_name+'.fits',flux_model_conv)


    def plot_filter(self,filter_name,figsize=(6,4),legend=False,title=None,figname=None):
        '''
        Plots a filter in the database. To see the available filter use SedFitter().print_default_filters
        
        Parameters
        ----------
        filter_name: array of str
                array of strings with the name of the filter.
                E.g. if one filter is given ['filter1'],
                if two filters (or more) are given ['filter1','filter2']
        legend: bool, optional
                set the legend with the name of the filter. Default is False.
        figname: str, optional
                A path, or a Python file-like object. Note that fname is used verbatim,
                and there is no attempt to make the extension. Default is None.
                Note that one can choose the format of the figure by changing the extension,
                e.g., figure.pdf would generate the figure in PDF format.
            
        Returns
        ----------
        A figure with the normalised response for the given filter_name
        '''
        
        master_dir = self.master_dir

        existing_filters = os.listdir(master_dir+'/Model_SEDs/parfiles/')

        plt.figure(figsize=figsize)
        
        for filt in filter_name:
            if filt+'.txt' in existing_filters:
                lambda_array,response_array = np.loadtxt(master_dir+'/Model_SEDs/parfiles/'+filt+'.txt',unpack=True)
                plt.step(lambda_array,response_array/np.max(response_array),label=filt)
            else:
                print('WARNING! The filter ' + filt + ' is not in the database')
                return()
                
        plt.xlabel('$\lambda\,(\mu\mathrm{m})$')
        plt.ylabel('Filter Response (arbitrary units)')
        if legend:
            plt.legend()
        if title is not None:
            plt.title(title)
        if figname is not None:
            plt.savefig(figname, dpi=300, bbox_inches="tight")
            print('Image saved in ',figname)
        plt.show()


    def table2latex(self,table,keys=['chisq','mcore','sigma','rcore','mstar','theta_view','av','massenv','theta_w_esc','mdotd','lbol_iso','lbol'],tablename=None):
        '''
        Outputs the astropy table into a latex table given the keys properly
        formatted and with the correct units
        
        Parameters
        ----------
        table: `astropy.table`
            Table to be transformed into latex

        keys: str array
            Array with the keys that wants to be used.
            Default is keys=['chisq','mcore','sigma','rcore','mstar','theta_view','av','massenv','theta_w_esc','mdotd','lbol_iso','lbol']
                    
        tablename: str
            A path, or a Python file-like object. Note that fname is used verbatim,
            and there is no attempt to make the extension. Default is None.
            Note that one can choose the format of the figure by changing the extension.
            
        Returns
        ----------
        Latex formatted table in the path with the name tablename
        '''
        
        master_dir = self.master_dir
        
        col_names = np.array(['SED_number', 'chisq', 'chisq_nonlim', 'mcore', 'sigma', 'mstar', 'theta_view', 'dist', 'av', 'rcore', 'massenv', 'theta_w_esc',
                       'rstar', 'lstar', 'tstar', 'mdisk', 'rdisk', 'mdotd', 'lbol', 'lbol_iso', 'lbol_av', 't_now'])

        col_latex_names = np.array(['SED_number','$\chi^2$','$\chi^2_\mathrm{nonlimit}$','$M_\mathrm{c}$','$\Sigma_\mathrm{cl}$','$m_*$','$\\theta_\mathrm{view}$','$d$','$A_V$',
                     '$R_\mathrm{core}$','$M_\mathrm{env}$','$\\theta_\mathrm{w,esc}$','$R_*$', '$L_*$', '$T_*$', '$m_\mathrm{disk}$', '$r_\mathrm{disk}$',
                     '$\dot{M}_\mathrm{disk}$','$L_\mathrm{bol}$','$L_\mathrm{bol,iso}$','$L_\mathrm{bol,av}$','$t_\mathrm{now}$'])

        latex_names_table = Table(data=col_latex_names,names=col_names)

        formats_table = np.array(['','%.2f','%.2f','%.0f','%.3f','%.0f','%.0f','%.0f','%.2f','%.2f',
                                  '%.2f','%.0f','%.2f','%.2f','%.2f','%.2f','%.2f','%.1e','%.1e','%.1e',
                                  '%.1e','%.1e'])

        latex_formats_table = Table(data=formats_table,names=col_names)

       
        formats_dict = {}
        for name_value, format_value in zip(latex_names_table[keys].as_array()[0],
                                            latex_formats_table[keys].as_array()[0]):
            formats_dict[name_value] = format_value

        if tablename is not None:
            ascii.write(table[keys],
                        format='latex',
                        formats=formats_dict,
                        names=list(latex_names_table[keys].as_array()[0]),
                        output=tablename)
        else:
            print('Please, set tablename to a str to save the table, including the name of the table and extension, e.g., table.txt')
        
        return


    def get_average_model(self,models,number_of_models=5,chisq_cut=None,core_radius_cut=None,method=None,tablename=None):
        '''
        Get the average model by means of calculating the geometric mean
        for all parameters except for av, theta_view, and theta_w_esc from which arithmetic mean is calculated

        models must be the astropy table obtained from the SED results that has column names defined.
        
        There are two methods for the average. The first just take the number of models specified (default 5)
        and the second takes into account a chisq cut and/or radius_cut. When set the chisq_cut, it will consider
        all models below that chisq_cut value. When set the core_radius_cut, it will consider all models with
        the model output core radius is smaller than the core_radius_cutvalue.
        
        Parameters
        ----------
        models: table, `astropy.table`
            models to be averaged over. It should be an astropy table got from the get_model_info() function.

        number_of_models: int
            Number of models to be used for the average for method 1.

        chisq_cut: float
            Chisq cut value to consider all models below
                        
        core_radius_cut: float
            Core radius cut value to consider all models below

        method: str, {'liu','moser'} optional
            Defines the method for backwards compatibilty with previous works.
            liu takes the best 5 models OR FEWER with the condition chi2<chis2min+5 and Rc<2Rap
            moser takes the best 10 models OR FEWER with the condition chis2<2*chi2min.
            The chis2 and core radius conditions have to be given as usual.
            The method only takes the best 5 (or 10) OR FEWER models satisfying the condition.

        tablename: str, optional
            A path, or a Python file-like object. Note that fname is used verbatim,
            and there is no attempt to make the extension. Default is None.
            Note that one can choose the format of the figure by changing the extension.
                    
        Returns
        ----------
        table: `astropy.table` with all the model information based on the given keys
        '''
        
        master_dir = self.master_dir
        
        average_model_table = models[0:int(number_of_models)]

        columns_names = ['method','number_of_models_used',
                         'mcore','Dmcore','sigma','Dsigma','mstar','Dmstar','theta_view','Dtheta_view',
                         'dist','av','Dav','rcore','Drcore','massenv','Dmassenv','theta_w_esc','Dtheta_w_esc',
                         'rstar','Drstar','lstar','Dlstar','tstar','Dtstar','mdisk','Dmdisk','rdisk','Drdisk',
                         'mdotd','Dmdotd','lbol','Dlbol','lbol_iso','Dlbol_iso','lbol_av','Dlbol_av','t_now','Dt_now']

        units = [u.dimensionless_unscaled,u.dimensionless_unscaled,
                 u.M_sun,u.dimensionless_unscaled,u.g*u.cm**-2,u.dimensionless_unscaled,u.M_sun,u.dimensionless_unscaled,u.deg,u.deg,
                 u.pc,u.mag,u.mag,u.pc,u.dimensionless_unscaled,u.M_sun,u.dimensionless_unscaled,u.deg,u.deg,
                 u.R_sun,u.dimensionless_unscaled,u.L_sun,u.dimensionless_unscaled,u.K,u.dimensionless_unscaled,u.M_sun,u.dimensionless_unscaled,u.au,u.dimensionless_unscaled,
                 u.M_sun/u.yr,u.dimensionless_unscaled,u.L_sun,u.dimensionless_unscaled,u.L_sun,u.dimensionless_unscaled,u.L_sun,u.dimensionless_unscaled,u.yr,u.dimensionless_unscaled]

        data=np.array(['average model method 1',number_of_models,
                       gmean(average_model_table['mcore']),gstd(average_model_table['mcore']),
                       gmean(average_model_table['sigma']),gstd(average_model_table['sigma']),
                       gmean(average_model_table['mstar']),gstd(average_model_table['mstar']),
                       np.mean(average_model_table['theta_view']),np.std(average_model_table['theta_view']),
                       np.mean(average_model_table['dist']),
                       np.mean(average_model_table['av']),np.std(average_model_table['av']),
                       gmean(average_model_table['rcore']),gstd(average_model_table['rcore']),
                       gmean(average_model_table['massenv']),gstd(average_model_table['massenv']),
                       np.mean(average_model_table['theta_w_esc']),np.std(average_model_table['theta_w_esc']),
                       gmean(average_model_table['rstar']),gstd(average_model_table['rstar']),
                       gmean(average_model_table['lstar']),gstd(average_model_table['lstar']),
                       gmean(average_model_table['tstar']),gstd(average_model_table['tstar']),
                       gmean(average_model_table['mdisk']),gstd(average_model_table['mdisk']),
                       gmean(average_model_table['rdisk']),gstd(average_model_table['rdisk']),
                       gmean(average_model_table['mdotd']),gstd(average_model_table['mdotd']),
                       gmean(average_model_table['lbol']),gstd(average_model_table['lbol']),
                       gmean(average_model_table['lbol_iso']),gstd(average_model_table['lbol_iso']),
                       gmean(average_model_table['lbol_av']),gstd(average_model_table['lbol_av']),
                       gmean(average_model_table['t_now']),gstd(average_model_table['t_now'])],dtype=object)

        dtype = [str,int,float,float,float,float,float,float,float,float,float,float,float,
                 float,float,float,float,float,float,float,float,float,float,float,float,float,
                 float,float,float,float,float,float,float,float,float,float,float,float,float]

        final_average_table = Table(data = data, names = columns_names, units=units, dtype=dtype)


        if chisq_cut is not None and core_radius_cut is None:
            average_model_table_chisq = models[models['chisq']<=chisq_cut]

            number_of_models = len(average_model_table_chisq)
            if number_of_models==0:
                raise ValueError('The considered constraints in chisq_cut produce an empty table. Please consider to relax them.')
            if number_of_models==1:
                raise ValueError('The considered constraints in chisq_cut produce a table with 1 row, no mean or dispersion makes sense. Please consider to relax them.')

            if method is not None:
                if method == 'liu':
                    if len(average_model_table_chisq)>=5:
                        average_model_table_chisq = average_model_table_chisq[0:5]
                    else:
                        average_model_table_chisq = average_model_table_chisq
                elif method == 'moser':
                    if len(average_model_table_chisq)>=10:
                        average_model_table_chisq = average_model_table_chisq[0:10]
                    else:
                        average_model_table_chisq = average_model_table_chisq
                else:
                    raise ValueError('Method must be either liu, moser or None')
            
            number_of_models = len(average_model_table_chisq)
            
            data=np.array(['average model method 2',number_of_models,
                           gmean(average_model_table_chisq['mcore']),gstd(average_model_table_chisq['mcore']),
                           gmean(average_model_table_chisq['sigma']),gstd(average_model_table_chisq['sigma']),
                           gmean(average_model_table_chisq['mstar']),gstd(average_model_table_chisq['mstar']),
                           np.mean(average_model_table_chisq['theta_view']),np.std(average_model_table_chisq['theta_view']),
                           np.mean(average_model_table_chisq['dist']),
                           np.mean(average_model_table_chisq['av']),np.std(average_model_table_chisq['av']),
                           gmean(average_model_table_chisq['rcore']),gstd(average_model_table_chisq['rcore']),
                           gmean(average_model_table_chisq['massenv']),gstd(average_model_table_chisq['massenv']),
                           np.mean(average_model_table_chisq['theta_w_esc']),np.std(average_model_table_chisq['theta_w_esc']),
                           gmean(average_model_table_chisq['rstar']),gstd(average_model_table_chisq['rstar']),
                           gmean(average_model_table_chisq['lstar']),gstd(average_model_table_chisq['lstar']),
                           gmean(average_model_table_chisq['tstar']),gstd(average_model_table_chisq['tstar']),
                           gmean(average_model_table_chisq['mdisk']),gstd(average_model_table_chisq['mdisk']),
                           gmean(average_model_table_chisq['rdisk']),gstd(average_model_table_chisq['rdisk']),
                           gmean(average_model_table_chisq['mdotd']),gstd(average_model_table_chisq['mdotd']),
                           gmean(average_model_table_chisq['lbol']),gstd(average_model_table_chisq['lbol']),
                           gmean(average_model_table_chisq['lbol_iso']),gstd(average_model_table_chisq['lbol_iso']),
                           gmean(average_model_table_chisq['lbol_av']),gstd(average_model_table_chisq['lbol_av']),
                           gmean(average_model_table_chisq['t_now']),gstd(average_model_table_chisq['t_now'])],dtype=object)
            
            final_average_table.add_row(vals=data)


        elif chisq_cut is None and core_radius_cut is not None:
            distance = models['dist'][0] #pc
            core_radius_cut_au = core_radius_cut * distance #arcsec x pc = au
            core_radius_cut_pc = core_radius_cut_au*u.au.to(u.pc) #from au to pc

            average_model_table_rcore = models[models['rcore']<=core_radius_cut_pc]

            number_of_models = len(average_model_table_rcore)
            if number_of_models==0:
                raise ValueError('The considered constraints in core_radius_cut produce an empty table. Please consider to relax them.')
            if number_of_models==1:
                raise ValueError('The considered constraints in chisq_cut produce a table with 1 row, no mean or dispersion makes sense. Please consider to relax them.')

            if method is not None:
                if method == 'liu':
                    if len(average_model_table_rcore)>=5:
                        average_model_table_rcore = average_model_table_rcore[0:5]
                    else:
                        average_model_table_rcore = average_model_table_rcore
                elif method == 'moser':
                    if len(average_model_table_rcore)>=10:
                        average_model_table_rcore = average_model_table_rcore[0:10]
                    else:
                        average_model_table_rcore = average_model_table_rcore
                else:
                    raise ValueError('Method must be either liu, moser or None')
            
            number_of_models = len(average_model_table_rcore)
            
            data=np.array(['average model method 2',number_of_models,
                           gmean(average_model_table_rcore['mcore']),gstd(average_model_table_rcore['mcore']),
                           gmean(average_model_table_rcore['sigma']),gstd(average_model_table_rcore['sigma']),
                           gmean(average_model_table_rcore['mstar']),gstd(average_model_table_rcore['mstar']),
                           np.mean(average_model_table_rcore['theta_view']),np.std(average_model_table_rcore['theta_view']),
                           np.mean(average_model_table_rcore['dist']),
                           np.mean(average_model_table_rcore['av']),np.std(average_model_table_rcore['av']),
                           gmean(average_model_table_rcore['rcore']),gstd(average_model_table_rcore['rcore']),
                           gmean(average_model_table_rcore['massenv']),gstd(average_model_table_rcore['massenv']),
                           np.mean(average_model_table_rcore['theta_w_esc']),np.std(average_model_table_rcore['theta_w_esc']),
                           gmean(average_model_table_rcore['rstar']),gstd(average_model_table_rcore['rstar']),
                           gmean(average_model_table_rcore['lstar']),gstd(average_model_table_rcore['lstar']),
                           gmean(average_model_table_rcore['tstar']),gstd(average_model_table_rcore['tstar']),
                           gmean(average_model_table_rcore['mdisk']),gstd(average_model_table_rcore['mdisk']),
                           gmean(average_model_table_rcore['rdisk']),gstd(average_model_table_rcore['rdisk']),
                           gmean(average_model_table_rcore['mdotd']),gstd(average_model_table_rcore['mdotd']),
                           gmean(average_model_table_rcore['lbol']),gstd(average_model_table_rcore['lbol']),
                           gmean(average_model_table_rcore['lbol_iso']),gstd(average_model_table_rcore['lbol_iso']),
                           gmean(average_model_table_rcore['lbol_av']),gstd(average_model_table_rcore['lbol_av']),
                           gmean(average_model_table_rcore['t_now']),gstd(average_model_table_rcore['t_now'])],dtype=object)
                    
            final_average_table.add_row(vals=data)


        elif chisq_cut is not None and core_radius_cut is not None:
            distance = models['dist'][0] #pc
            core_radius_cut_au = core_radius_cut * distance #arcsec x pc = au
            core_radius_cut_pc = core_radius_cut_au*u.au.to(u.pc) #from au to pc

            average_model_table_chisq_rcore = models[(models['chisq']<=chisq_cut) & (models['rcore']<=core_radius_cut_pc)]

            number_of_models = len(average_model_table_chisq_rcore)
            if number_of_models==0:
                raise ValueError('The considered constraints either in chisq_cut or core_radius_cut produce an empty table. Please consider to relax them.')
            if number_of_models==1:
                raise ValueError('The considered constraints in chisq_cut or core_radius_cut produce a table with 1 row, no mean or dispersion makes sense. Please consider to relax them.')

            if method is not None:
                if method == 'liu':
                    if len(average_model_table_chisq_rcore)>=5:
                        average_model_table_chisq_rcore = average_model_table_chisq_rcore[0:5]
                    else:
                        average_model_table_chisq_rcore = average_model_table_chisq_rcore
                elif method == 'moser':
                    if len(average_model_table_chisq_rcore)>=10:
                        average_model_table_chisq_rcore = average_model_table_chisq_rcore[0:10]
                    else:
                        average_model_table_chisq_rcore = average_model_table_chisq_rcore
                else:
                    raise ValueError('Method must be either liu, moser or None')
            
            number_of_models = len(average_model_table_chisq_rcore)
            
            data=np.array(['average model method 2',number_of_models,
                           gmean(average_model_table_chisq_rcore['mcore']),gstd(average_model_table_chisq_rcore['mcore']),
                           gmean(average_model_table_chisq_rcore['sigma']),gstd(average_model_table_chisq_rcore['sigma']),
                           gmean(average_model_table_chisq_rcore['mstar']),gstd(average_model_table_chisq_rcore['mstar']),
                           np.mean(average_model_table_chisq_rcore['theta_view']),np.std(average_model_table_chisq_rcore['theta_view']),
                           np.mean(average_model_table_chisq_rcore['dist']),
                           np.mean(average_model_table_chisq_rcore['av']),np.std(average_model_table_chisq_rcore['av']),
                           gmean(average_model_table_chisq_rcore['rcore']),gstd(average_model_table_chisq_rcore['rcore']),
                           gmean(average_model_table_chisq_rcore['massenv']),gstd(average_model_table_chisq_rcore['massenv']),
                           np.mean(average_model_table_chisq_rcore['theta_w_esc']),np.std(average_model_table_chisq_rcore['theta_w_esc']),
                           gmean(average_model_table_chisq_rcore['rstar']),gstd(average_model_table_chisq_rcore['rstar']),
                           gmean(average_model_table_chisq_rcore['lstar']),gstd(average_model_table_chisq_rcore['lstar']),
                           gmean(average_model_table_chisq_rcore['tstar']),gstd(average_model_table_chisq_rcore['tstar']),
                           gmean(average_model_table_chisq_rcore['mdisk']),gstd(average_model_table_chisq_rcore['mdisk']),
                           gmean(average_model_table_chisq_rcore['rdisk']),gstd(average_model_table_chisq_rcore['rdisk']),
                           gmean(average_model_table_chisq_rcore['mdotd']),gstd(average_model_table_chisq_rcore['mdotd']),
                           gmean(average_model_table_chisq_rcore['lbol']),gstd(average_model_table_chisq_rcore['lbol']),
                           gmean(average_model_table_chisq_rcore['lbol_iso']),gstd(average_model_table_chisq_rcore['lbol_iso']),
                           gmean(average_model_table_chisq_rcore['lbol_av']),gstd(average_model_table_chisq_rcore['lbol_av']),
                           gmean(average_model_table_chisq_rcore['t_now']),gstd(average_model_table_chisq_rcore['t_now'])],dtype=object)
                    
            final_average_table.add_row(vals=data)

        #For format consistency
        final_average_table['mcore'].info.format = '%.5f'
        final_average_table['Dmcore'].info.format = '%.5f'
        final_average_table['sigma'].info.format = '%.5f'
        final_average_table['Dsigma'].info.format = '%.5f'
        final_average_table['mstar'].info.format = '%.5f'
        final_average_table['Dmstar'].info.format = '%.5f'
        final_average_table['theta_view'].info.format = '%.5f'
        final_average_table['Dtheta_view'].info.format = '%.5f'
        final_average_table['av'].info.format = '%.5f'
        final_average_table['Dav'].info.format = '%.5f'
        final_average_table['rcore'].info.format = '%.5f'
        final_average_table['Drcore'].info.format = '%.5f'
        final_average_table['massenv'].info.format = '%.5f'
        final_average_table['Dmassenv'].info.format = '%.5f'
        final_average_table['theta_w_esc'].info.format = '%.5f'
        final_average_table['Dtheta_w_esc'].info.format = '%.5f'
        final_average_table['theta_w_esc'].info.format = '%.5f'
        final_average_table['Dtheta_w_esc'].info.format = '%.5f'
        final_average_table['rstar'].info.format = '%.5f'
        final_average_table['Drstar'].info.format = '%.5f'
        final_average_table['lstar'].info.format = '%.5e'
        final_average_table['Dlstar'].info.format = '%.5f'
        final_average_table['tstar'].info.format = '%.5e'
        final_average_table['Dtstar'].info.format = '%.5f'
        final_average_table['mdisk'].info.format = '%.5f'
        final_average_table['Dmdisk'].info.format = '%.5f'
        final_average_table['rdisk'].info.format = '%.5e'
        final_average_table['Drdisk'].info.format = '%.5f'
        final_average_table['mdotd'].info.format = '%.5e'
        final_average_table['Dmdotd'].info.format = '%.5f'
        final_average_table['theta_w_esc'].info.format = '%.5e'
        final_average_table['Dtheta_w_esc'].info.format = '%.5f'
        final_average_table['lbol'].info.format = '%.5e'
        final_average_table['Dlbol'].info.format = '%.5f'
        final_average_table['lbol_iso'].info.format = '%.5e'
        final_average_table['Dlbol_iso'].info.format = '%.5f'
        final_average_table['lbol_av'].info.format = '%.5e'
        final_average_table['Dlbol_av'].info.format = '%.5f'
        final_average_table['t_now'].info.format = '%.5e'
        final_average_table['Dt_now'].info.format = '%.5f'

        
        if tablename is not None:
            ascii.write(final_average_table,tablename)
            print('Table saved in ',tablename)
            
        return(final_average_table)

    #Indepent Plots
    #TODO: Put all the plots in this class
    
class PentagonPlot(object):
    '''
    A class used to plot a pentagon plot
    
    Attributes
    ----------
    plot: method
        Creates the pentagon plot
    '''
        
    # modified idea taken from
    # https://stackoverflow.com/questions/33028843/how-to-remove-polar-gridlines-and-add-major-axis-ticks
    def __init__(self, fig, titles, labels, rect=None):
        if rect is None:
            rect = [0.05, 0.05, 0.95, 0.95]

        self.n = len(titles)
        self.angles = [a if a <=360. else a - 360. for a in np.arange(90, 90+360, 360.0/self.n)]
        self.axes = [fig.add_axes(rect, projection="polar", label="axes%d" % i)
                         for i in range(self.n)]

        self.ax = self.axes[0]
        self.ax.set_thetagrids(self.angles, labels=titles, fontsize=10, color="black")

        for ax in self.axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)
            self.ax.yaxis.grid(False)

        for ax, angle, label in zip(self.axes, self.angles, labels):
            ax.set_rgrids(np.round(np.linspace(0.0,1.0,5),3), labels=label, angle=angle, fontsize=10)
            ax.spines["polar"].set_visible(False)
            ax.set_ylim(0, 1.0)
            ax.xaxis.grid(True,color='black',linestyle='-')
            pos=ax.get_rlabel_position()
            ax.set_rlabel_position(pos+7)

    def __prepare_plot(self, values, *args, **kw):
        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        self.ax.plot(angle, values, *args, **kw)
        self.ax.fill(angle, values, alpha=0.1)
        
    #TODO: find out why when putting self in pentagon_plot it does not work
    def plot(models,figname=None):

        
        models_phm = models['mcore','sigma','mstar','theta_view','av']
        models_chisq = models['chisq']

        fig = plt.figure(figsize=(5, 5))

        titles = [r'$M_\mathrm{c}$ ($M_\odot$)',
                  r'$\Sigma$ ($\mathrm{g\,cm^{-2}}$)        ',
                  r'$m_*$ ($M_\odot$)',
                  r'      $\theta_\mathrm{view}$ ($\mathrm{deg}$)',
                  r'        $A_V$ (mag)']


        labels = [np.round(np.linspace(0.0,1.0,5)*models_phm['mcore'].max(),3),
                  np.round(np.linspace(0.0,1.0,5)*models_phm['sigma'].max(),3),
                  np.round(np.linspace(0.0,1.0,5)*models_phm['mstar'].max(),3),
                  np.round(np.linspace(0.0,1.0,5)*models_phm['theta_view'].max(),2),
                  np.round(np.linspace(0.0,1.0,5)*models_phm['av'].max(),2)]

        norm_models = np.array([models_phm['mcore'].max(),models_phm['sigma'].max(),models_phm['mstar'].max(),
                                models_phm['theta_view'].max(),models_phm['av'].max()])

        penta_plot = PentagonPlot(fig, titles, labels)
        #TODO: get rid of this try except.
        #Temproary fix to avoid the function braking when only 1 row is provided
        try:
            for i in range(len(models)):
                penta_plot.__prepare_plot(np.array(list(models_phm[i]))/norm_models, "-", lw=1, alpha=.5, label=r'$\chi^2={0:.2f}$'.format(models_chisq[i]))
        except:
            penta_plot.__prepare_plot(np.array(list(models_phm))/norm_models, "-", lw=1, alpha=.5, label=r'$\chi^2={0:.2f}$'.format(models_chisq))
            

        penta_plot.ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10),
                             fancybox=True, shadow=True, ncol=5)

        fig = plt.gcf()
        fig.set_size_inches(5, 5, forward=True)
        if figname is not None:
            plt.savefig(figname, dpi=300, bbox_inches="tight", pad_inches=1)
            

class ModelPlotter(FitterContainer):
    '''
    A class used to plot the SED model results
    '''
        
    def __init__(self,fittercontainer):
        self.master_dir = fittercontainer.master_dir
        self.extc_law = fittercontainer.extc_law
        self.lambda_array = fittercontainer.lambda_array
        self.flux_array = fittercontainer.flux_array
        self.err_flux_array = fittercontainer.err_flux_array
        self.upper_limit_array = fittercontainer.upper_limit_array
        self.dist = fittercontainer.dist
        self.__best_model = None


    def plot_best_sed(self,models=None,figsize=(6,4),xlim=[1.0,1000.0],ylim=[1.0e-12,1.0e-5],marker='k*',title=None,figname=None):
        '''
        Plots the best SED model
        
        Parameters
        ----------
        models: `astropy.table`
            Table with the sed results.

        figsize: tuple
            specify the size of the figure. Default is (6,4)

        xlim: array
            array to specify the limits in the x axis. Default is [1.0,1000.0]

        ylim: array
            array to specify the limits in the y axis. Default is [1.0e-12,1.0e-5]
            
        marker: str
            defines the marker style and color like matplotlib. Default is 'k*'

        title: str, optional
            Defines the title of the plot. Default is None.

        figname: str, optional
                    A path, or a Python file-like object. Note that fname is used verbatim,
                    and there is no attempt to make the extension. Default is None.
                    Note that one can choose the format of the figure by changing the extension,
                    e.g., figure.pdf would generate the figure in PDF format.
            
        Returns
        ----------
        plot of the best SED
        '''
        
        master_dir = self.master_dir
        norm_extc_law = self.extc_law
        best_model = models[models['chisq']==models['chisq'].min()]
        
        plt.figure(figsize=figsize)
        sed_best = np.loadtxt(master_dir+'/Model_SEDs/sed/'+best_model['SED_number'].data[0]+'.dat',unpack=True)
        lambda_model = sed_best[0] #micron
        flux_model = sed_best[1]*Lsun2erg_s/(4.0*np.pi*(pc2cm*self.dist)**2.0) #from Lsun to erg s-1 cm-2
        flux_model_extincted = flux_model*10.0**(-0.4*best_model['av'].data[0]*norm_extc_law)
        plt.plot(lambda_model,flux_model_extincted,'k-',linewidth=1.0)
        source_nu_Fnu = c_micron_s/self.lambda_array*self.flux_array*Jy2erg_s_cm2 #erg s-1 cm-2
        #next two lines is to avoid lower limits with large errors
        #to go out of the figure define a false 0.5 error for those marked as upper limit
        error_flux_for_SED_plot = self.err_flux_array
        error_flux_for_SED_plot[self.upper_limit_array] = 0.5
        plt.errorbar(self.lambda_array,source_nu_Fnu,yerr=error_flux_for_SED_plot*source_nu_Fnu,fmt=marker,
                     uplims=self.upper_limit_array)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\lambda\,(\mathrm{\mu m})$')
        plt.ylabel(r'$\nu F_\nu\,(\mathrm{erg\,s^{-1}\,cm^{-2}})$')
        if title is not None:
            plt.title(title)
        if figname is not None:
            plt.savefig(figname, dpi=300, bbox_inches="tight")
            print('Image saved in ',figname)
        plt.show()

        
    def plot_multiple_seds(self,models=None,figsize=(6,4),xlim=[1.0,1000.0],ylim=[1.0e-12,1.0e-5],marker='k*',cmap='rainbow_r',colorbar=True,title=None,figname=None):
        '''
        Plots multiple SEDs.
        
        Parameters
        ----------
        models: `astropy.table`
            Table with the sed results.

        figsize: tuple
            specify the size of the figure. Default is (6,4)

        xlim: array
            array to specify the limits in the x axis. Default is [1.0,1000.0]

        ylim: array
            array to specify the limits in the y axis. Default is [1.0e-12,1.0e-5]
            
        marker: str
            defines the marker style and color like matplotlib. Default is 'k*'

        cmap: str
            defines the colormap like matplotlib of the SEDs. Default is 'rainbow_r'

        colorbar: bool
            Turn on or off the colorbar of the plot. Default is True
            
        title: str, optional
            Defines the title of the plot. Default is None.

        figname: str, optional
                    A path, or a Python file-like object. Note that fname is used verbatim,
                    and there is no attempt to make the extension. Default is None.
                    Note that one can choose the format of the figure by changing the extension,
                    e.g., figure.pdf would generate the figure in PDF format.
            
        Returns
        ----------
        plot of the multiple SEDs
        '''
        
        master_dir = self.master_dir
        norm_extc_law = self.extc_law
        
        cmap = plt.cm.ScalarMappable(cmap=cmap,
                                     norm=colors.LogNorm(vmin=models['chisq'].min(),
                                                         vmax=models['chisq'].max()))
        #cmap.set_array([])

        plt.figure(figsize=figsize)
        for ind_mod in models:
            sed = np.loadtxt(master_dir+'/Model_SEDs/sed/'+ind_mod['SED_number']+'.dat',unpack=True)
            lambda_model = sed[0] #micron
            flux_model = sed[1]*Lsun2erg_s/(4.0*np.pi*(pc2cm*self.dist)**2.0) #from Lsun to erg s-1 cm-2
            flux_model_extincted = flux_model*10.0**(-0.4*ind_mod['av']*norm_extc_law)
            plt.plot(lambda_model,flux_model_extincted,'-',linewidth=0.5,zorder=-ind_mod['chisq'],c=cmap.to_rgba(ind_mod['chisq']))
           
        best_model = models[models['chisq']==models['chisq'].min()]
        sed = np.loadtxt(master_dir+'/Model_SEDs/sed/'+best_model['SED_number'].data[0]+'.dat',unpack=True)
        lambda_model = sed[0] #micron
        flux_model = sed[1]*Lsun2erg_s/(4.0*np.pi*(pc2cm*self.dist)**2.0) #from Lsun to erg s-1 cm-2
        flux_model_extincted = flux_model*10.0**(-0.4*best_model['av']*norm_extc_law)
        plt.plot(lambda_model,flux_model_extincted,'k-',linewidth=1.0,zorder=-best_model['chisq'].data[0])
        source_nu_Fnu = c_micron_s/self.lambda_array*self.flux_array*Jy2erg_s_cm2 #erg s-1 cm-2
        #next two lines is to avoid lower limits with large errors
        #to go out of the figure define a false 0.5 error for those marked as upper limit
        error_flux_for_SED_plot = self.err_flux_array
        error_flux_for_SED_plot[self.upper_limit_array] = 0.5
        plt.errorbar(self.lambda_array,source_nu_Fnu,yerr=error_flux_for_SED_plot*source_nu_Fnu,fmt=marker,
                     uplims=self.upper_limit_array)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\lambda\,(\mathrm{\mu m})$')
        plt.ylabel(r'$\nu F_\nu\,(\mathrm{erg\,s^{-1}\,cm^{-2}})$')
        if title is not None:
            plt.title(title)
        if colorbar:
            cbar = plt.colorbar(cmap,label=r'$\chi^2$')
            cbar.set_ticks([models['chisq'].min(),5,10,50,models['chisq'].max()])
            cbar.set_ticklabels([np.around(models['chisq'].min(),decimals=1),5,10,50,np.around(models['chisq'].max(),decimals=1)])
        if figname is not None:
            plt.savefig(figname, dpi=300, bbox_inches="tight")
            print('Image saved in ',figname)
        plt.show()
    
    
    def plot2d(self,models,figsize=(10,5),marker='s',marker_size=50,cmap='rainbow_r',title=None,figname=None):
        '''
        2D Plots of the SEDs results.
        
        Parameters
        ----------
        models: `astropy.table`
            Table with the sed results.

        figsize: tuple
            specify the size of the figure. Default is (10,5)

        marker: str
            defines the marker style like matplotlib. Default is 's'

        marker_size: int
            defines the marker size. Default is 50

        cmap: str
            defines the colormap like matplotlib of the 2D plot. Default is 'rainbow_r'

        title: str, optional
            Defines the title of the plot. Default is None.

        figname: str, optional
                    A path, or a Python file-like object. Note that fname is used verbatim,
                    and there is no attempt to make the extension. Default is None.
                    Note that one can choose the format of the figure by changing the extension,
                    e.g., figure.pdf would generate the figure in PDF format.
            
        Returns
        ----------
        plot the 2D results
        '''

        fig, axs = plt.subplots(1, 3, constrained_layout=True,figsize=figsize)

        triple_MC_sigma = unique(models,keys=['mcore','sigma'])
        ax1 = axs[0]
        sct1 = ax1.scatter(triple_MC_sigma['mcore'],
                           triple_MC_sigma['sigma'],
                           linewidths=1, alpha=.8,
                           #edgecolor='k',
                           s = marker_size,
                           c=triple_MC_sigma['chisq'],
                           marker=marker,
                           norm=colors.LogNorm(),
                           cmap=cmap,
                           zorder=-1)

        sct1_best = ax1.scatter(triple_MC_sigma['mcore'][triple_MC_sigma['chisq']==triple_MC_sigma['chisq'].min()],
                                triple_MC_sigma['sigma'][triple_MC_sigma['chisq']==triple_MC_sigma['chisq'].min()],
                                linewidths=1,
                                alpha=1.0,
                                s = 100,
                                color='k',
                                marker='+',
                                zorder=1)

        ax1.set_xscale('log')
        ax1.set_yscale('log')

        ax1.plot([0,triple_MC_sigma['mcore'][triple_MC_sigma['chisq']==triple_MC_sigma['chisq'].min()]],
                 [triple_MC_sigma['sigma'][triple_MC_sigma['chisq']==triple_MC_sigma['chisq'].min()],
                  triple_MC_sigma['sigma'][triple_MC_sigma['chisq']==triple_MC_sigma['chisq'].min()]],
                 'k--',alpha=0.1)

        ax1.plot([triple_MC_sigma['mcore'][triple_MC_sigma['chisq']==triple_MC_sigma['chisq'].min()],
                  triple_MC_sigma['mcore'][triple_MC_sigma['chisq']==triple_MC_sigma['chisq'].min()]],
                 [0,triple_MC_sigma['sigma'][triple_MC_sigma['chisq']==triple_MC_sigma['chisq'].min()]],
                 'k--',alpha=0.1)


        ax1.set_xlim(8.0,600.0)
        ax1.set_ylim(0.05,4.0)

        ax1.set_xticks([10.0,40.0,120.0,480.0])
        ax1.set_yticks(np.unique(models['sigma']))

        ax1.set_xticklabels([10.0,40.0,120.0,480.0])
        ax1.set_yticklabels(np.unique(models['sigma']))

        ax1.minorticks_off()#turning off minor ticks
        
        ax1.set_xlabel(r'$M_\mathrm{c}\,(M_\odot)$')
        ax1.set_ylabel(r'$\Sigma_\mathrm{cl}\,(\mathrm{g\,cm^{-2}})$')

        ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')


        triple_MC_ms = unique(models,keys=['mcore','mstar'])
        ax2 = axs[1]
        sct2 = ax2.scatter(triple_MC_ms['mcore'],
                           triple_MC_ms['mstar'],
                           linewidths=1,
                           alpha=.8,
                           #edgecolor='k',
                           s = marker_size,
                           c=triple_MC_ms['chisq'],
                           marker=marker,
                           norm=colors.LogNorm(),
                           cmap=cmap,
                           zorder=-1)

        sct2_best = ax2.scatter(triple_MC_ms['mcore'][triple_MC_ms['chisq']==triple_MC_ms['chisq'].min()],
                                triple_MC_ms['mstar'][triple_MC_ms['chisq']==triple_MC_ms['chisq'].min()],
                                linewidths=1,
                                alpha=1.0,
                                s = 100,
                                color='k',
                                marker='+',
                                zorder=1)

        ax2.set_xscale('log')
        ax2.set_yscale('log')

        ax2.plot([0,triple_MC_ms['mcore'][triple_MC_ms['chisq']==triple_MC_ms['chisq'].min()]],
                 [triple_MC_ms['mstar'][triple_MC_ms['chisq']==triple_MC_ms['chisq'].min()],
                  triple_MC_ms['mstar'][triple_MC_ms['chisq']==triple_MC_ms['chisq'].min()]],
                 'k--',alpha=0.1)

        ax2.plot([triple_MC_ms['mcore'][triple_MC_ms['chisq']==triple_MC_ms['chisq'].min()],
                  triple_MC_ms['mcore'][triple_MC_ms['chisq']==triple_MC_ms['chisq'].min()]],
                 [0,triple_MC_ms['mstar'][triple_MC_ms['chisq']==triple_MC_ms['chisq'].min()]],
                 'k--',alpha=0.1)

        ax2.set_xlim(8.0,600.0)
        ax2.set_ylim(0.4,200.0)

        ax2.set_xticks([10.0,40.0,120.0,480.0])
        ax2.set_yticks([0.5,2.0,8.0,32.0,128.0])

        ax2.set_xticklabels([10.0,40.0,120.0,480.0])
        ax2.set_yticklabels([0.5,2.0,8.0,32.0,128.0])

        ax2.minorticks_off()#turning off minor ticks
        
        ax2.set_xlabel(r'$M_\mathrm{c}\,(M_\odot)$')
        ax2.set_ylabel(r'$m_*\,(M_\odot)$')

        ax2.set_aspect(1.0/ax2.get_data_ratio(), adjustable='box')


        triple_sigma_ms = unique(models,keys=['sigma','mstar'])
        ax3 = axs[2]
        sct3 = ax3.scatter(triple_sigma_ms['sigma'],
                           triple_sigma_ms['mstar'],
                           linewidths=1,
                           alpha=.8,
                           #edgecolor='k',
                           s = marker_size,
                           c=triple_sigma_ms['chisq'],
                           marker=marker,
                           norm=colors.LogNorm(),
                           cmap=cmap,
                           zorder=-1)

        sct3_best = ax3.scatter(triple_sigma_ms['sigma'][triple_sigma_ms['chisq']==triple_sigma_ms['chisq'].min()],
                                triple_sigma_ms['mstar'][triple_sigma_ms['chisq']==triple_sigma_ms['chisq'].min()],
                                linewidths=1,
                                alpha=1.0,
                                s = 100,
                                color='k',
                                marker='+',
                                zorder=1)

        ax3.set_xscale('log')
        ax3.set_yscale('log')

        ax3.plot([0.0,triple_sigma_ms['sigma'][triple_sigma_ms['chisq']==triple_sigma_ms['chisq'].min()]],
                 [triple_sigma_ms['mstar'][triple_sigma_ms['chisq']==triple_sigma_ms['chisq'].min()],
                  triple_sigma_ms['mstar'][triple_sigma_ms['chisq']==triple_sigma_ms['chisq'].min()]],
                 'k--',alpha=0.1)

        ax3.plot([triple_sigma_ms['sigma'][triple_sigma_ms['chisq']==triple_sigma_ms['chisq'].min()],
                  triple_sigma_ms['sigma'][triple_sigma_ms['chisq']==triple_sigma_ms['chisq'].min()]],
                 [0.0,triple_sigma_ms['mstar'][triple_sigma_ms['chisq']==triple_sigma_ms['chisq'].min()]],
                 'k--',alpha=0.1)


        ax3.set_xlim(0.05,4.0)
        ax3.set_ylim(0.4,200.0)

        ax3.set_xticks(np.unique(models['sigma']))
        ax3.set_yticks([0.5,2.0,8.0,32.0,128.0])

        ax3.set_xticklabels(np.unique(models['sigma']))
        ax3.set_yticklabels([0.5,2.0,8.0,32.0,128.0])
        
        ax3.minorticks_off()#turning off minor ticks

        ax3.set_xlabel(r'$\Sigma_\mathrm{cl}\,(\mathrm{g\,cm^{-2}})$')
        ax3.set_ylabel(r'$m_*\,(M_\odot)$')

        ax3.set_aspect(1.0/ax3.get_data_ratio(), adjustable='box')

        cbar = fig.colorbar(sct3,label=r'$\chi^2$',shrink=0.36,pad=0.01,orientation='vertical',ax=axs)
        cbar.set_ticks([triple_sigma_ms['chisq'].min(),5,10,50,triple_sigma_ms['chisq'].max()])
        cbar.set_ticklabels([np.around(triple_sigma_ms['chisq'].min(),decimals=1),5,10,50,np.around(triple_sigma_ms['chisq'].max(),decimals=1)])
        
        if title is not None:
            fig.suptitle(title,y=0.8)
            
        if figname is not None:
            plt.savefig(figname, dpi=300, bbox_inches="tight")
            print('Image saved in ',figname)
            
        plt.show()
    
    def plot3d(self,models,figsize=(8,8),marker_size=200,cmap='rainbow_r',title=None,figname=None):
        '''
        3D plot for the SED model results
        
        Parameters
        ----------
        models: `astropy.table`
            Table with the sed results.

        figsize: tuple
            specify the size of the figure. Default is (8,8)

        marker_size: int
            defines the marker size. Default is 200

        cmap: str
            defines the colormap like matplotlib of the 2D plot. Default is 'rainbow_r'

        title: str, optional
            Defines the title of the plot. Default is None.

        figname: str, optional
                    A path, or a Python file-like object. Note that fname is used verbatim,
                    and there is no attempt to make the extension. Default is None.
                    Note that one can choose the format of the figure by changing the extension,
                    e.g., figure.pdf would generate the figure in PDF format.
            
        Returns
        ----------
        plot the 3D results
        '''
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        sct = ax.scatter(np.log10(models['mcore']),
                         np.log10(models['mstar']),
                         np.log10(models['sigma']),
                         linewidths=1, alpha=.8,
                         edgecolor='k',
                         s = marker_size,
                         c=models['chisq'],
                         norm=colors.LogNorm(),
                         cmap=cmap,
                         zorder=-1)

        cbar = fig.colorbar(sct,label=r'$\chi^2$',shrink=0.5,pad=0.1)
        cbar.set_ticks([np.min(models['chisq']),5,10,50,
                        np.max(models['chisq'])])
        cbar.set_ticklabels([np.around(np.min(models['chisq']),
                                       decimals=1),5, 10,50,np.around(np.max(models['chisq']),decimals=1)])


        #plot the best 5 models
        sct = ax.scatter(np.log10(models['mcore'])[0:5],
                         np.log10(models['mstar'])[0:5],
                         np.log10(models['sigma'])[0:5],
                         linewidths=1, alpha=.8,
                         edgecolor='k',
                         s = 100,
                         color='k',
                         marker='.',
                         zorder=1)

        ax.set_xlabel(r'$M_c\,(M_\odot)$')
        ax.set_ylabel(r'$m_*\,(M_\odot)$')
        ax.set_zlabel(r'$\Sigma\,(\mathrm{g\,cm^{-2}})$')

        ax.set_xticks(np.log10([10.0,40.0,120.0,480.0]))
        ax.set_yticks(np.log10([0.5,2.0,8.0,32.0,128.0]))
        ax.set_zticks(np.log10([0.1,0.316,1.0,3.16]))

        ax.set_xticklabels([10.0,40.0,120.0,480.0])
        ax.set_yticklabels([0.5,2.0,8.0,32.0,128.0])
        ax.set_zticklabels([0.1,0.316,1.0,3.16])

        #setting the viewing defaults
        ax.azim = -60
        ax.dist = 10
        ax.elev = 15

        #plotting eye-guide planes
        x_mesh = np.arange(models['mcore'].min(),models['mcore'].max())
        y_mesh = np.arange(models['mstar'].min(),models['mstar'].max())

        xx, yy = np.meshgrid(x_mesh,y_mesh)

        zz = np.zeros((np.shape(xx)))

        ax.plot_surface(np.log10(xx),np.log10(yy),np.log10(zz+3.16),alpha=0.1)
        ax.plot_surface(np.log10(xx),np.log10(yy),np.log10(zz+1.0),alpha=0.1)
        ax.plot_surface(np.log10(xx),np.log10(yy),np.log10(zz+0.316),alpha=0.1)
        ax.plot_surface(np.log10(xx),np.log10(yy),np.log10(zz+0.1),alpha=0.1)
        
        if title is not None:
            plt.title(title)

        if figname is not None:
            plt.savefig(figname, dpi=300, bbox_inches="tight")
            print('Image saved in ',figname)

        plt.show()
