import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import astropy.units as u
from photutils.datasets import make_wcs
from photutils.aperture import SkyCircularAperture
from photutils.aperture import aperture_photometry


def mag_AB(sum_flux):
    '''
    Standard AB_magnitude function to convert flux (in jansky) to magnitude.
    '''
    if sum_flux == 0:
        return 99.0
    elif sum_flux < 0:
        return -99.0
    else:
        return -2.5*np.log10((sum_flux)/(3631 *u.Jy))


def magnitude_in_each_band(skycoords, data, header, size_arcsec = 0.1):
    '''
    Calculate the magnitude of the source. 
    Give skycoords of determined transients.
    Uses the data image and the header to find the pixel coordinates. 
    The size of the cutout is given in arcseconds, with standard being 0.1 arcsec.
    The pixel value is 10 nJy/pixel. This is being accounted for when summing the flux.
    '''
    size_deg = size_arcsec/3600

    magnitudes = []
    for i in range(len(skycoords)):
        magnitude = []
        for j in range(len(header)):
            wcs = WCS(header[j])
            x_center, y_center = wcs.world_to_pixel(skycoords[i])
            pixel_scale = np.abs(wcs.proj_plane_pixel_scales()[0])*3600 # Pixel scale in degrees/pixel   
            #print(pixel_scale.value)
            size_pixels = int(size_arcsec/pixel_scale.value) # Size of the cutout in pixels

            x_start = int(x_center - size_pixels / 2)
            x_end = int(x_center + size_pixels / 2)
            y_start = int(y_center - size_pixels / 2)
            y_end = int(y_center + size_pixels / 2)

            roi = data[j][y_start:y_end, x_start:x_end]
            finite_roi = roi[np.isfinite(roi)]
            roi_sum = (np.sum(finite_roi)*10 * 10**-(9) * u.Jy) # Sum of the flux in the cutout in Jy
            magnitude.append(mag_AB(roi_sum))
        magnitudes.append(magnitude)
    return magnitudes


def magnitude_in_each_band_with_aperture(skycoords, data, header, size_arcsec = 0.1):
    '''
    Calculate the magnitude of the source. 
    Give skycoords of determined transients.
    Uses the data image and the header to find the pixel coordinates. 
    The size of the cutout is given in arcseconds, with standard being 0.1 arcsec.
    The pixel value is 10 nJy/pixel. This is being accounted for when summing the flux.
    '''
    

    magnitudes = []
    for i in range(len(skycoords)):
        magnitude = []
        for j in range(len(header)):
            wcs = WCS(header[j])
            aperture = SkyCircularAperture(skycoords[i], r=size_arcsec * u.arcsec)
            pix_aperture = aperture.to_pixel(wcs)

            phot_table = aperture_photometry(data[j], pix_aperture, method='center')
            finite_roi = phot_table['aperture_sum'][np.isfinite(phot_table['aperture_sum'])]
            roi_sum = (np.sum(finite_roi.value)*10 * 10**-(9) * u.Jy)
            # Sum of the flux in the cutout in Jy
            #print(roi_sum)
            magnitude.append(mag_AB(roi_sum))
        #print(magnitude)
        magnitudes.append(magnitude)
    return magnitudes


def find_limiting_magnitude(filter, data, header, aperture_size = 0.1, n_apertures = 100, if_hist = 1):
    '''
    Takes data and header of a filter (write name of filter in the filter augment)
    It finds flux in random apertures and makes a histogram
    To not plot the histogram, set if_hist = 0
    The 1 sigma of the fluxses is set as the background flux level
    The limiting magnitude is calculated from the background flux level
    The limiting magnitude is returned in AB magnitudes
    '''


    #make 100 random apertures
    x = np.random.randint(1000, data.shape[1]-1000, n_apertures)
    y = np.random.randint(1000, data.shape[0]-1000, n_apertures)

    #convert to sky coordinates
    w = WCS(header)
    coords = w.pixel_to_world(x, y)
    #make a random aperture
    apertures = []
    for i in range(len(coords)):
        aperture = SkyCircularAperture(coords[i], r=aperture_size * u.arcsec)
        apertures.append(aperture)
    
    # find flux in each aperture in ADU
    fluxes_ADU = []
    for i in range(len(apertures)):
        flux = aperture_photometry(data, apertures[i] ,wcs = w)
        fluxes_ADU.append(flux['aperture_sum'][0])

    if if_hist == 1:
    # make and plot histogram of fluxes
        plt.hist(fluxes_ADU, bins=50)
        plt.xlabel('flux (ADU)')
        plt.ylabel('number of apertures')
        plt.title(filter)
        plt.show()
    
    # find the 1-sigma background flux level
    sigma = np.std(fluxes_ADU, ddof=1)
    bkg_ADU_1s = sigma
    bkg_ADU_3s = 3 * sigma
    bkg_ADU_5s = 5 * sigma
    bkg_ADU = np.array([bkg_ADU_1s, bkg_ADU_3s, bkg_ADU_5s])
    bkg = bkg_ADU * 10 * 10 ** -9 * u.Jy # convert to Jy
    #print('1-sigma background flux level:', bkg)
    #convert to magnitude
    bkg_mag = -2.5*np.log10((bkg)/(3631 *u.Jy))
    
    return bkg_mag


def flux_from_mag_AB(magnitudes):
    """
    Convert AB magnitude back to flux in Janskys.
    
    Parameters:
    magnitude (float or array-like): AB magnitude(s) to convert
    
    Returns:
    Quantity: flux in Janskys
    """
    return 3631 * u.Jy * 10**(-0.4 * magnitudes)