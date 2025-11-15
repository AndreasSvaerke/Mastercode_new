import numpy as np
from astropy.stats import SigmaClip
from photutils.background import StdBackgroundRMS
from photutils.detection import DAOStarFinder

def find_sources_with_annuli(image_data, fwhm=4.0, peak_max_value=63000, r_out=100):
    # Robust background noise estimation
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_rms = StdBackgroundRMS(sigma_clip)
    std = bkg_rms(image_data)

    # Detect sources
    daofind = DAOStarFinder(fwhm=fwhm, threshold=10.0 * std, peakmax=peak_max_value)
    sources = daofind(image_data)
    
    if sources is None:
        return None  # No sources found
    
    # Filter out saturated sources and check annulus bounds
    height, width = image_data.shape
    valid_sources = []
    
    for source in sources:
        x, y = source['xcentroid'], source['ycentroid']
        if (x - r_out >= 0 and x + r_out < width and
            y - r_out >= 0 and y + r_out < height):
            valid_sources.append(source)
    
    return valid_sources, daofind


def find_pos_from_source(valid_sources):
    positions = []
    for source in valid_sources:
        positions.append(np.transpose((source['xcentroid'],source['ycentroid'])))
    return positions

