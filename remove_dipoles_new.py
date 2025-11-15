import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from photutils.detection import DAOStarFinder
from photutils.background import StdBackgroundRMS


def initialize_dao_star_finder(image_data, sigma = 3.0, threshold=10.0, fwhm=4.0, peak_max_value=63000):
    '''
    Initalizes DAOStarFinder with sigma, threshold, fwhm and peak_value.
    Returns the daofind, that can then be used on images to find objects.
    '''
    sigma_clip = SigmaClip(sigma=sigma)
    bkg_rms = StdBackgroundRMS(sigma_clip)
    std = bkg_rms(image_data)

    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std, peakmax=peak_max_value)

    return daofind

def remove_dipoles_with_flux(skycoords, data):
    '''
    Finds objects in skycoord intervals on the data image.
    Finds objects for a positive and negative values.
    Checks that the difference between the values is above a threshold.
    If the flux is below the threshold the skycoords are considered a dipole.
    If the flux is above the threshold the skycoords are considered a non-dipole.
    Returns a list of the skycoords with the dipoles removed
    '''
    dipole_indices = []
    non_dipole_indices = []
    for i, (x, y) in enumerate(skycoords):
        x = int(x)
        y = int(y)
        inter = 25
        # Extract a 50x50 cutout
        cutout = data[y-inter:y+inter, x-inter:x+inter]
        neg_cutout = -cutout

        # Initialize DAOStarFinder for each cutout
        dao_finder = initialize_dao_star_finder(cutout, threshold=2.7, fwhm=3)
        pos_sources = dao_finder(cutout)
        neg_sources = dao_finder(neg_cutout)

        if pos_sources is not None:
            pos_max_flux = max(pos_sources['flux'])
            if neg_sources is not None:
                neg_max_flux = max(neg_sources['flux'])
                flux_diff = pos_max_flux - neg_max_flux
                flux_threshold = 0.5 * np.abs(pos_max_flux)
                #print(f"pos_flux: {pos_max_flux}, neg_flux: {neg_max_flux}, flux_diff: {flux_diff}, flux_threshold: {flux_threshold}")
                if flux_diff < flux_threshold:
                    #print(f'dipole_{i}')
                    dipole_indices.append(i)
                if flux_diff > flux_threshold:
                    #print(f'non-dipole_{i}')
                    non_dipole_indices.append(i)
            if neg_sources is None:
                # No negative sources found, so not a dipole
                #print(f'non-dipole_{i}')
                non_dipole_indices.append(i)

    return skycoords[non_dipole_indices], skycoords[dipole_indices]




def plot_non_dipoles(removed_dipoles, data, header):
    for i in range(len(removed_dipoles)):
        x,y = removed_dipoles[i]  # 4 miss in non-dipole # 3-4 miss in dipole
        x = int(x)
        y = int(y)
        cutout = data[2][y-25:y+25, x-25:x+25]
        neg_cutout = -cutout
        
        dao_finder = initialize_dao_star_finder(cutout, 2.7, 3)
        pos_sources = dao_finder(cutout)
        neg_sources = dao_finder(neg_cutout)
        mean = np.mean(cutout)
        std = np.std(cutout)
        neg_mean = np.mean(neg_cutout)
        neg_std = np.std(neg_cutout)
        neg_flux_min = 5 * neg_std #find 5 sigma background as minimum
        #find skycoords
        w = WCS(header[2])
        ra_dec = w.pixel_to_world(x, y)
        sky_coords = (ra_dec)
        import matplotlib.pyplot as plt
        print('RA:', sky_coords.ra.deg, 'DEC:', sky_coords.dec.deg)
        plt.subplot(1, 2, 1)
        plt.imshow(cutout, vmin = mean - 5 * std, vmax = mean + 5 * std)
        plt.title("Original Cutout")
        plt.subplot(1, 2, 2)
        plt.imshow(neg_cutout, vmin = neg_mean - 5 * neg_std, vmax = neg_mean + 5 * neg_std)
        plt.title("Negated Cutout")
        plt.tight_layout()
        plt.show()
        if pos_sources is not None and neg_sources is not None:
            #print(f'pos_flux: {pos_sources["flux"]}')
            #print(f'neg_flux: {neg_flux_min}')
            if max(neg_sources['flux']) > neg_flux_min:
                pos_flux = max(pos_sources['flux']) * 10 * u.nJy
                neg_flux = max(neg_sources['flux']) * 10 * u.nJy
                print(f'neg_flux_min: {neg_flux_min}')
                print(f'pos_flux: {pos_flux}, neg_flux: {neg_flux}')
                print(f'flux diff: {pos_flux-neg_flux}, flux threshold: {0.4 * np.abs(pos_flux)}')
                if ( pos_flux-neg_flux) < 0.4 * np.abs(pos_flux):
                    print('dipole')
                else:
                    print('not dipole')
            if max(neg_sources['flux']) < neg_flux_min:
                print('not dipole')
        if neg_sources is None:
            print('no neg sources found')
            print('not dipole')


def plot_dipoles(dipoles, data, header):
    for i in range(len(dipoles)):
        x,y = dipoles[i]  # 4 miss in non-dipole # 3-4 miss in dipole
        x = int(x)
        y = int(y)
        cutout = data[2][y-25:y+25, x-25:x+25]
        neg_cutout = -cutout
        
        dao_finder = initialize_dao_star_finder(cutout, 2.7, 3)
        pos_sources = dao_finder(cutout)
        neg_sources = dao_finder(neg_cutout)
        mean = np.mean(cutout)
        std = np.std(cutout)
        neg_mean = np.mean(neg_cutout)
        neg_std = np.std(neg_cutout)
        neg_flux_min = 5 * neg_std #find 5 sigma background as minimum
        #find skycoords
        w = WCS(header[2])
        ra_dec = w.pixel_to_world(x, y)
        sky_coords = (ra_dec)
        import matplotlib.pyplot as plt
        print('RA:', sky_coords.ra.deg, 'DEC:', sky_coords.dec.deg)
        plt.subplot(1, 2, 1)
        plt.imshow(cutout, vmin = mean - 5 * std, vmax = mean + 5 * std)
        plt.title("Original Cutout")
        plt.subplot(1, 2, 2)
        plt.imshow(neg_cutout, vmin = neg_mean - 5 * neg_std, vmax = neg_mean + 5 * neg_std)
        plt.title("Negated Cutout")
        plt.tight_layout()
        plt.show()
        if pos_sources is not None and neg_sources is not None:
            #print(f'pos_flux: {pos_sources["flux"]}')
            #print(f'neg_flux: {neg_flux_min}')
            if max(neg_sources['flux']) > neg_flux_min:
                pos_flux = max(pos_sources['flux']) * 10 * u.nJy
                neg_flux = max(neg_sources['flux']) * 10 * u.nJy
                print(f'neg_flux_min: {neg_flux_min}')
                print(f'pos_flux: {pos_flux}, neg_flux: {neg_flux}')
                print(f'flux diff: {pos_flux-neg_flux}, flux threshold: {0.4 * np.abs(pos_flux)}')
                if ( pos_flux-neg_flux) < 0.4 * np.abs(pos_flux):
                    print('dipole')
                else:
                    print('not dipole')
            if max(neg_sources['flux']) < neg_flux_min:
                print('not dipole')
        if neg_sources is None:
            print('no neg sources found')
            print('not dipole')