import numpy as np
from astropy.stats import SigmaClip
from photutils.detection import DAOStarFinder
from photutils.background import StdBackgroundRMS
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.spatial import cKDTree

def initialize_dao_star_finder(image_data, sigma = 3.0, threshold=10.0, fwhm=4.0, peak_max_value=63000):
    sigma_clip = SigmaClip(sigma=sigma)
    bkg_rms = StdBackgroundRMS(sigma_clip)
    std = bkg_rms(image_data)

    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std, peakmax=peak_max_value)

    return daofind

def remove_dipoles(image_data, header, coords, sigma, daofind_threshold, daofind_fwhm, arcsec_radius=2.0):
    # Initialize list to store non-dipole sources
    non_dipole_sources = []

    # Get pixel scale (arcsec/pixel) from WCS
    wcs = WCS(header)
    pixel_scale = np.mean(np.abs(wcs.pixel_scale_matrix.diagonal())) * 3600.0  # arcsec/pixel
    pixel_radius = arcsec_radius / pixel_scale

    for x, y in coords:
        x_int = int(x)
        y_int = int(y)

        # Skip sources too close to edge
        if x_int - 25 < 0 or y_int - 25 < 0 or x_int + 25 >= image_data.shape[1] or y_int + 25 >= image_data.shape[0]:
            continue

        # Cutouts
        cutout = image_data[y_int-25:y_int+25, x_int-25:x_int+25]
        neg_cutout = -cutout
        # Initialize DAOStarFinder with specified parameters
        daofind = initialize_dao_star_finder(cutout, sigma, threshold=daofind_threshold, fwhm=daofind_fwhm)
        # Run DAOStarFinder
        pos_sources = daofind(cutout)
        neg_sources = daofind(neg_cutout)

        is_dipole = False

        if pos_sources is not None and neg_sources is not None:
            # Shift positions back to full-image coordinates
            neg_x = neg_sources['xcentroid'] + (x_int - 25)
            neg_y = neg_sources['ycentroid'] + (y_int - 25)

            # Use KDTree for distance check
            neg_tree = cKDTree(np.vstack((neg_x, neg_y)).T)
            dist, _ = neg_tree.query([[x, y]], distance_upper_bound=pixel_radius)

            if dist[0] < pixel_radius:
                is_dipole = True

        if not is_dipole:
            non_dipole_sources.append((x, y))

    return non_dipole_sources

def remove_dipoles_with_flux(image_data, header, coords, sigma, daofind_threshold, daofind_fwhm, flux_threshold, arcsec_radius=2.0):
    # Initialize list to store non-dipole sources
    non_dipole_sources = []

    # Get pixel scale (arcsec/pixel) from WCS
    wcs = WCS(header)
    pixel_scale = np.mean(np.abs(wcs.pixel_scale_matrix.diagonal())) * 3600.0  # arcsec/pixel
    pixel_radius = arcsec_radius / pixel_scale

    for x, y in coords:
        x_int = int(x)
        y_int = int(y)

        # Skip sources too close to edge
        if x_int - 25 < 0 or y_int - 25 < 0 or x_int + 25 >= image_data.shape[1] or y_int + 25 >= image_data.shape[0]:
            continue

        # Cutouts
        cutout = image_data[y_int-25:y_int+25, x_int-25:x_int+25]
        neg_cutout = -cutout

        # Initialize DAOStarFinder with specified parameters
        daofind = initialize_dao_star_finder(cutout, sigma, threshold=daofind_threshold, fwhm=daofind_fwhm)
        # Run DAOStarFinder
        pos_sources = daofind(cutout)
        neg_sources = daofind(neg_cutout)

        is_dipole = False

        if pos_sources is not None and neg_sources is not None:
            # Shift positions back to full-image coordinates
            # pos_x = pos_sources['xcentroid'] + (x_int - 25)
            # pos_y = pos_sources['ycentroid'] + (y_int - 25)
            neg_x = neg_sources['xcentroid'] + (x_int - 25)
            neg_y = neg_sources['ycentroid'] + (y_int - 25)

            # Use KDTree for distance check
            neg_tree = cKDTree(np.vstack((neg_x, neg_y)).T)
            #print(pixel_radius)
            dist, idx = neg_tree.query([[x, y]], distance_upper_bound=pixel_radius)

            if dist[0] < pixel_radius:
                # Check if the flux of the negative source is 5 sigma above the background
                if 5 * np.std(-cutout) < abs(neg_sources['flux'][idx[0]]):
                    # Flux check
                    pos_flux = max(pos_sources['flux'])
                    neg_flux = -neg_sources['flux'][idx[0]]
                    flux_ratio = abs(pos_flux - abs(neg_flux)) / max(abs(pos_flux), abs(neg_flux))

                    # Find dipoles based on flux ratio
                    if flux_ratio < flux_threshold:
                        is_dipole = True

                    # Find dipoles based on symmetry in flux
                    # flux_symmetry = abs(pos_flux + neg_flux) / (abs(pos_flux) + abs(neg_flux))
                    # if flux_symmetry < flux_threshold:
                    #     is_dipole = True


        if not is_dipole:
            non_dipole_sources.append((x, y))

        if neg_sources is None:
            non_dipole_sources.append((x, y))

    return non_dipole_sources