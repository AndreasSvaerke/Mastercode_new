from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors
#from matplotlib.colors import LogNorm  # Optional for visualization
from PIL import Image
import os

from python_files.convert_sky_phy import find_coor_in_ICRS

def find_cutout_and_stack_and_save_koki(obs, data, header, position, path, cutout_size = 50):
    '''
    (obs takes in a string of what you want to be the name of the cutout)
    The function takes data and header of 3 images, the old, new, and difference images.
    It takes the positions as pixel coordinates of the objects.
    It creates cutouts at the specified positions with coutout_size x cutout_size dimensions.
    It places the cutouts bedide each other and saves them as a single image.
    The cutouts are saved in the path specified by the user, with names of the specified obs,
    together with the RA and DEC of the object.
    '''
    wcs_old = WCS(header[0])
    wcs_new = WCS(header[1])
    wcs_diff = WCS(header[2])
    mean_old, std_old = np.mean(data[0]), np.std(data[0])
    
    for i in range(len(position)):
        p = position[i]

        cutout_old = Cutout2D(data[0], p, size=cutout_size, wcs=wcs_old)
        cutout_new = Cutout2D(data[1], p, size=cutout_size, wcs=wcs_new)
        cutout_diff = Cutout2D(data[2], p, size=cutout_size, wcs=wcs_diff)
        vertical_bar = np.zeros((cutout_size,3))

        position_sky = find_coor_in_ICRS(header[2], p[0],p[1])

        ra = position_sky.ra.deg
        dec = position_sky.dec.deg

        output_filename = f"{obs}_{ra}_{dec}.png"
        output_path = os.path.join(path, output_filename)

        vertical_bar[:] = np.nan

        stacked = np.hstack((cutout_old.data, vertical_bar, cutout_new.data, vertical_bar, cutout_diff.data ))
        

        plt.imshow(stacked, cmap = 'viridis', vmin = mean_old - std_old , vmax = mean_old + std_old)
        plt.axis('off')
        plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
        plt.close()



def find_cutout_and_stack_and_save_christian(obs, data, header, position, path, cutout_size = 50):

    wcs_old = WCS(header[0])
    wcs_new = WCS(header[1])
    wcs_diff = WCS(header[2])
    
    
    for i in range(len(position)):
        p = position[i]

        cutout_old = Cutout2D(data[0], p, size=cutout_size, wcs=wcs_old)
        cutout_new = Cutout2D(data[1], p, size=cutout_size, wcs=wcs_new)
        cutout_diff = Cutout2D(data[2], p, size=cutout_size, wcs=wcs_diff)

        position_sky = find_coor_in_ICRS(header[2], p[0],p[1])

        ra = position_sky.ra.deg
        dec = position_sky.dec.deg

        output_filename = f"{obs}_{ra}_{dec}.png"
        output_path = os.path.join(path, output_filename)

        stacked = (np.stack((cutout_old.data, cutout_new.data ,cutout_diff.data ), axis = -1))
        
        stacked_normalized = (stacked - np.min(stacked)) / (np.max(stacked) - np.min(stacked))
        stacked_uint8 = (255 * stacked_normalized).astype(np.uint8)

        image = Image.fromarray(stacked_uint8)
        image.save(output_path)
