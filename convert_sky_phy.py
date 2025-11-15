from astropy.wcs import WCS
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

def find_coor_in_ICRS(header,x,y):
    '''
    takes the header corresponding to the image data and the x,y pixel coordinates of that image.
    It returns the ICRS coordinates of the object in the image.
    '''
    w = WCS(header)
    sky = w.pixel_to_world(x,y)   
    return sky

def find_skycoord_in_physical(header,sky):
    '''
    takes the header corresponding to the image data and the skycoordinates of that image.
    It returns the (x,y) pixel coordinates of the object in the image in a numpy array.
    '''
    x_obs = []
    y_obs = []
    w = WCS(header)
    for i in range(len(sky)):
        x,y = w.world_to_pixel(sky[i])
        x_obs.append(x)
        y_obs.append(y)
    
    x_y_obs = np.array([x_obs,y_obs])
    x_y_obs = x_y_obs.T

    return x_y_obs


def find_skycoord_in_physical_from_nonSkycoords(header,transients):
    '''
    takes the header corresponding to the image data and the skycoordinates of that image.
    It returns the (x,y) pixel coordinates of the object in the image in a numpy array.
    NonSkycoords refers to the Ra and Dec of the object if it is not in SkyCoord format, here you just use float values for RA and DEC.
    '''
    sky = convert_ra_dec_to_skycoords(transients)
    x_obs = []
    y_obs = []
    w = WCS(header)
    for i in range(len(sky)):
        x,y = w.world_to_pixel(sky[i])
        x_obs.append(x)
        y_obs.append(y)
    
    x_y_obs = np.array([x_obs,y_obs])
    x_y_obs = x_y_obs.T

    return x_y_obs

def find_x_y_coords_vector(position):
    '''
    Takes a list of postions in x,y pixel coordinates.
    Returns two arrays, one for x and one for y.
    '''
    position_arr = np.array(position)
    x_coords, y_coords = position_arr[:,0], position_arr[:,1]
    return x_coords, y_coords

def find_coor_in_ICRS_vector(header, position):
    '''
    Takes the header of the cooresponding image and a list of positions in x,y pixel coordinates.
    Returns the ICRS coordinates of the objects in the image using numpy for fast computing.
    '''
    x_coords, y_coords = find_x_y_coords_vector(position)
    sky_coords = find_coor_in_ICRS(header, x_coords, y_coords)
    return sky_coords


def convert_ra_dec_to_skycoords(transients):
    """
    Convert RA and Dec to SkyCoord object.
    Parameters of transients should be ra in [1] and dec in [2]
    """  
    skycoords = []
    for i in range(len(transients)):
        ra = float(transients[i][1])
        dec = float(transients[i][2])

        skycoords.append(SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs'))

    return skycoords