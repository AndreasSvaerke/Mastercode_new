import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

def find_skycoords(skycoords):
    '''
    Takes a list of (id, ra, dec, prediction) tuples and returns the sky coordinates as a list of SkyCoord objects.
    '''
    ra = []
    dec = []
    skycoordinates = []    
    for skycoord in skycoords:
        ra.append(float(skycoord[1]))
        dec.append(float(skycoord[2]))
    
    for i in range(len(ra)):
        skycoordinates.append((SkyCoord(ra=ra[i]*u.deg, dec=dec[i]*u.deg , frame='icrs')))


    return skycoordinates

def find_skycoords_2(transients):
    '''
    This function takes a list with the ID, RA, and DEC of the transients.
    It then returns two lists of SkyCoord objects, one for each observation (obs17 and obs18).
    '''
    skycoords_obs17 = []
    skycoords_obs18 = []

    for obj in transients:
        if obj[0].startswith('obs17'):
            RA = (float(obj[1]))
            DEC = (float(obj[2]))
            skycoords_obs17.append(SkyCoord(ra=RA * u.deg, dec=DEC * u.deg, frame='icrs'))

        if obj[0].startswith('obs18'):
            RA = (float(obj[1]))
            DEC = (float(obj[2]))
            skycoords_obs18.append(SkyCoord(ra=RA * u.deg, dec=DEC * u.deg, frame='icrs'))

    return skycoords_obs17, skycoords_obs18