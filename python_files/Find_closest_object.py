import numpy as np

def from_sky_to_ra_dec(skycoords):
    ra = []
    dec = []
    for i in range(len(skycoords)):
        ra.append(skycoords[i].ra.deg)
        dec.append(skycoords[i].dec.deg)
    return ra, dec

def circle_distance(ra1, dec1, ra2, dec2):
    ra1_rad = np.radians(ra1)
    dec1_rad = np.radians(dec1)
    ra2_rad = np.radians(ra2)
    dec2_rad = np.radians(dec2)
    return np.degrees(np.arccos(np.sin(dec1_rad)*np.sin(dec2_rad) + np.cos(dec1_rad)*np.cos(dec2_rad)*np.cos(ra1_rad - ra2_rad)))

def find_closest_object(rows, cosmos_ra, cosmos_dec, skycoords):
    '''
    Uses double loop resulting in O(N x M) complexity and is therefore very slow.
    Takes rows from Cosmos calalog and the RA and DEC and finds the closest object to the given skycoords.
    Returns the ID of the closest object, the distance to the object and the redshift of the object.
    The skycoords is a tuple with the value of (RA, DEC).
    '''
    determined_objects = []
    determined_distance = []
    determined_redshift = []
    
    for i in range(len(skycoords)):
        transient_ra, transient_dec = skycoords[i]
        min_distance = 1e10
        closest_object = None
        final_redshift = None
        for j in range(len(cosmos_ra)):
            distance = circle_distance(transient_ra, transient_dec, cosmos_ra[j], cosmos_dec[j])
            if distance < min_distance:
                min_distance = distance
                closest_object = rows[j]['ID']
                final_redshift = rows[j]['LP_zfinal']
        determined_objects.append(closest_object)
        determined_distance.append(min_distance)
        determined_redshift.append(final_redshift)
    return determined_objects, determined_distance, determined_redshift


from scipy.spatial import cKDTree

def radec_to_xyz(ra, dec):
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)
    return np.vstack((x, y, z)).T

def find_closest_objects_kdtree(rows, cosmos_ra, cosmos_dec, skycoords):
    '''
    Uses cKDTree and is therefore much faster.
    Takes rows from Cosmos calalog and the RA and DEC and finds the closest object to the given skycoords.
    Returns the ID of the closest object, the distance to the object and the redshift of the object.
    The skycoords is a tuple with the value of (RA, DEC).
    '''
    cosmos_xyz = radec_to_xyz(np.array(cosmos_ra), np.array(cosmos_dec))
    tree = cKDTree(cosmos_xyz)

    transient_ra = np.array([ra for ra, dec in skycoords])
    transient_dec = np.array([dec for ra, dec in skycoords])
    transients_xyz = radec_to_xyz(transient_ra, transient_dec)

    # Query the nearest neighbors
    distances, indices = tree.query(transients_xyz, k=1)

    determined_objects = [rows[i]['ID'] for i in indices]
    determined_zfinal = [rows[i]['LP_zfinal'] for i in indices]
    determiend_zPDF = [rows[i]['LP_zPDF'] for i in indices]
    determined_zPDF_u68 = [rows[i]['LP_zPDF_u68'] for i in indices]
    determined_zPDF_l68 = [rows[i]['LP_zPDF_l68'] for i in indices]
    

    # Convert distance from chord length back to angular distance
    angular_distance = np.degrees(2 * np.arcsin(np.clip(distances / 2, 0, 1)))

    return determined_objects, angular_distance, determined_zfinal, determiend_zPDF, determined_zPDF_u68, determined_zPDF_l68