#Making a polygon of the COSMOS field to generate random points within it
import numpy as np
from shapely.geometry import Point, Polygon

def generate_random_catalog(entire_reference_catalog, sample_size=10000, seed=42):
    '''
    Generates a random catalog within the created COSMOS field polygon.
    Is therefore important to use the entire COSMOS-Web reference catalog.
    '''
    min_ra = min(entire_reference_catalog['ra'])
    min_dec = min(entire_reference_catalog['dec'])
    max_ra = max(entire_reference_catalog['ra'])
    max_dec = max(entire_reference_catalog['dec'])

    cosmos_polygon = Polygon([
        (150.325, min_dec), (min_ra, 1.952), (149.92, max_dec), (max_ra, 2.465)  # polygon edges
    ])

    def generate_valid_points(num_points, seed=seed):
        ra_dec_list = []
        np.random.seed(seed)
        while len(ra_dec_list) < num_points:
            # Generate candidates within broad bounds
            ra = np.random.uniform(min_ra, max_ra)
            dec = np.random.uniform(min_dec, max_dec)
            # Check if point is inside the polygon
            if cosmos_polygon.contains(Point(ra, dec)):
                ra_dec_list.append((ra, dec))
        return np.array(ra_dec_list).T  # Returns (ra_array, dec_array)

    rand_ra, rand_dec = generate_valid_points(sample_size)

    random_catalog = {
        'ra': rand_ra,
        'dec': rand_dec,
        'id': np.arange(1, sample_size + 1),  # Unique IDs
    }

    return random_catalog