import numpy as np
import matplotlib.pyplot as plt

def calculate_angular_correlation(real_angular, fake_angular, bin_edges=None):
    """
    Calculate angular correlation function with error handling.
    Input the angular distances found from cKDtree from the
    two lists you want to find the correlation between
    """
    if bin_edges is None:
        bin_edges = np.logspace(-2, -1, 20)  # Default bins: 0.01° to 0.1°
    
    # Count pairs in each bin
    pair_counts, bin_edges = np.histogram(real_angular, bins=bin_edges)
    random_counts, _ = np.histogram(fake_angular, bins=bin_edges)
    
    # Calculate bin centers and areas
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_areas = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)
    pair_density = pair_counts / bin_areas
    
    # Calculate correlation function with error handling
    with np.errstate(divide='ignore', invalid='ignore'):
        norm_real = pair_counts / len(real_angular)
        norm_random = random_counts / len(fake_angular)
        correlation = (norm_real / norm_random) - 1
    
    return {
        'pair_counts': pair_counts,
        'random_counts': random_counts,
        'pair_density': pair_density,
        'correlation': correlation,
        'bin_edges': bin_edges,
        'bin_centers': bin_centers
    }

def plot_angular_correlation(results, figsize=(16, 6)):
    """Plot the angular correlation analysis."""
    plt.figure(figsize=figsize)
    
    # Raw pair counts
    plt.subplot(1, 3, 1)
    plt.stairs(results['pair_counts'], results['bin_edges'], fill=True, color='blue')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Angular separation (degrees)')
    plt.ylabel('Number of pairs')
    plt.title('Raw pair counts')
    plt.grid(True, which="both", ls="--", alpha=0.3)
    
    # Normalized density
    plt.subplot(1, 3, 2)
    plt.stairs(results['pair_density'], results['bin_edges'], fill=True, color='green')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Angular separation (degrees)')
    plt.ylabel('Pairs per square degree')
    plt.title('Surface density')
    plt.grid(True, which="both", ls="--", alpha=0.3)
    
    # Correlation function
    from scipy.optimize import curve_fit
    def power_law(theta, A, gamma):
        return A * theta**(-gamma)

    valid = np.isfinite(results['correlation'])
    print(valid)
    popt, pcov = curve_fit(power_law, 
                        results['bin_centers'][valid], 
                        results['correlation'][valid])
    print(popt)
    plt.subplot(1, 3, 3)
    plt.stairs(results['correlation'], results['bin_edges'], fill=True, color='lightcoral')
    plt.plot(results['bin_centers'], power_law(results['bin_centers'], popt[0], popt[1]), 
            'r--', label=f'Fit: Aθ$^{{-γ}}$, A = {popt[0]:.2f}, γ={popt[1]:.2f}')
    plt.xscale('log')
    plt.yscale('linear')  # Try 'log' if power law
    plt.xlabel('Angular separation (degrees)')
    plt.ylabel('Correlation function ω(θ)')
    plt.title('Correlation with power-law fit')
    plt.axhline(0, color='k', linestyle='--')
    plt.legend()
    plt.grid(True)


#HOW TO USE:

# # Convert the transients' RA/Dec to Cartesian coordinates
# transients_xyz = radec_to_xyz(ra_deg, dec_deg)

# fake_transient_ra = np.array([float(ra) for id, ra, dec, z in fake_transients['05_1']])
# fake_transient_dec = np.array([float(dec) for id, ra, dec, z in fake_transients['05_1']])
# fake_transient_id = np.array([float(id) for id, ra, dec, z in fake_transients['05_1']])
# sorted_indices = np.argsort(fake_transient_id)
# fake_transient_id = fake_transient_id[sorted_indices]

# fake_transients_xyz = radec_to_xyz(fake_transient_ra, fake_transient_dec)


# from scipy.spatial import cKDTree
# from sklearn.neighbors import BallTree
# galaxy_datasets = {
#     '05_1': galaxy_positions_05_1,
#     '1_15': galaxy_positions_1_15,
#     '15_2': galaxy_positions_15_2,
#     '2_25': galaxy_positions_2_25,
#     '25_3': galaxy_positions_25_3,
#     '3_35': galaxy_positions_3_35,
#     '35_4': galaxy_positions_35_4,
#     '4_45': galaxy_positions_4_45,
#     '45_5': galaxy_positions_45_5,
#     'fake': fake_galaxies,
#     'entire': np.array([ref_ra, ref_dec]).T
# }

# radius_deg = 1.0 / 60 # 1 arcmin
# radius_rad = np.deg2rad(radius_deg)  # Convert to radians for chord length

# # The maximum Euclidean distance between two unit vectors separated by angle θ is 2*sin(θ/2)
# max_euclidean_dist = 2 * np.sin(radius_rad / 2)

# # Initialize dictionaries to store results
# trees = {}
# Btrees = {}

# distances = {}
# counts = {}
# indices = {}

# clustering_strength = {}

# # Process each dataset
# for name, positions in galaxy_datasets.items():
#     #print(positions)
#     ra = [ra for ra, dec in positions]
#     dec = [dec for ra, dec in positions]
#     #print(ra)
#     #print(dec)
#     # Create KDTree for this galaxy set
#     trees[name] = cKDTree(radec_to_xyz(ra, dec))
#     Btrees[name] = BallTree(radec_to_xyz(ra, dec))
    
#     # Query with fake transients (k=1 finds the single nearest neighbor)
#     distances[name], indices[name] = trees[name].query(fake_transients_xyz, k=1)
    


#     counts[name], _ = Btrees[name].query_radius(fake_transients_xyz, r=max_euclidean_dist, return_distance=True)

#     # Convert to angular separation in degrees
#     angular_distances = {name: 2 * np.degrees(np.arcsin(dist/2)) for name, dist in distances.items()}
#     #angular_distance = {name: np.degrees(2 * np.arcsin(np.clip(dist / 2, 0, 1))) for name, dist in distances.items()}

#     clustering_strength[name] = np.mean([len(c) for c in counts[name]])


# from angular_correlation_try import calculate_angular_correlation, plot_angular_correlation
# # Main analysis
# #bin_edges = np.logspace(-3, 0, 20) # from 0.001° to 1° (logarithmic spacing)
# bin_edges = np.logspace(-2, -1, 20) # from 0.01° to 0.1° (logarithmic spacing)
# real_angular = angular_distances['05_1']
# fake_angular = angular_distances['fake']

# # Calculate and plot
# results = calculate_angular_correlation(real_angular, fake_angular, bin_edges)
# plot_angular_correlation(results)