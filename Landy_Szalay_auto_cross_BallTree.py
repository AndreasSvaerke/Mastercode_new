import numpy as np
import astropy.units as u
from sklearn.neighbors import BallTree

def angular_bins(min_bin,max_bin,bins,type='linear-bin'):
        #theta_bins,dtheta_bin=np.linspace(min_bin,max_bin,bins,retstep=True)
        #return theta_bins, dtheta_bin
        if type=='linear-bin':
            theta_edges=np.linspace(min_bin,max_bin,bins+1)
        if type=='log-bin':
            theta_edges=np.logspace(np.log10(min_bin.to('arcmin').value),
                                    np.log10(max_bin.to('arcmin').value),
                                    bins+1)
            theta_edges=theta_edges*u.arcmin
        theta_bins=np.zeros(bins)*u.arcmin
        for i in range(bins):
            theta_bins[i]=(theta_edges[i]+theta_edges[i+1])/2
        return theta_bins, theta_edges

import numpy as np
import astropy.units as u
from sklearn.neighbors import BallTree

def Landy_Szalay_estimator_tree(data_objects, random_objects,
                                min_bin=1*u.arcmin, max_bin=30*u.arcmin,
                                bins=10, type='linear-bin', estimator_type='Natural'):
    """
    Landy-Szalay estimator using BallTree pair counting.
    two_point_correlation returns cumulative counts for each radius supplied,
    so we take np.diff to get per-bin counts.
    """

    # angular bins in radians (BallTree expects radians)
    theta_bins, theta_edges = angular_bins(min_bin, max_bin, bins, type=type)
    theta_edges_rad = theta_edges.to(u.radian).value  # len = bins + 1

    N_data = len(data_objects)
    N_random = len(random_objects)

    # Convert to radians for BallTree: shape (N, 2) with (lon, lat)
    data_radec = np.vstack([data_objects.ra.rad, data_objects.dec.rad]).T
    random_radec = np.vstack([random_objects.ra.rad, random_objects.dec.rad]).T

    # Build trees
    tree_data = BallTree(data_radec, metric='haversine')
    tree_rand = BallTree(random_radec, metric='haversine')

    # --- DD: cumulative counts at each edge, then per-bin = diff ---
    cum_DD = tree_data.two_point_correlation(data_radec, theta_edges_rad)  # len = bins+1
    DD = np.diff(cum_DD)  # len = bins

    # --- DR: cross cumulative counts, then per-bin ---
    # note: tree_rand.two_point_correlation(data_radec, ...) gives cumulative cross-counts
    cum_DR = tree_rand.two_point_correlation(data_radec, theta_edges_rad)
    DR = np.diff(cum_DR)

    # scale DR to match your original normalization (if desired)
    # your previous code multiplied DR by (N_data / N_random)
    DR = (N_data / N_random) * DR

    # --- RR: from random catalog ---
    cum_RR = tree_rand.two_point_correlation(random_radec, theta_edges_rad)
    RR = np.diff(cum_RR)
    RR = (N_data * (N_data - 1)) / (N_random * (N_random - 1)) * RR

    # --- Landy-Szalay estimator ---
    # Use full Landy-Szalay form (recommended)
    with np.errstate(divide='ignore', invalid='ignore'):
        if estimator_type == 'Landy-Szalay':
            ACF = (DD - 2 * DR + RR) / RR
        elif estimator_type == 'Natural':
            ACF = DD / RR - 1
    # fallback / alternative (yours was: ACF = DD / RR - 1)

    # --- ERROR: avoid div-by-zero but don't artificially inflate small non-zero counts ---
    DD_safe = DD.clip(min=1)
    ERROR = 1.0 / np.sqrt(DD_safe)

    return ACF, ERROR, theta_bins


def Landy_Szalay_estimator_cross_tree(
        data_objects_1, data_objects_2,
        random_objects_1, random_objects_2,
        min_bin=1*u.arcmin, max_bin=30*u.arcmin,
        bins=10, type='linear-bin'):
    """
    Landy-Szalay cross-correlation estimator using BallTree pair counting.
    Returns arrays of length = bins.
    """

    # angular bins (in radians for BallTree)
    theta_bins, theta_edges = angular_bins(min_bin, max_bin, bins, type=type)
    theta_edges_rad = theta_edges.to(u.radian).value  # len = bins+1

    # number counts
    N_data_1 = len(data_objects_1)
    N_data_2 = len(data_objects_2)
    N_random_1 = len(random_objects_1)
    N_random_2 = len(random_objects_2)

    # convert to radians (lon, lat format)
    def radec_array(cat):
        return np.vstack([cat.ra.rad, cat.dec.rad]).T

    data1_radec = radec_array(data_objects_1)
    data2_radec = radec_array(data_objects_2)
    rand1_radec = radec_array(random_objects_1)
    rand2_radec = radec_array(random_objects_2)

    # build BallTrees
    tree_data1 = BallTree(data1_radec, metric='haversine')
    tree_data2 = BallTree(data2_radec, metric='haversine')
    tree_rand1 = BallTree(rand1_radec, metric='haversine')
    tree_rand2 = BallTree(rand2_radec, metric='haversine')

    # --- pair counts (cumulative -> per-bin via np.diff) ---
    cum_D1D2 = tree_data1.two_point_correlation(data2_radec, theta_edges_rad)
    D1D2 = np.diff(cum_D1D2) / (N_data_1 * N_data_2)

    cum_D1R2 = tree_data1.two_point_correlation(rand2_radec, theta_edges_rad)
    D1R2 = np.diff(cum_D1R2) / (N_data_1 * N_random_2)

    cum_D2R1 = tree_data2.two_point_correlation(rand1_radec, theta_edges_rad)
    D2R1 = np.diff(cum_D2R1) / (N_data_2 * N_random_1)

    cum_R1R2 = tree_rand1.two_point_correlation(rand2_radec, theta_edges_rad)
    R1R2 = np.diff(cum_R1R2) / (N_random_1 * N_random_2)

    # --- Landy-Szalay cross-correlation ---
    with np.errstate(divide='ignore', invalid='ignore'):
        CCF = (D1D2 - D1R2 - D2R1 + R1R2) / R1R2

    # --- Poisson error ---
    D1D2_safe = (D1D2 * N_data_1 * N_data_2).clip(min=1)
    ERROR = 1.0 / np.sqrt(D1D2_safe)

    return CCF, ERROR, theta_bins



def Davis_Peebles_estimator_cross_tree(
        data_objects_1, data_objects_2,
        random_objects_1,
        min_bin=1*u.arcmin, max_bin=30*u.arcmin,
        bins=10, type='linear-bin'):
    """
    Davis-Peebles estimator for cross-correlation using BallTree pair counting:
        w(theta) = D1D2 / D2R1 - 1
    Returns arrays of length = bins.
    """

    # angular bins (in radians for BallTree)
    theta_bins, theta_edges = angular_bins(min_bin, max_bin, bins, type=type)
    theta_edges_rad = theta_edges.to(u.radian).value  # len = bins+1

    # number counts
    N_data_1 = len(data_objects_1)
    N_data_2 = len(data_objects_2)
    N_random_1 = len(random_objects_1)

    # convert to radians (lon, lat format)
    def radec_array(cat):
        return np.vstack([cat.ra.rad, cat.dec.rad]).T

    data1_radec = radec_array(data_objects_1)
    data2_radec = radec_array(data_objects_2)
    rand1_radec = radec_array(random_objects_1)

    # build BallTrees
    tree_data1 = BallTree(data1_radec, metric='haversine')
    tree_data2 = BallTree(data2_radec, metric='haversine')
    tree_rand1 = BallTree(rand1_radec, metric='haversine')

    # --- D1D2 counts (per bin) ---
    cum_D1D2 = tree_data1.two_point_correlation(data2_radec, theta_edges_rad)
    D1D2 = np.diff(cum_D1D2) / (N_data_1 * N_data_2)

    # --- D2R1 counts (per bin) ---
    cum_D2R1 = tree_data2.two_point_correlation(rand1_radec, theta_edges_rad)
    D2R1 = np.diff(cum_D2R1) / (N_data_2 * N_random_1)

    # --- Davis-Peebles estimator ---
    with np.errstate(divide='ignore', invalid='ignore'):
        CCF = D1D2 / D2R1 - 1

    # --- Poisson error ---
    D1D2_safe = (D1D2 * N_data_1 * N_data_2).clip(min=1)
    ERROR = 1.0 / np.sqrt(D1D2_safe)

    return CCF, ERROR, theta_bins
