from scipy.spatial import cKDTree
import numpy as np

def remove_nearby_sources(sources, threshold=500, max_neighbors=5):
    '''
    Takes a list of sources with positions in pixel coordinates.
    It checks wheter there are more neighbors than max_neighbors, within the threshold distance.
    It then removes both the source and the neighbors that are too close to it.
    The function returns the sources that are not too close to each other.
    '''
    sources = np.array(sources)

    # Build cKDTree for fast spatial search
    tree = cKDTree(sources)

    # Find all points within 'threshold' distance
    all_neighbors = tree.query_ball_tree(tree, r=threshold)

    sources_to_remove = set()

    for i, neighbors in enumerate(all_neighbors):
        neighbors = [j for j in neighbors if j != i]  # Exclude self

        if len(neighbors) > max_neighbors:
            sources_to_remove.add(i)

            # Sort the neighbors by distance
            dists = np.linalg.norm(sources[neighbors] - sources[i], axis=1)
            sorted_neighbors = [j for _, j in sorted(zip(dists, neighbors))]

            # Add the "extra" nearby ones
            sources_to_remove.update(sorted_neighbors[max_neighbors:])

    # Filter out the marked sources
    mask = np.ones(len(sources), dtype=bool)
    mask[list(sources_to_remove)] = False

    return sources[mask]


def remove_nearby_sources_leaves_central(sources, threshold=500, max_neighbors=5):
    '''
    Takes a list of sources with positions in pixel coordinates.
    It checks wheter there are more neighbors than the max_neighbors, within the threshold distance.
    It then removes the neighbors that are too close.

    The function returns the sources with postions x,y in pixel.
    '''
    sources = np.array(sources)

    # Build cKDTree for fast spatial search
    tree = cKDTree(sources)

    # Find all points within 'threshold' distance
    all_neighbors = tree.query_ball_tree(tree, r=threshold)

    sources_to_remove = set()

    for i, neighbors in enumerate(all_neighbors):
        neighbors = [j for j in neighbors if j != i]  # Exclude self

        if len(neighbors) > max_neighbors:
            # Sort the neighbors by distance to the current source
            dists = np.linalg.norm(sources[neighbors] - sources[i], axis=1)
            sorted_neighbors = [j for _, j in sorted(zip(dists, neighbors))]

            # Keep the current source (i), and max_neighbors closest neighbors
            sources_to_remove.update(sorted_neighbors[max_neighbors:])

    # Filter out the marked sources
    mask = np.ones(len(sources), dtype=bool)
    mask[list(sources_to_remove)] = False

    return sources[mask]