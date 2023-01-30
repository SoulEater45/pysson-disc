import numpy as np


def __hypersphere_point(n: int, m: int=1):
    positions = np.random.normal(size=(m, n))
    return positions / np.sum(positions**2, axis=1)

def poisson_disc_samples(bounds, radii, repetitions=30, distance=np.hypot):
    """Generate samples in [0, width] x [0, height] using a Poisson-disc
    algorithm.

    Parameters
    ----------
    width, height : float
        The width and height of the box.
    r : float
        The minimum distance between samples.
    k : int, optional
        The number of samples to try before marking a point as inactive.
    distance : callable, optional
        A function that takes two (x, y) points and returns the distance
        between them.

    Returns
    -------
    samples : list of (x, y) pairs
        The generated samples.

    Notes
    -----
    Based on algorithm 4.3 from Bridson, "Fast Poisson Disk Sampling in
    Arbitrary Dimensions".

    """

    if isinstance(radii, np.ndarray):
        radii_sampling = lambda: np.random.choice(radii)
        rmax, rmin = [max(radii), min(radii)]
    elif callable(radii):
        radii_sampling = lambda: radii()
    else:
        raise("radii must be an numpy.ndarray or a function")

    dimensions = len(bounds)
    cell_size = 2 * rmin / np.sqrt(dimensions)
    grid_radius = lambda r: np.ceil((rmax + r) / cell_size)

    # Step 0
    cell_bounds = np.ceil(bounds / cell_size).astype(int)
    background_grid = np.zeros(tuple(cell_bounds) + (dimensions + 1, )) * np.nan

    # Step 1
    position = np.random.uniform(size=dimensions) * bounds
    position_index = tuple(np.ceil(position / cell_size).astype(int))

    background_grid[position_index] = np.append(position, radii_sampling())
    active = [background_grid[position_index]]

    # Step 2
    safe_stop = 10
    safe_counter = 0
    while ~len(active) and safe_counter < safe_stop:
        active_index = np.random.randint(len(active))
        active_position = active[active_index][:-1]
        active_radius = active[active_index][-1]

        found = False
        for j in range(repetitions):
            test_radius = radii_sampling()
            test_position = active_position + test_radius * np.squeeze(__hypersphere_point(dimensions, 1) * (test_radius + active_radius + np.min((test_radius, active_radius)) * np.random.rand()))
            
            test_grid_radius = grid_radius(test_radius)
            position_index = tuple(np.ceil(test_position / cell_size).astype(int))
            

        # Delete after testing
        safe_counter += 1

    # extract coresponding points and radii from grid
    samples_mask = np.any(~np.isnan(background_grid), axis=-1)
    samples = background_grid[samples_mask]
    points = samples[:, :-1]
    radii = samples[:, -1]

    return points, radii