import numpy as np
from numba import jit

def __hypersphere_point(n: int, m: int=1):
    positions = np.random.normal(size=(m, n))
    return positions / np.sqrt(np.sum(positions**2, axis=1))

@jit(nopython=False)
def poisson_disc_samples(bounds, radii, repetitions=30):
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
        raise(f"""
            radii: 
                expected: numpy.ndarray or lambda function
                received: {type(radii)}""")

    dimensions = len(bounds)
    cell_size = 2 * rmin / np.sqrt(dimensions)
    grid_radius = lambda r: np.ceil((rmax + r) / cell_size).astype(int)
    # grid_radius = lambda r: int((rmax + r) // cell_size)

    # Step 0
    cell_bounds = np.ceil(bounds / cell_size).astype(int)
    background_grid = np.zeros(tuple(cell_bounds) + (dimensions + 1, )) * np.nan

    # Step 1
    radius = radii_sampling()
    position = np.random.uniform(size=dimensions) * (bounds - 2 * radius) + radius
    position_index = tuple(np.floor(position / cell_size).astype(int))

    background_grid[position_index] = np.append(position, radius)
    active_list = [background_grid[position_index]]

    # Step 2
    while len(active_list) > 0:
        active_index = np.random.randint(len(active_list))
        active_position = active_list[active_index][:-1]
        active_radius = active_list[active_index][-1]

        found = False
        for j in range(repetitions):
            test_radius = radii_sampling()
            test_position = active_position + np.squeeze(__hypersphere_point(dimensions, 1) * (test_radius + active_radius + np.min((test_radius, active_radius)) * np.random.rand()))

            if np.any(test_position < test_radius) or np.any(test_position > bounds - test_radius):
                continue
            
            test_grid_radius = grid_radius(test_radius)
            position_index = tuple(np.floor(test_position / cell_size).astype(int))

            # selecting negibouring cells:
            # https://stackoverflow.com/a/56103114/4434071
            indexing = np.ix_([*(np.r_[np.max([0, position_index_coordinate - test_grid_radius]):np.min([shape_coordinate, position_index_coordinate + test_grid_radius + 1])] for (position_index_coordinate, shape_coordinate) in zip(position_index, background_grid.shape[:-1]))])
            neighbourhood = background_grid[indexing]
            if neighbourhood.size == 0 or np.all(np.isnan(neighbourhood)):
                continue
            neighbourhood = neighbourhood.reshape(-1, neighbourhood.shape[-1])

            # neighbourhood_positions = neighbourhood[..., :-1]
            neighbourhood_positions = neighbourhood[:, :-1]
            neighbourhood_radii = neighbourhood[:, -1]
            collision = np.any(np.sum((test_position[np.newaxis, :] - neighbourhood_positions)**2, axis=-1) < (test_radius + neighbourhood_radii)**2)
            if collision:
                continue

            data = np.append(test_position, test_radius)
            active_list.append(data)
            background_grid[position_index] = data
            found = True
            break

        if not found:
            active_list.pop(active_index)
    
    print(len(active_list))

    # extract coresponding points and radii from grid
    samples_mask = np.any(~np.isnan(background_grid), axis=-1)
    samples = background_grid[samples_mask]
    points = samples[:, :-1]
    radii = samples[:, -1]

    return points, radii