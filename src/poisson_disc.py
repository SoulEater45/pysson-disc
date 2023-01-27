import numpy as np

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

    if isinstance(radii, np.typing.ArrayLike):
        radii_sampling = lambda: np.random.choice(radii)
    elif callable(radii):
        radii_sampling = lambda: radii()
    else:
        raise("radii must be an array or a function")

    # cell_size = r / np.sqrt(2)

    # grid_width = int(np.ceil(width / cell_size))
    # grid_height = int(np.ceil(height / cell_size))
    # grid = [None] * (grid_width * grid_height)

    # def grid_coords(p):
    #     return int(p[0] // cell_size), int(p[1] // cell_size)

    # def fits(p):
    #     x, y = grid_coords(p)
    #     for i in range(max(x - 2, 0), min(x + 3, grid_width)):
    #         for j in range(max(y - 2, 0), min(y + 3, grid_height)):
    #             g = grid[i + j * grid_width]
    #             if g is not None and distance(p, g) < r:
    #                 return False
    #     return True

    # active = []

    # def add_to_grid(p):
    #     x, y = grid_coords(p)
    #     grid[x + y * grid_width] = p
    #     active.append(p)

    # p0 = np.random.uniform(0, width), np.random.uniform(0, height)
    # add_to_grid(p0)

    # while active:
    #     i = np.random.randint(len(active))
    #     p = active[i]
    #     found = False
    #     for j in range(k):
    #         a = 2 * np.pi * np.random.uniform(0, 1)
    #         d = np.random.uniform(r, 2 * r)
    #         x = p[0] + d * np.cos(a)
    #         y = p[1] + d * np