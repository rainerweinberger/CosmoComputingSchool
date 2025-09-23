import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import h5py

import numpy as np
from scipy.spatial import cKDTree

def two_point_correlation(positions, r_bins, box_size=None, seed=None):
    """
    Simple 2-point correlation function using DD/RR - 1.

    Parameters
    ----------
    positions : ndarray, shape (N, 3)
        Particle positions.
    r_bins : array_like
        Radial separation bin edges.
    box_size : float, optional
        Periodic cube size. If None, open boundaries.
    seed : int, optional
        RNG seed for reproducibility (for random catalog).

    Returns
    -------
    r_centers : ndarray
        Bin centers.
    xi : ndarray
        Two-point correlation function values.
    """
    rng = np.random.default_rng(seed)
    N = len(positions)

    # KDTree for data
    tree_data = cKDTree(positions, boxsize=box_size)

    # Data-data counts
    DD = np.diff([tree_data.count_neighbors(tree_data, r) for r in r_bins])

    # Random catalog
    if box_size is None:
        mins, maxs = positions.min(0), positions.max(0)
        randoms = rng.uniform(mins, maxs, size=(N, 3))
    else:
        randoms = rng.uniform(0, box_size, size=(N, 3))

    tree_rand = cKDTree(randoms, boxsize=box_size)

    # Random-random counts
    RR = np.diff([tree_rand.count_neighbors(tree_rand, r) for r in r_bins])

    # Normalize pair counts
    DD = DD / (N * (N - 1) / 2)
    RR = RR / (N * (N - 1) / 2)

    xi = DD / RR - 1
    r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
    return r_centers, xi



simulation_directory = str(sys.argv[1])

if len(sys.argv) > 2:
    snapshot=np.int32(sys.argv[2])
else:
    snapshot=62

# simulation parameters
Boxsize = 50 # Mpc/h
HubbleParam = 0.6774 # h
UnitMass = 1.0e10
Volume = Boxsize * Boxsize * Boxsize 



halo_filename = simulation_directory + "/fof_tab_%03d.hdf5" % snapshot 
particle_filename = simulation_directory + "/snapshot_%03d.hdf5" % snapshot 

## large scale structure visualization
try:
    data = h5py.File(particle_filename, "r")
except:
    print("could not open "+particle_filename)
pos = np.array(data["PartType1"]["Coordinates"], dtype=np.float64)

L = 50.0

r_bins = np.linspace(0, 25, 21)

r, xi = two_point_correlation(pos, r_bins, box_size=L)

fig = plt.figure(figsize=(5.0,4.0))
ax = plt.axes((0.15,0.15,0.8,0.8))

ax.plot(r, xi)

ax.set_xscale("log")
ax.set_yscale("log")

if not os.path.exists( simulation_directory+"/plots" ):
    os.mkdir( simulation_directory+"/plots" )
fig.savefig(simulation_directory + "/plots/twopoint_%03d.png"%snapshot)
