
""" load libraries """
import sys    ## system calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    ## plot stuff
from scipy.interpolate import interp1d    ## inetrpolation
plt.rcParams['text.usetex'] = True

createReferenceSolution = False
makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.argv[1])

FloatType = np.float64  # double precision: np.float64, for single use np.float32
Gcgs = 6.6738e-8

""" compare to reference output """
for i_snap in np.arange(65):
    filename = '/snapshot_%03d.hdf5'%i_snap
    print(simulation_directory+filename)
    try:
        data = h5py.File(simulation_directory+filename, "r")
    except:
        break
    ## figure
    pos = np.array( data['PartType2']['Coordinates'], dtype=np.float64 )
    
    fig = plt.figure(figsize=(5.,4.))
    ax = plt.axes((0.15,0.15,0.83,0.83))
    
    ax.scatter(pos[:, 0], pos[:, 1], marker='.', c='k', s=0.01, alpha=1.0, rasterized=True, zorder=2)
    ax.set_aspect('equal')
    ax.set_xlim([-160,160])
    ax.set_ylim([-120,120])

    ax.set_ylabel("[kpc]")
    ax.set_ylabel("[kpc]")
    ax.set_xlabel("[kpc]")
    ax.set_xlabel("[kpc]")
    ax.set_xlabel("[kpc]")

    if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )
    fig.savefig(simulation_directory+'/plots/stars_position_%03d.png'%i_snap)
    plt.close(fig)
    
#create movie
rate = 16 #frames per s
pix_fmt = 'yuv420p'
from subprocess import call
call(["ffmpeg",  "-r",  str(rate),  "-i",  simulation_directory+'/plots/stars_position_%03d.png',  '-pix_fmt', pix_fmt, simulation_directory + "/movie_galaxymerger.mp4"])
    
""" if everything is ok """
sys.exit(0)