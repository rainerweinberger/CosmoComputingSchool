import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import h5py

simulation_directory = str(sys.argv[1])

# simulation parameters
Boxsize = 50 # Mpc/h
HubbleParam = 0.6774 # h
UnitMass = 1.0e10
Volume = Boxsize * Boxsize * Boxsize 


trajectories = np.zeros((32**3, 3, 64), dtype=np.float64)

mass = 1.0


for snapshot in np.arange(63):
    halo_filename = simulation_directory + "/fof_tab_%03d.hdf5" % snapshot 
    particle_filename = simulation_directory + "/snapshot_%03d.hdf5" % snapshot 

    ## large scale structure visualization
    try:
        data = h5py.File(particle_filename, "r")
    except:
        print("could not open "+particle_filename)
    pos = np.array(data["PartType1"]["Coordinates"], dtype=np.float64)
    
    # make use of the fact that ids range from 1 to NumPart (included) -> index is ParticleIDs - 1 
    ids = np.array(data["PartType1"]["ParticleIDs"], dtype=np.int32) - 1


    trajectories[ids,:,snapshot] = pos[:,:]

    ### scatter plot
    fig = plt.figure(figsize=(6.9,6.9))
    ax = plt.axes([0.1,0.1,0.87,0.87])

    if(pos.shape[0] > 32**3):
        i_select = np.random.uniform(low=0.0, high=pos.shape[0], size=32**3).astype(np.int)
    else:
        i_select = np.arange(pos.shape[0])
    
    ax.scatter(pos[i_select, 0], pos[i_select, 1], marker='.', s=0.2, alpha=1.0, rasterized=True)
    # for i in np.arange(32**3):
    #     if i % 10 == 0:
    #         ax.plot(trajectories[i,0,:snapshot], trajectories[i,1,:snapshot], c='k', lw=0.1)


    ax.set_xlim([0,Boxsize])
    ax.set_ylim([0,Boxsize])
    ax.set_xlabel('[Mpc/h]')
    ax.set_ylabel('[Mpc/h]')

    if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )
    fig.savefig(simulation_directory+'/plots/largeScaleStructure_%03d.png'%snapshot, dpi=300)
    plt.close(fig)

#create movie
rate = 16 #frames per s
pix_fmt = 'yuv420p'
from subprocess import call
call(["ffmpeg",  "-r",  str(rate),  "-i",  simulation_directory+'/plots/largeScaleStructure_%03d.png',  '-pix_fmt', pix_fmt, simulation_directory + "/movie_cosmobox.mp4"])