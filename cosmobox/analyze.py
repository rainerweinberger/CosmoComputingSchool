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
Volume = Boxsize * Boxsize * Boxsize / HubbleParam / HubbleParam / HubbleParam


for snapshot in np.arange(63):
    halo_filename = simulation_directory + "/fof_subhalo_tab_%03d.hdf5" % snapshot 
    particle_filename = simulation_directory + "/snap_%03d.hdf5" % snapshot 

    ## large scale structure visualization
    try:
        data = h5py.File(particle_filename, "r")
    except:
        print("could not open "+particle_filename)
    pos = np.array(data["PartType1"]["Coordinates"], dtype=np.float64) / HubbleParam / 1000.

    fig = plt.figure(figsize=(6.9,6.9))
    ax = plt.axes([0.1,0.1,0.87,0.87])

    if(pos.shape[0] > 32**3):
        i_select = np.random.uniform(low=0.0, high=pos.shape[0], size=32**3).astype(np.int)
    else:
        i_select = np.arange(pos.shape[0])
    ax.scatter(pos[i_select, 0], pos[i_select, 1], marker='.', s=0.05, alpha=0.5, rasterized=True)

    ## halos
    try:
        data = h5py.File(halo_filename, "r")
        """ get simulation data """
        ## simulation data
        GrpPos = np.array(data["Group"]["GroupPos"], dtype=np.float64) / HubbleParam / 1000.
        GrpR200c = np.array(data["Group"]["Group_R_Crit200"], dtype=np.float64) / HubbleParam / 1000.
        M200c = np.array(data["Group"]["GroupMass"], dtype=np.float64) * UnitMass
        M200c = np.sort(M200c)[::-1]
        CumMassFunction = np.cumsum(np.ones(M200c.shape) ) / Volume
        for i in np.arange(GrpR200c.shape[0]):
            ax.add_artist(plt.Circle((GrpPos[i,0], GrpPos[i,1]), GrpR200c[i], color='k', fill=False))
        bx = plt.axes([0.70,0.74,0.26,0.22])
        bx.plot(M200c, CumMassFunction)
        bx.set_xscale('log')
        bx.set_yscale('log')
        bx.set_xlim([9e12,5e14])
        bx.set_ylim([9e-7,2e-4])
        bx.set_xlabel(r"$M_{200,c}\ \mathrm{[M_\odot]}$")
        bx.set_ylabel(r"$n(>M)\ \mathrm{[Mpc^{-3}]}$")
    except:
        print("could not open "+halo_filename)

    ax.set_xlim([0,Boxsize/HubbleParam])
    ax.set_ylim([0,Boxsize/HubbleParam])
    ax.set_xlabel('[Mpc]')
    ax.set_ylabel('[Mpc]')

    if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )
    fig.savefig(simulation_directory+'/plots/largeScaleStructure_%03d.png'%snapshot, dpi=300)
    plt.close(fig)

#create movie
rate = 16 #frames per s
pix_fmt = 'yuv420p'
from subprocess import call
call(["ffmpeg",  "-r",  str(rate),  "-i",  simulation_directory+'/plots/largeScaleStructure_%03d.png',  '-pix_fmt', pix_fmt, simulation_directory + "/movie_cosmobox.mp4"])