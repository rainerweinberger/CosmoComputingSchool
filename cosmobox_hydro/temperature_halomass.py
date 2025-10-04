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

snapshot = 62


halo_filename = simulation_directory + "/fof_subhalo_tab_%03d.hdf5" % snapshot 
particle_filename = simulation_directory + "/snap_%03d.hdf5" % snapshot 

## large scale structure visualization
try:
    data = h5py.File(particle_filename, "r")
    halodata = h5py.File(halo_filename, "r")
except:
    print("could not open "+particle_filename)

h = np.float64(data['Parameters'].attrs['HubbleParam']) 
UnitVel = np.float64(data['Parameters'].attrs['UnitVelocity_in_cm_per_s']) 
UnitLength = np.float64(data['Parameters'].attrs['UnitLength_in_cm'])
mu = 0.6
SOLAR_MASS = 1.989e33
BOLTZMANN = 1.38065e-16
PROTONMASS = 1.67262178e-24
pos = np.array(data["PartType0"]["Coordinates"], dtype=np.float64)
utherm = np.array(data["PartType0"]["InternalEnergy"], dtype=np.float64)
density = np.array(data["PartType0"]["Density"], dtype=np.float64) * UnitMass / UnitLength / UnitLength / UnitLength * SOLAR_MASS * h * h 
mass = np.array(data["PartType0"]["Masses"], dtype=np.float64) * UnitMass / h
temperature = (2./3.) * utherm * UnitVel * UnitVel / BOLTZMANN * PROTONMASS * mu

GrpPos = np.array(halodata["Group"]["GroupPos"], dtype=np.float64) 
GrpR200c = np.array(halodata["Group"]["Group_R_Crit200"], dtype=np.float64) 
M200c = np.array(halodata["Group"]["GroupMass"], dtype=np.float64) * UnitMass / h
GrpTemp = np.zeros(M200c.shape, dtype=np.float64)

for i_halo in np.arange(len(M200c)):
    dx = pos[:,0]-GrpPos[i_halo,0]
    dy = pos[:,1]-GrpPos[i_halo,1]
    dz = pos[:,2]-GrpPos[i_halo,2]

    rr = np.sqrt(dx*dx + dy*dy + dz*dz)

    i_select = np.where((rr< GrpR200c[i_halo]))[0]

    if len(i_select) > 0:
        GrpTemp[i_halo] = np.max(temperature[i_select])

fig = plt.figure(figsize=(3.32,2.5))
# main axes element
ax = plt.axes((0.18,0.18,0.8,0.8))
ax.scatter(M200c, GrpTemp, s=1, c='r')


ax.set_xlim([8e12,9e14])
ax.set_ylim([8e5,2e8])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'halo mass [M$_\odot$]')
ax.set_ylabel(r'temperature [K]')

if not os.path.exists( simulation_directory+"/plots" ):
    os.mkdir( simulation_directory+"/plots" )
fig.savefig(simulation_directory+'/plots/temperature_halomass_%03d.png'%snapshot, dpi=300)
plt.close(fig)
