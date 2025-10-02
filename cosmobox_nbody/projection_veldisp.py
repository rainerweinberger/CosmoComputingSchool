import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import h5py

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
vel = np.array(data["PartType1"]["Velocities"], dtype=np.float64)
h = np.float64(data['Parameters'].attrs['HubbleParam']) 
mass = np.float64(data['Header'].attrs['MassTable'][1]) * 1e10 / h

### density map
fig = plt.figure(figsize=(6.9,6.9))
ax = plt.axes([0.1,0.1,0.87,0.87])
cax = plt.axes([0.12, 0.12, 0.02, 0.3])

from scipy.spatial import KDTree
from matplotlib.colors import LogNorm
searchtree = KDTree(pos, boxsize=Boxsize)

# calculate "smoothing length"
k = 8
dist, i_part = searchtree.query(pos, k=k)

hsml = dist[:,-1]

# project on a canvas
Nplot=256
Edges1d = np.linspace(0.0, Boxsize, Nplot+1, endpoint=True, dtype=np.float64)
Grid1d = 0.5 * (Edges1d[1:] + Edges1d[:-1])
dx = Boxsize/np.float64(Nplot)
Canvas = np.zeros((Nplot, Nplot), dtype=np.float64)
Norm = np.zeros((Nplot, Nplot), dtype=np.float64)

for i in np.arange(32**3):

    sigma2 = np.var(vel[i_part[i], 2])

    jx = np.int32((pos[i,0]-hsml[i])/dx - 1)
    mx = np.int32((pos[i,0]+hsml[i])/dx + 1)


    while jx < mx:
        
        jy = np.int32((pos[i,1]-hsml[i])/dx - 1)
        my = np.int32((pos[i,1]+hsml[i])/dx + 1)
        while jy < my:
            wk = dx * dx / (np.pi * hsml[i] * hsml[i])

            if jx>=0 and jy>=0 and jx<Nplot and jy<Nplot: # note: no periodic wrapping of canvas yet
                Canvas[jx, jy] += sigma2 * mass * wk / dx / dx
                Norm[jx, jy] += mass * wk / dx / dx
            jy += 1
        jx += 1

Canvas /= Norm

Canvas[...] = np.sqrt(Canvas[...])

pc = ax.imshow(Canvas.T, cmap=plt.get_cmap("cubehelix_r"), extent=(0, Boxsize, 0, Boxsize), origin="lower", norm=LogNorm(vmin=80,vmax=1500))

plt.colorbar(pc, cax=cax)
cax.set_ylabel("velocity dispersion [km/s]")

ax.set_xlim([0,Boxsize])
ax.set_ylim([0,Boxsize])
ax.set_xlabel('[Mpc/h]')
ax.set_ylabel('[Mpc/h]')

if not os.path.exists( simulation_directory+"/plots" ):
    os.mkdir( simulation_directory+"/plots" )
fig.savefig(simulation_directory+'/plots/VelDispProjection_%03d.png'%snapshot, dpi=300)
plt.close(fig)
