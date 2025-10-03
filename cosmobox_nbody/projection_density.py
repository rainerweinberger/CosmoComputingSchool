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

# coefficients for 2d normalized kernel
KERNEL_COEFF_1 = (5.0 / 7 * 2.546479089470)
KERNEL_COEFF_2 = (5.0 / 7 * 15.278874536822)
KERNEL_COEFF_5 = (5.0 / 7 * 5.092958178941)

halo_filename = simulation_directory + "/fof_tab_%03d.hdf5" % snapshot 
particle_filename = simulation_directory + "/snapshot_%03d.hdf5" % snapshot 

## large scale structure visualization
try:
    data = h5py.File(particle_filename, "r")
except:
    print("could not open "+particle_filename)
pos = np.array(data["PartType1"]["Coordinates"], dtype=np.float64)
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
k=16
dist, _ = searchtree.query(pos, k=k)
hsml = dist[:,-1]

# project on a canvas
Nplot=512
Edges1d = np.linspace(0.0, Boxsize, Nplot+1, endpoint=True, dtype=np.float64)
Grid1d = 0.5 * (Edges1d[1:] + Edges1d[:-1])
dx = Boxsize/np.float64(Nplot)
Canvas = np.zeros((Nplot, Nplot), dtype=np.float64)
Norm = 0.0

for i in np.arange(32**3):

    jx = np.int32((pos[i,0]-hsml[i])/dx - 1)
    mx = np.int32((pos[i,0]+hsml[i])/dx + 1)

    while jx < mx:
        
        jy = np.int32((pos[i,1]-hsml[i])/dx - 1)
        my = np.int32((pos[i,1]+hsml[i])/dx + 1)
        while jy < my:
            xx = pos[i,0] - jx*dx
            yy = pos[i,1] - jy*dx
            r = np.sqrt(xx*xx + yy*yy)
            u = r/hsml[i]
            hinv3 = 1./hsml[i]/hsml[i]/hsml[i]
            wk = 0
            if u < 0.5:
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u)
            else:
                if u < 1.:
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u)

            if jx>=0 and jy>=0 and jx<Nplot and jy<Nplot: # note: no periodic wrapping of canvas yet
                Canvas[jx, jy] += mass * wk
                Norm += wk
            jy += 1
        jx += 1

Canvas[...] *= 32**3 / Norm / dx * h / dx * h

pc = ax.imshow(Canvas.T, cmap=plt.get_cmap("YlOrBr"), extent=(0, Boxsize, 0, Boxsize), origin="lower", norm=LogNorm(vmin=1e-5*np.max(Canvas), vmax=np.max(Canvas)))

plt.colorbar(pc, cax=cax)
cax.set_ylabel(r"projected density [M$_\odot$ Mpc$^{-2}$]")

ax.set_xlim([0,Boxsize])
ax.set_ylim([0,Boxsize])
ax.set_xlabel('[Mpc/h]')
ax.set_ylabel('[Mpc/h]')

if not os.path.exists( simulation_directory+"/plots" ):
    os.mkdir( simulation_directory+"/plots" )
fig.savefig(simulation_directory+'/plots/DensityProjection_%03d.png'%snapshot, dpi=300)
plt.close(fig)
