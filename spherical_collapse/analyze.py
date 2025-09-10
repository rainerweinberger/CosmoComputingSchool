import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib.animation as animation
import sys
import os

FloatType = np.float64

simulation_directory = str(sys.argv[1])

if not os.path.exists(simulation_directory + '/plots'):
    os.makedirs(simulation_directory + '/plots')

N_snap = 31

for i_snap in np.arange(N_snap):
    with h5py.File(simulation_directory + "/output/snapshot_%03d.hdf5"%i_snap) as snap:
        pos = np.array(snap["PartType1"]["Coordinates"], dtype=FloatType)
        vel = np.array(snap["PartType1"]["Velocities"], dtype=FloatType)

        fig = plt.figure(figsize=(5.0,5.0))
        ax = plt.axes((0.2,0.2,0.78,0.78))

        ax.scatter(pos[:,0], pos[:,1], s=2, c='k')

        ax.set_xlim((45, 55))
        ax.set_ylim((45, 55))
        #ax.set_ylim((-1500, 1500))

        fig.savefig(simulation_directory + "/plots/pos_xy_%03d.png"%i_snap)

        plt.close(fig)


# animation
with h5py.File(simulation_directory + "/output/snapshot_%03d.hdf5"%0) as snap:
    a = snap["Header"].attrs["Time"]
    boxhalf = 0.5 * np.float32(snap["Parameters"].attrs["BoxSize"])
    pos = np.array(snap["PartType1"]["Coordinates"], dtype=FloatType) - boxhalf
    #print(a)
    pos *= a
    vel = np.array(snap["PartType1"]["Velocities"], dtype=FloatType)

    fig = plt.figure(figsize=(5.0,5.0))
    ax = plt.axes((0.2,0.2,0.78,0.78))

    scatter = ax.scatter(pos[:,0], pos[:,1], s=1, c='k')

    ax.set_xlim((-5, 5))
    ax.set_ylim((-5, 5))


def update(i_snap):
    with h5py.File(simulation_directory + "/output/snapshot_%03d.hdf5"%i_snap) as snap:
        a = snap["Header"].attrs["Time"]
        boxhalf = 0.5 * np.float32(snap["Parameters"].attrs["BoxSize"])
        pos = np.array(snap["PartType1"]["Coordinates"], dtype=FloatType) - boxhalf
        #print(a)
        pos *= a
        vel = np.array(snap["PartType1"]["Velocities"], dtype=FloatType)


    scatter.set_offsets(pos[:,0:2])
    return 


ani = animation.FuncAnimation(fig=fig, func=update, frames=N_snap, interval=100)

ani.save(simulation_directory + "/plots/collapse.gif")
plt.show()
