# create_ics.py
# 
# initial condition creation for a cosmological spherical collapse problem

import numpy as np
import matplotlib.pyplot as plt
import h5py

FloatType = np.float64
IntType = np.int32

def sample_background_grid(n_bg_grid, simulation_box):
    pos = []
    for ix in np.arange(n_bg_grid):
        posx = ( (ix+0.5)/n_bg_grid - 0.5 ) * simulation_box 

        for iy in np.arange(n_bg_grid):
            posy = ( (iy+0.5)/n_bg_grid - 0.5 ) * simulation_box

            for iz in np.arange(n_bg_grid):
                posz = ( (iz+0.5)/n_bg_grid - 0.5 ) * simulation_box

                pos.append([posx, posy, posz])
    return np.array(pos)

def sample_sphere(Mass, Radius, m_particle):
    #list of particle positions in spherical coordinates
    r_particle = []
    phi_particle = []
    costheta_particle = []
    cum_mass_particles = 0

    radii = np.linspace(0,Radius,100)
    menc = (radii/Radius)**3 * Mass # top-hat density profile

    for i, r in enumerate(radii):
        while cum_mass_particles < menc[i]:
            r_particle.append(r)
            phi_particle.append(np.random.uniform(0,2*np.pi))
            costheta_particle.append(np.random.uniform(-1,1))
            cum_mass_particles += m_particle
    r_particle = np.array(r_particle)
    phi_particle = np.array(phi_particle)
    theta_particle = np.arccos(np.array(costheta_particle))

    return np.array([r * np.sin(phi_particle) * np.cos(theta_particle), r * np.cos(phi_particle) * np.cos(theta_particle), r * np.sin(theta_particle)], dtype=np.float64)


#constants and units
KM_PER_S = FloatType(1e5) # km/s in cm/s
SOLAR_MASS = FloatType(1.989e+33)
MEGAPARSEC = FloatType(3.085678e+24) # Mpc in cm
GRAVITY = FloatType(6.6738e-08) # gravitational constant in cgs

Omega0 = 0.308
overdensity = FloatType(0.67) # 0.67 # 1
h = FloatType(0.7)
hubble = FloatType(100) * h * KM_PER_S / MEGAPARSEC
rho_crit = 3. * hubble**2 / 8. / np.pi / GRAVITY
a_start = 0.0625 # z=15 

simulation_box = FloatType(50.) #100 Mpc/h in units of Mpc/h
simulation_volume = (simulation_box * MEGAPARSEC / h)**3
mass_in_box = Omega0 * rho_crit / (1e10*SOLAR_MASS/h) * simulation_volume

print("mass in box in Msun/h: %g"%(mass_in_box * 1e10))

radius_sphere = FloatType(15.) # in units of comoving Mpc/h
sphere_volume = 4./3. * np.pi * (radius_sphere * MEGAPARSEC / h)**3
H_start = hubble * np.sqrt(Omega0/a_start**3 + (1-Omega0))
rho_crit_start = 3. * H_start**2 / 8. / np.pi / GRAVITY * a_start**3 # comoving density
mass_sphere = Omega0 * overdensity * rho_crit_start / (1e10*SOLAR_MASS/h) * sphere_volume 

print("mass in sphere in Msun/h %g"%(mass_sphere * 1e10))

# set up a background grid
n_bg_grid = 32
n_bg_part = n_bg_grid**3
m_bg_part = (mass_in_box-mass_sphere)/n_bg_part

pos2 = sample_background_grid(n_bg_grid, simulation_box)
pos2 += 0.5 * simulation_box
vel2 = np.zeros(pos2.shape, dtype=FloatType)

# set up a sphere with desired overdensity
pos = sample_sphere(mass_sphere, radius_sphere, m_bg_part).T
pos += 0.5 * simulation_box
n_part = pos.shape[0]
if n_part > 0:
    m_part = mass_sphere/n_part
else:
    m_part = 0.0

print("mass part type 1", m_part, "number of particles ", n_part)
print("mass part type 2", m_bg_part, "number of particles ", n_bg_part)

mass_table = np.array([0, m_part, m_bg_part], dtype=FloatType)
vel = np.zeros(pos.shape, dtype=FloatType)

# perturbation
vel += np.random.normal(loc=0, scale = 30., size=pos.shape)

IC = h5py.File('./ICs_spherical_collapse.hdf5', 'w')

## create hdf5 groups
header = IC.create_group("Header")
part1 = IC.create_group("PartType1")
part2 = IC.create_group("PartType2")

## header entries
NumPart = np.array([0, n_part, n_bg_part], dtype = IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(3, dtype = IntType) )
header.attrs.create("MassTable", mass_table )
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", simulation_box)
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
header.attrs.create("Flag_DoublePrecision", 1)

## part type 1
part1.create_dataset("ParticleIDs", data = np.arange(1, n_part+1) )
part1.create_dataset("Coordinates", data = pos)
part1.create_dataset("Velocities", data = vel)
## part type 2
part2.create_dataset("ParticleIDs", data = np.arange(n_part+1, n_part+1+n_bg_part) )
part2.create_dataset("Coordinates", data = pos2)
part2.create_dataset("Velocities", data = vel2)

## close file
IC.close()

