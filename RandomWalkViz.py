import numpy as np
import matplotlib.pyplot as plt
from RandomWalkViz1 import get_hole_xys, set_axis


#=== Inputs and outputs from RandomWalkSim ===#

# Set random seed for reproducibility
np.random.seed(2727)

# Specify size of box in Angstroms
size = np.array([0.1, 0.1, 0.05])

# Specify location of hole centers
hole_xys = np.array([
	[0.01, 0.01],
	[0.01, 0.03],
	])

# Specify half hole width
hole_half_width = 0.005

# Specify number of particles
n_particles = int(1e6)

# Specify mass of particle (unit: kg)
# In this case it's dioxygen
m_particle = 15.999*1.66053907e-27

# Specify temperature as a list (unit: K)  
temp = 300

# Flux data from RandomWalkSim with above specifications

PHI = np.array([50040100.0, 60469876.543209866, 74285156.24999999, 93019591.83673468, 118307222.22222224, 153114000.0, 199537500.0, 256761111.11111107, 318817500.0, 354960000.0])


# The effective volume of reactants which is being used per time unit
# dVdt = PHI * A_cat / rho
alpha = (2*hole_half_width)**2 / (n_particles / np.prod(size))
dVdt = np.multiply(PHI, alpha).astype(float)

print(dVdt)

# Initiate list for the following data:
# Radii for sphere when n_holes = 1 (rs_sph) and cylinders when n_holes != 1 (rs_cyl)
# Heights for cylinder (hs_cyl) and spherical cap/ellipsoid (hs_sc) where rs_sph = hs_cyl + hs_sc
rs_sph, rs_cyl, hs_cyl, hs_sc = [], [], [], []


# Get radius for half sphere when only one catalyst is present
r_sph = float((3*dVdt[-1]/(2*np.pi)) ** (1/3))


# Make figure
fig = plt.figure(figsize = (14, 6))

# Iterate through hole spacings
for spacing_idx, hole_xys in enumerate(get_hole_xys(hole_half_width, size)):
    
    # Get the number of holes
    n_holes = hole_xys.shape[0]

    if n_holes == 0:
            continue

    # Make a half sphere as representation of effective volume if only one hole
    elif n_holes == 1:
        # Make axis object
        ax = fig.add_subplot(2, 5, spacing_idx + 1, projection='3d')

        # Set axis with presets
        set_axis(ax, size, hole_xys, hole_half_width, r_sph, 0, 0, 0)
        plt.title(f'{n_holes} hole', fontsize = 8)

    # Make cut off half spheres as representation of effective volume if more than one hole
    else:
        # Find the radius of the cylinder given the effective volume 
        # and the radius of the half sphere when only one catalyst is present
        # at that temperature
        r_cyl = np.sqrt(- 4**(2/3) * ((2 * np.pi * r_sph**3 - 3 * dVdt[spacing_idx]) * np.pi**2)**(2/3) + 4 * r_sph**2 * np.pi**2) / (2 * np.pi)
        rs_cyl.append(r_cyl)

        # Find the height of the cylinder via pythagorean theorem
        h_cyl = np.sqrt(r_sph**2 - r_cyl**2)
        hs_cyl.append(h_cyl)

        # Find height of ellipse/spherical cap
        h_sc = r_sph - h_cyl
        hs_sc.append(h_sc)

        # Make axis object
        ax = fig.add_subplot(2, 5, spacing_idx + 1, projection='3d')

        # Set axis with presets
        set_axis(ax, size, hole_xys, hole_half_width, r_sph, r_cyl, h_sc, h_cyl)
        plt.title(f'{n_holes} holes', fontsize = 8)


plt.suptitle(f'Visualization of effective volume of reactants used per time unit')
plt.show()


# ==== Collected data === # 
print(rs_cyl)
print(hs_cyl)
print(hs_sc)