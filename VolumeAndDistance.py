import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import scipy.integrate as integrate

def get_cubic_volume(d, r):
	'''
	Return the volume of a (possibly) truncated hemisphere on a cubicly arranged surface,
	i.e. where the hemisphere has four nearest neighbors
	'''
	
	# Get volume of isolated hemisphere
	V0 = 2*np.pi*r**3 / 3
	
	# If no overlap between hemispheres, i.e. the hemispheres are isolated
	if d >= 2*r:
		return V0
	
	# If only two spheres are overlapping at a time
	# (i.e. the distance between the hemispheres are greater than sqrt(2)*r;
	# at smaller distances also the neighboring hemispheres starts to overlapm and the calculation,
	# becomes different)
	elif d >= (2**0.5*r):
	
		# Get half the volumetric overlap between two hemispheres
		S = ((np.pi/24)*(d**3 - 12*d*r**2) + V0) / 2
		
		# Subtract the half overlaps from the four neighboring hemispheres
		# from the isolated hemisphere volume
		return V0 - 4*S
	
	# If more than two spheres are overlapping at a time
	# (i.e. each hemisphere has a square-shaped base area)
	else:
	
		# Do a double numerical integration to get the volume of the truncated hemisphere
		result = integrate.dblquad(lambda y, x:(r**2 - x**2 - y**2)**0.5, -d/2, d/2, -d/2, d/2)
		return result[0]

# Set random seed for reproducibility
np.random.seed(2727)

# Specify radius of hemisphere
r = 1.

# Specify list of distances between hemispheres
step = 0.1
ds = np.arange(0.1, 3.0+step, step)*r

# Initiate array of volumes
Vs = np.zeros_like(ds)

# Initiate array of density of particles on a surface which is 10x10 (j_d(d/r))
rhos = np.zeros_like(ds)

# Iterate through distances between hemispheres
for d_idx, d in enumerate(ds):

	# Get volume of (possibly truncated) hemispheres
	Vs[d_idx] = get_cubic_volume(d, r)

    # Get volume densities
	rhos[d_idx] = Vs[d_idx] / (d**2)


# Make plot of volume
fig, ax = plt.subplots()

# Plot volumes vs. distances
ax.plot(ds, Vs, marker='o', lw=0.5, color='black', ms=3.)

# Get volume of isolated hemisphere
V0 = 2*np.pi*r**3 / 3

# Set x axis limits
xlim = ax.get_xlim()
xlim = (0., xlim[1])
ax.set_xlim(xlim)

# Set y axis minimum at zero
ylim = ax.get_ylim()
ylim = (0., 2.5)
ax.set_ylim(ylim)

# Plot volume of isolated hemisphere
ax.plot(xlim, [V0]*2, ls='dashed', lw=0.7, color='pink', zorder=0)
ax.set_xlim(xlim)

# Plot lines where the function expression changes
ylim = ax.get_ylim()
ax.plot([2*r]*2, [0., ylim[1]], ls='dashed', lw=0.7, color='orange')
ax.plot([2**0.5*r]*2, [0., ylim[1]], ls='dashed', lw=0.7, color='orange')
ax.set_ylim(ylim)

# Set axis labels
ax.set_xlabel('d / r')
ax.set_ylabel('Volume of a hemisphere')
ax.set_title('Volume of a Hemisphere vs d / r')

plt.show()

# Make plot of volume density
fig, ax = plt.subplots()

# Plot volumes vs. distances
ax.plot(ds, rhos, marker='o', lw=0.5, color='black', ms=3.)

# Get volume of isolated hemisphere
V0 = 2*np.pi*r**3 / 3

# Set x axis limits
xlim = ax.get_xlim()
xlim = (0., xlim[1])
ax.set_xlim(xlim)

# Set y axis minimum at zero
ylim = ax.get_ylim()
ylim = (0., 1.05)
ax.set_ylim(ylim)

# # Plot volume of isolated hemisphere
# ax.plot(xlim, [V0]*2, ls='dashed', lw=0.7, color='pink', zorder=0)
# ax.set_xlim(xlim)

# Plot lines where the function expression changes
ylim = ax.get_ylim()
ax.plot([2*r]*2, [0., ylim[1]], ls='dashed', lw=0.7, color='orange')
ax.plot([2**0.5*r]*2, [0., ylim[1]], ls='dashed', lw=0.7, color='orange')
ax.set_ylim(ylim)

# Set axis labels
ax.set_xlabel('d / r')
ax.set_ylabel('Volume density of a hemisphere')
ax.set_title('Volume Density of a Hemisphere vs d / r')

plt.show()