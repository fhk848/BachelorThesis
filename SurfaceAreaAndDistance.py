import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

def get_cubic_area(d, r):
	'''
	Return the area of a (possibly) truncated hemisphere on a cubicly arranged surface,
	i.e. where the hemisphere has four nearest neighbors
	'''
	
	# If no overlap between hemispheres, i.e. the hemispheres are isolated
	if d >= 2*r:
		return 2*np.pi*r**2
	
	# If spheres are overlapping
	# Use a sampling algorithm to approximate the surface area by picking random points
	# from the surface of a sphere, and counting the ratio of points within the desired
	# square area
	n_rand = int(1e6)
	rand1 = np.random.random_sample(size=n_rand)
	rand2 = np.random.random_sample(size=n_rand)
	
	# Get uniformly distributed spherical polar angles
	# see https://mathworld.wolfram.com/SpherePointPicking.html
	phis = 2*np.pi*rand1 # 0...2*pi
	thetas = np.arccos(2*rand2 - 1) # 0...pi
	
	# Get corresponding x and y Cartesian coordinates
	xs = r*np.sin(thetas)*np.cos(phis)
	ys = r*np.sin(thetas)*np.sin(phis)
	
	# Count the ratio of x and y coordinates that are within the desired square in the xy-plane,
	# i.e., -d/2 < x < d/2, and -d/2 < y < d/2, 
	n_within = sum((-d/2 < xs) * (xs < d/2) * (-d/2 < ys) * (ys < d/2))
	ratio = n_within / n_rand
	return 0.5*(ratio * 4*np.pi*r**2)

# Set random seed for reproducibility
np.random.seed(2727)

# Specify radius of hemisphere
r = 1.

# Specify list of distances between hemispheres
step = 0.05
ds = np.arange(0.1, 3.00, step)*r

# Initiate array of areas
As = np.zeros_like(ds)

# Initiate list of density 
rhos = np.zeros_like(ds)

# Iterate through distances between hemispheres
for d_idx, d in enumerate(ds):

	# Get area of (possibly truncated) hemispheres
	As[d_idx] = get_cubic_area(d, r)

    # Get area densities
	rhos[d_idx] = As[d_idx] / (d**2)

# Make plot for surface area
fig, ax = plt.subplots()

# Plot volumes vs. distances
ax.plot(ds, As, marker='o', lw=0.5, color='black', ms=3.)

# Get area of isolated hemisphere
A0 = 2*np.pi*r**2

# Set x axis limits
xlim = ax.get_xlim()
xlim = (0., xlim[1])
ax.set_xlim(xlim)

# Set y axis minimum at zero
ylim = ax.get_ylim()
ylim = (0., 8.0)
ax.set_ylim(ylim)

# Plot area of isolated hemisphere
ax.plot(xlim, [A0]*2, ls='dashed', lw=0.7, color='pink', zorder=0)
ax.set_xlim(xlim)

# Plot lines where the function expression changes
ylim = ax.get_ylim()
ax.plot([2*r]*2, [0., ylim[1]], ls='dashed', lw=0.7, color='orange')
ax.plot([2**0.5*r]*2, [0., ylim[1]], ls='dashed', lw=0.7, color='orange')
ax.set_ylim(ylim)

# Set axis labels
ax.set_xlabel('d / r')
ax.set_ylabel('Surface area of a hemisphere')
ax.set_title('Surface Area of a Hemisphere vs d / r')

plt.show()


# Make plot for surface area densities
fig, ax = plt.subplots()

# Plot volumes vs. distances
ax.plot(ds, rhos, marker='o', lw=0.5, color='black', ms=3.)

# Get area of isolated hemisphere
A0 = 2*np.pi*r**2

# Set x axis limits
xlim = ax.get_xlim()
xlim = (0., xlim[1])
ax.set_xlim(xlim)

# Set y axis minimum at zero
ylim = ax.get_ylim()
ylim = (0., 1.8)
ax.set_ylim(ylim)

# # Plot area of isolated hemisphere
# ax.plot(xlim, [A0]*2, ls='dashed', lw=0.7, color='darkorange', zorder=0)
# ax.set_xlim(xlim)

# Plot lines where the function expression changes
ylim = ax.get_ylim()
ax.plot([2*r]*2, [0., ylim[1]], ls='dashed', lw=0.7, color='orange')
ax.plot([2**0.5*r]*2, [0., ylim[1]], ls='dashed', lw=0.7, color='orange')
ax.set_ylim(ylim)

# Set axis labels
ax.set_xlabel('d / r')
ax.set_ylabel('Surface area density of a hemisphere')
ax.set_title('Surface Area Density of a Hemisphere vs d / r')

plt.show()