import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import math
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


def get_j_k(delta_G_bind, U, temp):
	# Define Boltzman constant in eV/K
	k_B = 8.617333262 * 10 ** (-5)

	# Define elementary charge in C
	# Not multiply with e as 1 V is equal eV
	e = 1.602176634 * 10 ** (-19)

	# Calculate j_k 
	j_k = float(np.exp(-(np.absolute(delta_G_bind)-(U))/(k_B * temp)))

	return j_k


def get_j_tot(j_k, j_d, a):
	
	# Calculate j_tot
	j_tot = (1/((1/(a*j_k)) + (1/j_d)))

	return j_tot

# Set random seed for reproducibility
np.random.seed(2727)

# Specify factor a (a*j_k) so j_k and j_d are comparable
a = 1

# Specify temperature in K
temp = 300

# Specify Gibb's free energy when particle binds to catalyst in eV
delta_G_binds = np.arange(-0.5, 0.51, 0.01)

# Specify potential in V
Us = np.arange(-0.2, 0.22, step = 0.02)

# Specify radius of hemisphere
r = 1

# Specify list of distances between hemispheres
step = 0.2
ds = np.arange(0.1, 3.0+step, step)*r

# Initiate array of volumes
Vs = np.zeros_like(ds)

# ==== j_tots as a function of ∆G_bind and V(d/r) ==== #

# Initiate array of j_k
j_ks = np.zeros_like(delta_G_binds)

# Initiate array of j_total
j_tots = []

# Iterate through distances between hemispheres
for d_idx, d in enumerate(ds):

	# Get volume of (possibly truncated) hemispheres
	Vs[d_idx] = get_cubic_volume(d, r)

	# Iterate through potientials
	for G_idx, g in enumerate(delta_G_binds):
		
		# Get j_k given the binding energy ∆G_bind
		# In this case U = 0 V and temp = 300 K
		j_ks[G_idx] = get_j_k(g, 0, temp)

		# Get j_total for given ∆G_bind and d/r
		j_tots = np.append(j_tots, get_j_tot(j_ks[G_idx], Vs[d_idx], a))

# Reshape j_tots
j_tots = np.reshape(j_tots, (len(ds), len(delta_G_binds)))

# Make 3D plot

# Repeat and reshape binding energies so compatible with PolyCollection
delta_G_bindsVert = np.repeat(delta_G_binds, len(ds), axis=0).reshape((len(delta_G_binds), len(ds)))

# Vectorize the coordinates so compatible with PolyCollection
verts = []
for d_idx in range(len(ds)):
	xs = np.concatenate([[delta_G_bindsVert[0, d_idx]], delta_G_bindsVert[:, d_idx], [delta_G_bindsVert[-1, d_idx]]])
	ys = np.concatenate([[0],j_tots[d_idx, :],[0]])
    
	verts.append(list(zip(xs, ys)))

# Set axes and plot
poly = PolyCollection(verts, facecolors = plt.colormaps['viridis_r'](np.linspace(0, 1, len(verts))), alpha = 0.7)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.add_collection3d(poly, zs=ds, zdir='y')

ax.set_xlim3d(delta_G_bindsVert.min(), delta_G_bindsVert.max())
ax.set_xlabel(r'$\Delta G_{i} (eV)$')
ax.set_ylim3d(ds.min(), ds.max())
ax.set_ylabel('d/r')
ax.set_zlim3d(j_tots.min(), j_tots.max())
ax.set_zlabel('Current density')
ax.set_title('Volcano plots')

plt.show()



# ==== j_tots as a function of U and V(d/r) ==== #

# Initiate array of j_total
j_tots = []

# Initiate array of j_k
j_ks = np.zeros_like(Us)

# Iterate through distances between hemispheres
for d_idx, d in enumerate(ds):

	# Get volume of (possibly truncated) hemispheres
	Vs[d_idx] = get_cubic_volume(d, r)

	# Iterate through potientials
	for U_idx, u in enumerate(Us):
		
		# Get j_k given the potential U
		# In this case ∆G_bind = 0 and temp = 300 K
		j_ks[U_idx] = get_j_k(0., u, temp)

		# Get j_total for given U and d/r
		j_tots = np.append(j_tots, get_j_tot(j_ks[U_idx], Vs[d_idx], a))


# Reshape j_tots
j_tots = np.reshape(j_tots, (len(ds), len(Us)))

# Make 2D plot
for d_idx, d in enumerate(ds):
	plt.plot(Us, j_tots[d_idx,:], label = (f'd/r = {d:.1f}'))

plt.title(f'Polarization curves')
plt.legend(loc = 'center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('Potential (V)')
plt.ylabel('Current density')
plt.show()


# Make 3D plot 

# Repeat and reshape binding energies so compatible with PolyCollection
UsVert = np.repeat(Us, len(ds), axis=0).reshape((len(Us), len(ds)))

# Vectorize the coordinates so compatible with PolyCollection
verts = []
for d_idx in range(len(ds)):
	xs = np.concatenate([[UsVert[0, d_idx]], UsVert[:, d_idx], [UsVert[-1, d_idx]]])
	ys = np.concatenate([[0],j_tots[d_idx, :],[0]])
    
	verts.append(list(zip(xs, ys)))

# Set axes and plot
poly = PolyCollection(verts, facecolors = plt.colormaps['viridis_r'](np.linspace(0, 1, len(verts))), alpha = 0.7)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.add_collection3d(poly, zs=ds, zdir='y')

ax.set_xlim3d(UsVert.min(), UsVert.max())
ax.set_xlabel('Potential (V)')
ax.set_ylim3d(ds.min(), ds.max())
ax.set_ylabel('d/r')
ax.set_zlim3d(j_tots.min(), j_tots.max()+1)
ax.set_zlabel('Current density')
ax.set_title('Polarization Curves')

plt.show()