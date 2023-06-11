import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def get_j_d(r):
    # j_d is proportional to reaction volume
    # Calculate reaction volume hemisphere given the reaction radius

    j_d = (2*np.pi*r**3)/3

    return j_d

def get_j_k(delta_G_i, U, temp, a):
	# Define Boltzman constant in eV/K
	k_B = 8.617333262 * 10 ** (-5)

	# Define elementary charge in C
	# Not multiply with e as 1 V is equal to 1 eV
	e = 1.602176634 * 10 ** (-19)

	# Calculate j_k 
	j_k = float(a*np.exp(-(np.absolute(delta_G_i)-(U))/(k_B * temp)))

	return j_k

def get_J_i(j_k, j_d):
	
	# Calculate j_tot
	J_i = (1/((1/(j_k)) + (1/j_d)))

	return J_i

# Specify the length of teh quadratic surface
# as number of atoms
l_surface = 100

# Specify temperature in K
temp = 300

# Specify list of reaction radius of active site compared to the atomic radius 
# For now the atomic radius is set to 1 unit 
r = np.arange(1., 11., step = 1)

# Specify a so j_k and j_d are comparable
a = 1

# Specify list of potentials in V
Us = np.arange(0.0, 0.6, step = 0.01)

# Set random seed for reproducibility
np.random.seed(2727)

# Get random array of ∆G_i for a surface of 100x100 atoms 
# chosen by an uniform distrubution within [-1.0 eV, 1.0 eV)
G_i = np.random.uniform(-0.25, 0.25, size = (l_surface, l_surface))

# Initiate list of J_tot
J_tots = []

# Iterate through reaction radii
for r_i in r:

	# Calculate j_d
	j_d = get_j_d(r_i)

	# Iterate through potentials
	for u in Us:

		# Initiate list of J_i
		J_is = []

		# Iterate through each catalytic site	
		for i in range(l_surface):
			for j in range(l_surface):
				# Calculate J_i for each atom
				J_is = np.append(J_is, get_J_i(get_j_k(G_i[i, j], u, temp, a), j_d))

		J_tots = np.append(J_tots, np.sum(J_is))

# Reshape j_tots
J_tots = np.reshape(J_tots, (len(r), len(Us)))

# Make 2D plot
for r_idx, r_i in enumerate(r):
	plt.plot(Us, J_tots[r_idx,:], label = (f'r = {r_i}'))

plt.title(f'Polarization Curves for a Random HEA Surface (Uniform Distribution)')
plt.legend()
plt.xlabel('Potential (V)')
plt.ylabel('Current')
plt.show()

# # Make 3D plot

# # Repeat and reshape binding energies so compatible with PolyCollection
# UsVert = np.repeat(Us, len(r), axis=0).reshape((len(Us), len(r)))

# # Vectorize the coordinates so compatible with PolyCollection
# verts = []
# for r_idx in range(len(r)):
# 	xs = np.concatenate([[UsVert[0, r_idx]], UsVert[:, r_idx], [UsVert[-1, r_idx]]])
# 	ys = np.concatenate([[0],J_tots[r_idx, :],[0]])
    
# 	verts.append(list(zip(xs, ys)))

# # Set axes and plot
# poly = PolyCollection(verts, facecolors = plt.colormaps['viridis_r'](np.linspace(0, 1, len(verts))), alpha = 0.7)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.add_collection3d(poly, zs=r, zdir='y')

# ax.set_xlim3d(UsVert.min(), UsVert.max())
# ax.set_xlabel(r'$U (V)$')
# ax.set_ylim3d(r.min(), r.max())
# ax.set_ylabel('r')
# ax.set_zlim3d(J_tots.min(), J_tots.max())
# ax.set_zlabel(r'$J$')

# plt.show()

