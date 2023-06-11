import numpy as np
import matplotlib.pyplot as plt
from RandomWalkSim1 import RandomWalk, get_hole_xys, get_std, set_axis
from copy import deepcopy

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

for spacing_idx, hole_xys in enumerate(get_hole_xys(hole_half_width, size)):
	pass
n_spacings = spacing_idx + 1

# Initiate list of fluxes through holes
flux_holes = np.zeros(n_spacings-1)

# Iterate through hole spacings
for spacing_idx, hole_xys in enumerate(get_hole_xys(hole_half_width, size)):
	
	# Get the number of holes
	n_holes = hole_xys.shape[0]
	
	if n_holes == 0:
		continue
	#print(f'no. of holes: {n_holes}')
	
	# Initiate counter of the number of particles that hit each hole
	# and the top as the last item
	n_each = np.zeros(n_holes+1, dtype=int)
	
	# Do random walk of particle
	rand_walk = RandomWalk(size, hole_xys, hole_half_width, temp, m_particle)

	# Initiate history of positions
	rs = []

	# Iterate through initial positions
	for particle_idx, r_init in enumerate(np.random.random_sample(size=(n_particles, 3)) * size):

		# Set initial position of particle
		rand_walk.set_initial_position(r_init)
	
		# Set initial point of random walk
		rs.append([deepcopy(rand_walk.r)])
	
		# Do random walk until stopping
		while True:

			# Make particle take a step			
			status = rand_walk.take_step(get_std(temp, m_particle))

			# Append position to history
			rs[particle_idx].append(deepcopy(rand_walk.r))
		
			# Stop random walk if None is not returned
			if status is not None:
				n_each[status] += 1
				break
	
	# Get areas of each hole and the top as the last item
	areas = np.zeros(n_holes+1)
	areas[:-1] = (2*hole_half_width)**2
	areas[-1] = np.prod(size[:2])
		
	# Get flux through each hole/catalyst
	# i.e. the number of particles per time per area
	# the time is the same for all holes/catalysts, so here is compared the number of particles
	# per area
	#flux_each = n_each / areas
		
	# Get the combined flux through the holes/catalysts
	flux_holes[spacing_idx] = np.sum(n_each[:-1]) / np.sum(areas[:-1])
	#flux_holes[spacing_idx] = n_each[0] / areas[0]
	print(f'Flux w/ {n_holes} holes at {temp} K = {flux_holes[spacing_idx]}')

	# # Make figure that shows each particles random walk
	# fig = plt.figure()

	# # Make axis object
	# ax = fig.add_subplot(projection='3d')

	# # Set axis with presets
	# set_axis(ax, size, hole_xys, hole_half_width)

	# # Iterate through particle trajectories
	# for rs0 in rs:
	
	# 	# Make into numpy
	# 	rs0 = np.array(rs0)
	
	# 	# Plot random walk
	# 	ax.plot(*rs0.T, marker='o', markersize=1, lw=1.)

	# plt.show()

dist_index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

plt.bar(dist_index, flux_holes, color ='pink',
        width = 0.4)
 
plt.xlabel("Distance index")
plt.ylabel("Flux (number of particles per time per area)")
plt.title("Flux vs distance index")
plt.show()
