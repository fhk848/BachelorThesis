import numpy as np

class RandomWalk():

	def __init__(self, size, hole_xys, hole_half_width, temp, m_particle):
		'''	
		Parameters
		----------
			size				(dx, dy, dz) size of box
			hole_xys			[(x0, y0), (x1, y1), ...] xy-coordinates of hole centers
			hole_half_width		half width of the holes
			temp				temperature (unit: K)
			m_particle			mass of particle (unit: kg)
		'''
		self.size = size
		self.hole_xys = hole_xys
		self.hole_half_width = hole_half_width
		self.temp = temp
		self.m_particle = m_particle
		

	def set_initial_position(self, r):
		'''
	    Parameters
		----------
			r_init		(x, y, z) initial position of particle
		'''
		self.r = r


	def take_step(self, std):
		'''
        Walk randomly in 3D until stopping criterion is met

		Return
		------
			-1			if the particle diffuses to a position above the box
			0..N		the index of one of the N holes that the particle diffuses to
		'''
		# Take random step in x, y, z with step sizes
		# picked from a Maxwell-Boltzmann distribution	
		r_step = np.random.normal(loc = 0, scale = std, size = 3)
	
		# Update position
		self.r += r_step
		
		# Stop if particle is above box
		if self.r[2] > self.size[2]:
			return -1
		
		# Make position periodic in x and y	
		self.r[:2] = self.r[:2] % self.size[:2]
		
		# === If particle hit a hole === #
		# If below box
		if self.r[2] < 0.:
	
			# Iterate through hole coordinates
			for hole_idx, (x, y) in enumerate(self.hole_xys):
		
				# If particle is within x coordinates
				if x - self.hole_half_width < self.r[0] < x + self.hole_half_width:
			
					# If particle is within y coordinates
					if y - self.hole_half_width < self.r[1] < y + self.hole_half_width:
					
                        # Return hole index
						return hole_idx
		
			# If not hole was hit, then bounce off the bottom
			# by mirroring the z coordinate to above the bottom
			self.r[2] *= -1.
		
			return None

def set_axis(ax, size, hole_xys, hole_half_width):
	'''
	Prepare axis with preset parameters
	'''	
	# Clear axis of its content
	ax.clear()
	
	# Show box edges
	xs = np.array([0., size[0]])
	ys = np.array([0., size[1]])
	zs = np.array([0., size[2]])
	
	# Set shared keywords
	surf_kw = dict(color='green', alpha=0.1, edgecolors='black')
	
	for x in xs:
		xs_ = np.zeros((2,2)) + x
		ys_ = np.repeat(ys, 2).reshape(2,2)
		zs_ = np.tile(zs, 2).reshape(2,2)
		ax.plot_surface(xs_, ys_, zs_, **surf_kw)
	
	for y in ys:
		xs_ = np.repeat(xs, 2).reshape(2,2)
		ys_ = np.zeros((2,2)) + y
		zs_ = np.tile(zs, 2).reshape(2,2)
		ax.plot_surface(xs_, ys_, zs_, **surf_kw)
	
	for z in zs:
		xs_ = np.repeat(xs, 2).reshape(2,2)
		ys_ = np.tile(ys, 2).reshape(2,2)
		zs_ = np.zeros((2,2)) + z
		ax.plot_surface(xs_, ys_, zs_, **surf_kw)
	
	# Show holes
	for x, y in hole_xys:
		hole_xs = x + np.array([1, -1]) * hole_half_width
		hole_ys = y + np.array([1, -1]) * hole_half_width
		
		xs_ = np.repeat(hole_xs, 2).reshape((2,2))
		ys_ = np.tile(hole_ys, 2).reshape((2,2))
		zs_ = np.zeros((2,2))
		ax.plot_surface(xs_, ys_, zs_, color='tomato', alpha=0.5)
	
	# Remove spines, ticks, etc.
	ax.set_axis_off()
	
	# Set equal aspect ratio
	ax.set_box_aspect((size[0], size[1], size[2]))
	ax.set_xlim3d([0., size[0]])
	ax.set_ylim3d([0., size[1]])
	ax.set_zlim3d([0., size[2]])
	
	# Set camera view
	ax.view_init(elev=10., azim=-70.)

def get_hole_xys(half_width, box_size):
	'''
	Return hole center coordinates for the possible regular hole spacings
	given the box dimensions and hole half width
	'''
	width = 2*half_width
	max_holes_x = box_size[0] / width
	max_holes_y = box_size[1] / width
	
	for n_holes_x in np.arange(max_holes_x, -1, step=-1).astype(int):
		
		if n_holes_x == 0:
			yield np.array([])
		
		elif n_holes_x == 1:
			yield np.array([[box_size[0]/2, box_size[1]/2]])
			
		else:
			
			# Get spacing between holes in each dimension
			spacing_x = box_size[0] / n_holes_x
			spacing_y = box_size[1] / n_holes_x
			
			# Get coordinates
			xs = np.arange(spacing_x/2, box_size[0]+spacing_x/2, spacing_x)
			ys = np.arange(spacing_y/2, box_size[1]+spacing_y/2, spacing_y)
			
			# Repeat to get coordinates on a grid
			xs = np.repeat(xs, n_holes_x).reshape(-1, 1)
			ys = np.tile(ys, n_holes_x).reshape(-1, 1)
			
			# Yield grid of hole coordinates
			yield np.concatenate((xs, ys), axis=1)

def get_std(temp, m_particle):
	'''
	Return standard deviations given the temperature of the system
	and the mass of the particles (unit: Å/fs)
	'''

	# Boltzmann constant (unit: J/K)
	k_B = 1.380649e-23

	alpha = temp / m_particle
	std = float(((k_B * alpha)**(1/2))/1e5)
		

	return std
