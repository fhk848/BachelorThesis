from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations


def set_axis(ax, size, hole_xys, hole_half_width, radius, radiusC, radiusSC, heightC):
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
	
	# Show holes and effective volume of reactants used
	for x, y in hole_xys:
		hole_xs = x + np.array([1, -1]) * hole_half_width
		hole_ys = y + np.array([1, -1]) * hole_half_width
		
		xs_ = np.repeat(hole_xs, 2).reshape((2,2))
		ys_ = np.tile(hole_ys, 2).reshape((2,2))
		zs_ = np.zeros((2,2))
		ax.plot_surface(xs_, ys_, zs_, color='tomato', alpha = 0.5)


		if heightC == 0:
			# Show spheres as a representation of effective volume
			# if only one hole (no interference with other catalysts)
			Xsp, Ysp, Zsp = data_for_sphere(x, y, radius)
			ax.plot_surface(Xsp, Ysp, Zsp, color="blue", alpha = 0.2)

		else:
			# Show cylinders/"trunkerede" spheres as a representation 
			# of effective volume if more than one hole
			Xc, Yc, Zc = data_for_cylinder(x, y, radiusC, heightC)
			ax.plot_surface(Xc, Yc, Zc, color='blue', alpha = 0.2)

			# Spherical cap that starts where the cylinder ends in the z-direction
			Xsc, Ysc, Zsc = data_for_spherical_cap(x, y, radiusC, radiusSC, heightC)
			ax.plot_surface(Xsc, Ysc, Zsc, color = 'blue', alpha = 0.2)
			
	
	# Remove spines, ticks, etc.
	ax.set_axis_off()
	
	# Set equal aspect ratio
	ax.set_box_aspect((size[0], size[1], size[2]))
	ax.set_xlim3d([0., size[0]])
	ax.set_ylim3d([0., size[1]])
	ax.set_zlim3d([0., size[2]])
	
	# Set camera view
	ax.view_init(elev=35., azim=-70.)

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

def data_for_sphere(center_x, center_y, radius):
	"""
	Get data to draw half spheres
	"""
	u, v = np.mgrid[0:2*np.pi:20j, 0:0.5*np.pi:20j]
	x_grid = center_x + radius*np.cos(u)*np.sin(v)
	y_grid = center_y + radius*np.sin(u)*np.sin(v)
	z_grid = radius*np.cos(v)
	
	return x_grid, y_grid, z_grid

def data_for_spherical_cap(center_x, center_y, radiusC, radiusSC, heightC):
	"""
	Get data to draw the spherical cap on top of the cylinder 
	to make up the cut of sphere
	"""
	u, v = np.mgrid[0:2*np.pi:20j, 0:0.5*np.pi:20j]
	x_grid = center_x + radiusC*np.cos(u)*np.sin(v)
	y_grid = center_y + radiusC*np.sin(u)*np.sin(v)
	z_grid = heightC + radiusSC*np.cos(v)
	
	return x_grid, y_grid, z_grid

def data_for_cylinder(center_x, center_y, radiusC, heightC):
    """
	Get data to draw cylinders
	"""
    z = np.linspace(0, heightC, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radiusC*np.cos(theta_grid) + center_x
    y_grid = radiusC*np.sin(theta_grid) + center_y
    
    return x_grid, y_grid, z_grid
