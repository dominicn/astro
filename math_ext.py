import math
import numpy as np

def magnitude(vec):
	return np.sqrt(vec.dot(vec))
	
def normalise(vec):
	return vec / magnitude(vec)

def rotate_about_x_axis(vec, alpha):
	return np.array([
		[1, 0, 0],
		[0, math.cos(alpha), math.sin(alpha)],
		[0, -math.sin(alpha), math.cos(alpha)]
		]).dot(vec)
		
def rotate_about_y_axis(vec, beta):
  return np.array([
		[math.cos(beta), 0, -math.sin(beta)],
		[0, 1, 0],
		[math.sin(beta), 0, math.cos(beta)]
		]).dot(vec)
		
def rotate_about_z_axis(vec, gamma):
	return np.array([
		[math.cos(gamma), math.sin(gamma), 0],
		[-math.sin(gamma), math.cos(gamma), 0],
		[0, 0, 1]
		]).dot(vec)
		
def rotate_about(vec, axis_vec, theta):
	# Decompose vec into a component parallel to axis_vec and a
	# component orthogonal to axis_vec.
	par_vec = (vec.dot(axis_vec) / axis_vec.dot(axis_vec)) * axis_vec
	orth_vec = vec - par_vec
	
	# Find a vector orthogonal to both par_vec and orth_vec.
	# This is effectively the third vector in a coordinate set one
	# axis of which is the axis we want to rotate around.
	orth2_vec = np.cross(axis_vec, orth_vec)
	
	# The problem now looks like the 2d rotation problem.
	# We're going to rotate the component of vec which is orthogonal
	# to axis in the plane which is orthogonal to axis.
	rotated_orth_vec = magnitude(orth_vec) * (
		math.cos(theta) * orth_vec / magnitude(orth_vec)
		+
		math.sin(theta) * orth2_vec / magnitude(orth2_vec))
		
	# We can now add back in the component of vec which is parallel
	# to axis - it has been unchanged by the rotation.
	return rotated_orth_vec + par_vec

