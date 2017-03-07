# Functions for computing properties of orbits or
# the effects of delta-v on orbits using the Kepler model.
# Class representing an orbit (irrespective of time) using the classical
# orbital elements.

import math
import consts
import math_ext
import orbit_intersections
import numpy as np

def get_period_elliptical_orbit(a, M):
	mu = consts.G * M
	return (2 * math.pi) / math.sqrt(mu) * math.pow(a, 1.5)
	
def get_semi_major_axis_elliptical_orbit(tp, M):
	mu = consts.G * M
	return math.pow((math.sqrt(mu) * tp) / (2 * math.pi), 2/3)
	
def get_speed_circular_orbit(r, M):
	mu = consts.G * M
	return math.sqrt(mu / r)

def get_eccentricy_vec(r, v, mu):
	v_sq_minus_mu_over_r_times_r = (v.dot(v) - mu / math_ext.magnitude(r)) * r
	r_dot_v_times_v = r.dot(v) * v
	return (v_sq_minus_mu_over_r_times_r - r_dot_v_times_v) / mu
	
class Orbit:
	def __init__(self, mu, p, e, i, OMEGA, PI):
		self.mu = float(mu)
		self.p = float(p)
		self.e = float(e)
		self.i = float(i)
		self.OMEGA = float(OMEGA)
		self.PI = float(PI)
		
		# Unit vectors in the perifocal coordinate system can be constructed
		# from the above. Firstly, W is a unit vector in the direction of h,
	  # the angular momentum vector. It is perpendicular to the orbital plane.
		W_vec = np.array([0, 0, 1])
		W_vec = math_ext.rotate_about_x_axis(W_vec, -i)
		if not math.isnan(OMEGA):
			# If inclination is zero or pi then the orbital plane is the equatorial
			# plane, W = K and we have nothing to do here, which is handy because
			# OMEGA will be NaN.
			W_vec = math_ext.rotate_about_z_axis(W_vec, -OMEGA)
		
		# P is a unit vector in the direction of periapsis, i.e. it is a unit
		# vector in the direction of e, unless the orbit is circular, in which
		# case it is a unit vector in the direction of the line of nodes, unless
		# the orbit is also equatorial in which case it is I.
		P_vec = np.array([1, 0, 0])
		if not math.isnan(PI):
			if not math.isnan(OMEGA):
				P_vec = math_ext.rotate_about_z_axis(P_vec, -OMEGA)
				P_vec = math_ext.rotate_about(P_vec, W_vec, PI - OMEGA)
			else:
				P_vec = math_ext.rotate_about_z_axis(P_vec, -PI)
		else:
			if not math.isnan(OMEGA):
				P_vec = math_ext.rotate_about_z_axis(P_vec, -OMEGA)
				
		# Q is perpendicular to both P and W.
		Q_vec = np.cross(P_vec, W_vec)
		
		self.P_vec = P_vec
		self.Q_vec = Q_vec
		self.W_vec = W_vec
		
	def __repr__(self):
		return "Orbit({0}, {1}, {2}, {3}, {4}, {5})".format(
			self.mu, self.p, self.e, self.i, self.OMEGA, self.omega)
			
	def __str__(self):
		return "Orbit(mu={0:.4g}, p={1:.4g}, e={2:.4g}, i={3:.4g}, OMEGA={4:.4g}, PI={5:.4g})".format(
			self.mu, self.p, self.e, self.i, self.OMEGA, self.PI)
	
	def is_coplanar(self, other):
		# Inclination must be the same.
		if math.fabs(self.i - other.i) > 0.0000001:
			return False
		
		# If both orbits are equatorial then we don't need to compare the
		# longitudes of ascending nodes (both would be nan anyway).
		# Better to do this check first though, because if one is equatorial
		# and the other near as dammit we'll count them as coplanar, though
		# a comparisons of the OMEGA values would find one nan and the other not.
		if min(math.fabs(self.i), math.fabs(other.i)) < 0.0000001:
			return True
			
		# Inclination is the same, and it's far enough from equatorial that
		# comparing intersections with the equator is worth doing. OMEGA is the
		# longitude of the *ascending* node. If the two OMEGAs are on opposite
		# sides of the orbited body (i.e. separated by pi) then the orbits are
		# coplanar, just in opposing directions.
		if (math.fabs(self.OMEGA - other.OMEGA) % math.pi) > 0.0000001:
			return False
			
		return True
	
	def get_intersections(self, other):
		if not self.is_coplanar(other):
			# TODO support non-coplanar orbits, it's actually much easier
			return []
		
		if math.isnan(self.PI) and math.isnan(other.PI):
			# Two coplanar circular orbits either intersect everywhere (because
			# they are the same orbit) or nowhere. Finding points of intersection
			# isn't something we can do either way.
			return []
			
		if math.isnan(other.PI):
			# If the other orbit is circular then any value of k would give us
			# the same answer. Might as well choose zero.
			k = 0
		elif math.isnan(self.PI):
			# If this orbit is circular, we still want to know the rotation of
			# the other orbit relative to the point where we're measuring angles
			# from so that the intersection angles we return (which are relative
			# to this orbit) are correct.
			k = other.PI
		else:
			k = other.PI - self.PI
			
		return orbit_intersections.get_intersections(
			self.p,
			self.e,
			other.p,
			other.e,
			k)
		
	# nu is the angle from periapsis
	def get_r(self, nu):
		return self.p / (1 + self.e * math.cos(nu))
		
	# nu is the angle from periapsis
	def get_r_vec(self, nu):
		r = self.get_r(nu)
		return (
			r * math.cos(nu) * self.P_vec
			+ r * math.sin(nu) * self.Q_vec
			)
			
	def get_r_vecs(self, from_nu=0.0, to_nu=math.pi*2, samples=100):
		result = []
		for i in range(0, samples):
			nu = (from_nu - to_nu) / (samples - 1) * i
			r_vec = self.get_r_vec(nu)
			result.append(r_vec)
		return result
		
	@classmethod
	def from_r_v_mu(cls, r_vec, v_vec, mu):
		h_vec = np.cross(r_vec, v_vec)
		h = math_ext.magnitude(h_vec)
		n_vec = np.cross(np.array([0, 0, 1]), h_vec)
		n = math_ext.magnitude(n_vec)
		e_vec = get_eccentricy_vec(r_vec, v_vec, mu)
		p = np.square(h) / mu
		e = math_ext.magnitude(e_vec)
		i = math.acos(h_vec[2] / h)
		if n == 0:
		  # n_vec, the line of nodes, is the cross product of h_vec (which is
		  # perpendicular to the orbital plane) and K, which is perpendicular
		  # to the equatorial plane. If n is zero that means these two vectors
		  # are parallel, which means the orbital plane IS the equatorial plane
		  # and therefore there is no line of nodes.
		  OMEGA = math.nan
		  if e != 0:
		  	# If e is non-zero then PI is the angle from I to e_vec, which
		  	# points to periapsis.
		  	PI = math.acos(e_vec[0] / e)
		  else:
		  	# BUT e can also be zero. The magnitude of e_vec, e, is equal to the
		  	# eccentricity of the orbit. If the orbit is circular that's zero,
		  	# which makes sense because there is no periapsis to point to.
		  	# In the case where both n and e are zero the orbit is circular and
		  	# equatorial and PI is undefined.
		  	PI = math.nan
		else:
		 	OMEGA = math.acos(n_vec[0] / n)
		 	if n_vec[1] < 0:
		 		OMEGA = 2 * math.pi - OMEGA
		 	if e != 0:
		 		# The most general case, where neither e nor n are zero:
		 		# PI is then defined somewhat oddly as the angle from I to the line
		 		# of nodes (in the equatorial plane) plus the angle from the line of
		 		# nodes to e (which is in the orbital plane).
		 		omega = math.acos(n_vec.dot(e_vec) / (n * e))
		 		if e_vec[2] < 0:
		 			omega = 2 * math.pi - omega
		 		PI = omega + OMEGA
		 	else:
		 		PI = math.nan

		return cls(mu, p, e, i, OMEGA, PI)
