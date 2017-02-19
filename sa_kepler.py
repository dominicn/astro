# Simulation model using the universal variable formulation of
# the analytical solution to the Kepler problem (which is to
# predict the position of an orbiting body given time, using
# the simplifying assumption that the only acceleration we
# need to worry about is the gravitational force attracting
# our body to the one large mass it is orbiting).

import numpy as np
import math

def magnitude(vec):
	return np.sqrt(vec.dot(vec))

class SingleAttractor_Kepler:
	
	def __init__(self, attractor_mass_kg, r0, v0):
		self.attractor_mass_kg = attractor_mass_kg
		self.r = r0
		self.t = 0
		self.set_trajectory(v0)
	
	def set_trajectory(self, v0):
		G = 6.67408e-11
		mu = G * self.attractor_mass_kg
		t0 = self.t
		r0 = self.r.copy()
		v0 = v0.copy()
		r0_mag = magnitude(r0)
		v0_mag = magnitude(v0)
		energy = math.pow(v0_mag, 2) / 2 - mu / r0_mag
		a = - (mu / (2 * energy))
		r0_dot_v0 = r0.dot(v0)
		root_mu = math.sqrt(mu)
		r0_dot_v0_over_root_mu = r0_dot_v0 / root_mu
		one_minus_r0_mag_over_a = 1 - (r0_mag / a)

		def get_C_or_S_series(z, C_not_S):
			epsilon = 0.0000000000000001
			exponent = 0
			fact = 2 if C_not_S else 3
			sign = 1
			result = 0
			prev_result = -1
			while (math.fabs(result - prev_result) > epsilon):
				prev_result = result
				result += math.pow(z, exponent) / math.factorial(fact) * sign
				sign = -sign
				fact += 2
				exponent += 1
			return result	
		
		def get_C(z):
			if z > 0.00001:
				return (1 - math.cos(math.sqrt(z))) / z
			elif z < -0.00001:
				return (1 - math.cosh(math.sqrt(-z))) / z
			else:
				return get_C_or_S_series(z, True)
				
		def get_S(z):
			if z > 0.00001:
				return (math.sqrt(z) - math.sin(math.sqrt(z))) / math.pow(z, 1.5)
			elif z < -0.00001:
				return (math.sinh(math.sqrt(-z)) - math.sqrt(-z)) / math.pow(-z, 1.5)
			else:
				return get_C_or_S_series(z, False)
		
		def get_t(x, z, C, S):
			return ((
				one_minus_r0_mag_over_a * S * math.pow(x, 3)
				+		
				r0_dot_v0_over_root_mu * C * math.pow(x, 2)
				+
				r0_mag * x
				) / root_mu
				)
			
		def get_dt_dx(x, z, C, S):
			return ((
				C * math.pow(x, 2)
				+
				r0_dot_v0_over_root_mu * (1 - (z * S)) * x
				+
				r0_mag * (1 - (z * C))
				) / root_mu
				)
			
		def get_x_C_S(t):
			epsilon = 0.000001
			x_n = 1.58
			while (True):
				z_n = math.pow(x_n, 2) / a
				C_n = get_C(z_n)
				S_n = get_S(z_n)
				t_n = get_t(x_n, z_n, C_n, S_n)
				print("x_n = {0}, t_n = {1}".format(x_n, t_n))
				if math.fabs(t_n - t) < epsilon:
					return (x_n, C_n, S_n)
				dt_dx_n = get_dt_dx(x_n, z_n, C_n, S_n)
				x_n += (t - t_n) / dt_dx_n
				
		def get_f(x, C):
			return 1 - (math.pow(x, 2) / r0_mag) * C
			
		def get_g(t, x, S):
			return t - (math.pow(x, 3) / root_mu) * S
			
		def get_r(t):
			x, C, S = get_x_C_S(t - t0)
			f = get_f(x, C)
			g = get_g(t - t0, x, S)
			return f * r0 + g * v0
			
		def update(t):
			self.r = get_r(t)
			self.t = t
	
		self.update = update

