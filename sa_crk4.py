# Simulation model using Cowell's method (which is just adding
# up all the accelerations in the obvious way) with Runge-Kutta
# fourth order numerical integration (RK4).

import numpy as np
import math
import consts

class SingleAttractor_CowellRK4:
	
	def __init__(self, attractor_mass_kg, r0, v0):
		self.r = r0
		self.t = 0
		self.v = v0
		self.mu = consts.G * attractor_mass_kg
	
	# Runge-Kutta 4th order method applied to
	# two first order differential equations:
	# dr/dt = v
	# dv/dt = a = get_a(r)
	def update(self, t):
		h = t - self.t
		k0 = h * self.v
		l0 = h * self.get_a(self.r)
		k1 = h * (self.v + l0/2)
		l1 = h * self.get_a(self.r + k0/2)
		k2 = h * (self.v + l1/2)
		l2 = h * self.get_a(self.r + k1/2)
		k3 = h * (self.v + l2)
		l3 = h * self.get_a(self.r + k2)
		self.r += 1/6 * (k0 + 2*k1 + 2*k2 + k3)
		self.v += 1/6 * (l0 + 2*l1 + 2*l2 + l3)
		self.t = t
	
	# cowell's method is to compute all the
	# accelerations from all bodies and sum
	# them - follows directly from the
	# Newtonian equations of motion
	def get_a(self, r):
		r_mag_cu = math.pow(r.dot(r), 1.5)
		return self.mu / r_mag_cu * -r

