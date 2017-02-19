# Simulation model I came up with before getting into the
# Fundamentals of Astrodynamics book. I now know that the
# bleeding obvious "just add up all the accelerations"
# approach is called Cowell's method and is widely used,
# but my approach to solving the differential equations
# by holding acceleration constant over a timeslice and
# then using the Newtonian equations to solve precisely
# preforms very poorly compared to more general methods
# of numerical integration, e.g. Runge-Kutta.

import consts
import math
import numpy as np

class SingleAttractor_Schoolboy:
	
	def __init__(self, attractor_mass_kg, r0, v0):
		self.r = r0
		self.v = v0
		self.t = 0.0
		self.mu = consts.G * attractor_mass_kg
		
	def update(self, t):
		r_mag_cu = math.pow(self.r.dot(self.r), 1.5)
		a = self.mu / r_mag_cu * -self.r
		
		duration = t - self.t
		self.r += (
			duration * self.v +
			0.5 * a * np.square(duration)
			)
		self.v += duration * a
		self.t = t

