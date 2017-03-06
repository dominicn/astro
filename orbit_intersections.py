# Function to compute the intersections between coplanar orbits
# embedded in a script which draws graphs to demonstrate its correctness.

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math
import numpy as np

# This uses the equation of a conic section in polar coordinates
# with the addition of a rotation around the origin supplied by
# using nu = theta + k as the angle to feed into the standard
# equation.
def get_r(p, e, k, theta):
	nu = theta + k
	return p / (1 + e * math.cos(nu))

def draw_conic_section(p, e, k, steps, colour):
	xs = []
	ys = []
	for step in range(0, steps + 1):
			theta = 2 * math.pi / steps * step
			r = get_r(p, e, k, theta)
			x = r * math.cos(theta)
			y = r * math.sin(theta)
			xs.append(x)
			ys.append(y)
	plt.plot(xs, ys, colour)
	
def label_point(r, nu):
	x = r * math.cos(nu)
	y = r * math.sin(nu)
	plt.plot(x, y, 'ko')
	plt.annotate(
		'r={0:0.1f}, theta={1:0.1f}'.format(r, nu),
		xy=(x, y),
		textcoords='offset points',
		xytext=(4, 4)
		)

# A single orbit can be described with two varibles, p and e.
# However not only do those variables say nothing of the plane of
# the orbit, they also do not specify a rotation within the plane.
# Instead, theta = 0 will always point to periapsis.
# To compute the intersections between two coplanar orbits we need
# to consider the parameters of both orbits, p1, e1, p2 and e2, and
# also the angle between them, k.
# This function was derived as follows:
# 1. We're looking for angles from the origin (theta) where both
#    orbits have the same distance from the origin (r):
#    p1 / (1 + e1 cos theta) = p2 / (1 + e2 cos (theta + k))
#    (Note the theta + k on the right to give a rotation to one
#    orbit relative to the other.)
# 2. Use trig identities to get rid of cos (theta + k) and rearrange,
#    until we have a set of terms in cos squared theta, a set of terms
#    in cos theta and a set of terms which don't contain theta.
# 3. That's a quadratic formula ax^2 + bx + c where x = cos theta so
#    solve using standard methods to get zero, one or two real roots.
# 4. Each of those is a possible cos theta. One possible cos theta
#    gives two possible thetas (because cos theta = cos -theta). Check
#    each result by computing r with each equation and keep the ones
#    that match to within some degree of tolerance.
def get_intersections(p1, e1, p2, e2, k):
	cos_thetas = []
	for sign in [1, -1]:
		cos_sq_theta_term = (
			sign * np.square(p1 * e2 * math.sin(k))
			+ np.square(p2 * e1)
			- 2 * p1 * e2 * p2 * e1 * math.cos(k)
			+ np.square(p1 * e2 * math.cos(k))
			)
		cos_theta_term = 2 * (
			np.square(p2) * e1
			- p2 * e1 * p1
			- p1 * e2 * p2 * math.cos(k)
			+ np.square(p1) * e2 * math.cos(k)
			)
		const_term = (
			-sign * np.square(p1 * e2 * math.sin(k))
			+ np.square(p2)
			- 2 * p1 * p2
			+ np.square(p1)
			)
		
		# numpy's roots function can solve quadratic
		# equations but it returns complex numbers and
		# therefore gives all the roots whereas we just
		# want the real ones.
		roots = np.roots([
			cos_sq_theta_term,
			cos_theta_term,
			const_term
			])
		for root in roots:
			cos_thetas.append(root.real)

		# acos always returns a positive number, but
		# cos(-theta) = cos(theta)
		# so both are possible and we'll try both to see
		# which work. We sometimes get values outside the
		# acceptable range for acos (-1 to 1). This can be
		# due to mathemtical error or because the intersections
		# aren't real, so we clamp to the desired range if we're
		# a little over and otherwise discard the point.
		intersections = []
		epsilon = 0.0000000001
		for cos_theta in cos_thetas:
			if math.fabs(cos_theta) > 1.0:
				if math.fabs(cos_theta) < 1.0 + epsilon:
					cos_theta = 1.0 if cos_theta > 0 else -1.0
				else:
					continue
			positive_theta = math.acos(cos_theta)
			for theta in [positive_theta, -positive_theta]:
				r1 = get_r(p1, e1, 0, theta)
				r2 = get_r(p2, e2, k, theta)
				relative_error = math.fabs((r1 - r2) / p1)
				if relative_error < epsilon:
					intersections.append((theta, r1))

	return intersections
		
if __name__ == "__main__":
	fig = plt.figure()
	plt.axes().set_aspect('equal', 'datalim')
	draw_conic_section(2, 0.5, 0, 100, 'r')
	draw_conic_section(3.5, 0.9, math.pi / 3, 100, 'g')
	#label_point(0, 0)

	intersections = get_intersections(2, 0.5, 3.5, 0.9, math.pi / 3)
	for (theta, r) in intersections:
		print("Intersection at r={0}, theta={1}".format(r, theta))
		label_point(r, theta)
	
	plt.show()

