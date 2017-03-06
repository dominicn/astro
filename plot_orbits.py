import numpy as np
import math
from triplots import new_triplot, plot_path, label_point, show
from orbits import Orbit

orbits = [
	Orbit.from_r_v_mu(
		np.array([4, 0, 0]),
		np.array([0, 0.65, 0]),
		1),
	Orbit.from_r_v_mu(
		np.array([0, 4, 0]),
		np.array([-.65, 0, 0]),
		1)	
	]

new_triplot()	

colours = ['r', 'g', 'b']
for i in range(0, len(orbits)):
	orbit = orbits[i]
	print("orbits[{0}] = {1}".format(i, orbit))
	path = orbit.get_r_vecs()
	colour = colours[i % len(colours)]
	plot_path(path, colour)

for nu, r in orbits[0].get_intersections(orbits[1]):
	r_vec = orbits[0].get_r_vec(nu)
	print("intersection at nu={0}, r={1}, r_vec={2}".format(nu, r, r_vec))
	label_point(r_vec, "intersection", 'k')

show()

