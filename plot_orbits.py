import numpy as np
import math
from triplots import new_triplot, plot_curve, show
from orbits import Orbit

orbit = Orbit.from_r_v_mu(
	np.array([0, math.sqrt(8), math.sqrt(8)]),
	np.array([-0.71, 0, 0]),
	1)
	
#orbit = Orbit.from_r_v_mu(
#	np.array([4, 0, 0]),
#	np.array([0, 0, 0.5]),
#	1)
	
print("P_vec = " + str(orbit.P_vec))
print("Q_vec = " + str(orbit.Q_vec))
print("W_vec = " + str(orbit.W_vec))

print(__name__)
print(str(orbit))

path = orbit.get_r_vecs()
print(str(path))

new_triplot()
plot_curve(path, 'r')
show()

