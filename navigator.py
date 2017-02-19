import numpy as np
import math
import orbits

get_t_and_v_to_reach_orbit(r0, v0, mu, dest_orbit):
	curr_orbit = orbits.Orbit.from_r_v_mu(r0, v0, mu)
	
