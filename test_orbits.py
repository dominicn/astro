import unittest
import numpy as np
import orbits
import math
import math_ext

class TestOrbits(unittest.TestCase):
	
	def from_r_v_mu_test(self, r_vec, v_vec, mu, expected_str):
		orbit = orbits.Orbit.from_r_v_mu(r_vec, v_vec, mu)
		
		# Little bit of code duplicated from part of the method under test
		# but it allows us to test assertions about the P and W vectors
		# via a different route to the one taken by the code under test.
		# W is always a unit vector in the direction of h, and P is a unit
		# vector in the direction of e if e has magnitude (i.e. the orbit
		# isn't circular) otherwise n if n has magnitude (i.e. the orbit
		# isn't equatorial) otherwise it's I.
		h_vec = np.cross(r_vec, v_vec)
		h = math_ext.magnitude(h_vec)
		expected_W_vec = h_vec / h
		e_vec = orbits.get_eccentricy_vec(r_vec, v_vec, mu)
		e = math_ext.magnitude(e_vec)
		if e > 0:
			expected_P_vec = e_vec / e
		else:
			n_vec = np.cross(np.array([0, 0, 1]), h_vec)
			n = math_ext.magnitude(n_vec)
			if n > 0:
				expected_P_vec = n_vec / n
			else:
				expected_P_vec = np.array([1, 0, 0])

		np.testing.assert_allclose(
			orbit.W_vec,
			expected_W_vec,
			err_msg="W_vec mismatch, orbit=" + str(orbit),
			atol=1e-7)
		np.testing.assert_allclose(
			orbit.P_vec,
			expected_P_vec,
			err_msg="P_vec mismatch, orbit=" + str(orbit),
			atol=1e-7)
		self.assertEqual(expected_str, str(orbit))
			
	def test_equatorial_circular_direct_orbit(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, 0.5, 0]),
			1,
			"Orbit(mu=1, p=4, e=0, i=0, OMEGA=nan, PI=nan)")
	
	def test_equatorial_circular_retrograde_orbit(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, -0.5, 0]),
			1,
			"Orbit(mu=1, p=4, e=0, i=3.142, OMEGA=nan, PI=nan)")
	
	def test_polar_circular_orbit_ascending_at_i(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, 0, 0.5]),
			1,
			"Orbit(mu=1, p=4, e=0, i=1.571, OMEGA=0, PI=nan)")
			
	def test_polar_circular_orbit_descending_at_i(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, 0, -0.5]),
			1,
			"Orbit(mu=1, p=4, e=0, i=1.571, OMEGA=3.142, PI=nan)")
			
	def test_polar_circular_orbit_ascending_at_j(self):
		self.from_r_v_mu_test(
			np.array([0, 4, 0]),
			np.array([0, 0, 0.5]),
			1,
			"Orbit(mu=1, p=4, e=0, i=1.571, OMEGA=1.571, PI=nan)")
	
	def test_polar_circular_orbit_descending_at_j(self):
		self.from_r_v_mu_test(
			np.array([0, 4, 0]),
			np.array([0, 0, -0.5]),
			1,
			"Orbit(mu=1, p=4, e=0, i=1.571, OMEGA=4.712, PI=nan)")
			
	def test_inclined_circular_orbit_ascending_at_i(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, math.sqrt(1/8), math.sqrt(1/8)]),
			1,
			"Orbit(mu=1, p=4, e=2.22e-16, i=0.7854, OMEGA=0, PI=0)")
			
	def test_equatorial_eccentric_direct_orbit_periapsis_at_i(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, 1, 0]),
			1,
			"Orbit(mu=1, p=16, e=3, i=0, OMEGA=nan, PI=0)")
			
	def test_equatorial_eccentric_direct_orbit_periapsis_at_j(self):
		self.from_r_v_mu_test(
			np.array([0, 4, 0]),
			np.array([-1, 0, 0]),
			1,
			"Orbit(mu=1, p=16, e=3, i=0, OMEGA=nan, PI=1.571)")
			
	def test_equatorial_eccentric_retrograde_orbit_periapsis_at_i(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, -1, 0]),
			1,
			"Orbit(mu=1, p=16, e=3, i=3.142, OMEGA=nan, PI=0)")
			
	def test_equatorial_eccentric_retrograde_orbit_periapsis_at_j(self):
		self.from_r_v_mu_test(
			np.array([0, 4, 0]),
			np.array([1, 0, 0]),
			1,
			"Orbit(mu=1, p=16, e=3, i=3.142, OMEGA=nan, PI=1.571)")
			
	def test_polar_eccentric_orbit_periapsis_at_i_ascending_at_i(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, 0, 1]),
			1,
			"Orbit(mu=1, p=16, e=3, i=1.571, OMEGA=0, PI=0)")
			
	def test_polar_eccentric_orbit_periapsis_at_k_ascending_at_i(self):
		self.from_r_v_mu_test(
			np.array([0, 0, 4]),
			np.array([-1, 0, 0]),
			1,
			"Orbit(mu=1, p=16, e=3, i=1.571, OMEGA=0, PI=1.571)")

	def test_polar_eccentric_orbit_periapsis_at_i_descending_at_i(self):
		self.from_r_v_mu_test(
			np.array([4, 0, 0]),
			np.array([0, 0, -1]),
			1,
			"Orbit(mu=1, p=16, e=3, i=1.571, OMEGA=3.142, PI=6.283)")
			
	def test_polar_eccentric_orbit_periapsis_at_k_descending_at_i(self):
		self.from_r_v_mu_test(
			np.array([0, 0, 4]),
			np.array([1, 0, 0]),
			1,
			"Orbit(mu=1, p=16, e=3, i=1.571, OMEGA=3.142, PI=4.712)")

	def test_inclined_eccentric_orbit_periapsis_at_kj_ascending_at_i(self):
		self.from_r_v_mu_test(
			# apoapsis at radius 4 but in kj direction
			np.array([0, math.sqrt(8), math.sqrt(8)]),
			np.array([-1, 0, 0]),
			1,
			"Orbit(mu=1, p=16, e=3, i=0.7854, OMEGA=0, PI=1.571)")
			
	def test_inclined_eccentric_orbit_periapsis_at_kmi_ascending_at_j(self):
		self.from_r_v_mu_test(
			# apoapsis at radius 4 but in k, -i direction
			np.array([-math.sqrt(8), 0, math.sqrt(8)]),
			np.array([0, -1, 0]),
			1,
			"Orbit(mu=1, p=16, e=3, i=0.7854, OMEGA=1.571, PI=3.142)")

	def test_inclined_eccentric_orbit(self):
		self.from_r_v_mu_test(
			np.array([0.5, 1, 0]),
			np.array([1, 0, 2]),
			1,
			"Orbit(mu=1, p=6, e=4.502, i=1.991, OMEGA=1.107, PI=7.145)")
			
	def test_inclined_eccentric_orbit_2(self):
		self.from_r_v_mu_test(
			np.array([-0.2, -3.1, -8.7]),
			np.array([-3, 0.1, 4.0]),
			2,
			"Orbit(mu=2, p=471.7, e=76.14, i=1.879, OMEGA=3.547, PI=9.265)")
	
unittest.main(verbosity=2)
