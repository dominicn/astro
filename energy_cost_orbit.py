import numpy as np

altitude_geo_orbit_m = 35786000
speed_geo_orbit_ms = 3070

altitude_leo_m = 200000
speed_leo_ms = 7780

def print_energy_orbit(altitude_m, speed_ms):
	earth_mass_kg = 5.972e+24
	G = 6.67408e-11
	u = earth_mass_kg * G

	radius_earth_m = 12742000 / 2

	# required to make potential energy at earth's surface zero
	potential_energy_const = u / radius_earth_m

	potential_energy_orbit = (
		potential_energy_const -
		u / (altitude_m + radius_earth_m)
		)

	kinetic_energy_orbit = np.square(speed_ms) * .5

	print("Kinetic energy = " + str(kinetic_energy_orbit))
	print("Potential energy = " + str(potential_energy_orbit))

print("geo")
print_energy_orbit(altitude_geo_orbit_m, speed_geo_orbit_ms)

print("leo")
print_energy_orbit(altitude_leo_m, speed_leo_ms)
