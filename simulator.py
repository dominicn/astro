import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math
import numpy as np
import time
import consts
import orbits
import sa_schoolboy
import sa_kepler
import sa_crk4

orbit_radius_m = consts.earth_radius_m + 200000

orbit_speed_ms = orbits.get_speed_circular_orbit(
	orbit_radius_m,
	consts.earth_mass_kg)

orbit_period_s = orbits.get_period_elliptical_orbit(
	orbit_radius_m,
	consts.earth_mass_kg)

print("orbit period s = " + str(orbit_period_s))
		
def run_sim(duration):
	start_time = time.time()
	
	bodies = [
		sa_schoolboy.SingleAttractor_Schoolboy(
			consts.earth_mass_kg,
			np.array([orbit_radius_m, 0.0, 0.0]),
			np.array([0.0, orbit_speed_ms, 0.0])
			),
		sa_kepler.SingleAttractor_Kepler(
			consts.earth_mass_kg,
			np.array([orbit_radius_m, 0.0, 0.0]),
			np.array([0.0, orbit_speed_ms, 0.0])
			),
		sa_crk4.SingleAttractor_CowellRK4(
			consts.earth_mass_kg,
			np.array([orbit_radius_m, 0.0, 0.0]),
			np.array([0.0, orbit_speed_ms, 0.0])
			)
		]
	
	times = []
	positions = []
	
	t = 0.0
	times.append(t)
	positions_this_ts = []
	for body in bodies:
		positions_this_ts.append(body.r.copy())
	positions.append(positions_this_ts)
	
	while (t < duration):
		
		t += 5
		if (t > duration):
			t = duration
		times.append(t)
		
		positions_this_ts = []
		for body in bodies:
			body.update(t)
			positions_this_ts.append(body.r.copy())
		positions.append(positions_this_ts)
			
	end_time = time.time()
	print("Simulated {0} time slices in {1}s".format(
		len(times) - 1,
		end_time - start_time
		))
	return times, positions

times, positions = run_sim(orbit_period_s)

print(len(positions))
print(len(positions[0]))
print(len(positions[-1]))

for i in range(0, 3):
	print("{0} Drift = {1}".format(
		i,
		positions[0][i] - positions[-1][i]
		)
		)

positions_matrix = np.array(positions)

def samples(ds, n):
	if len(ds) <= n:
		return ds.copy()
		
	result = []
	for i in range(0, n - 1):
		idx_to_append = int(i * len(ds) / (n - 1))
		result.append(ds[idx_to_append])
	result.append(ds[len(ds)-1])
	return result
	
sampled_times = samples(times, 5)
sampled_positions = samples(positions, 5)
sampled_positions_matrix = np.array(sampled_positions)

def draw_trajectory(ax, xs, ys, ts, colour, num_t_anns):
	ax.plot(xs, ys, colour)
	xs_samples = samples(xs, num_t_anns)
	ys_samples = samples(ys, num_t_anns)
	ts_samples = samples(ts, num_t_anns)
	ax.plot(xs_samples, ys_samples, colour + 'o')
	for i in range(0, len(ts_samples)):
		ax.annotate(
			't={0:0.1f}'.format(ts_samples[i]),
			xy=(xs_samples[i], ys_samples[i]),
			textcoords='offset points',
			xytext=(4, 4)
			)
			
colours = ['b', 'g', 'r']

def draw_trajectories(subplot, x_idx, y_idx):
	for i in range(0, 3):
		draw_trajectory(
			subplot,
			positions_matrix[:, i, x_idx],
			positions_matrix[:, i, y_idx],
			times,
			colours[i],
			5)

fig = plt.figure()

xy = plt.subplot(223)
xy.set_xlabel('x')
xy.set_ylabel('y', rotation=0)
draw_trajectories(xy, 0, 1)

zy = plt.subplot(224, sharey=xy)
zy.get_yaxis().set_visible(False)
zy.set_xlabel('z')
draw_trajectories(zy, 2, 1)

xz = plt.subplot(221, sharex=xy)
xz.get_xaxis().set_visible(False)
xz.set_ylabel('z', rotation=0, labelpad=15)
draw_trajectories(xz, 0, 2)

fig.autofmt_xdate()

plt.show()

