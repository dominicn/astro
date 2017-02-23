# Code to draw trajectories, labels etc. on a three view setup where
# each view looks directly down one of the axes.

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.pyplot import show
import numpy as np

def new_triplot():
	fig = plt.figure()

	xy = plt.subplot(223, aspect=1)
	xy.set_xlabel('x')
	xy.set_ylabel('y', rotation=0, labelpad=15)
	
	zy = plt.subplot(224, sharey=xy, aspect=1)
	zy.get_yaxis().set_visible(False)
	zy.set_xlabel('z')

	xz = plt.subplot(221, sharex=xy, aspect=1)
	xz.get_xaxis().set_visible(False)
	xz.set_ylabel('z', rotation=0, labelpad=15)

	fig.autofmt_xdate()
	
	global _curr_axes
	_curr_axes = (
		(xy, 0, 1),
		(zy, 2, 1),
		(xz, 0, 2)
		)

def plot_curve(trajectory, colour, time_samples=4):
	global _curr_axes
	def slice_trajectory(idx):
		print("slice_trajectory(%d) =" % idx)
		result = list(map(lambda p: p[idx], trajectory))
		print(str(result))
		return result
	for axes, x_idx, y_idx in _curr_axes:
		axes.plot(slice_trajectory(x_idx), slice_trajectory(y_idx), colour)

if __name__ == "__main__":
	new_triplot()
	plot_curve([
		np.array([0.0, 0, 0]),
		np.array([0.1, 1, 0]),
		np.array([0.3, 2, 0]),
		np.array([0.2, 3, 0]),
		np.array([0.25, 4, 0]),
		np.array([0.6, 5, 0]),
		np.array([1.0, 6, 0])
		], 'b')
	show()

