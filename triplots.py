# Code to draw trajectories, labels etc. on a three view setup where
# each view looks directly down one of the axes.

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.pyplot import show
import numpy as np

def new_triplot():
	fig = plt.figure()

	xy = plt.subplot(223)
	xy.set_xlabel('x')
	xy.set_ylabel('y', rotation=0, labelpad=15)
	
	zy = plt.subplot(224, sharey=xy)
	zy.get_yaxis().set_visible(False)
	zy.set_xlabel('z')

	xz = plt.subplot(221, sharex=xy)
	xz.get_xaxis().set_visible(False)
	xz.set_ylabel('z', rotation=0, labelpad=15)

	fig.autofmt_xdate()
	
	global _curr_axes
	_curr_axes = (
		(xy, 0, 1),
		(zy, 2, 1),
		(xz, 0, 2)
		)

def plot_path(path, colour, time_samples=4):
	global _curr_axes
	def slice(idx):
		return list(map(lambda p: p[idx], path))
	for axes, x_idx, y_idx in _curr_axes:
		axes.plot(slice(x_idx), slice(y_idx), colour)
		
def label_point(pos, label, colour):
	global _curr_axes
	for axes, x_idx, y_idx in _curr_axes:
		axes.plot(pos[x_idx], pos[y_idx], colour + 'o')
		axes.annotate(
			label,
			xy=(pos[x_idx], pos[y_idx]),
			textcoords='offset points',
			xytext=(4, 4)
			)

if __name__ == "__main__":
	new_triplot()
	plot_path([
		np.array([0.0, 0, 0]),
		np.array([0.1, 1, 0]),
		np.array([0.3, 2, 0]),
		np.array([0.2, 3, 0]),
		np.array([0.25, 4, 0]),
		np.array([0.6, 5, 0]),
		np.array([1.0, 6, 0])
		], 'b')
	label_point(np.array([1.0, 6, 0]), 'end', 'k')
	show()

