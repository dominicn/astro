import numpy as np
import math

def rk4_step(h, r0, v0, get_a):
	k0 = h * v0
	l0 = h * get_a(r0)
	k1 = h * (v0 + l0/2)
	l1 = h * get_a(r0 + k0/2)
	k2 = h * (v0 + l1/2)
	l2 = h * get_a(r0 + k1/2)
	k3 = h * (v0 + l2)
	l3 = h * get_a(r0 + k2)
	r1 = r0 + 1/6 * (k0 + 2*k1 + 2*k2 + k3)
	v1 = v0 + 1/6 * (l0 + 2*l1 + 2*l2 + l3)
	return (r1, v1)
	
def rk4_steps(n, h, r0, v0, get_a):
	r = r0
	v = v0
	for i in range(0, n):
		r, v = rk4_step(h, r, v, get_a)
	return r, v
	
def get_a(r):
	return np.array([-10, 0, 0])

for steps_exp in range(0, 7):
	steps = 10 ** steps_exp
	h = 2 / steps
	r, v = rk4_steps(steps, h, np.array([0, 0, 0]), np.array([10, 0, 0]), get_a)
	print("h={0}, r={1}, v={2}".format(h, r, v))


