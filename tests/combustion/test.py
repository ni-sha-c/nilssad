import numpy as np
def solve_fd(u0, steps):
	u1= u0
	dt = 1.e-4
	dx = 1.e-2
	for k in range(steps):
		u1[0] = 1.0
		for i in range(1,101):
			u1[i] = u0[i]*(1.0 - dt/dx) + u0[i-1]*dt/dx
		u0 = u1
	return u1
def solve_sp(u0, steps):
	u1= u0
	dt = 1.e-4
	dx = 1.e-2
	for k in range(steps):
		u1[0] = 1.0
		for i in range(1,101):
			u1[i] = u0[i]*(1.0 - dt/dx) + u0[i-1]*dt/dx
		u0 = u1
	return u1
def cheb_pts(k,n):
	return np.cos((2.0*k - 1.0)*np.pi/2.0/n)
	
