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
	n = len(u0)
	dt = 1.e-4
	D = cheb_diff_matrix(n)
	D = D[0:n-1,0:n-1]
	u1 = np.copy(u0)
	for k in range(steps):
		u1[-1] = 1.0
		u1[0:n-1] = np.linalg.solve(np.eye(n-1,n-1)-dt*D,u0[0:n-1])
		u0 = u1
	return u1
def cheb_pts(k,n):
	return np.cos(k*np.pi/n)
def cheb_diff_matrix(n):
	D = np.eye(n+1,n+1)
	x = [cheb_pts(k,n) for k in range(n+1)]
	x = np.array(x)
	c = np.ones(n+1)
	c[0] = c[n] = 2.0
	if n>1 :
		for i in range(1,n):
			D[i,i] = -1.e0*x[i]/(1.e0 - x[i]**2.0)/2.e0 
			for j in range(i+1,n):
				D[i,j] = (-1.e0)**(i+j)/(x[i]-x[j])
				D[j,i] = -1.e0*D[i,j]
	D[0,0] = (2.0*(n**2.0)+1.e0)/6.e0
	D[n,n] = -1.e0*D[0,0]
	D[0,1:n+1] = [ 2.e0*(-1)**j/(x[0]-x[j]) for j in range(1,n+1)]
	D[0,n] = D[0,n]/2.e0
	D[1:n+1,0] = -1.e0*D[0,1:n+1]/4.e0
	D[n,0] = D[n,0]*4.e0
	D[n,1:n] = [ 2.e0*(-1)**(n+j)/(x[n]-x[j]) for j in range(1,n)]
	D[1:n,n] = -1.e0*D[n,1:n]/4.e0
	return D


