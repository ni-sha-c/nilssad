module cheb
	implicit none
	integer, parameter :: Ncheb = 10 ! this means that the interpolating 
!Lagrange polynomials have a maximum degree of Ncheb - 1
contains
	subroutine solve_advection(u0, steps)
		implicit none
		real(kind=8), dimension(Ncheb) :: u0, u1
		real(kind=8) :: dt
		integer : steps, k
	
		dt = 1.e-4
		D = cheb_diff_matrix(Ncheb)
		do k = 1, Ncheb, 1
			u1(k) = u0(k)
		end do
		do k = 1, steps, 1
			u1(Ncheb) = 1.0
			u1(1:Ncheb) = np.linalg.solve(np.eye(n-1,n-1)-dt*D,u0[0:n-1])
			u0 = u1
		end do
	end subroutine solve_advection
	function cheb_pts(k,n):
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


