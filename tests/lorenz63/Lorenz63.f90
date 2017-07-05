! Lorenz ' 63 system

module Lorenz63
	implicit none

	double precision, parameter :: dt = 0.001d0
	integer, parameter :: d = 3	

contains

subroutine Xnp1(X,Xnp1_res,r)

	implicit none
	double precision, dimension(d):: X, Xnp1_res
    double precision, dimension(d):: k1, k2, k3, k4
	double precision, dimension(d):: ddt
	integer :: i, imax
	double precision :: r
	
	imax = size(X)
	
    call dXdt(X,ddt,r)
    do i = 1, imax, 1
		k1(i) = dt*ddt(i)
	   	Xnp1_res(i) = X(i) + 0.5d0*k1(i) 
	end do
	call dXdt(Xnp1_res,ddt,r)
    do i = 1, imax, 1
		k2(i) = dt*ddt(i)
		Xnp1_res(i) = X(i) + 0.5d0*k2(i) 
	end do
	call dXdt(Xnp1_res,ddt,r)
    do i = 1, imax, 1
		k3(i) = dt*ddt(i)
		Xnp1_res(i) = X(i) + k3(i) 
	end do
	call dXdt(Xnp1_res,ddt,r)
 	do i = 1, imax, 1
		k4(i) = dt*ddt(i) 
	end do
  
	do i = 1, imax, 1
    	Xnp1_res(i) = X(i) + 1.d0/6.d0*k1(i) + &
               1.d0/3.d0*k2(i) + 1.d0/3.d0*k3(i) + &
                1.d0/6.d0*k4(i)   

	end do

end subroutine Xnp1
subroutine sys_params(sigma, b)
	
	implicit none
	double precision, intent(out) :: sigma, b
	sigma = 10.d0
	b = 8.d0/3.d0

end subroutine sys_params
subroutine dXdt(X,dXdt_res,r)
	implicit none
	double precision, dimension(d) :: X
	double precision :: r
	double precision, intent(out), dimension(d) :: dXdt_res
	double precision :: sigma, b
	integer :: i

	call sys_params(sigma,b) 	
	dXdt_res(1) = -sigma*X(1) + sigma*X(2)
	dXdt_res(2) = X(1)*(r + 28.d0 -X(3)) - X(2)
	dXdt_res(3) = X(1)*X(2) - b*X(3)  
        
end subroutine dXdt
subroutine Objective(X,J,s)
	implicit none
	double precision, intent(in), dimension(d) :: X
	double precision, intent(out) :: J
	double precision, intent(in) :: s
	
	J = (X(3) - 28.d0)**2.0

end subroutine Objective
subroutine dfdX(X,dfdX_res)

	implicit none
	double precision, dimension(d) :: X
	double precision, intent(out), dimension(d,d) :: dfdX_res 	
	integer :: i,j
	double precision:: sigma,b,r

	dfdX_res = 0.d0
	call sys_params(sigma,b)
	r = 28.d0
 
	dfdX_res(1,1) = -1.d0*sigma
	dfdX_res(1,2) =	sigma 
	
	dfdX_res(2,1) = r - X(3)
	dfdX_res(2,2) = -1.d0
	dfdX_res(2,3) = -1.d0*X(1)	
	
	dfdX_res(3,1) = X(2)
	dfdX_res(3,2) = X(1)
	dfdX_res(3,3) = -1.d0*b
    
end subroutine dfdX  

subroutine dvdt(X,v1,dvdt_res)

		implicit none
		double precision, dimension(d) :: X
		double precision, intent(out), dimension(d,1) :: dvdt_res
		double precision, dimension(d,d) :: dfdX_res
		double precision, dimension(d,1) :: v1
		double precision, dimension(d,1):: dfdX_times_v
		integer :: i		
	
		call dfdX(X,dfdX_res)
		!Manual matrix-vector product
		do i = 1,d,1
			dvdt_res(i,1) = dfdX_res(i,1)*v1(1,1) + dfdX_res(i,2)*v1(2,1) + dfdX_res(i,3)*v1(3,1)
		end do
		dvdt_res(1,1) = dvdt_res(1,1) + 0.d0
		dvdt_res(2,1) = dvdt_res(2,1) + X(1)
		dvdt_res(3,1) = dvdt_res(3,1) + 0.d0	

end subroutine dvdt
subroutine rk45_full(X,v,vnp1)
!Assumes full perturbation vector.
	implicit none
	double precision , dimension(d) :: X
	double precision , intent(out), dimension(d,1) :: vnp1
	double precision , dimension(d,1) :: v
	double precision , dimension(d,1) :: v1,k1,k2,k3,k4
	double precision , dimension(d,1):: dvdt_res

	v1 = v
	call dvdt(X,v1,dvdt_res)	
	k1 = dt*dvdt_res
	call dvdt(X,v1 + 0.5d0*k1,dvdt_res)
	k2 = dt*dvdt_res
	call dvdt(X,v1 + 0.5d0*k2,dvdt_res)
	k3 = dt*dvdt_res
	call dvdt(X,v1 + k3,dvdt_res)
	k4 = dt*dvdt_res
	vnp1 = v1 + 1.d0/6.d0*k1 + &
		 1.d0/3.d0*k2 + 1.d0/3.d0*k3 + &
		1.d0/6.d0*k4
    

end subroutine rk45_full 

end module Lorenz63
