! Lorenz ' 63 system

module Lorenz63
	implicit none

	double precision, parameter :: dt = 0.001d0
	double precision, parameter :: sigma = 10.d0, b = 8.d0/3.d0	

contains

subroutine Xnp1(X,Xnp1_res,r,d)

	implicit none
	integer :: d 
	double precision, dimension(d):: X
	double precision, dimension(d), intent(out) :: Xnp1_res
    double precision, dimension(d):: k1, k2, k3, k4
	double precision, dimension(d):: ddt
	integer :: i
	double precision :: r
	
    call dXdt(X,ddt,r,d)
    do i = 1, d, 1
		k1(i) = dt*ddt(i)
	   	Xnp1_res(i) = X(i) + 0.5d0*k1(i) 
	end do
	call dXdt(Xnp1_res,ddt,r,d)
    do i = 1, d, 1
		k2(i) = dt*ddt(i)
		Xnp1_res(i) = X(i) + 0.5d0*k2(i) 
	end do
	call dXdt(Xnp1_res,ddt,r,d)
    do i = 1, d, 1
		k3(i) = dt*ddt(i)
		Xnp1_res(i) = X(i) + k3(i) 
	end do
	call dXdt(Xnp1_res,ddt,r,d)
 	do i = 1, d, 1
		k4(i) = dt*ddt(i) 
	end do
  
	do i = 1, d, 1
    	Xnp1_res(i) = X(i) + 1.d0/6.d0*k1(i) + &
               1.d0/3.d0*k2(i) + 1.d0/3.d0*k3(i) + &
                1.d0/6.d0*k4(i)   

	end do

end subroutine Xnp1
subroutine dXdt(X,dXdt_res,r,d)
	implicit none
	integer :: d
	double precision, dimension(d) :: X
	double precision :: r
	double precision, intent(out), dimension(d) :: dXdt_res
	double precision :: sigma, b
	integer :: i
	
	dXdt_res(1) = -sigma*X(1) + sigma*X(2)
	dXdt_res(2) = X(1)*(r + 28.d0 -X(3)) - X(2)
	dXdt_res(3) = X(1)*X(2) - b*X(3)  
        
end subroutine dXdt
subroutine Objective(X,J,s,d)
	implicit none
	integer :: d
	double precision, intent(in), dimension(d) :: X
	double precision, intent(out) :: J
	double precision, intent(in) :: s
	
	J = (X(3) - 28.d0)**2.0

end subroutine Objective
end module Lorenz63
