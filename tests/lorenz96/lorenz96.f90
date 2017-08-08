! Lorenz ' 96 system

module lorenz96
    implicit none
    
    integer, parameter :: d = 32 

contains

subroutine Xnp1(Xorig,Xnp1_res,M)

        
    implicit none
	real(kind=8) :: M
	real(kind=8), dimension(D) :: Xorig
    real(kind=8), dimension(D+3):: k1, k2, k3, k4, Xnp1_res_temp
	real(kind=8), dimension(D) :: Xnp1_res
    real(kind=8), dimension(D+3):: ddt
	real(kind=8), dimension(D+3) :: X
	integer :: i, Dext
	real(kind=8) :: dt
	
	!if(present(M) .eqv. .false.) then
	!	M = 8.1d0
	!end if
	dt = 0.01d0
	do i = 1, d, 1
		X(i+2) = Xorig(i)
	end do
	X(d+3) = X(3)
	X(2) = X(d+2)
	X(1) = X(d+1)
		
    call dXdt(X,ddt,M)
    do i = 3, d+2, 1
        k1(i) = dt*ddt(i)
        Xnp1_res_temp(i) = X(i) + 0.5d0*k1(i)
    end do
    k1(1) = k1(d+1)
    k1(2) = k1(d+2)
    k1(d+3) = k1(3)
    Xnp1_res_temp(1) = Xnp1_res_temp(d+1)
    Xnp1_res_temp(2) = Xnp1_res_temp(d+2)
    Xnp1_res_temp(d) = Xnp1_res_temp(3)


    call dXdt(Xnp1_res_temp,ddt,M)
    do i = 3, d+2, 1
        k2(i) = dt*ddt(i)
        Xnp1_res_temp(i) = X(i) + 0.5d0*k2(i)
    end do

    k2(1) = k2(d+1)
    k2(2) = k2(d+2)
    k2(d+3) = k2(3)
    Xnp1_res_temp(1) = Xnp1_res_temp(d+1)
    Xnp1_res_temp(2) = Xnp1_res_temp(d+2)
    Xnp1_res_temp(d+3) = Xnp1_res_temp(3)


    
    call dXdt(Xnp1_res_temp,ddt,M)
    do i = 3, d+2, 1
        k3(i) = dt*ddt(i)
        Xnp1_res_temp(i) = X(i) + k3(i)
    end do

    k3(1) = k3(d+1)
    k3(2) = k3(d+2)
    k3(d+3) = k3(3)
    Xnp1_res_temp(1) = Xnp1_res_temp(d+1)
    Xnp1_res_temp(2) = Xnp1_res_temp(d+2)
    Xnp1_res_temp(d+3) = Xnp1_res_temp(3)


    
    call dXdt(Xnp1_res_temp,ddt,M)
    do i = 3, d+2, 1
        k4(i) = dt*ddt(i)
    end do

   
    do i = 3, d+2, 1 
        Xnp1_res(i-2) = X(i) + 1.d0/6.d0*k1(i) + &
               1.d0/3.d0*k2(i) + 1.d0/3.d0*k3(i) + &
                1.d0/6.d0*k4(i)    

    end do

end subroutine Xnp1
subroutine dXdt(X,ddt_res,M)
    implicit none
	real(kind=8) :: M
	real(kind=8), dimension(D+3) :: X
	real(kind=8), dimension(D+3) :: ddt_res
	integer :: i
	real(kind=8) :: dt
    
    dt = 0.01d0
    !if(present(M) .eqv. .false.) then
	!	M = 8.1d0
	!end if
    do i= 3,d+2,1
		ddt_res(i) = (-X(i-2) + X(i+1))*X(i-1) - X(i) + M	
	end do

    
end subroutine dXdt
subroutine Objective(X,J,s)
	implicit none
	double precision, dimension(d) :: X
	double precision :: J
	double precision :: s
    integer :: i 
    
    J = 0.d0
    do i = 1, d, 1
        J = J + X(i)
    end do
	J = J/d


end subroutine Objective

end module lorenz96
