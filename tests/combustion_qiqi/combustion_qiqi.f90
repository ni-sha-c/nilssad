!Simple model of combustion in gas turbines

module combustion_qiqi
	use cheb
	implicit none
	
	
	double precision, parameter :: dt = 0.005d0
	integer, parameter :: N = 10
	integer, parameter :: d = 2*N + Ncheb + 1
	integer, parameter :: N_p = 5	
	double precision, parameter :: sigma = 10., b = 8./3.
	double precision, dimension(Ncheb+1,Ncheb+1) :: D_cheb

contains
subroutine Xnp1(X,Xnp1_res,param_active,params_passive)

	implicit none
	double precision :: c2, beta, xf, tau
	double precision, dimension(d):: X, Xnp1_res
    double precision, dimension(d):: k1, k2, k3, k4
	double precision, dimension(d):: ddt
	integer :: i, imax 
	double precision, dimension(N_p-1) :: params_passive
	double precision :: param_active 
	

	c2 = params_passive(1)
	beta = params_passive(2)
	xf = params_passive(3)
	tau = params_passive(4)
		
    call dXdt(X,ddt,param_active,c2,beta,xf,tau)
    do i = 1, d, 1
		k1(i) = dt*ddt(i)
	   	Xnp1_res(i) = X(i) + 0.5d0*k1(i) 
	end do
	call dXdt(Xnp1_res,ddt,param_active,c2,beta,xf,tau)
    do i = 1, d, 1
		k2(i) = dt*ddt(i)
		Xnp1_res(i) = X(i) + 0.5d0*k2(i) 
	end do
	call dXdt(Xnp1_res,ddt,param_active,c2,beta,xf,tau)
    do i = 1, d, 1
		k3(i) = dt*ddt(i)
		Xnp1_res(i) = X(i) + k3(i) 
	end do
	call dXdt(Xnp1_res,ddt,param_active,c2,beta,xf,tau)
 	do i = 1, d, 1
		k4(i) = dt*ddt(i) 
	end do
  
	do i = 1, d, 1
    	Xnp1_res(i) = X(i) + 1.d0/6.d0*k1(i) + &
               1.d0/3.d0*k2(i) + 1.d0/3.d0*k3(i) + &
                1.d0/6.d0*k4(i)   

	end do
	Xnp1_res(d) = 0.d0
	!Xnp1_res(d) = uf(Xnp1_res,xf) 

end subroutine Xnp1
double precision function uf(X,xf)
	
	implicit none
	double precision, dimension(N) :: X
	double precision :: uf0, xf
	integer :: i

	uf0 = 0.d0
	do i = 1, N, 1
	
		uf0 = uf0 + X(i)*cos(i*pi*xf)	

	end do 			
	uf = uf0
end function uf 
double precision function zeta(i,c1,c2)
	
	implicit none
	integer :: i
	double precision :: c1, c2
	zeta = c1*(i**2.0) + c2*(i**0.5d0)

end function zeta
subroutine dXdt(X,dXdt_res,c1,c2,beta,xf,tau)
	implicit none
	double precision :: c1, c2, beta
	double precision :: xf, tau
	double precision, dimension(d) :: X
	double precision, intent(out), dimension(d) :: dXdt_res
	integer :: i,j
	
	do i = 1, N, 1
		dXdt_res(i) = i*pi*X(N+i)
		dXdt_res(N+i) = -1.d0*i*pi*X(i) - zeta(i,c1,c2)*X(N+i) &
						- 2.d0*beta*(((1.d0/3.d0 + X(2*N+1))**2.0d0 &
						+ 0.001d0)**0.25d0 &
						  - (1.d0/3.d0)**0.5d0)* &
						 sin(i*pi*xf)
	end do

	
	do i = 1, Ncheb + 1, 1
	
		dXdt_res(2*N+i) = 0.d0
		do j = 1, Ncheb + 1, 1
			dXdt_res(2*N+i) = dXdt_res(2*N+i) - &
							(1.d0/tau)*D_cheb(i,j)*X(2*N+j)
		end do   	


	end do

end subroutine dXdt
subroutine Objective(X,J,param_active,params_passive)
	implicit none
	double precision, intent(in), dimension(d) :: X
	double precision, intent(out) :: J
	double precision :: xf
	double precision, dimension(N_p-1):: params_passive
	double precision :: param_active
	integer :: t

	J = 0.d0
	xf = params_passive(N_p-2)
	do t = 1, N, 1
		J = J - X(N+t)*sin(t*pi*xf)
	end do

end subroutine Objective
subroutine FlameOutput(X,Xtmtau,pf,ufvar,heat,param_active,params_passive)
	implicit none
	double precision, intent(in), dimension(d) :: X, Xtmtau
	double precision, intent(out) :: pf,ufvar,heat
	double precision :: xf
	double precision, dimension(N_p-1):: params_passive
	double precision :: param_active
	integer :: t

	pf = 0.d0
    ufvar = 0.d0
    heat = 0.d0
	xf = params_passive(N_p-2)
	do t = 1, N, 1
		pf = pf - X(N+t)*sin(t*pi*xf)
        ufvar = ufvar + X(t)*cos(t*pi*xf)
        heat = heat + Xtmtau(t)*cos(t*pi*xf)
    end do
    heat = ((1.d0/3.d0 + heat)**2.d0 + 0.001d0)**0.25d0 - (1.d0/3.d0)**0.5d0
end subroutine FlameOutput

subroutine dfdX(X,dfdX_res)

	implicit none
	double precision, dimension(d) :: X
	double precision, intent(out), dimension(d,d) :: dfdX_res 	
	integer :: i,j
	double precision:: r

	dfdX_res = 0.d0
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

end module combustion_qiqi
