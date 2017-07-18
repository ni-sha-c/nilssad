subroutine head_inhomogeneous(c1,c1prime,params_passive,y,X0,v0,nSteps)
	use combustion_juniper
	implicit none 
	double precision, dimension(d) :: X
	double precision, dimension(d) :: Xnp1_res
	double precision, intent(in) :: c1
	double precision :: c1prime
	double precision, dimension(N_p-1) :: params_passive
	double precision, dimension(N_p) :: params
	double precision, dimension(d), intent(in) :: X0, v0
	double precision, dimension(d), intent(out) :: y
	integer :: t, t1, t2, inttau, counter
	integer, intent(in) :: nSteps
	double precision :: sprime, tau
	double precision, dimension(d) :: zeroarray
	double precision, dimension(:,:), allocatable :: Xtmtau

!$openad INDEPENDENT(c1)
	tau = params_passive(N_p-1)
	do t = 1, d, 1 
		X(t) = X0(t) + v0(t)*(c1-c1prime)
		zeroarray(t) = 0.d0
	end do

	inttau = int(tau/dt)
	allocate(Xtmtau(d,inttau))
	
	do t = 1, inttau, 1
		call Xnp1(X,Xnp1_res,zeroarray,c1,params_passive)
		do t1 = 1, d, 1
			X(t1) = Xnp1_res(t1)
			Xtmtau(t1,t) = Xnp1_res(t1)		
		end do
	end do
	counter = 0	
	do t = 1, nSteps, 1
		counter = counter + 1
		call Xnp1(X,Xnp1_res,Xtmtau(:,counter),c1,params_passive)
		do t1 = 1, d, 1
			X(t1) = Xnp1_res(t1)
			do t2 = 1, inttau-1, 1
				Xtmtau(t1,t2) = Xtmtau(t1,t2+1)		
			end do
			Xtmtau(t1,inttau) = Xnp1_res(t1)
		end do
		if(counter .eq. inttau - 1) then 
			counter = 0
		end if 
	end do

	do t = 1, d, 1
		y(t) = X(t)
	end do
	!y(1) = 1.d0
	!y(2) = s
	!y(3) = 1.d0
!$openad DEPENDENT(y)

end subroutine
