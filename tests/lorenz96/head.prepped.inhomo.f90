subroutine head_inhomogeneous(s,sprime,y,X0,v0,nSteps)
	use lorenz96
	implicit none 
	double precision, dimension(d) :: X
	double precision, dimension(d) :: Xnp1_res
	double precision, intent(in) :: s
	double precision, dimension(d) :: X0, v0
	double precision, dimension(d) :: y
	integer :: t, t1
	integer :: nSteps
	double precision :: sprime


!$openad INDEPENDENT(s)
	do t = 1, d, 1
		X(t) = X0(t) + v0(t)*(s-sprime)
	end do
	do t = 1, nSteps, 1
		call Xnp1(X,Xnp1_res,s)
		do t1 = 1, d, 1
			X(t1) = Xnp1_res(t1)
		end do
	end do
	do t = 1, d, 1
		y(t) = X(t)
	end do
!$openad DEPENDENT(y)

end subroutine


