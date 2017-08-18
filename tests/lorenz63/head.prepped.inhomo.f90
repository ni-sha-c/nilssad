subroutine head_inhomogeneous(s,sprime,y,X0,v0,nSteps,d)
	use Lorenz63
	implicit none 
	double precision, dimension(d) :: X
	double precision, dimension(d) :: Xnp1_res
	double precision, intent(in) :: s
	double precision, dimension(d), intent(in) :: X0, v0
	double precision, dimension(d), intent(out) :: y
	integer :: t, t1
	integer, intent(in) :: nSteps, d
	double precision :: sprime


!$openad INDEPENDENT(s)
	do t = 1, d, 1 
		X(t) = X0(t) + v0(t)*(s-sprime)
	end do
	do t = 1, nSteps, 1
		call Xnp1(X,Xnp1_res,s,d)
		do t1 = 1, d, 1 
			X(t1) = X0(t1) + v0(t1)*(s-sprime)
		end do
	end do
	do t = 1, d, 1
		y(t) = X(t)
	end do
!$openad DEPENDENT(y)

end subroutine


