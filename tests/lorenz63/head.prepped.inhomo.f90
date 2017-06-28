subroutine head_inhomogeneous(s,y,X0,v0,nSteps)
	use Lorenz63
	implicit none 
	double precision, dimension(3) :: X
	double precision, dimension(3) :: Xnp1_res
	double precision, intent(in) :: s
	double precision, dimension(3), intent(in) :: X0, v0
	double precision, dimension(3), intent(out) :: y
	integer :: t
	integer, intent(in) :: nSteps
	double precision :: ds



!$openad INDEPENDENT(s) 
	ds = s - 2.8005e1
	X(1) = X0(1) + v0(1)*ds
	X(2) = X0(2) + v0(2)*ds
	X(3) = X0(3) + v0(3)*ds
	do t = 1, nSteps, 1
		call Xnp1(X,Xnp1_res,s)
		X = Xnp1_res
	end do
	y = X 
	!Let's differentiate v with respect to time

!$openad DEPENDENT(y)
end subroutine


