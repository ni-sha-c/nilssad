subroutine head(s,zres)
	use Lorenz63
	implicit none 
	double precision, dimension(3) :: X
	double precision, dimension(3) :: Xnp1_res
	double precision, intent(in) :: s
	double precision, intent(out) :: zres
	integer :: t



!$openad INDEPENDENT(s) 
	X(1) = -2.4d0
	X(2) = -3.7d0
	X(3) = 14.98d0
	do t = 1, 1000, 1
		call Xnp1(X,Xnp1_res,s)
		X = Xnp1_res
	end do

	zres = Xnp1_res(3) 
	!Let's differentiate v with respect to time

!$openad DEPENDENT(zres)
end subroutine


