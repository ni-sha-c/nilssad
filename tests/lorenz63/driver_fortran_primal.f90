subroutine driver_fortran_primal(X0, s, nSteps, X1)
	use lorenz63_passive
	implicit none 
	double precision, dimension(d) :: X0, X1
	double precision, dimension(:), allocatable :: J
	integer :: t, nSteps
	double precision:: s

	allocate(J(nSteps))

	
	do t = 1, nSteps, 1
		call Xnp1(X0, X1, s)
		call Objective(X1,J(t),s)
		X0 = X1
	end do

end subroutine driver_fortran_primal
