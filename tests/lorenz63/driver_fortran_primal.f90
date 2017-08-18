subroutine driver_fortran_primal(X0, s, nSteps, X1)
	use lorenz63_passive
	implicit none 
	double precision, dimension(:) :: X0
	double precision, dimension(:), intent(out) :: X1
	double precision, dimension(:), allocatable :: J
	integer :: t, nSteps, d
	double precision:: s

	allocate(J(nSteps))
	d = size(X0)
	
	do t = 1, nSteps, 1
		call Xnp1(X0, X1, s, d)
		call Objective(X1,J(t),s, d)
		X0 = X1
	end do
	print *, "X at the end of ", nSteps, " steps is: ", X0
end subroutine driver_fortran_primal
