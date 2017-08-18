subroutine driver_fortran_ht(subspace_dimension, s, X0, v0, vt, nSteps)
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head_homogeneous
	integer :: d
	type(active), dimension(:), allocatable :: x
    type(active), dimension(:,:), allocatable :: y
	double precision, dimension(:) :: X0
	double precision, dimension(:, :) :: v0
	double precision, dimension(:, :), allocatable, intent(out) :: vt
	integer :: t, nSteps, subspace_dimension, t1, t2
	double precision :: s, eps
	character(len=128) :: arg1, arg2


	our_rev_mode%tape=.TRUE.
	

	d = size(X0)	
	eps = 1.d-4
	allocate(y(d,subspace_dimension))
	allocate(x(subspace_dimension))
	allocate(vt(d,subspace_dimension))
	do t = 1, subspace_dimension, 1
		do t2 = 1, d, 1
			y(t2,t)%d = 0.d0
			y(t2,t)%d(t2) = 1.d0
		end do
		x(t)%v = eps
		call head_homogeneous(x(t)%v,s,y(:,t),X0,v0(:,t),nSteps,d)
		do t2 = 1, d, 1
			vt(t2,t) = x(t)%d(t2)
		end do
	
	end do
   	 
end subroutine driver_fortran_ht
