subroutine driver_fortran_ht(subspace_dimension, s, X0, v0, vt, nSteps)
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head_homogeneous
	type(active), dimension(:), allocatable :: x
    type(active), dimension(:,:), allocatable :: y
	double precision, dimension(d) :: X0, X1, Xorig, Xtemp1, Xtemp2
	double precision, dimension(d,subspace_dimension) :: v0,vt
	integer :: t, nSteps, subspace_dimension, t1, t2
	double precision :: s, eps
	character(len=128) :: arg1, arg2


	our_rev_mode%tape=.TRUE.
	

	print *, "subspace dimension is:", subspace_dimension
	eps = 1.d-4
	Xorig = X0
	allocate(y(d,subspace_dimension))
	allocate(x(subspace_dimension))
	do t = 1, subspace_dimension, 1
		do t2 = 1, d, 1
			y(t2,t)%d = 0.d0
			y(t2,t)%d(t2) = 1.d0
		end do
		x(t)%v = eps
		call head_homogeneous(x(t)%v,s,y(:,t),X0,v0(:,t),nSteps)
		do t2 = 1, d, 1
			vt(t2,t) = x(t)%d(t2)
		end do
	
	end do
   	 
end subroutine driver_fortran_ht
