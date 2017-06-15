program driver
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head
	type(active) :: x, y
	double precision, dimension(3) :: X0
	double precision, dimension(3,1) :: v0,v
	double precision :: dzds_fd
	integer :: t
	double precision :: zs, z0
	x%v=2.8D1
	y%d=1.0D0
	our_rev_mode%tape=.TRUE.
	call head(x,y)
	print *, 'driver running for x =',x%v
	print *, '            yields y =',y%v,' dy/dx =',x%d

	!check using tangent equation
	X0(1) = x%v
	X0(2) = x%v
	X0(3) = x%v
	v0(1,1) = 1.d0
	v0(2,1) = 1.d0
	v0(3,1) = 1.d0
	do t = 1, 1000, 1
		call rk45_full(X0,v0,v)
		call Xnp1(X0,X0,x%v)
		v0 = v
	end do	

	z0 = X0(3)
	dzds_fd = v(3,1)
	print *,'dy/dx from the tangent equation is: ', dzds_fd

	!Check using finite difference	
	X0(1) = x%v
	X0(2) = x%v
	X0(3) = x%v

	do t = 1, 1000, 1
		call Xnp1(X0,X0,x%v+0.0005)
	end do	
	zs = X0(3)
	print *,'dy/dx from finite difference is: ', (zs-z0)/0.0005
	print *, '    dz^1000/dz = ',1.D0*dzds_fd-x%d
end program driver
