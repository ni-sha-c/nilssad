program driver
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head
	type(active) :: x, y
	double precision, dimension(3) :: X0
	double precision, dimension(3,1) :: v0,v
	double precision :: dzds_t
	integer :: t
	double precision :: zs, z0, ds
	ds = 0.005d0
	x%v=2.8D1
	y%d=1.0D0
	our_rev_mode%tape=.TRUE.
	X0(1) = x%v
	X0(2) = x%v
	X0(3) = x%v
	!Get on the attractor.
	do t = 1, 10000, 1
		call Xnp1(X0,X0,x%v)
	end do	

	v0(1,1) = 2.d0
	v0(2,1) = 0.d0
	v0(3,1) = 1.d0


	call head(x,y,X0,v0)
	print *, 's =',x%v
	print *, 'From AD, dzbar/ds = ',x%d

	!check using tangent equation
	do t = 1, 100, 1
		call rk45_full(X0,v0,v)
		call Xnp1(X0,X0,x%v)
		v0 = v
	end do	

	z0 = X0(3)
	dzds_t = v(3,1)
	print *,'From the tangent equation, dzbar/dx = ', dzds_t

end program driver
