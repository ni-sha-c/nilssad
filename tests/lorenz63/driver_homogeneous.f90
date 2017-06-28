program driver
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head_homogeneous
	type(active) :: x
    type(active), dimension(3):: y
	double precision, dimension(3) :: X0
	double precision, dimension(3,1) :: v0,v
	double precision :: dzds_t
	integer :: t, nSteps
	double precision :: zs, z0, ds
	character(len=128) :: arg


	x%v = 0.005d0
	y%d=1.0D0
	our_rev_mode%tape=.TRUE.


	
	if (command_argument_count() .ne. 1) then
        print *, "Need number of time steps"
        call exit(-1)
    end if

	call get_command_argument(1, arg)
    Read(arg, '(i10)') nSteps


	Open(1, file="input_primal.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) X0
    Close(1)


	Open(1, file="input_tangent.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) v0(:,1)
    Close(1)


	Open(1, file="param.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) s
    Close(1)


	
	call head_homogeneous(s,x,X0,v0,y,nSteps)
	print *, 'epsilon =',x%v
	print *, 'From AD, dX/depsilon = ',x%d

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
