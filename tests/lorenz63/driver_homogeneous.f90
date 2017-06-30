program driver
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head_homogeneous
	type(active) :: x
    type(active), dimension(3) :: y
	double precision, dimension(3) :: X0, X1, Xorig
	double precision, dimension(3,1) :: v0,v, vorig
	integer :: t, nSteps
	double precision :: s
	character(len=128) :: arg


	ds = 0.005d0
	y(1)%d = [1.d0, 0.d0, 0.d0]
	y(2)%d = [0.d0, 1.d0,0.d0]
	y(3)%d = [0.d0, 0.d0,1.d0]
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
	Xorig = X0

	Open(1, file="input_htangents.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) v0(:,1)
    Close(1)
	vorig = v0

	Open(1, file="param.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) s
    Close(1)


	x%v = 1.d-4	
	call head_homogeneous(x,s,y,X0,v0,nSteps)
		
	print *, 's =',s
	do t = 1, 3, 1
		print *, 'dy/dx at position ', t, 'is ',  x%d(t)
	end do

	!check using tangent equation
	X0 = X0 + v0(:,1)*5.d-3
	v0 = 0.d0
	do t = 1, nSteps, 1
		call rk45_full(X0,v0,v)
		call Xnp1(X0,X1,x%v)
		X0 = X1
		v0 = v
	end do	

	print *,'From the tangent equation = ', v

	!check using finite difference
	Xorig = Xorig + vorig(:,1)*5.d-3
	v0 = 0.d0
	do t = 1, nSteps, 1
		call rk45_full(X0,v0,v)
		call Xnp1(Xorig,X1,x%v+0.005d0)
		Xorig = X1
		v0 = v
	end do	

	print *,'From finite difference = ', -1.d0*(X0-Xorig)/0.005d0
end program driver
