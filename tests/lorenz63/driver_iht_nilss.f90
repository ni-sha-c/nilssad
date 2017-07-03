program driver
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head_inhomogeneous
	integer, parameter:: d=3
	type(active) :: x
    type(active), dimension(d) :: y
	double precision, dimension(d) :: X0, X1, Xorig
	double precision, dimension(d,1) :: v0,v, vorig
	double precision :: dzds_t, xprime
	integer :: t, nSteps
	double precision :: ds
	character(len=128) :: arg


	ds = 0.005d0
	do t = 1, d, 1
		y(t)%d = 0.d0
		y(t)%d(t) = 1.d0
	end do
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

	Open(1, file="input_ihtangent.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) v0(:,1)
    Close(1)
	vorig = v0

	Open(1, file="param.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) xprime
    Close(1)
	x%v = xprime
	xprime = xprime - ds

	
	call head_inhomogeneous(x,xprime,y,X0,v0,nSteps)
		
	Open(1, file="output_ihtangent.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')
    Write(1) x%d
    Close(1)



end program driver
