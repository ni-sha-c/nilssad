program driver
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head_homogeneous
	type(active) :: x
    type(active), dimension(3) :: y
	double precision, dimension(3) :: X0, X1, Xorig
	double precision, dimension(3,1) :: v0,v,vorig
	integer :: t, nSteps
	double precision :: s
	character(len=128) :: arg


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
	vorig = v	

	Open(1, file="param.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) s
    Close(1)


	x%v = 1.d-4	
	call head_homogeneous(x,s,y,X0,v0,nSteps)

	Open(1, file="output_htangents.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')
    Write(1) x%d
    Close(1)
		

end program driver
