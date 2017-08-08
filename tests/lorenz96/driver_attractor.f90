program flow
	use lorenz96_passive
	implicit none 
	double precision, dimension(d) :: X0, X1
	double precision, dimension(:), allocatable :: J
	integer :: t, nSteps
	double precision:: s
	character(len=128) :: arg




	
	if (command_argument_count() .ne. 1) then
        print *, "Need number of time steps"
        call exit(-1)
    end if

	call get_command_argument(1, arg)
    Read(arg, '(i10)') nSteps

	allocate(J(nSteps))

	Open(1, file="input_attractor.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) X0
    Close(1)

	
	Open(1, file="param.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) s
    Close(1)
	
	Open(1, file="output_attractor.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')

	
	do t = 1, nSteps, 1
		call Xnp1(X0, X1, s)
		
		call Objective(X1,J(t),s)
		X0 = X1
		
	end do
	print *, X0
	Write(1) X1    
    Close(1)

    Open(1, file="J_history.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')
    Write(1) J
    Close(1)


end program flow
