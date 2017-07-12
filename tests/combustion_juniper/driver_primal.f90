program flow
	use combustion_juniper_passive
	implicit none 
	double precision, dimension(d) :: X0, X1, zeroarray
	double precision, dimension(:), allocatable :: J
	integer :: t, nSteps, inttau
	double precision:: c1,c2,beta,xf,tau
	double precision, dimension(:,:), allocatable :: Xtmtau
	character(len=128) :: arg




	
	if (command_argument_count() .ne. 1) then
        print *, "Need number of time steps"
        call exit(-1)
    end if

	call get_command_argument(1, arg)
    Read(arg, '(i10)') nSteps

	allocate(J(nSteps))

	c2 = 0.01d0
	xf = 0.3d0
	beta = 0.75d0
	tau = 0.02d0


	Open(1, file="input_primal.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) X0
    Close(1)

	
	Open(1, file="param.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) c1
    Close(1)

	inttau = int(tau/dt)
	allocate(Xtmtau(d,inttau))
	do t = 1, d, 1
		X1(t) = X0(t) 
		zeroarray(t) = 0.d0
	end do

	do t = 1, inttau, 1
		call Xnp1(X1,Xnp1_res,zeroarray,c1,c2,beta,xf)
		do t1 = 1, d, 1
			X1(t1) = Xnp1_res(t1)
			Xtmtau(t1,t) = Xnp1_res(t1)		
		end do
	end do

		
	do t = 1, nSteps, 1
		call Xnp1(X1,Xnp1_res,Xtmtau(:,mod(t,inttau)),c1,c2,beta,xf)
		do t1 = 1, d, 1
			X1(t1) = Xnp1_res(t1)
			do t2 = 1, inttau-1, 1
				Xtmtau(t1,t2) = Xtmtau(t1,t2+1)		
			end do
			Xtmtau(t1,inttau) = Xnp1_res(t1)
		end do
		call Objective(Xnp1_res,J(t),xf)
	end do

	

	Open(1, file="output_primal.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')
    Write(1) X1
    Close(1)

    Open(1, file="objective.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')
    Write(1) J
    Close(1)


end program flow
