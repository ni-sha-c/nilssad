program flow
	use combustion_juniper_passive
	implicit none 
	double precision, dimension(d) :: X0, X1, zeroarray, Xnp1_res
	double precision, dimension(:), allocatable :: J
	integer :: t, nSteps, inttau, t1, t2, counter
	double precision:: param_active, tau
	double precision, dimension(N_p) :: params
	double precision, dimension(N_p-1) :: params_passive
	double precision, dimension(:,:), allocatable :: Xtmtau
	character(len=128) :: arg




	
	if (command_argument_count() .ne. 1) then
        print *, "Need number of time steps"
        call exit(-1)
    end if

	call get_command_argument(1, arg)
    Read(arg, '(i10)') nSteps

	allocate(J(nSteps))


	Open(1, file="input_primal.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) X0
    Close(1)

	Open(1, file="param.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
	do t = 1, N_p, 1
    	Read(1) params(t)
	end do
    Close(1)
	param_active = params(1)
	params_passive = params(2:N_p)
	tau = params_passive(N_p-1)
	inttau = int(tau/dt)
	allocate(Xtmtau(d,inttau))
	do t = 1, d, 1
		X1(t) = X0(t) 
		zeroarray(t) = 0.d0
	end do

	do t = 1, inttau, 1
		call Xnp1(X1,Xnp1_res,zeroarray,param_active,params_passive)
		do t1 = 1, d, 1
			X1(t1) = Xnp1_res(t1)
			Xtmtau(t1,t) = Xnp1_res(t1)		
		end do
	end do
	counter = 0		
	do t = 1, nSteps, 1
		counter = counter + 1
		call Xnp1(X1,Xnp1_res,Xtmtau(:,counter),param_active,params_passive)
		do t1 = 1, d, 1
			X1(t1) = Xnp1_res(t1)
			do t2 = 1, inttau-1, 1
				Xtmtau(t1,t2) = Xtmtau(t1,t2+1)		
			end do
			Xtmtau(t1,inttau) = Xnp1_res(t1)
		end do
		call Objective(Xnp1_res,J(t),param_active,params_passive)
		if(counter .eq. inttau) then
			counter = 0
		end if
	end do

	
	!print *, X1
	Open(1, file="output_primal.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')
    Write(1) X1
    Close(1)

    Open(1, file="objective.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')
    Write(1) J
    Close(1)


end program flow
