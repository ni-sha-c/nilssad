program driver
	use OAD_active
	use OAD_rev
	use lorenz63_passive
	implicit none 
	external head_homogeneous
	integer, parameter:: d=3
	type(active), dimension(:), allocatable :: x
    type(active), dimension(:,:), allocatable :: y
	double precision, dimension(d) :: X0, X1, Xorig, Xtemp1, Xtemp2
	double precision, dimension(:,:), allocatable :: v0,v,vorig
	integer :: t, nSteps, subspace_dimension, t1, t2
	double precision :: s, eps
	character(len=128) :: arg1, arg2


	our_rev_mode%tape=.TRUE.
	eps = 1.d-4

	
	if (command_argument_count() .ne. 2) then
        print *, "Need number of time steps and subspace dimension"
        call exit(-1)
    end if

	call get_command_argument(1, arg1)
    Read(arg1, '(i10)') nSteps

	call get_command_argument(2, arg2)
    Read(arg2, '(i10)') subspace_dimension


	Open(1, file="input_primal.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) X0
    Close(1)
	Xorig = X0


	
		
	allocate(v0(d,subspace_dimension))
	allocate(vorig(d,subspace_dimension))
	allocate(v(d,subspace_dimension))

	allocate(y(d,subspace_dimension))
	allocate(x(subspace_dimension))
	
	Open(1, file="input_htangents.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    do t = 1, subspace_dimension, 1
		Read(1) v0(:,t)		
    end do
	Close(1)
	vorig = v0


	Open(1, file="param.bin", form="unformatted", access="stream", &
            status="old", convert='big_endian')
    Read(1) s
    Close(1)


		
	Open(1, file="output_htangents.bin", form="unformatted", access="stream", &
         status='replace', convert='big_endian')
	
	Xtemp1 = Xorig
	do t = 1, nSteps, 1
		call Xnp1(Xtemp1,X1,s)
		Xtemp1 = X1
	end do	
 
	
	do t = 1, subspace_dimension, 1
		do t2 = 1, d, 1
			y(t2,t)%d = 0.d0
			y(t2,t)%d(t2) = 1.d0
		end do
		x(t)%v = eps
		

		call head_homogeneous(x(t)%v,s,y(:,t),X0,v0(:,t),nSteps)

		Write(1) x(t)%d		

		Xtemp2 = Xorig + eps*v0(:,t)
		do t2 = 1, nSteps, 1
			call Xnp1(Xtemp2,X1,s)
			Xtemp2 = X1
		end do	

		print *, "From AD: ", x(t)%d
		print *, "From FD: ", (Xtemp2 - Xtemp1)/eps



	end do
    Close(1)
		

end program driver
