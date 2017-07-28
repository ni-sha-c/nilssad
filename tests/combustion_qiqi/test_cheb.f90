program test_cheb
	use cheb
	implicit none
	real(kind=8), dimension(Ncheb+1,Ncheb+1) :: D
	real(kind=8), dimension(Ncheb+1) :: u0, u1, dudt
	integer :: steps,i
	u0 = 0.d0
	u0(1:Ncheb/2) = 1.d0
	call cheb_diff_matrix(D)
	steps = 1000
	open(1, file="advection_solution.bin", form="unformatted", access="stream", status="replace", convert="big_endian")	
	do i = 1, steps, 1
		dudt = -1.d0*matmul(D,u0)
		call step_forward_euler(u0,dudt,u1)
		u0 = u1
		u0(Ncheb+1) = 1.d0
		write(1) u0
	end do
	print *, "differentiation matrix for N = 2"


end program test_cheb
