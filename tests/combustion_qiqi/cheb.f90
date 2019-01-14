module cheb
	implicit none
	integer, parameter :: Ncheb = 5 ! this means that the interpolating 
!Lagrange polynomials have a maximum degree of Ncheb 
	REAL, PARAMETER :: Pi = 3.1415927
contains
	real(kind=8) function cheb_pts(k,n)
		implicit none
		integer:: k,n
		cheb_pts = cos(k*pi/n)
	end function cheb_pts

	subroutine step_forward_euler(u0, dudt, u1)
		implicit none
		real(kind=8), dimension(Ncheb + 1) :: u0, dudt, u1
		real(kind=8) :: dt
		integer :: k
	
		dt = 1.e-4
		do k = 1, Ncheb + 1, 1
			u1(k) = u0(k) + dt*dudt(k)
		end do
	end subroutine step_forward_euler
		
	
	function cheb_diff_matrix()
		implicit none
		real(kind=8), dimension(Ncheb + 1, Ncheb + 1) :: D
		real(kind=8), dimension(Ncheb + 1) :: x
		real(kind=8), dimension(Ncheb + 1, Ncheb + 1):: cheb_diff_matrix
		integer :: i,j,k
		D = 0.d0
	
		do k = 1, Ncheb + 1, 1
			x(k) = cheb_pts(k-1,Ncheb) 
		end do
		if(Ncheb > 1) then
			do i = 2, Ncheb, 1
				D(i,i) = -1.d0*x(i)/(1.d0 - x(i)**2.d0)/2.d0
				do j  = i+1, Ncheb, 1
					D(i,j) = (-1.d0)**(i+j)/(x(i)-x(j))
					D(j,i) = -1.d0*D(i,j)
				end do
			end do
		end if
		D(1,1) = (2.0*(Ncheb**2.0)+1.d0)/6.d0
		D(Ncheb+1,Ncheb+1) = -1.d0*D(1,1)
		do j = 2, Ncheb + 1, 1
			D(1,j) = 2.d0*(-1)**(j-1)/(x(1)-x(j))
		end do
		D(1,Ncheb+1) = D(1,Ncheb+1)/2.d0
		D(2:Ncheb+1,1) = -0.25d0*D(1,2:Ncheb+1)
		D(Ncheb+1,1) = 4.d0*D(Ncheb+1,1)
		do j = 2, Ncheb, 1
			D(Ncheb+1,j) = 2.d0*(-1)**(Ncheb+j-1)/(x(Ncheb+1)-x(j))
		end do
		D(2:Ncheb,Ncheb+1) = -1.d0*D(Ncheb+1,2:Ncheb)/4.d0 
		do j = 1, Ncheb + 1, 1
			do i = 1, Ncheb + 1, 1	
				cheb_diff_matrix(j,i) = D(Ncheb - j + 2,Ncheb-i + 2)
			end do
		end do
	end function cheb_diff_matrix 
end module cheb
