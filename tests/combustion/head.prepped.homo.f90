subroutine head_homogeneous(dtau,inttau,y,X0,v0,nSteps)
	use combustion_juniper
	implicit none 
	double precision, dimension(d) :: X
	double precision, dimension(d) :: Xnp1_res
	integer, intent(in) :: inttau
	double precision, dimension(d), intent(in) :: X0, v0
	double precision, dimension(d), intent(out) :: y
	integer :: t, t1, t2
	integer, intent(in) :: nSteps
	double precision, intent(in) :: ds
	double precision, dimension(d,inttau) :: Xtmtau = 0.d0


!$openad INDEPENDENT(ds)

	
	do t = 1, d, 1 
		X(t) = X0(t) + v0(t)*dtau
	end do

	do t = 1, inttau, 1
		call Xnp1(X,Xnp1_res,0.d0)
		do t1 = 1, d, 1
			X(t1) = Xnp1_res(t1)
			Xtmtau(t1,t) = Xnp1_res(t1)		
		end do
	end do
	
	do t = 1, nSteps - inttau, 1
		call Xnp1(X,Xnp1_res,Xtmtau(:,mod(t,inttau)))
		do t1 = 1, d, 1
			X(t1) = Xnp1_res(t1)
			do t2 = 1, inttau-1, 1
				Xtmtau(t1,t2) = Xtmtau(t1,t2+1)		
			end do
			Xtmtau(t1,inttau) = Xnp1_res(t1)
		end do
	end do
	do t = 1, d, 1
		y(t) = X(t)
	end do
	!y(1) = 1.d0
	!y(2) = s
	!y(3) = 1.d0
!$openad DEPENDENT(y)

end subroutine


