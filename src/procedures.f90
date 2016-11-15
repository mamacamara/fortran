!module procedures
module constants
	implicit none 
	integer, parameter :: precI=kind(1), precR=kind(1.0)
	real ( kind = precR ) , save :: pi = 4*atan(1.0) !Declaration statique de Pi 
end module constants
!	use constants 
	! Densite de proba de la loi normale
	! De moyenne mu et decart type sigma
	function loi_normale(xx,x_mu,sigma) result(j)
		use constants 
		implicit none
		real(kind=precR)  :: xx, sigma,x_mu,j
		j = (1/(sigma*sqrt(2*pi)))*exp((-1/2)*((xx-x_mu)/sigma)**2)
	end function loi_normale
	! Integration par la methode des rectangles	
	! de loi_normale(,mu,sigma) sur [x_min;x_max] avec le pas h_step
	! le resultat est mis sur result
	subroutine methode_rect(loi_normale,sigma,x_mu,x_min,x_max,h_step,resultat)
		use constants 
	 	implicit none  
		real (kind = precR) ::  f0,fN,loi_normale,xi
		real ( kind = precR), intent(in) :: x_min
		real ( kind = precR), intent(in) :: x_max
		real ( kind = precR), intent(in) :: h_step
		real ( kind = precR), intent(out) :: resultat
	        integer ( kind = precI) ::  N   ,i=1                                 
        	real ( kind = precR), intent(in) :: x_mu, sigma
		N = (x_max-x_min)*h_step 
		f0 = loi_normale(x_min,x_mu,sigma)
		fN = loi_normale(x_max,x_mu,sigma)
		!i = 1
		resultat = loi_normale(x_min,x_mu,sigma) ! h*f(x_min)
		do while (i<N)
			xi = x_min + h_step * i 
			resultat = resultat +  loi_normale(xi,x_mu,sigma)
			i = i + 1 
		end do 
		resultat = resultat * h_step 
	end subroutine methode_rect   
program mon_main 
	use constants 
!	implicit none
	implicit none
	real(kind = precR):: fonction,loi_normale
	fonction = loi_normale(1.0,0.0,1.0)
	write(*,*) "La valeur de la loi normale centrée réduite au point 0. est"
	write(*,*) fonction
end program mon_main 

			
		
		
	

		
	
!end mocleardule procedures 
