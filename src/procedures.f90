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
		N = (x_max-x_min)*h_step ! Le nombre de points discret 
		f0 = loi_normale(x_min,x_mu,sigma) !f(a)
		fN = loi_normale(x_max,x_mu,sigma) ! f(b)
		!i = 1
		resultat = loi_normale(x_min,x_mu,sigma) ! on commence par f(x_min)=f(x_0) pour la methode du point de gauche
		do while (i<N)
			xi = x_min + h_step * i! retrouve le prochain point discret 
			resultat = resultat +  loi_normale(xi,x_mu,sigma)
			i = i + 1 
		end do 
		resultat = resultat * h_step 
	end subroutine methode_rect   
	subroutine methode_trap(loi_normale,sigma,x_mu,x_min,x_max,h_step,resultat)
		use constants 
	 	implicit none  
		real (kind = precR) ::  f0,fN,loi_normale,xi,somme
		real ( kind = precR), intent(in) :: x_min
		real ( kind = precR), intent(in) :: x_max
		real ( kind = precR), intent(in) :: h_step
		real ( kind = precR), intent(out) :: resultat
	        integer ( kind = precI) ::  N   ,i = 0                                 
        	real ( kind = precR), intent(in) :: x_mu, sigma
		N = (x_max-x_min)*h_step ! Le nombre de points discret 
		f0 = loi_normale(x_min,x_mu,sigma) !f(a)
		fN = loi_normale(x_max,x_mu,sigma) ! f(b)
		!i = 1
		
		somme = 0 ! on commence par f(x_min)=f(x_0) pour la methode du point de gauche
		do while (i<N)
			xi = x_min + h_step * i! retrouve le prochain point discret 
			somme = somme +  loi_normale(xi,x_mu,sigma)
			i = i + 1 
		end do 
		somme = 2*somme
		resultat = (h_step/2)*(f0+fN+somme)
		 
	end subroutine methode_trap 
program mon_main 
	use constants 
!	implicit none
	implicit none
	real(kind = precR):: fonction, loi_normale, valeur, x1, x2, moy, ecart, pas,maval
	fonction = loi_normale(1.0,0.0,1.0)
	write(*,*) "La valeur de la loi normale centrée réduite au point 0. est"
	write(*,*) fonction
	x1 = -1.0
	x2 = 1.0 
	moy = 1.0 
	ecart = 2 
	pas = 0.01
	open (unit=10,file=donnee.dat,action="read",status="old") 
  	call methode_trap(loi_normale,ecart,moy,x1,x2,pas,maval)
  	call methode_rect(loi_normale,ecart,moy,x1,x2,pas,valeur)
	write(*,*) "La valeur de l'intégrale de la loi normale par la methode des rectangles est: ",valeur 
	write(*,*) "La valeur de l'intégrale de la loi normale par la methode des trapèzes est: ",maval 
		
end program mon_main 

			
		
	

		
	
!end mocleardule procedures 
