!ir y regresar fourier 
!saco primera y segunda derivada

module bffourier
	contains
		
		!subroutine FFTPS(N,k0,f, bp)
		subroutine FFTPS(N,k0,f, bp,bs)		!f= funcion, bp=primera der, bs=segunda deri			
			implicit none
			
			integer, intent(in) ::n
			real*8 :: k0
			real*8, dimension(n) :: f
			real*8, dimension(n) :: bp, bs
			
			integer :: i, j
			real*8, dimension(n) :: a, as, ap 	!vector de trabajo
			
			real*8, dimension(2*n+15) :: wsave
	
			call drffti(n,wsave)		!caracteristicas de wsave 
			
			a=f				!para conservar funcion original
				
			call drfftf(n,a,wsave)		!llevo fun a fourier
			
			a=a/n				!def de coeficientes de fourier 
			
			!call drfftb(n,a,wsave)		!regreso espacio fisico	(prueba pa ver ida y venida de fourier)
			!bp=a
			
!----------------------------------------------------------------------------------------------------------------------------------
!primera rerivada
!-----------------------------------------------------------------------------------------------------------------------------------
			
			ap(1) = 0.d0			!primera derivada que era mean value (cte)
			ap(n) = 0.d0			!ultima derivada por razones numericas (estabilidad)
			
			!even odd build
			
				do i=2, n-1, 2
					j=i+1
					
					ap(i)=-(i/2)*k0*a(j)
					ap(j)=-(i/2)*k0*a(i)
				enddo
				
			bp=ap
			
			call drfftb(n,bp, wsave)	
			
!-----------------------------------------------------------------------------------------------------------------------------------
!segunda derivada
!-----------------------------------------------------------------------------------------------------------------------------------

			as(1) = 0.d0			!primera derivada que era mean value (cte)
			as(n) = 0.d0			!ultima derivada por razones numericas (estabilidad)
			
			!even odd build
			
				do i=2, n-1, 2
					j=i+1
					
					as(i)=-((i/2)*k0)**2*a(i)
					as(j)=-((i/2)*k0)**2*a(j)
				enddo
				
			bs=as
			
			call drfftb(n,bs, wsave)	




			
		end subroutine

end module
