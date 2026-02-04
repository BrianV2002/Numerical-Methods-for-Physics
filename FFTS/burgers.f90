!para compilar 
!gfortran -c -fallow-argument-mismatch fft1d.f modulo.f90 burgers.f90
!gfortran -fallow-argument-mismatch fft1d.f modulo.f90 burgers.f90
program burgers

	use bffourier
	implicit none
	
	integer, parameter :: m = 2, nx = 128
	integer :: i
	real*8, dimension(nx) :: x, f, faux, ap, as
	real*8 :: k0, L, pi, dx, nu, kmax, tmax, dt, t
	integer :: nt, itprint, it, cont

	L = 10.d0
	pi = dacos(-1.d0)
	k0 = 2.d0*pi/l
	dx = L/float(nx)
	nu = 0.1d0
	
	!Condición inicial
	do i = 1, nx
		x(i) = dfloat(i-1)*dx
		f(i) = dsin(m*k0*x(i))
		write(100,*) x(i), f(i)
	enddo
	close(100)
	
	!GDL
	kmax = (dfloat(nx)/2.d0)*k0
	dt = 0.5d0*(2.d0*nu/((maxval(f))*2+nu*2*kmax*2))
	tmax = 30.d0
	nt = int(tmax/dt)
	itprint = int(nt/300)
	
	cont = 1
	t = 0.d0
	
	do it = 1, nt
	
		t = t + dt
		
		call FFTPS(nx,k0,f ,ap,as)
		
		do i = 1, nx
			faux(i) = f(i) - dt*(f(i)*ap(i)-nu*as(i))
		enddo
		
		f = faux
		
		if((it/itprint)*itprint.eq.it) then
			do i = 1, nx
				write(100+cont,*) x(i), f(i), ap(i), as(i)
			enddo
			close(100+cont)
			cont = cont + 1
		endif
		
	enddo
	
	print*, 'nt =', nt
	print*, 'itprint =', itprint
	print*, 'dt =', dt
	print*, 'cont =', cont
	print*, 'FIN DEL PROGRAMA'
	

endprogram
