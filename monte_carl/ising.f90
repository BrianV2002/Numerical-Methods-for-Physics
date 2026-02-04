program Ising

   use rand


   implicit none

   integer, parameter :: N = 100
   real*8,  parameter :: JT = 1.5d0        !J/T -- facilidad de los espines al cambiar
   
   integer :: itprint, sem
   integer :: k1, k2, i, j
   integer, dimension(N,N) :: s             ! matriz de configuracion de spin
   
   
   real*8 :: z, m, A                        ! contadores y observables
   real*8 :: de, omega

   integer, parameter :: con = 1.d6       ! pasos Monte Carlo
   integer, parameter :: nterm = 2.d4  
   
   itprint = 100                            ! cada cuantos pasos se imprime
   sem = -123                          
 
 !  s = 1
  ! m = dfloat(N*N)                         
   z = 0.d0    !mediciones                   
   A = 0.d0     ! suma de m
   
   
   
	do i = 1, N
	 do j = 1, N
	     if (ran2(sem) .gt. 0.5d0) then
		 s(i,j) = 1
	     else
		 s(i,j) = -1
	     endif
	 enddo
	enddo
   
   m = sum(s) 
   
   
   open(20, file='S_inicial.dat', status='unknown')
      do i = 1, N
         write(20,*) (s(i,j), j=1,N)
      enddo
   close(20)
 

   open(10, file='datos.dat', status='replace')

   do j = 1, con

      ! seleccion aleatoria del spin
      
      k1 = int(ran2(sem)*N + 1)
      k2 = int(ran2(sem)*N + 1)

      de = 0.d0                             ! cambio de energia

      


      i = k1 - 1
      if (i .gt. 0) then
         de = de + s(i,k2)
      endif

      i = k1 + 1
      if (i .le. N) then
         de = de + s(i,k2)
      endif

      i = k2 - 1
      if (i .gt. 0) then
         de = de + s(k1,i)
      endif

      i = k2 + 1
      if (i .le. N) then
         de = de + s(k1,i)
      endif

      de = de * JT * 2.d0*dfloat(s(k1,k2))

      omega = exp(-de)                      !tasa de aceptacion


      if (omega .gt. 1.d0 .or. ran2(sem) .lt. omega) then
         m = m - 2.d0*dfloat(s(k1,k2))
         s(k1,k2) = -s(k1,k2)
      endif

      if (j .gt. nterm) then
         z = z + 1.d0
         A = A + abs(m)
      endif

  
      if (mod(j,itprint) .eq. 0 .and. j .gt. nterm) then
      
            write(10,*) j, A/z
            print*, '<A/z> = ', A/z
            
         endif


   enddo

   close(10)
   
   
   open(21, file='S_final.dat', status='replace')

  do i = 1, N
   write(21,*) (s(i,j), j=1,N)
  enddo

 close(21)


end program

