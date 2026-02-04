program poisson

       implicit none
       
       integer,parameter :: m = 1, nx = 124
       
       real*8 :: l, pi, dx, k0, csi ,media
       
       real*8, dimension(nx) :: x, phi, phith , b
       real*8, dimension(nx) :: rho,rhoth
       
       real*8, dimension(nx) :: A1, A2, A3 ! Elementos de las diagonales
       real*8, dimension(nx-1) :: work2, work3, work4
       
       integer :: i , job, info
       
       
       l = 10.d0
       pi = acos(-1.d0)
       k0 = 2.d0*pi/l
       dx = l/dfloat(nx)
       
       
       do i = 1 , nx
            x(i) = (i-1)*dx
       enddo
       
       rho(1:nx) = -(m*k0)**2*dsin(m*k0*x)             !dado rho encuentro phi ,  A=1 (amplitud)
       phith(1:nx) = dsin(m*k0*x)                     ! en condisiones ideales
       
       A1 = 1.d0  ! todos los elemtos del vector se llena con 1
       A2 = -2.d0
       A3 = 1.d0
       
       !directamente ya factorizamos
       
       call s3p_fa(nx,A1,A2,A3,INFO,work2,work3,work4) !factorizacion
       
       !por periocidad
       
       b(1) = (5.d0/6.d0)*dx**2*(rho(2)/10.d0 + rho(1) + rho(nx)/10.d0)  
       b(nx) = (5.d0/6.d0)*dx**2*(rho(1)/10.d0 + rho(nx) + rho(nx-1)/10.d0)
       
        do i = 2 , nx -1
           b(i) = (5.d0/6.d0)*dx**2*(rho(i+1)/10.d0 + rho(i) + rho(i-1)/10.d0)
        enddo
        
        job = 0 
        
       call s3p_sl(nx,A1,A2,A3,b, phi, job, work2, work3, work4)
       
       media = sum(phi)/nx
       phi = phi - media
       
        do i = 1 , nx
             write(101,*) x(i), phith(i) , phi(i)
        enddo
        
       close(101)
       
       csi = (sum(phith - phi)**2/dfloat(nx))**0.5
       
       print*, 'error respecto a la teórica (M.pade):' , csi
       
       
       ! pilas que no tiene contains 
       
       

end program
