program prueba

    use bffourier

    implicit none

    integer, parameter :: m = 1, nx = 2048
    integer :: i, info
    real*8, dimension(nx) :: x, phi
    real*8, dimension(nx) :: rho, rhose, rhopr, rhop, rhos
    real*8 :: k0, l, pi, dx, csi, csis

    l  = 10.d0
    pi = dacos(-1.d0)
    k0 = 2.d0 * pi / l
    dx = l / float(nx)

    do i = 1, nx
        x(i) = (i - 1) * dx
    end do

    rho(1:nx)   =  dsin(m * k0 * x)
    rhopr(1:nx) =  m*k0 * dcos(m * k0 * x)        ! primera derivada exacta
    rhose(1:nx) = -(m*k0)**2 * dsin(m*k0*x)       ! segunda derivada exacta
  
    call fftps(nx,k0,rho, rhop,rhos)
    
    csi  = sqrt( sum((rhopr - rhop)**2)  / float(nx) )
    csis = sqrt( sum((rhose - rhos)**2) / float(nx) )
    
    print *, "csi primera derivada  = ", csi
    print *, "csi segunda derivada  = ", csis

end program

