program monte_carlo

    use rand
    
    implicit none

    integer, parameter :: N = 10000000
    real*8 :: error, e, y, x
    integer :: i
    integer :: sem
    
    

    e = 0.0d0
    error = 0.0d0
    sem = -1234


    do i = 1, N      
        x = ran1(sem)
        y = 4.0d0/(1.0d0 + x*x)
        e = e + y
       error = error + y*y
    end do

    e = e/N
    error = dsqrt(error/N-(e**2))/dsqrt(dfloat(N))
    
    print*,e
    print*, ' '
    print*, error
    
   
    open(15, file='valor7.dat', status='unknown')
    
    write(15,*) N, e, error
    
    
    close(15)

end program monte_carlo

