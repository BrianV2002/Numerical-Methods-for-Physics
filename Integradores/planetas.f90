program planetas

    use sm
    use timestep
    
    implicit none
    
    integer, parameter :: neq=8
    real :: dt, tend, tprint, t
    real :: x01, vx01, y01, vy01, x02, vx02, y02, vy02
    real :: ek, ep, en0, en, r
    real, dimension(neq) :: y, dy
    integer :: it, nend, nprint
    

    dt = 0.001 
    tend   = 1000.
    tprint = 0.01

    m1 = 5.
    m2 = 10.

    x01 = -1.66
    vx01 = -0.347
    y01 = 0
    vy01 = -0.179

    x02 = 13.33
    vx02 = 0.174
    y02 = -13
    vy02 = 0.089
    

    nend = nint(tend/dt)
    nprint = nint(tprint/dt)
    
    t = 0.

    y(1) = x01
    y(2) = vx01
    y(3) = y01
    y(4) = vy01
    y(5) = x02
    y(6) = vx02
    y(7) = y02
    y(8) = vy02

    r = sqrt((y(1)-y(5))**2 + (y(3)-y(7))**2)
    
    ep = -(m1*m2)/r
    ek = 0.5*(m1*(y(2)**2 + y(4)**2) + m2*(y(6)**2 + y(8)**2))
    
    en0 = ek + ep
    

    open(15, file='orbitas.dat', status='unknown')
    
    write(15,*) t, y(1:neq), 0.
    

    do it = 1, nend
    
        t = t + dt
    
        call rk2(neq, t, dt, y, dy)
    
        r = sqrt( (y(1)-y(5))**2 + (y(3)-y(7))**2 )
        ep = -(m1*m2)/r
        ek = 0.5*( m1*(y(2)**2 + y(4)**2) + m2*(y(6)**2 + y(8)**2) )
        en = ek + ep
        if (mod(it, nprint) .eq. 0) then
            write(15,*) t, y(1:neq), en-en0
        endif
    
    enddo
    
    close(15)
    
end program

