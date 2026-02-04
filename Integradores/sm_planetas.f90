module sm
implicit none

real :: m1, m2  

contains

subroutine dery(neq, t, y, dy)
    implicit none
    integer, intent(in) :: neq
    real, intent(in) :: t
    real, dimension(neq), intent(in) :: y
    real, dimension(neq), intent(out) :: dy

    real :: rx, ry, r, r3

    ! posiciones relativas 
    rx = y(5) - y(1)
    ry = y(7) - y(3)

    r  = sqrt( rx*rx + ry*ry )
    r3 = r**3


    ! cuerpo 1
    dy(1) = y(2)                       ! dx1/dt = vx1
    dy(2) = m2 * rx / r3              ! dvx1/dt
    dy(3) = y(4)
    dy(4) = m2 * ry / r3

    ! cuerpo 2
    dy(5) = y(6)
    dy(6) = -m1 * rx / r3               
    dy(7) = y(8)
    dy(8) = -m1 * ry / r3

end subroutine dery

end module sm

