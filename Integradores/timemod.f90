module timestep

    use sm
    
    contains    
    
    subroutine eulerf(neq,t,dt,y,dy)
    
        implicit none
        
        integer, intent(in) :: neq
        real, intent(in) :: t, dt
        real, dimension(neq), intent(inout) :: y
        real, dimension(neq), intent(out)   :: dy
        
        call dery(neq, t, y, dy)
        y = y + dt * dy
        
    end subroutine eulerf
    
    
    subroutine rk2(neq,t,dt,y,dy)
    
        implicit none
        
        integer, intent(in) :: neq
        real , intent(in) :: t, dt
        real , dimension(neq), intent(inout) :: y
        real , dimension(neq), intent(out)   :: dy 
        real , dimension(neq) :: ys
        real :: ts
        
        ! Paso 1
        call dery(neq, t, y, dy)
        
        ys = y + 0.5*dt*dy
        ts = t + 0.5*dt
        
        ! Paso 2
        call dery(neq, ts, ys, dy)
        
        y = y + dt*dy
        
    end subroutine rk2

end module

