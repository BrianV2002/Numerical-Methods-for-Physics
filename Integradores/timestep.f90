module timestep
 
    use segundo_miembro
    
    contains
    
    subroutine rk2(neq,t,dt,y,dy)
    
        implicit none
        
         
         integer, intent(in) :: neq
         real , intent(in) :: t, dt
         real , dimension(neq) :: y,dy 
         real , dimension(neq) :: ys
         real :: ts
         
         call dery(Neq, t, y, dy)
         
         
         ys(1:neq) = y(1:neq) + 0.5*dt*dy(1:neq)
         ts = t + 0.5*dt
         
         call dery(neq,ts,ys,dy)
         y(1:neq) = y(1:neq) + dt*dy(1:neq)
         
     end subroutine rk2
     
end module


          
      
