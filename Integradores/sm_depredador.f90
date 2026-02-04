module segundo_miembro 

!implicit none

      real :: a, b, c


    
    contains 
    
             subroutine derY(neq,t,y,dy)
             
             !implicit none
             
              
             integer , intent(in) :: neq
             real :: t
             real, dimension(neq), intent(in) :: y
             real, dimension(neq), intent(out) :: dy
             
             !coordenada --> derivada --> coordenada
             
             dy(1) = a*y(1) - b*y(2)*y(1)
             dy(2) = b*y(1)*y(2) - c*y(2)
             
             end subroutine derY
             
end module
