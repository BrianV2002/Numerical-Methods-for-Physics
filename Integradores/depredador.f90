program depredadores

     use segundo_miembro
     use timestep
     
     implicit none
     
     
     
     
     
     
     integer,parameter :: neq = 2
     real,dimension(neq) :: y ,dy
     
     real  :: c1
     real :: tend, tprint, dt ,t 
     integer :: nend, nprint ,it
     
     
       tend = 500. ! TIEMPO TOTAL DE LA SIMULACION
       tprint = 0.01 ! INTERVALO DE TIEMPO DE IMPRESION
       dt = 0.001  ! TIEMPO PARA EL METODO
       
       
     nend = nint(tend/dt) 
     nprint = nint(tprint/dt)
     
     
    
      !----
       a = 0.012
       b = 0.0001
       c1 = 0.15 ! parametro de mortalidad-- solo coloco para q no muera antes de comer
       
      !---
     
     
     t=0.
     
     y(1)= 500        !condiciones iniciales
     y(2)= 10
        
        
     open(16, file='pandemia.dat', status='unknown')   
     
     write(16,*) t, y(1:neq)
     
     do it = 1 , nend
     
     
         t = it*dt
             
             
        call rk2(neq,t,dt, y, dy)
        
        
        if (t .le. 60.) then !!cambiar el parametro aleatorio
             c = 0.
        
        else 
           
            c = c1
            
        endif
            
     
       if(mod(it,nprint) .eq. 0) then
         write(16,*) t ,y(1:neq)
       endif
     
    enddo
    
    close(16)
    
endprogram 
        
        

      

     
     
