close, /all

nt = 50001

t = dblarr(nt)
y1 = dblarr(nt)
y2 = dblarr(nt)

device, decomposed=0		;hacer pantalla blanca
loadct, 5			;load color table


work = dblarr(3)

openr , 1, 'pandemia.dat' 

 for it = 0 , nt-1 do begin
 
    
      readf, 1, work
      t(it) = work(0)
      y1(it) = work(1)
      y2(it) = work(2)


 endfor
 

  !p.multi = [0, 1, 2]   

  plot, t, y1, color=0, background=255, chars=2, xs = 1,$
        xtitle='t', ytitle='y1(t)'

  oplot, t, y2, color=125, $
         thick=2
 
  plot,  y1 ,y2 , color = 0 , background= 255 , chars=2
 
  !p.multi = 0

 

 end

 
 
 
  
   
