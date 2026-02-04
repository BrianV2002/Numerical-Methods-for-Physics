close, /all

nx = 9800

A = dblarr(nx)
itprint = dblarr(nx)

device, decomposed=0		;hacer pantalla blanca
loadct, 5			;load color table

name = 'datos.dat'

openr, 1, name

for ix = 0, nx-1 do begin
   readf, 1, step, mag
   itprint(ix) = step
   A(ix) = mag
endfor

close, 1


plot, itprint, A, col=0, back=255, chars=2, $
      xtitle='MC step', ytitle='<|M|>', ystyle=1




end
