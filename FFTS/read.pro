close, /all

nx=128
nt=5425	

x = dblarr(nx)			; dimension de cada columna en fort.
fnum = dblarr(nx,nt+1)		; + cond inicial

name = string('fort.')

device, decomposed=0		;pantalla blanca
loadct, 5			;load color table

check = ' '			; teclado 

for it =0, nt do begin			; no nt-1 pq incluyo condicion inicial

	titulo = string('Timestep =', it)	;s
	if it le 899 then begin
		ind = string(format='(I03)', 100+it)
	endif
	
	if it gt 899 then begin
		ind = string(format='(IO4)', 100+it)
	endif
	
	namef = name + ind
	
	openr, 1, namef
	
	for ix = 0, nx-1 do begin
		readf, 1, a, b
		x(ix) = a
		fnum(ix,it) = b
	endfor
	
	plot, x, fnum(*,it), col=0, back=255, chars=2, psym=3, $
        xtitle='x', ytitle='f(x)', title=titulo, $
        xrange=[0, 10], yrange=[-1,1]


	
	read, check						;darle al teclado
	
	close,1

endfor

end
