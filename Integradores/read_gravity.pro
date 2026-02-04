close, /all


nt = 100001 

tiempo = dblarr(nt)
x1 = dblarr(nt)
vx1 = dblarr(nt)
y1 = dblarr(nt)
vy1 = dblarr(nt)
x2 = dblarr(nt)
vx2 = dblarr(nt)
y2 = dblarr(nt)
vy2 = dblarr(nt)
Etot = dblarr(nt)

work = dblarr(10)

openr, 1, 'orbitas.dat'

for it = 0, nt-1 do begin
    readf, 1, work
	tiempo(it) = work(0)
	x1(it) = work(1)
	vx1(it) = work(2)
	y1(it) = work(3)
	vy1(it) = work(4)
	x2(it) = work(5)
	vx2(it) = work(6)
	y2(it) = work(7)
	vy2(it) = work(8)
	Etot(it) = work(9)
   
endfor

close, 1


set_plot, 'ps'
device, file='grafica.eps', /color
loadct, 5

!p.multi = [0,1,2]    ; 1 columna, 2 filas

plot, x1, y1, col=0, back=255, /isotropic, $
      xtitle='x', ytitle='y'

oplot, x2, y2, col=120

xyouts, min(x1), max(y1), 'm1', col=0
xyouts, min(x2), max(y2), 'm2', col=120


plot, tiempo, Etot, col=0, back=255, $
      xtitle='t', ytitle='E - E0'


!p.multi = 0
device, /close
set_plot, 'x'
end

