close, /all

nx = 124

x = dblarr(nx)
phith = dblarr(nx)
phi = dblarr(nx)


name = 'fort.101'

device, decomposed = 0
loadct, 5

openr, 1, name

for ix = 0, nx-1 do begin
    readf, 1, a, b, c
    x(ix) = a
    phith(ix) = b
    phi(ix) = c
endfor

plot, x, phi, col=0, back=255, chars=2, $
    xtitle='x', ytitle='f(x)'

oplot, x, phith, col=125, ps=2

close, 1

end
