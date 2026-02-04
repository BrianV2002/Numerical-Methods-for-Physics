close, /all


N = 100                
scale_factor = 4

si = dblarr(N, N)
sf = dblarr(N, N)


openr, 1, 'S_inicial.dat'
readf, 1, si
close, 1

openr, 2, 'S_final.dat'
readf, 2, sf
close, 2


window, 0, xsize=900, ysize=500
device, decomposed=0       
loadct, 0                 


; rebin: Resizes the 8x8 matrix to 400x400 so it looks big. 
; /sample: Keeps the "blocky" pixel look (nearest neighbor interpolation)
img_ini = rebin(si, N*scale_factor, N*scale_factor, /sample)
img_fin = rebin(sf, N*scale_factor, N*scale_factor, /sample)


tv, bytscl(img_ini, min=-1, max=1), 50, 50

tv, bytscl(img_fin, min=-1, max=1), 500, 50


; /device uses pixel coordinates
xyouts, 150, 20, 'E.Inicial', chars=2, color=255, /device
xyouts, 600, 20, 'E.Final', chars=2, color=255, /device

end
