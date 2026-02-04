program main

  use vegas_mod
  use rand   
  
  implicit none

!! por ser varaibles locales

  integer :: ndim, ncall, itmx, nprn, init, sem
  real*8 :: tgral, sd, chi2a
  
  real*8 :: region(20) 
  

  real*8, external :: fxn


  ndim = 1        !Dimensiones de la integral
  ncall = 10000    !Puntos por iteracion
  itmx = 1000      !Numero de iteraciones
  nprn = 1       !Mostrar output
  init = -1     ! -1 = Inicializacion completa




  
  ! Integral de 0.0 a 1.0
  region(1) = 0.d0       !Minimo x
  region(1+ndim) = 1.d0  !Maximo x
  
  call vegas(region, ndim, fxn, init, ncall, itmx, nprn, tgral, sd, chi2a)


  write(*,*) 'Valor de la integral: ', tgral
  write(*,*) 'Error estimado:', sd
  write(*,*) 'Chi2', chi2a

end program main

!!afuera del main porque es funcion externa
!mejor hacerlo en modulo

function fxn(x, wgt)
  implicit none
  
  real*8 :: fxn
  real*8 :: x(*)   ! Vector de coordenadas
  real*8  :: wgt    

  ! Como es 1D, usamos x(1)
  fxn = 4.d0/(1.d0 + x(1)**2)
  
end function fxn

