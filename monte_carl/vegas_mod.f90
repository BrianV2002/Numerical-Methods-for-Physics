module vegas_mod
 
 use rand
 
 contains 
 
 subroutine vegas(region,ndim,fxn,init,ncall,itmx,nprn,tgral,sd,chi2a)
 
 
  implicit none

  integer :: ndim, init, ncall, itmx, nprn !in
  real*8  :: region(2*ndim)   !in
  real*8 :: tgral, sd, chi2a  !out

   real*8, external :: fxn !! uan funcion que esta en el main, vive afuera.

  integer, parameter :: ndmx=50, mxdim=10
  real*8, parameter :: ALPH=1.5d0, TINY=1.d-30

  integer :: i, it, j, k, mds, nd, ndo, ng, npg, idum
  integer :: ia(mxdim), kg(mxdim)
  
  real*8 :: calls, dv2g, dxg, f, f2, f2b, fb, rc, ti, tsi, wgt, xjac, xn, xnd, xo, schi, si, swgt
  real*8 :: d(ndmx, mxdim), di(ndmx, mxdim), dt(mxdim), dx(mxdim)
  real*8 :: r(ndmx), x(mxdim), xin(ndmx)
    

  real*8 :: xi(ndmx, mxdim)
  save xi !!por lo que el algoritmo vegas es adaptativo


   idum = 123
   
  if (init .le. 0) then
     mds = 1
     ndo = 1
     
    do j = 1, ndim 
    
     xi(1,j) = 1.d0
     
    enddo
   
  end if

  if (init .le. 1) then
     si = 0.d0
     swgt = 0.d0
     schi = 0.d0
  end if

  if (init .le. 2) then
     nd = ndmx
     ng = 1
     
     if (mds .ne. 0) then
     
        ng = int((ncall/2.d0 + 0.25d0)**(1.d0/ndim))
        mds = -1
            
           if (2*ng - ndmx .ge. 0) then
           npg = ng/ndmx + 1
           nd = ng/npg
           ng = npg*nd
        end if
     end if

     k = ng**ndim
     npg = max(ncall/k,2)
     calls = dfloat(npg)*dfloat(k) !! asegurar que guarde
     dxg = 1.d0/ng
     dv2g = (calls*dxg**ndim)**2/(npg/npg/(npg-1))
     xnd = nd !dflaot
     dxg = dxg*xnd
     xjac = 1.d0/calls

     do j = 1, ndim
        dx(j) = region(j+ndim) - region(j)
        xjac = xjac*dx(j)
     enddo

     if (nd .ne. ndo) then
      do i = 1 , max(nd,ndo)
      
       r(i) = 1.d0
      enddo
      
        do j = 1, ndim
           call rebin(dfloat(nd)/xnd, nd, r, xin, xi(:,j))
        enddo
        ndo = nd
     end if
     
      if (nprn .ge. 0) then
      
          write(*,*) ndim, calls, itmx, nprn , alph, mds, nd
          
          do j = 1, ndim
             write(*,*) ' xl(',j,')= ', region(j), ' xu(',j,')= ', region(j+ndim)
          end do
       end if
    end if

  do it = 1, itmx
  
 
     ti = 0.d0
     tsi = 0.d0
     
     do j = 1 , ndim 
     kg(j) = 1
       do i = 1 , nd 
       
     d(i,j)  = 0.d0
     di(i,j) = 0.d0
       enddo
       
     enddo
        
          

!Se hace así porque es la única forma eficiente de recorrer una matriz 
!de dimensiones variables sin usar recursividad compleja. 
!El iterate envuelve el proceso de 
!"Visitar una celda -> Calcular -> Buscar la siguiente".

       iterate: do

          fb = 0.d0
          f2b = 0.d0
    
          do k = 1, npg
             wgt = xjac
             do j = 1, ndim

                xn = (dfloat(kg(j)) - ran2(idum))*dxg + 1.d0
                ia(j) = max(min(int(xn),ndmx),1)
    
                if (ia(j) .gt. 1) then
                   xo = xi(ia(j),j) - xi(ia(j)-1,j)
                   rc = xi(ia(j)-1,j) + (xn-dfloat(ia(j)))*xo
                else
                   xo = xi(ia(j),j)
                   rc = (xn-dfloat(ia(j)))*xo
                end if
    
                x(j) = region(j) + rc*dx(j)
                wgt = wgt*xo*xnd
                
             enddo
    
             f = wgt*fxn(x,wgt)
             f2 = f*f
             fb = fb + f
             f2b = f2b + f2
    
             do j = 1, ndim
                di(ia(j),j) = di(ia(j),j) + f
                if (mds .ge. 0) d(ia(j),j) = d(ia(j),j) + f2
             end do
          end do
    
          f2b = sqrt(f2b*dfloat(npg))
          f2b = (f2b-fb)*(f2b+fb)
          
          if (f2b .le. 0.d0) f2b = TINY
    
          ti = ti +fb
          tsi = tsi + f2b
          
          if (mds .lt. 0) then
               do j = 1, ndim
                  d(ia(j),j) = d(ia(j),j) + f2b
               end do
          end if
    

          do k = ndim,1,-1
             kg(k) = mod(kg(k),ng) + 1
             
             ! Si el contador avanza sin reiniciar, volvemos al inicio del grid_loop
             if (kg(k) .ne. 1) cycle iterate
          end do
          
          
          exit iterate
    
       enddo iterate
      
     tsi = tsi*dv2g
     wgt = 1.d0/tsi
     si = si + wgt*ti
     schi = schi + wgt*ti*ti
     swgt = swgt + wgt

     tgral = si/swgt
     chi2a = max((schi-si*tgral)/(it-0.99d0),0.d0)
     sd = sqrt(1.d0/swgt)
     
     tsi = sqrt(tsi)
     
     if (nprn .ge. 0) then
          write(*,*)   it ,' integral =', ti, '+-', tsi
          write(*,*)  tgral, '+/- ', sd, ' chi2a =', chi2a
       end if
     
      do j = 1, ndim
      !Refine the grid. Consult references to understand the subtlety
       xo=d(1,j)

       xn=d(2,j)

       d(1,j)=(xo+xn)/2.

       dt(j)=d(1,j)
       
          do  i=2,nd-1
     rc=xo+xn
     xo=xn
     xn=d(i+1,j)
     d(i,j)=(rc+xn)/3.
     dt(j)=dt(j)+d(i,j)
          enddo 
          
    d(nd,j)=(xo+xn)/2.
    dt(j)=dt(j)+d(nd,j)
    
       enddo 
       
         do  j=1,ndim
         
    rc=0.d0
    
    do  i=1,nd
    
    if(d(i,j).lt.TINY) d(i,j)=TINY
    
     r(i)=((1.d0 -d(i,j)/dt(j))/(log(dt(j))-log(d(i,j))))**ALPH
     rc = rc +r(i)
    
    enddo 
    
      call rebin(rc/xnd,nd,r,xin,xi(:,j))
    
      enddo 
    
      enddo 
       
    
end subroutine vegas




  subroutine rebin(rc, nd, r, xin, xi)
    implicit none
    integer :: nd
    real*8:: rc
    real*8, dimension(:) :: r
    real*8,dimension(:) :: xin
    real*8 ,dimension(:) :: xi
    integer :: i, k
    real*8 :: dr, xn, xo

    k = 0
    xo = 0.0d0
    dr = 0.0d0

    do i = 1, nd-1
       do while (rc .gt. dr)
          k = k + 1
          dr = dr + r(k)
       end do
       if (k .gt. 1) xo = xi(k-1)
       xn = xi(k)
       dr = dr - rc
       xin(i) = xn - (xn-xo)*dr/r(k)
    end do

    do i = 1, nd-1
       xi(i) = xin(i)
    end do
    xi(nd) = 1.0d0
  end subroutine rebin

end module
