module rand


  contains 
  

  real function ran0(idum)
  
    implicit none

    integer, intent(inout) :: idum
    integer :: k
    
    !integer , parameter :: ia, im, iq, ir, mask
    !real, integer :: am
   

    integer, parameter :: ia = 16807, im = 2147483647,  iq = 127773, ir = 2836 , mask = 123459876
    real, parameter :: am = 1./im
    
    
    
        !“Minimal” random number generator of Park and Miller. Returns a uniform    random deviate
!between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
!to initialize the sequence; idum must not be altered between calls for successive deviates
!in a sequence.


    idum = ieor(idum, mask)
    k = idum/iq
    idum = ia*(idum - k*iq) - ir*k


    if (idum .eq. 0) then
    idum = idum + im
    
    endif

    ran0 = am*idum
    idum = ieor(idum,mask)

  end function ran0
  
  
  
  !---------------------------------------
  
  
  real function ran1(idum)
  
    implicit none
    
    




    
    integer, intent(inout) :: idum

    integer, parameter :: ia = 16807, im = 2147483647 ,iq = 127773
    integer, parameter :: ir = 2836 , ntab = 32

    integer, parameter :: ndiv = 1 +(im-1)/ntab

    real,parameter :: am = 1.0/im , eps = 1.2e-7, rnmx = 1.- eps
    
    integer, save :: iv(ntab) = 0
    integer, save :: iy = 0

    integer :: j, k 
    !integer :: ,iv(ntab), iy
    
    
    !“Minimal” random number generator of Park and Miller with Bays-Durham shuffle and
!added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
!the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
!alter idum between successive deviates in a sequence. RNMX should approximate the largest
!floating value that is less than 1.

   
   !iv = 0
   !iy = 0


    !save iv, iy
    !data iv /ntab*0/, iy /0/


    if (idum .le. 0 .or. iy .eq. 0) then !Initialize.
       idum = max(-idum, 1)  !Be sure to prevent idum = 0.

       do j = ntab+8, 1, -1 !Load the shuffle table (after 8 warm-ups).
          k = idum/iq
          idum = ia*(idum - k*iq) - ir*k
          if (idum .lt. 0) then
          idum = idum + im
          endif
          if (j .le. ntab)then
           iv(j) = idum
           endif
       enddo
       iy = iv(1)
    end if


    k = idum/iq  !Start here when not initializing.
    idum = ia*(idum - k*iq) - ir*k  ! Compute idum=mod(IA*idum,IM) without overflows by Schrage’s method.
    if (idum .lt. 0) idum = idum + im !Will be in the range 1:NTAB.

    j  = 1 + iy/ndiv
    iy = iv(j)  !Output previously stored value and refill the shuffle table.
    iv(j) = idum

    ran1 = min(am*iy, rnmx)

  end function ran1
  
  
  

  real*8 function ran2(idum)
    implicit none
    


    integer:: idum


    integer, parameter :: im1 = 2147483563, im2 = 2147483399, ia1 = 40014, ia2 = 40692 
    integer, parameter :: iq1 = 53668, iq2 = 52774, ir1 = 12211, ir2 = 3791
    integer, parameter :: ntab = 32
    integer, parameter :: imm1 = im1 - 1
    integer, parameter :: ndiv = 1 + imm1/ntab


    real*8, parameter :: am = 1.0d0/2147483563.0d0
    real*8, parameter :: eps = 1.2d-16
    real*8, parameter :: rnmx = 1.0d0 - eps

    integer :: j, k
    integer, save :: idum2 = 123456789
    integer, save :: iy = 0
    integer, save :: iv(ntab) = 0


    if (idum .le. 0) then
       idum = max(-idum, 1)
       idum2 = idum

       do j = ntab + 8,1, -1
          k = idum/iq1
          idum = ia1 * (idum - k*iq1) - k*ir1
          if (idum .lt. 0) idum = idum + im1
          if (j .le. ntab) iv(j) = idum
       end do

       iy = iv(1)
    end if


    k = idum/iq1
    idum = ia1*(idum - k*iq1) - k*ir1
    
    if (idum .lt. 0) idum = idum + im1

    k = idum2/iq2
    idum2 = ia2*(idum2 - k*iq2) - k*ir2
    if (idum2 .lt. 0) idum2 = idum2 + im2

    j = 1 + iy/ndiv
    iy = iv(j) - idum2
    iv(j) = idum

    if (iy .lt. 1) iy = iy + imm1

    ran2 = min(am*dfloat(iy), rnmx)

  end function ran2





    real function expdev(idum)
    
    implicit none
    
    integer, intent(inout) :: idum
    !integer :: idum
    real :: dum
    
    

    if (dum .eq. 0) then
         dum = ran1(idum)    
    endif
    
  expdev = -log(dum)

    end function expdev
 

real function gasdev(idum)

  implicit none

  integer, intent(inout) :: idum
  real :: fac,rsq, v1, v2
  integer , save :: gset
  integer, save :: iset = 0


  if (idum .lt. 0) iset = 0

  if (iset .eq. 0) then
     do
        v1 = 2.0*ran1(idum) - 1.0
        v2 = 2.0*ran1(idum) - 1.0
        rsq = v1*v1 + v2*v2
        if (rsq .lt. 1.0 .and. rsq .ne. 0.0) exit
     end do

     fac = sqrt(-2.0*log(rsq)/rsq)
     gset = v1*fac
     gasdev = v2*fac
     iset = 1
  else
     gasdev = gset
     iset = 0
  end if

end function gasdev


 real function gamdev(ia, idum)
  implicit none

  integer, intent(in)    :: ia
  integer, intent(inout) :: idum
  integer :: j
  real :: am, e, s, v1, v2, x, y

  if (ia .lt. 1) stop "bad argument in gamdev"

  if (ia .lt. 6) then
     x = 1.0
     do j = 1, ia
        x = x * ran1(idum)
     end do
     gamdev = -log(x)
  else
  
     am = ia - 1.0
     s  = sqrt(2.0*am + 1.0)

     do
        v1 = ran1(idum)
        v2 = 2.0*ran1(idum) - 1.0
        if (v1*v1 + v2*v2 .le. 1.0) then
        
        y = v2 / v1
        x = s*y + am
        
        if (x .gt. 0.0) then

        e = (1.0 + y*y) * exp(am*log(x/am) - s*y)
        
        if (ran1(idum) .le. e) then
        exit
        endif
        
       endif
      endif 
     end do

     gamdev = x
  end if

end function gamdev



real function gammln(xx)
    implicit none
    real, intent(in) :: xx
    integer :: j
    real :: ser, stp, tmp, x, y
    real, dimension(6) :: cof

    cof(1)=76.18009172947146
    cof(2)=-86.50532032941677
    cof(3)=24.01409824083091
    cof(4)=-1.231739572450155
    cof(5)=0.001208650973866179
    cof(6)=-0.000005395239384953

    stp=2.5066282746310005

    x = xx
    y = x
    tmp = x + 5.5
    tmp = (x + 0.5)*log(tmp) - tmp
    ser = 1.000000000190015

    do j = 1, 6
        y = y + 1.0
        ser = ser + cof(j)/y
    end do

    gammln = tmp + log(stp*ser/x)
end function gammln



real function poidev(xm, idum)
  implicit none

  real, intent(in) :: xm
  integer, intent(inout) :: idum
  real :: pi  
  real :: alxm, em, g, sq, t, y
  real, save :: oldm = -1.0
  
  pi = acos(-1.)

  if (xm .lt. 12.) then

     if (xm .ne. oldm) then
        oldm = xm
        g = exp(-xm)
     end if

     em = -1.
     t  = 1.
     
     do
        em = em + 1.0
        t  = t * ran1(idum)
        if (t .le. g) exit
        
     end do

  else

     if (xm .ne. oldm) then
        oldm = xm
        sq = sqrt(2.0*xm)
        alxm = log(xm)
        g = xm*alxm - gammln(xm + 1.0)
     end if

     do
        y = tan(pi*ran1(idum))
        em = sq*y + xm
        if (em .ge. 0.) then

        em = int(em)
        t  = 0.9*(1.0 + y*y) * exp(em*alxm - gammln(em + 1.0) - g)
        endif
        
        if (ran1(idum) .le. t) exit
     enddo
  endif

  poidev = em

end function poidev


  real function bnldev(pp,n, idum)


  implicit none

  ! Arguments
  !real,intent(in) :: pp
  !integer, intent(in) :: n
  !integer, intent(inout) :: idum
  
  real :: pp
  integer :: n
  integer :: idum


  integer :: j
  integer, save :: nold = -1
  real, save :: pold = -1.0
  real, save :: pc, plog, pclog, en, oldg
  real :: am, em, g, p, sq, t, y
  real :: pi

  pi = acos(-1.0)

  if (pp .le. 0.5) then
     p = pp
  else
     p = 1. - pp
  end if

  am = n*p   

  if (n .lt. 25) then
     bnldev = 0.
     do j = 1, n
        if (ran1(idum) .lt. p)then
         bnldev = bnldev + 1.0
        endif
     end do
     
 
  else if (am .lt. 1.0) then
     g = exp(-am)
     t = 1.
     !do j = 0, n
     
      j = 0
      
     do while (t .ge. g .or. j .le. n)
        t = t*ran1(idum)
        j = j + 1
      enddo
      
        !if (t .lt. g) exit
        !bndevl = float(j)
     !end do
     j = n 
     bnldev = float(j)

  else
     if (n .ne. nold) then
        en= float(n)
        oldg = gammln(en +1.)
        nold = n
     end if
     
   if (n .ne. nold) then
        en = n
        oldg = gammln( en + 1)
        nold = n
    endif
    
    if (p .ne. pold) then
         pc = 1. -p
         plog = log(p)
         pclog = log(pc)
         pold = p
    endif
    
    sq = sqrt(2.*am*pc)
    !y = tan(pi*ran1(idum))
    em = sq*y + am
    if (em .lt. 0 .or. em .ge. en+1.) then
    
    y = tan(pi*ran1(idum))
    
    endif
    
    em = int(em) 
    t = 1.2*sq*(1. + y**2)*exp(oldg - gammln(em +1.)) - gammln(en - em+1) &
    + em *plog + (en- em)*pclog
    
    if (ran1(idum) .gt. t) then
      y = tan(pi*ran1(idum))
    endif
    
    bnldev = em
    
    
    endif
    
    if (p .ne. pp) then
    bnldev = n - bnldev
    endif

   endfunction bnldev

 end module rand
