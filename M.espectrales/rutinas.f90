 subroutine s3p_sl ( n, a1, a2, a3, b, x, job, work2, work3, work4 )
!
!*******************************************************************************
!
!! S3P_SL solves a tridiagonal periodic system factored by S3P_FA.
!
!
!  Modified:
!
!    03 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real A1(N), A2(N), A3(N), factor data from S3P_FA.
!
!    Input, real B(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution to the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
!    Input, real WORK2(N-1), WORK3(N-1), WORK4, factor data from S3P_FA.
!
  integer n
!
  real*8 a1(n)
  real*8 a2(n)
  real*8 a3(n)
  real*8 b(n)
  integer i
  integer ierror
  integer job
  real*8 work2(n-1)
  real*8 work3(n-1)
  real*8 work4
  real*8 x(n)
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  x(1:n) = b(1:n)

  if ( job == 0 ) then
!
!  Solve A1 * X1 = B1.
!
    call s3_np_sl ( n-1, a1(2), a2, a3, x, job )
!
!  X2 = B2 - A3 * X1
!
    x(n) = x(n) - a3(n) * x(1) - a1(n) * x(n-1)
!
!  Solve A4 * X2 = X2
!
    x(n) = x(n) / work4
!
!  X1 := X1 - inverse ( A1 ) * A2 * X2.
!
    x(1:n-1) = x(1:n-1) - work2(1:n-1) * x(n)

  else
!
!  Solve transpose ( A1 ) * X1 = B1.
!
    call s3_np_sl ( n-1, a1(2), a2, a3, x, job )
!
!  X2 := X2 - transpose ( A2 ) * B1
!
    x(n) = x(n) - a1(1) * x(1) - a3(n-1) * x(n-1)
!
!  Solve A4 * X2 = X2.
!
    x(n) = x(n) / work4
!
!  X1 := X1 - transpose ( inverse ( A1 ) * A3 ) * X2.
!
    x(1:n-1) = x(1:n-1) - work3(1:n-1) * x(n)

  end if

  return
end



!-------------------------------------------------------------------------------



  subroutine s3_np_sl ( n, a1, a2, a3, b, job )
!
!*******************************************************************************
!
!! S3_NP_SL solves a tridiagonal system factored by S3_NP_FA.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1), the factor information
!    returned by S3_NP_FA.
!
!    Input/output, real B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer n
!
  real*8 a1(2:n)
  real*8 a2(1:n)
  real*8 a3(1:n-1)
  real*8 b(n)
  integer i
  integer ierror
  integer job
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    do i = 2, n
      b(i) = b(i) - a1(i) * b(i-1)
    end do
!
!  Solve U * X = Y.
!
    do i = n, 1, -1
      b(i) = b(i) / a2(i)
      if ( i > 1 ) then
        b(i-1) = b(i-1) - a3(i-1) * b(i)
      end if
    end do

  else
!
!  Solve tranpose ( U ) * Y = B
!
    do i = 1, n
      b(i) = b(i) / a2(i)
      if ( i < n ) then
        b(i+1) = b(i+1) - a3(i) * b(i)
      end if
    end do
!
!  Solve transpose ( L ) * X = Y.
!
    do i = n-1, 1, -1
      b(i) = b(i) - a1(i+1) * b(i+1)
    end do

  end if

  return
end



!-------------------------------------------------------------------------------



  subroutine s3p_check ( n, ierror )
!
!*******************************************************************************
!
!! S3P_CHECK checks the dimensions of a tridiagonal periodic matrix.
!
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Output, integer IERROR, error flag.
!    0, the dimensions are legal.
!    1, N is less than 3.
!
  integer ierror
  integer n
!
  ierror = 0

  if ( n < 3 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'S3P_CHECK - Fatal error!'
    write ( *, * ) '  N must be at least 3.'
    write ( *, * ) '  The input value is N = ', n
  end if

  return
end



!-------------------------------------------------------------------------------


  
  subroutine s3_check ( n, ierror )
!
!*******************************************************************************
!
!! S3_CHECK checks the dimensions of a real tridiagonal matrix.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Output, integer IERROR, error flag.
!    0, no errors detected.
!    1, N was less than 2.
!
  integer ierror
  integer n
!
  ierror = 0

  if ( n < 2 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'S3_CHECK - Fatal error!'
    write ( *, * ) '  N must be at least 2.'
    write ( *, * ) '  The input N was ', n
  end if

  return
end




!-------------------------------------------------------------------------------

  
  
  
  subroutine s3p_fa ( n, a1, a2, a3, info, work2, work3, work4 )
!
!*******************************************************************************
!
!! S3P_FA factors a tridiagonal periodic matrix.
!
!
!  Discussion:
!
!    Once the matrix has been factored by S3P_FA, S3P_SL may be called
!    to solve linear systems involving the matrix.
!
!    The logical matrix has a form which is suggested by this diagram:
!
!      D1 U1          L1
!      L2 D2 U2
!         L3 D3 U3
!            L4 D4 U4
!               L5 D5 U5
!      U6          L6 D6
!
!    The algorithm treats the matrix as a border banded matrix:
!
!      ( A1  A2 )
!      ( A3  A4 )
!
!    where:
!
!      D1 U1          | L1
!      L2 D2 U2       |  0
!         L3 D3 U3    |  0
!            L4 D4 U4 |  0
!               L5 D5 | U5
!      ---------------+---
!      U6  0  0  0 L6 | D6
!
!  Method:
!
!    The algorithm rewrites the system as:
!
!         X1 + inverse(A1) A2 X2 = inverse(A1) B1
!
!      A3 X1 +             A4 X2 = B2
!
!    The first equation can be "solved" for X1 in terms of X2:
!
!         X1 = - inverse(A1) A2 X2 + inverse(A1) B1
!
!    allowing us to rewrite the second equation for X2 explicitly:
!
!      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
!
!  Modified:
!
!    03 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input/output, real A1(N), A2(N), A3(N).
!    On input, these arrays contain the subdiagonal, diagonal, and
!    superdiagonal entries of the coefficient matrix.  The special
!    cases are that A1(1) is the coefficient of X(N), and A3(N)
!    is the coefficient of X(1).
!
!    On output, the arrays have been modified to hold information
!    defining the border-banded factorization of submatrices A1
!    and A3.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
!    Output, real WORK2(N-1), WORK3(N-1), WORK4, factorization information.
!
  integer n
!
  real*8 a1(n)
  real*8 a2(n)
  real*8 a3(n)
  integer i
  integer ierror
  integer info
  integer job
  real*8 work2(n-1)
  real*8 work3(n-1)
  real*8 work4
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Compute inverse(A1):
!
  call s3_np_fa ( n-1, a1(2), a2, a3, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_FA - Fatal error!'
    write ( *, * ) '  S3_NP_FA returned INFO = ', info
    write ( *, * ) '  Factoring failed for column INFO.'
    write ( *, * ) '  The tridiagonal matrix A1 is singular.'
    write ( *, * ) '  This algorithm cannot continue!'
    write ( *, * ) ' '
    return
  end if
!
!  WORK2 := inverse(A1) * A2.
!
  work2(1) = a1(1)
  work2(2:n-2) = 0.0E+00
  work2(n-1) = a3(n-1)

  job = 0
  call s3_np_sl ( n-1, a1(2), a2, a3, work2, job )
!
!  WORK3 := inverse ( transpose ( A1 ) ) * tranpose ( A3 ).
!
  work3(1) = a3(n)
  work3(2:n-2) = 0.0E+00
  work3(n-1) = a1(n)

  job = 1
  call s3_np_sl ( n-1, a1(2), a2, a3, work3, job )
!
!  A4 := ( A4 - A3 * inverse(A1) * A2 )
!
  work4 = a2(n) - a3(n) * work2(1) - a1(n) * work2(n-1)

  if ( work4 == 0.0E+00 ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'S3P_FA - Fatal error!'
    write ( *, * ) '  The factored A4 submatrix is zero.'
    write ( *, * ) '  This algorithm cannot continue!'
    write ( *, * ) ' '
    return
  end if

  return
end
  




!-------------------------------------------------------------------------------

  
  
 
 
 subroutine s3_np_fa ( n, a1, a2, a3, info )
!
!*******************************************************************************
!
!! S3_NP_FA factors a tridiagonal system without pivoting.
!
!
!  Discussion:
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!    S3_NP_FA and S3_NP_SL may be preferable to the corresponding
!    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
!    in one step, and does not save the factorization.
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output real A1(2:N), A2(1:N), A3(1:N-1), the subdiagonal,
!    diagonal, and superdiagonal of the matrix.  On output, these are
!    overwritten by factorization information.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer n
!
  real*8 a1(2:n)
  real*8 a2(1:n)
  real*8 a3(1:n-1)
  integer i
  integer ierror
  integer info
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  info = 0

  do i = 1, n-1

    if ( a2(i) == 0.0E+00 ) then
      info = i
      write ( *, * ) ' '
      write ( *, * ) 'S3_NP_FA - Fatal error!'
      write ( *, * ) '  Zero pivot on step ', info
      return
    end if

    a1(i+1) = a1(i+1) / a2(i)
    a2(i+1) = a2(i+1) - a1(i+1) * a3(i)

  end do

  if ( a2(n) == 0.0E+00 ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_FA - Fatal error!'
    write ( *, * ) '  Zero pivot on step ', info
    return
  end if

  return
end 

!subroutine s3_check ( n, ierror )
!
!*******************************************************************************
!
!! S3_CHECK checks the dimensions of a real tridiagonal matrix.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Output, integer IERROR, error flag.
!    0, no errors detected.
!    1, N was less than 2.
!
!  implicit none
!
!  integer ierror
!  integer n
!
! ierror = 0

!  if ( n < 2 ) then
!    ierror = ierror + 1
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'S3_CHECK - Fatal error!'
!    write ( *, '(a)' ) '  N must be at least 2.'
!    write ( *, '(a,i6)' ) '  The input N was ', n
!  end if
!
!  return
!end

