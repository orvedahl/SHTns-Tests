!
! very basic, verification type tests
!
program test

  implicit none

  integer :: i
  real*8 :: xreal(1:10), ximag(1:10)
  complex*16 :: ar, ai, values(1:10)

  ai = (0.0d0, 1.0d0) ! imaginary unit
  ar = (1.0d0, 0.0d0) ! regular unit

  ! set real/imag arrays
  do i=1,10
    xreal(i) = 1.0d0*i
    ximag(i) = 1.0d0*i*i
  enddo

  ! try initializing complex arrays
  values(:) = ar*xreal(:) + ai*ximag(:)

  do i=1,10
    write(*,*) xreal(i), ximag(i), values(i)
  enddo

end program test
