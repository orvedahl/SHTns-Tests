!
! Module to hold variables
!
module variables

   use data_types

   implicit none

   ! constants
   real(kind=dpt), parameter :: zero    = 0.0_dpt
   real(kind=dpt), parameter :: small   = 1.0e-12_dpt
   real(kind=dpt), parameter :: half    = 0.5_dpt
   real(kind=dpt), parameter :: sixth   = 1.0_dpt/6.0_dpt
   real(kind=dpt), parameter :: third   = 1.0_dpt/3.0_dpt
   real(kind=dpt), parameter :: twothrd = 2.0_dpt/3.0_dpt
   real(kind=dpt), parameter :: two     = 2.0_dpt
   real(kind=dpt), parameter :: one     = 1.0_dpt
   real(kind=dpt), parameter :: pi      = acos(-1.0_dpt)
   real(kind=dpt), parameter :: twopi   = two*pi

end module variables

