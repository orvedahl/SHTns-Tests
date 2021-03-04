module data_types

  implicit none

  ! reals
  integer, parameter, public :: dpt = selected_real_kind(15,307)
  integer, parameter, public :: spt = selected_real_kind(6,37)

  ! integers
  integer, parameter, public :: i8_t = selected_int_kind(15)
  integer, parameter, public :: i4_t = selected_int_kind(9)

end module data_types
