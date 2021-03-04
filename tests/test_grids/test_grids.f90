program test_grids
  use data_types
  use grids

  implicit none

  integer :: io = 17

  call initialize_grids(N_radius, N_theta)

  ! write data to a formatted file for reading/comparison to python


  call cleanup_grids()

end program test_grids

