program test_grids
  use data_types
  use grids
  use input_params
  use job_info

  implicit none

  integer :: io = 17, i, N_theta

  call runtime_init()
  call write_job_info(".")

  N_theta = 3*(ell_max + 1)/2
  write(*,*) "Initializing grid ..."
  call initialize_grids(N_radius, N_theta)

  ! write data to a formatted file for reading/comparison to python
  write(*,*) "writing chebyshev grid file ..."
  open(unit=io, file="test_grid_chebyshev.txt", form='formatted')
  write(io,*) '#', 'index', "  ", "x(i)", "  ", "theta(i)"
  do i=1,Nr
     write(io,*) i, "  ", x_cheb(i), "  ", theta_cheb(i)
  enddo
  close(io)

  write(*,*) "writing legendre grid file ..."
  open(unit=io, file="test_grid_legendre.txt", form='formatted')
  write(io,*) '#', 'index', "  ", "x(i)", "  ", "theta(i)", "  ", "weights(i)"
  do i=1,Nth
     write(io,*) i, "  ", x_leg(i), "  ", theta(i), "  ", weights_leg(i)
  enddo
  close(io)

  write(*,*) "writing Plm file ..."
  open(unit=io, file="test_grid_Plm.txt", form='formatted')
  write(io,*) '#', 'index', "  ", "Plm(i,10,5)", "  ", &
                    "Plm(i,5,10)", "  ", "Plm(i,1,0)", "  ", "Plm(i,1,1)"
  do i=1,Nth
     write(io,*) i, "  ", Plm(i,10,5), "  ", Plm(i,5,10), "  ", Plm(i,1,0), "  ", Plm(i,1,1)
  enddo
  close(io)

  write(*,*) "writing fourier grid file ..."
  open(unit=io, file="test_grid_fourier.txt", form='formatted')
  write(io,*) '#', 'index', "  ", "phi(i)"
  do i=1,Nphi
     write(io,*) i, "  ", phi(i)
  enddo
  close(io)

  write(*,*) "clean up ..."
  call cleanup_grids()

end program test_grids

