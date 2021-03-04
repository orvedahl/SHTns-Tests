program test_grids
  use data_types
  use grids
  use input_params
  use spherical_harmonics

  implicit none

  integer :: io = 17, N_theta
  integer :: ell, m, i, j
  real(kind=dpt) :: ylm_r, ylm_i

  call runtime_init()
  call write_job_info(".")

  N_theta = 3*(ell_max + 1)/2
  write(*,*) "Initializing grid ..."
  call initialize_grids(N_radius, N_theta)

  !-----------------------------------------------------------------------

  open(unit=io, file="test_ylm.txt", form='formatted')
  write(io,*) '#','ell',"  ","m","  ","theta","  ","phi","  ","Re(Ylm)","  ","Im(Ylm)"

  do ell=0,4
     do m=0,ell
        do i=1,Nth,7
           do j=1,Nphi,15
              call Ylm(ell,m,theta(i), phi(j), ylm_r, ylm_i)
              write(io,*) ell,"  ",m,"  ",theta(i),"  ",phi(j),"  ",ylm_r, "  ", ylm_i
           enddo
        enddo
     enddo
  enddo

  close(io)

  write(*,*) "clean up ..."
  call cleanup_grids()

end program test_grids

