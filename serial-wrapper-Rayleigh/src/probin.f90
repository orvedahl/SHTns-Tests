
module input_params

   implicit none

   public :: runtime_init

   ! declarations/set default values
   integer, save, public :: n_r = 16
   integer, save, public :: ell_max = -1
   integer, save, public :: n_threads = 1
   real*8, save, public :: eps_polar = 1.d-10
   logical, save, public :: on_the_fly = .false.
   integer, save, public :: verbose = 2
   integer, save, public :: ntest = 1

   contains

   !====================================================================
   ! Read the namelist
   !====================================================================
   subroutine read_namelist(namelist_file)
   
      character(len=256), intent(in) :: namelist_file
      integer :: io=34, ierr, ierr2
   
      namelist /input/ n_r
      namelist /input/ ell_max
      namelist /input/ n_threads
      namelist /input/ eps_polar
      namelist /input/ on_the_fly
      namelist /input/ verbose
      namelist /input/ ntest

      open(unit=io, file=trim(namelist_file), status='old', action='read', &
           delim='apostrophe', iostat=ierr)

      if (ierr == 0) then ! open successful
         read(unit=io, nml=input, iostat=ierr2)
         close(unit=io)
         if (ierr2 /= 0) then ! error reading namelist
            write(*,*)
            write(*,*) "ERROR: Could not Read namelist: "//trim(namelist_file)
            write(*,*)
            stop 6
         endif
      else ! open unsuccessful
         write(*,*)
         write(*,*) "ERROR: Could not Open namelist: "//trim(namelist_file)
         write(*,*)
         close(unit=io)
         stop 6
      endif

   end subroutine read_namelist

   !====================================================================
   ! initialize the runtime parameters including commandline args/namelists
   !====================================================================
   subroutine runtime_init()

      integer :: narg, farg, ind, ind2
      character(len=256) :: fname, namelist_name
      logical :: found_inputs = .false.

      ! read namelist
      namelist_name = "input"
      call read_namelist(namelist_name)

      ! parse command line

   end subroutine runtime_init

end module input_params
