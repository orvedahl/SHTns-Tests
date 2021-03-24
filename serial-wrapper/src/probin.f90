
module input_params

   implicit none

   private
   public :: runtime_init, runtime_print_params

   ! declarations/set default values
   character(len=256), save, public :: inputs_file_used = "not supplied"
   character(len=256), save, public :: job_name = "not supplied"
   integer, save, public :: n_r = 16
   integer, save, public :: ell_max = -1
   integer, save, public :: n_threads = 1
   real*8, save, public :: eps_polar = 1.d-10
   logical, save, public :: on_the_fly = .false.
   integer, save, public :: verbose = 2
   integer, save, public :: ntest = 1

   contains

   !====================================================================
   ! initialize the runtime parameters including commandline args/namelists
   !====================================================================
   subroutine runtime_init()

      integer :: narg, farg, ind, ind2
      character(len=256) :: fname, namelist_name
      logical :: found_inputs = .false.

      ! get number of args
      narg = command_argument_count()

      ! accepts command line arguments of the form:
      !   ./a.out <inputs-file> --tol 3.5 --title "My Title" --tstop 4.5
      !   or
      !   ./a.out --tol 3.5 --title "My Title" --tstop 4.5 <inputs-file>
      !
      ! args of the form: 
      !   ./a.out --eps
      ! are not yet supported
      !
      ! all arguments are read in as character variables, 
      ! it is our job to convert them to integer or real
      !
      ! -->if the inputs file is the first arg, it should not have the "--"
      ! -->if narg == 1, the arg is the input file 
      !      because --eps is not supported
      !      but it should still not have the "--"
      ! -->if it is the last arg, the second to last arg should not have 
      !      the "--" (b/c --eps is not supported)

      !---------------------------------------------------------------
      ! Search for Inputs file
      !---------------------------------------------------------------
      if (narg > 1) then

         ! see if inputs file is first
         farg = 1
         call get_command_argument(farg, value=fname)

         ! search for the "--" string
         ind = index(fname, "--", back = .false.)
         if (ind /= 1) then ! this is not a cmd line keyword

             namelist_name = fname
             call read_namelist(namelist_name)

             ! adjust farg so while loop (see below) skips the inputs file
             farg = 2
             found_inputs = .true.

         endif

         if (.not. found_inputs) then ! keep looking

            ! see if inputs file is last so try the third to
            ! last character and see if it has the "--"
            ! so: ./a.out --xcol 2 --tol 2.4 inputfile
            ! this means the third to last (--tol) should have "--"
            if (narg > 2) then ! must be more than 2 args for this version

               call get_command_argument(narg - 2, value=fname)

               ind2 = index(fname, "--", back = .false.)
               if (ind2 == 1) then ! this is a cmd line keyword

                  call get_command_argument(narg, value=fname)
                  namelist_name = fname
                  call read_namelist(namelist_name)

                  ! adjust narg so while loop does not include inputs file
                  narg = narg - 1
                  found_inputs = .true.

               endif

            endif

         endif 

      ! only 1 cmd line arg => it is assumed to be the namelist
      elseif (narg == 1) then

         call get_command_argument(narg, value=fname)
         namelist_name = fname
         call read_namelist(namelist_name)

         return

      ! narg < 1 => there are no cmd line args
      else

         return

      endif    

      !---------------------------------------------------------------
      ! parse the cmd line args
      !---------------------------------------------------------------
      do while (farg <= narg)

         call get_command_argument(farg, value=fname)

         select case (fname)

            case ('--n_r')
               farg = farg + 1
               call get_command_argument(farg, value = fname)

               read(fname, *) n_r

            case ('--ell_max')
               farg = farg + 1
               call get_command_argument(farg, value = fname)

               read(fname, *) ell_max

            case ('--n_threads')
               farg = farg + 1
               call get_command_argument(farg, value = fname)

               read(fname, *) n_threads

            case ('--eps_polar')
               farg = farg + 1
               call get_command_argument(farg, value = fname)

               read(fname, *) eps_polar

            case ('--on_the_fly')
               farg = farg + 1
               call get_command_argument(farg, value = fname)

               read(fname, *) on_the_fly

            case ('--verbose')
               farg = farg + 1
               call get_command_argument(farg, value = fname)

               read(fname, *) verbose

            case ('--ntest')
               farg = farg + 1
               call get_command_argument(farg, value = fname)

               read(fname, *) ntest

            case ('--job_name')
               farg = farg + 1
               call get_command_argument(farg, value = job_name)

            case ("--")
               farg = farg + 1
               exit

            case default
               write(*,*)
               write(*,*) 'UNKNOWN Option: '//trim(fname)
               write(*,*)
               write(*,*) 'Reminder:'
               write(*,*) '   Inputs file must be first or last arg'
               write(*,*)
               stop

         end select

         farg = farg + 1

      enddo

      return

   end subroutine runtime_init

   !====================================================================
   ! Read the namelist
   !====================================================================
   subroutine read_namelist(namelist_file)
   
      ! INPUT
      character(len=256), intent(in) :: namelist_file

      ! LOCAL
      integer :: io=34, ierr, ierr2
   
      namelist /input/ n_r
      namelist /input/ ell_max
      namelist /input/ n_threads
      namelist /input/ eps_polar
      namelist /input/ on_the_fly
      namelist /input/ verbose
      namelist /input/ ntest

      inputs_file_used = namelist_file

      ! open namelist
      open(unit=io, file=trim(namelist_file), status='old', action='read', &
           delim='apostrophe', iostat=ierr)

      ! open successful
      if (ierr == 0) then

         read(unit=io, nml=input, iostat=ierr2)
         close(unit=io)

         ! error reading namelist
         if (ierr2 /= 0) then
            write(*,*)
            write(*,*) "ERROR: Could not Read namelist: "//trim(namelist_file)
            write(*,*)
            stop 6
         endif

      ! open unsuccessful
      else

         write(*,*)
         write(*,*) "ERROR: Could not Open namelist: "//trim(namelist_file)
         write(*,*)
         close(unit=io)
         stop 6

      endif

      return
   
   end subroutine read_namelist

   !====================================================================
   ! pretty printing routine
   !====================================================================
   subroutine runtime_print_params(io)

      integer, intent(in) :: io

      100 format (1x, a32, 1x, "=", 1x, a)
      101 format (1x, a32, 1x, "=", 1x, i10)
      102 format (1x, a32, 1x, "=", 1x, g20.10)
      103 format (1x, a32, 1x, "=", 1x, l)

      write (io,101) "n_r", n_r
      write (io,101) "ell_max", ell_max
      write (io,101) "n_threads", n_threads
      write (io,102) "eps_polar", eps_polar
      write (io,103) "on_the_fly", on_the_fly
      write (io,101) "verbose", verbose
      write (io,101) "ntest", ntest
      write (io,100) "job_name", trim(job_name)

   end subroutine runtime_print_params

end module input_params
