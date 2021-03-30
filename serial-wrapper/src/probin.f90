
module input_params

   implicit none

   private
   public :: runtime_init

   ! To set a new namelist element:
   !    1) declare the variable and give it a default value
   !    2) add the variable to the "namelist" line
   !    3) optional: add variable to parse_command_line_options subroutine
   !
   integer, save, public :: n_r = 16
   integer, save, public :: ell_max = -1

   integer, save, public :: n_threads = 1

   real*8, save, public :: eps_polar = 1.d-10
   logical, save, public :: on_the_fly = .false.
   integer, save, public :: verbose = 2

   integer, save, public :: ntest = 1

   integer, save, public :: nfields = 100
   logical, save, public :: run_timing = .false.
   integer, save, public :: nloops = 100
   character(len=1024), save, public :: timing_file = "timing.out"

   namelist /input/ n_r, ell_max, n_threads, &
                    eps_polar, on_the_fly, verbose, &
                    ntest, &
                    nfields, run_timing, nloops, timing_file

   interface Read_CMD_Line
      module procedure Read_CMD_Integer, Read_CMD_Double
      module procedure Read_CMD_Logical, Read_CMD_String
   end interface Read_CMD_Line

   contains

   !====================================================================
   ! Parse the command line for suitable values
   !====================================================================
   subroutine parse_command_line_options()
      call Read_CMD_Line("--nr", n_r)
      call Read_CMD_Line("--lmax", ell_max)
      call Read_CMD_Line("--n-threads", n_threads)
      call Read_CMD_Line("--eps-polar", eps_polar)
      call Read_CMD_Line("--on-the-fly", on_the_fly)
      call Read_CMD_Line("--verbose", verbose)
      call Read_CMD_Line("--ntest", ntest)
      call Read_CMD_Line("--run-timing", run_timing)
      call Read_CMD_Line("--nloops", nloops)
      call Read_CMD_Line("--nfields", nfields)
      call Read_CMD_Line("--output", timing_file)
   end subroutine parse_command_line_options

   !====================================================================
   ! Read the namelist
   !====================================================================
   subroutine read_namelist(namelist_file)
   
      character(len=256), intent(in) :: namelist_file
      integer :: io=34, ierr, ierr2
   
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

      character(len=256) :: namelist_filename

      namelist_filename = "input"

      call read_namelist(namelist_filename)

      call parse_command_line_options()

   end subroutine runtime_init

   !====================================================================
   ! Command line interface routines all named Read_CMD_<type>
   ! arguments are assumed to be of the format:
   !     ./a.out vname1 value1 vname2 value2 vname3 1 vname4 0 ...
   ! every variable must be assigned a value, even logicals (which use 0,1)
   !
   ! The calling sequence is the same for each routine:
   !
   !     Read_CMD_Line(search_string, variable)
   !
   !     where search_string is what will appear on the command line
   !     and variable is the actual variable that will take on the
   !     given value. for example:
   !          ...
   !          real(kind=dpt) :: tolerance
   !          integer :: magnetic_init_type
   !          ...
   !          call Read_CMD_Line("--tol", tolerance)
   !          call Read_CMD_Line("--mag-init", magnetic_init_type)
   !          ...
   !====================================================================
   subroutine Read_CMD_String(search_string, variable)
      character(*), intent(in) :: search_string
      character(*), intent(inout) :: variable
      integer :: i,n
      character(len=1024) :: argname, argval

      n = command_argument_count()

      do i=1,n,2 ! loop over all arguments, count by 2

         call get_command_argument(i, argname)  ! name of specified variable
         call get_command_argument(i+1, argval) ! value of variable

         if (trim(search_string) .eq. trim(argname)) then
             variable = trim(adjustl(argval)) ! remove space infront and trailing
         endif

      enddo
   end subroutine Read_CMD_String

   subroutine Read_CMD_Logical(search_string, variable)
      character(*), intent(in) :: search_string
      logical, intent(inout) :: variable
      integer :: i,n,itemp
      character(len=1024) :: argname, argval, argshift
      n = command_argument_count()

      do i=1,n,2 ! loop over all arguments, count by 2

         call get_command_argument(i, argname)  ! name of specified variable
         call get_command_argument(i+1, argval) ! value of variable

         if (search_string .eq. argname) then
             argshift = trim(adjustl(argval)) ! remove space infront and trailing

             read(argshift,*) itemp ! read from argshift, store value in itemp

             if (itemp .eq. 1) then
                 variable = .true.
             else
                 variable = .false.
             endif
         endif

      enddo
   end subroutine Read_CMD_Logical

   subroutine Read_CMD_Integer(search_string, variable)
      character(*), intent(in) :: search_string
      integer, intent(inout) :: variable
      integer :: i,n
      character(len=1024) :: argname, argval, argshift
      n = command_argument_count()

      do i=1,n,2 ! loop over all arguments, count by 2

         call get_command_argument(i, argname)  ! name of specified variable
         call get_command_argument(i+1, argval) ! value of variable

         if (search_string .eq. argname) then
             argshift = trim(adjustl(argval)) ! remove space infront and trailing

             read(argshift,*) variable ! read from argshift, store value in variable
         endif

      enddo
   end subroutine Read_CMD_Integer

   subroutine Read_CMD_Double(search_string, variable)
      character(*), intent(in) :: search_string
      real*8, intent(inout) :: variable
      integer :: i,n
      character(len=1024) :: argname, argval, argshift
      n = command_argument_count()

      do i=1,n,2 ! loop over all arguments, count by 2

         call get_command_argument(i, argname)  ! name of specified variable
         call get_command_argument(i+1, argval) ! value of variable

         if (search_string .eq. argname) then
             argshift = trim(adjustl(argval)) ! remove space infront and trailing

             read(argshift,*) variable ! read from argshift, store value in variable
         endif

      enddo
   end subroutine Read_CMD_Double

end module input_params
