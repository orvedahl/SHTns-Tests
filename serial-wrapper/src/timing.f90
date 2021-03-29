
module Timing

   use data_types

   implicit none

   ! a single timer
   type timer
      real(kind=dpt) :: time    ! current time, i.e., "start" the stopwatch
      real(kind=dpt) :: delta   ! elapsed time from previous "start", i.e., "stop" the stopwatch and record it
      real(kind=dpt) :: elapsed ! total elapsed time, including multiple start/stops

      contains

      procedure :: init => initialize_timer
      procedure :: startclock ! "start" the clock
      procedure :: stopclock  ! "stop" the clock and record the elapsed time
      procedure :: increment  ! add the elapsed time to the total elapsed time
   end type timer

   ! specify number of timers and indices for each one
   integer, parameter :: ntimes = 3
   integer, parameter :: initialization_time = 1
   integer, parameter :: loop_time = 2
   integer, parameter :: output_time = 3
   type(timer), allocatable :: stopwatch(:) ! a collection of times

   private
   public :: timer

   contains

   ! generic routine to get current time in seconds
   subroutine get_current_time(t)
      real(kind=dpt), intent(inout) :: t
      integer(kind=i8_t) :: rate, cnts

      call system_clock(cnts, rate)
      t = real(cnts, kind=dpt)/real(rate, kind=dpt)
   end subroutine get_current_time

   ! allocate a collection of timers
   subroutine initialize_timers()
      integer :: i
      allocate(stopwatch(1:ntimers))
      do i=1,ntimers
         stopwatch(i)%init()
      enddo
   end subroutine initialize_timers

   subroutine finalize_timers()
      deallocate(stopwatch)
   endsubroutine finalize_timers

   ! initialize components to zero
   subroutine initialize_timer(self)
      class(timer) :: self
      self%elapsed = 0.0_dpt
      self%time = 0.0_dpt
      self%delta = 0.0_dpt
   end subroutine initialize_timer

   ! call current time and store it
   subroutine startclock(self)
      class(timer) :: self
      call get_current_time(self%time)
   end subroutine startclock

   ! call current time and store the elapsed time from previous call to startclock
   subroutine stopclock(self)
      class(timer) :: self
      real(kind=dpt) :: t
      call get_current_time(t)
      self%delta = t - self%time
   end subroutine stopclock

   ! add elapsed time to cumulative record
   subroutine increment(self)
      class(timer) :: self
      call self%stopclock()
      self%elapsed = self%elapsed + self%delta
   end subroutine increment

end module Timing

