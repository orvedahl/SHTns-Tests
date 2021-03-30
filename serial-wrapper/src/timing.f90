
module Timing

   implicit none

   ! a single timer
   type timer
      real*8 :: time    ! current time, i.e., "start" the stopwatch
      real*8 :: delta   ! elapsed time from previous "start", i.e., "stop" the stopwatch and record it
      real*8 :: elapsed ! total elapsed time, including multiple start/stops

      contains

      procedure :: init => initialize_timer
      procedure :: startclock ! "start" the clock
      procedure :: stopclock  ! "stop" the clock and record the elapsed time
      procedure :: increment  ! add the elapsed time to the total elapsed time
   end type timer

   ! specify number of timers and indices for each one
   integer, parameter, public :: ntimers = 5
   integer, parameter, public :: wall_time = 1
   integer, parameter, public :: initialization_time = 2
   integer, parameter, public :: loop_time = 3
   integer, parameter, public :: to_physical = 4
   integer, parameter, public :: to_spectral = 5
   type(timer), allocatable :: stopwatch(:) ! a collection of timers

   private
   public :: timer, stopwatch, initialize_timers, finalize_timers

   contains

   ! generic routine to get current time in seconds
   subroutine get_current_time(t)
      real*8, intent(inout) :: t
      integer :: rate, cnts

      call system_clock(cnts, rate)
      t = real(cnts)/real(rate)
   end subroutine get_current_time

   ! allocate a collection of timers
   subroutine initialize_timers()
      integer :: i
      allocate(stopwatch(1:ntimers))
      do i=1,ntimers
         call stopwatch(i)%init()
      enddo
   end subroutine initialize_timers

   subroutine finalize_timers()
      deallocate(stopwatch)
   endsubroutine finalize_timers

   ! initialize components to zero
   subroutine initialize_timer(self)
      class(timer) :: self
      self%elapsed = 0.0d0
      self%time = 0.0d0
      self%delta = 0.0d0
   end subroutine initialize_timer

   ! call current time and store it
   subroutine startclock(self)
      class(timer) :: self
      call get_current_time(self%time)
   end subroutine startclock

   ! call current time and store the elapsed time from previous call to startclock
   subroutine stopclock(self)
      class(timer) :: self
      real*8 :: t
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

