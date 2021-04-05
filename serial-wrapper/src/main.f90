program sht_example

   use cheb_grids, only: initialize_grids, cleanup_grids
   use input_params ! namelist quantities
   use timing

   ! Rayleigh modules
#ifdef USE_SHTns
   use test_suite, only: test_SHTns, test_timing
   use Legendre_Transforms_SHTns, only: SHTns_Initialize, SHTns_Finalize
#else
   use test_suite, only: test_LT, test_timing
#endif
   use Legendre_Polynomials, only: Initialize_Legendre, Finalize_Legendre

   implicit none

   integer :: ntheta, n_phi, nm, m_max, i
   integer, allocatable :: mval(:)

   call initialize_timers() ! initialize timers
   call stopwatch(wall_time)%startclock()

   !=================================================================
   call stopwatch(initialization_time)%startclock()
   write(*,*)
   write(*,*) 'Begin initialization ...'
   write(*,*)
   call runtime_init() ! read namelists

   ! set some local variables
   ntheta = 3*(ell_max+1)/2
   nm = ell_max + 1
   allocate(mval(1:nm))
   do i=1,nm
      mval(i) = i-1
   enddo
   n_phi = 2*ntheta

   m_max = ell_max ! not stored in Legendre_Polynomials. ProblemSize, I think

   call initialize_grids(n_r) ! intialize my chebyshev grids

   ! initialize Legendre polynomials with parity
   call Initialize_Legendre(ntheta, ell_max, mval, .true.)

   ! intialize Legendre transforms
#ifdef USE_SHTns
   call SHTns_Initialize(n_threads, on_the_fly, verbose, eps_polar, &
                         n_phi, m_max, theta_contiguous)
#endif

   write(*,*)
   write(*,*) 'Initialization complete'
   call stopwatch(initialization_time)%increment()

   !=================================================================
   write(*,*)
   write(*,*) 'Setup/running tests ...'
   if (run_timing) then
      call test_timing(nloops, nfields, timing_file, max_walltime, min_walltime)
   else
#ifdef USE_SHTns
      call test_SHTns(ntest)
#else
      call test_LT(ntest)
#endif
   endif

   !=================================================================
   write(*,*)
   write(*,*) 'Cleanup ...'

   ! local stuff
   deallocate(mval)

   call cleanup_grids()
   call Finalize_Legendre()
#ifdef USE_SHTns
   call SHTns_Finalize()
#endif

   call finalize_timers()

   write(*,*)
   write(*,*) '--- Complete ---'
   write(*,*)

end program sht_example
