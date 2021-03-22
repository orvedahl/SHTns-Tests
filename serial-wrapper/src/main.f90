program sht_example

   use grids, only: initialize_grids, cleanup_grids
   use input_params ! namelist quantities

   ! Rayleigh modules
#ifdef USE_SHTns
   use test_suite, only: test_SHTns
   use Legendre_Transforms_SHTns, only: SHTns_Initialize, SHTns_Finalize
#else
   use test_suite, only: test_LT
#endif
   use Legendre_Polynomials, only: Initialize_Legendre, Finalize_Legendre

   implicit none

   integer :: ntheta, n_phi, nm, m_max, i
   integer, allocatable :: mval(:)

   !=================================================================
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

   call initialize_grids(n_r, ntheta) ! intialize my chebyshev/legendre grids

   ! initialize Legendre polynomials with parity
   call Initialize_Legendre(ntheta, ell_max, mval, .true.)

   ! intialize Legendre transforms
#ifdef USE_SHTns
   call SHTns_Initialize(n_threads, on_the_fly, verbose, eps_polar, n_phi, m_max)
#endif

   write(*,*)
   write(*,*) 'Initialization complete'

   !=================================================================
   write(*,*)
   write(*,*) 'Setup/running tests ...'
#ifdef USE_SHTns
   call test_SHTns(ntest)
#else
   call test_LT(ntest)
#endif

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

   write(*,*)
   write(*,*) '--- Complete ---'
   write(*,*)

end program sht_example
