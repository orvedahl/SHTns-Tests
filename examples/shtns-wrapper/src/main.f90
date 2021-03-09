program sht_example

   use data_types
   use variables
   use grids, only: initialize_grids, cleanup_grids

   use input_params, only: eps_polar, n_threads, on_the_fly, verbose
   use problemsize, only: l_max, n_theta, init_problemsize
   use legendre_polynomials, only: initialize_legendre, finalize_legendre
   use legendre_transforms

   implicit none

   integer :: lm, m
   real(kind=dpt) :: err

   ! physical space is complex, because it is fixed m, i.e., its been FFT-ed
   complex(kind=dpt),allocatable :: spectral1(:), spectral2(:)
   complex(kind=dpt),allocatable :: physical_m_fixed_1(:), physical_m_fixed_2(:)

   !=================================================================
   write(*,*)
   write(*,*) 'Begin initialization ...'
   write(*,*)
   ! read namelists and initialize Legendre polynomials
   call init_problemsize(eps_polar, n_threads, on_the_fly, verbose)

   ! this is new ...
   call initialize_legendre_transforms(eps_polar, n_threads, on_the_fly, verbose)

   ! my own intialization routines, independent of the above
   call initialize_grids(32, n_theta)
  
   write(*,*)
   write(*,*) 'Initialization complete'

   write(*,*) 'num lm modes',n_lm

   !=================================================================
   write(*,*)
   write(*,*) 'Setup tests ...'
   if (l_max .lt. 12) then
      write(*,*)
      write(*,*) 'ERROR: tests require l_max >= 12, l_max = ',l_max
      write(*,*)
      stop
   endif
   allocate(physical_m_fixed_1(1:n_theta), spectral1(1:n_lm))
   allocate(physical_m_fixed_2(1:n_theta), spectral2(1:n_lm))
   spectral1(:) = (zero,zero)
   spectral2(:) = (zero,zero)
   physical_m_fixed_1(:) = (zero,zero)
   physical_m_fixed_2(:) = (zero,zero)

   write(*,*)
   write(*,*) 'Test 1: f(l,m) = (2+3j)*Y_2^m + (1-4j)*Y_1^m + (-6+1j)*Y_12^m'
   m = 7
   lm = lm_index(2,m)
   spectral1(lm) = (2.0_dpt,3.0_dpt)
   lm = lm_index(1,m)
   spectral1(lm) = (1.0_dpt,-4.0_dpt)
   lm = lm_index(12,m)
   spectral1(lm) = (-6.0_dpt,1.0_dpt)

   call LT_ToPhysical_single(spectral1, physical_m_fixed_1, m) ! to Physical

   call LT_ToSpectral_single(physical_m_fixed_1, spectral2, m) ! to Spectral

   err = maxval(abs(spectral1 - spectral2))
   write(*,*)
   write(*,*) '   m = ',m
   write(*,*) '   max error in IC - LTinv(LT(IC)) = ',err
   write(*,*)

   !=================================================================
   write(*,*)
   write(*,*) 'Cleanup ...'
   deallocate(physical_m_fixed_1, spectral1)
   deallocate(physical_m_fixed_2, spectral2)

   call cleanup_grids() ! my routines
   call finalize_legendre()
   call finalize_legendre_transforms()

   write(*,*)
   write(*,*) '--- Complete ---'
   write(*,*)

end program sht_example
