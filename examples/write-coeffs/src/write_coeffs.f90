module utils
   use iso_fortran_env, only: real64

   implicit none

   integer, parameter :: dp=real64
   real(dp), parameter :: pi=acos(-1.0_dp)
  
   contains

    subroutine Y_32(costh, Yout)
        real(dp), intent(in) :: costh(:)
        complex(dp),intent(inout) :: Yout(:)

        real(dp) :: A

        A = sqrt(7.0_dp / (4.0_dp*pi*120.0_dp))
        Yout(:) = 15.0_dp*A*costh(:)*(1.0_dp-costh(:)*costh(:))
    end subroutine Y_32

    subroutine Y_2m(costh, Yout, m)
        integer, intent(in) :: m
        real(dp), intent(in) :: costh(:)
        complex(dp),intent(inout) :: Yout(:)

        real(dp) :: A

        if (m .eq. 0) then
           A = sqrt(5.0_dp / (4.0_dp*pi))
           Yout(:) = A*0.5_dp*(3.0_dp*costh(:)*costh(:)-1.0_dp)

        elseif (m .eq. 1) Then
           A = sqrt(5.0_dp / (4.0_dp*pi*6.0_dp))
           Yout(:) = -3.0_dp*A*costh(:)*sqrt(1.0_dp-costh(:)*costh(:))

        elseif (m .eq. 2) Then
           A = sqrt(5.0_dp / (4.0_dp*pi*24.0_dp))
           Yout(:) = 3.0_dp*A*(1.0_dp-costh(:)*costh(:))

        else
           Yout(:) = -1.0_dp
        endif
    end subroutine Y_2m

end module utils
  
program write_coeffs
  
    use iso_c_binding

    use utils

    implicit none
  
    include 'shtns.f03'
  
    type(shtns_info), pointer :: shtns
    real(dp), pointer :: cosTheta(:), sinTheta(:)
    type(c_ptr) :: shtns_c
  
    integer :: lmax, mmax, mres, nthreads
    integer :: nlat, nphi, lm, l, m, lmStart, lmStop
    integer :: norm, layout, i
    real(dp) :: eps_polar, maxerr, diff, temp
    logical :: do_real
    complex(dp), allocatable :: spectral(:), physical(:), spec_true(:), phys_true(:)
  
    ! set size of transforms, lmax is set in the input file
    lmax = 31
    mmax = lmax              ! max azimuthal wave number
    mres = 1                 ! azimuthal periodicity is 2*pi/mres
    nlat = int(3*(lmax+1)/2) ! number of points in theta
    nphi = 2*nlat            ! number of points in phi
  
    ! set polar optimisation threshold
    eps_polar = 1.0e-10_dp
  
    !-- Normalization and layout for SHTs
    norm = sht_orthonormal
    !layout = sht_gauss + sht_phi_contiguous ! use Legendre grid/quadrature
    layout = sht_gauss_fly + sht_phi_contiguous ! use on-the-fly Legendre grid/quadrature

    layout = layout + sht_south_pole_first + sht_scalar_only
  
    ! set output level: 0=none, 1=some, 2=debug
    call shtns_verbose(2)

    ! return number of OpenMP threads that will be used, try to use passed value
    !   send_threads > 0  : set max number of threads to use
    !   send_threads <= 0 : max number of threads is set to number of processors
    !   send_threads = 1  : disable OpenMP
    nthreads = shtns_use_threads(1)
  
    ! initialize/allocate transforms and build useful arrays
    !   lmax : max SH degree
    !   mmax : max azimuthal order
    !   mres : 2*pi/mres is the azimuthal periodicity, mmax*mres is max azimuthal order
    !   norm : describe normalization, Condon-Shortley is default
    ! returns a pointer (basically?) to a shtns_info structure with attributes (not all):
    !   nlm - number of (l,m) entries
    !   lmax, mmax, mres, nlat, nphi
    !   li - size(nlm) degree l for the given mode index
    !   mi - size(nlm) order m for the given mode index
    !   ct - cos(theta) array of size nlat
    !   st - sin(theta) array of size nlat
    shtns_c = shtns_create(lmax, mmax, mres, norm)

    ! attach a grid to the created SHT
    ! must occur *after* transforms are initialized and *before* any SHT call
    !   eps - polar optimization, threshold below which Plm is treated as 0 for high-m
    !           1e-14 = VERY safe
    !           1e-10 = safe
    !           1e-6  = agressive, but still good accuracy
    call shtns_set_grid(shtns_c, layout, eps_polar, nlat, nphi)
  
    ! map Fortran pointers to the generated/existing C pointers
    !   shtns is now an alias of shtns_c
    !   costheta is now an alias of shtns_c%ct or equivalently shtns%ct
    !   sintheta is now an alias of shtns_c%st or equivalently shtns%st
    call c_f_pointer(cptr=shtns_c, fptr=shtns)
    call c_f_pointer(cptr=shtns%ct, fptr=costheta, shape=[shtns%nlat])
    call c_f_pointer(cptr=shtns%st, fptr=sintheta, shape=[shtns%nlat])

    !do i=1,nlat
    !   write(*,*) costheta(i), acos(costheta(i))
    !enddo
  
    ! access various quantities stored in the structure
    write(*,*) 'Total number of modes NLM=', shtns%nlm
    write(*,*) 'l_max =',shtns%lmax
    write(*,*) 'Nth =',shtns%nlat
    write(*,*) 'Nphi =',shtns%nphi
  
    ! allocate space for some stuff to be used later
    allocate(spectral(1:shtns%nlm), physical(1:shtns%nlat), &
             spec_true(1:shtns%nlm), phys_true(1:shtns%nlat) )
    spec_true(:) = (0.0_dp,0.0_dp)
    phys_true(:) = (0.0_dp,0.0_dp)
  
    !=====================
    do_real = .true.
    l = 2 ; m = 1
    !=====================

    lm = shtns_lmidx(shtns_c,l,m)
    if (do_real) then
       spec_true(lm) = (1.0_dp,0.0_dp)
    else
       spec_true(lm) = (0.0_dp,1.0_dp)
    endif

    ! intialize physical space
    Call Y_2m(costheta, phys_true, m)
    !do i=1,nlat
    !   write(*,*) costheta(i), phys_true(i)
    !enddo

    lmstart = shtns_lmidx(shtns_c,0,m)   ! get lm index of (l,m)=(0,0)
    lmstop = shtns_lmidx(shtns_c,lmax,m) ! get lm index of (l,m)=(lmax,0)

    ! reset
    spectral(:) = (0.0_dp,0.0_dp)
    physical(:) = (0.0_dp,0.0_dp)
    call spat_to_sh_ml(shtns_c, m, phys_true, spectral(lmstart:lmstop), lmax) ! to spectral
    call sh_to_spat_ml(shtns_c, m, spectral(lmstart:lmstop), physical, lmax) ! to physical
    maxerr = maxval(abs(physical - phys_true))
    write(*,*) 'LegendreTransform: Phys-->Spec-->Phys max error', maxerr

    ! reset
    physical(:) = 0.0_dp
    spectral(:) = 0.0_dp
    call sh_to_spat_ml(shtns_c, m, spec_true(lmstart:lmstop), physical, lmax) ! to physical
    call spat_to_sh_ml(shtns_c, m, physical, spectral(lmstart:lmstop), lmax) ! to spectral
    maxerr = maxval(abs(spectral - spec_true))
    write(*,*) 'LegendreTransform: Spec-->Phys-->Spec max error', maxerr

    ! reset
    physical(:) = 0.0_dp
    spectral(:) = 0.0_dp
    call sh_to_spat_ml(shtns_c, m, spec_true(lmstart:lmstop), physical, lmax) ! to physical
    maxerr = maxval(abs(physical - phys_true))
    write(*,*) 'LegendreTransform: Spectral-->Physical max error', maxerr

    do i=1,nlat
       diff = abs(physical(i) - phys_true(i))
       write(*,*) real(phys_true(i)),real(physical(i)), real(physical(i))/real(phys_true(i))
    enddo

    ! reset
!    physical(:) = 0.0_dp
!    spectral(:) = 0.0_dp
!    call spat_to_sh_ml(shtns_c, m, phys_true, spectral(lmstart:lmstop), lmax) ! to spectral
!    maxerr = maxval(abs(spectral(lmstart:lmstop) - spec_true(lmstart:lmstop)))
!    write(*,*) 'LegendreTransform: Physical-->Spectral max error', maxerr

    ! clean up stuff we made
    deallocate(spectral, physical)

    ! clean up SHTns stuff
    call shtns_unset_grid(shtns_c)
    call shtns_destroy(shtns_c)

end program write_coeffs
