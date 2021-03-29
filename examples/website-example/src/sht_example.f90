program sht_example
  
    use iso_c_binding
    use iso_fortran_env, only: real64

    implicit none
  
    include 'shtns.f03'
  
    integer, parameter :: dp=real64
    real(dp), parameter :: pi=acos(-1.0_dp)
  
    type(shtns_info), pointer :: shtns
    real(dp), pointer :: cosTheta(:), sinTheta(:)
    type(c_ptr) :: shtns_c
  
    integer :: lmax, mmax, mres, nthreads
    integer :: nlat, nphi, lm, l, m, lmStart, lmStop
    integer :: norm, layout, i
    real(dp) :: eps_polar, err, err1, err2, err3
    complex(dp), allocatable :: spec(:), phys(:)
  
    ! set size of transforms
    lmax = 31
    mmax = lmax              ! max azimuthal wave number
    mres = 1                 ! azimuthal periodicity is 2*pi/mres
    nlat = int(3*(lmax+1)/2) ! number of points in theta
    nphi = 2*nlat            ! number of points in phi
  
    eps_polar = 1.0e-10_dp   ! set polar optimisation threshold
  
    !-- Normalization and layout for SHTs
    norm = sht_orthonormal
    layout = sht_gauss + sht_phi_contiguous + sht_scalar_only + sht_south_pole_first

    call shtns_verbose(2)
    nthreads = shtns_use_threads(0)
  
    ! initialize/allocate transforms and build useful arrays
    shtns_c = shtns_create(lmax, mmax, mres, norm)

    ! attach a grid to the created SHT
    call shtns_set_grid(shtns_c, layout, eps_polar, nlat, nphi)
  
    call c_f_pointer(cptr=shtns_c, fptr=shtns)
    call c_f_pointer(cptr=shtns%ct, fptr=costheta, shape=[shtns%nlat])
    call c_f_pointer(cptr=shtns%st, fptr=sintheta, shape=[shtns%nlat])
  
    ! allocate space for some stuff to be used later
    allocate( spec(1:shtns%nlm), phys(1:shtns%nlat))
    spec(:) = 0.0_dp

    ! set physical space data: pure Y_2^1
    l = 2
    m = 1
    phys(:) = -3.0_dp*sqrt(5.0_dp/(4.0_dp*pi*6.0_dp))*sintheta(:)*costheta(:)
  
    ! Legrendre Transform, i.e., no FFT at a given m value: Physical-->Spectral
    lmstart = shtns_lmidx(shtns_c,m,m)
    lmstop = shtns_lmidx(shtns_c,lmax,m)
    call spat_to_sh_ml(shtns_c, m, phys, spec(lmstart:lmstop), lmax)

    write(*,*)
    write(*,*)
    write(*,*) 'Only expect power in l=',l,'m=',m
    write(*,*) 'l    m    coefficient'
    do i=1,shtns%nlm
       if (abs(spec(i)) .gt. 1e-8) then
          write(*,*) shtns_lm2l(shtns_c,i), shtns_lm2m(shtns_c,i), spec(i)
       endif
    enddo

    lm = shtns_lmidx(shtns_c,l,m)

    ! compute error (assume lm is not first or last)
    err1 = maxval(abs(spec(1:lm-1)))
    err2 = abs(spec(lm) - 1.0_dp)
    err3 = maxval(abs(spec(lm+1:)))
    err = max(err1,err2,err3)
    write(*,*) 'LegendreTransform: Physical-->Spectral max error', err
    write(*,*)

    ! clean up stuff we made
    deallocate(spec, phys)

    ! clean up SHTns stuff
    call shtns_unset_grid(shtns_c)
    call shtns_destroy(shtns_c)
  
end program sht_example
