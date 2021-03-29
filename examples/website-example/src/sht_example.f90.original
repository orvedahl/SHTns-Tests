program sht_example
    !
    ! This is a modern Fortran example that shows how to set-up the spherical
    ! harmonic transforms using SHTns.
    !
  
    use iso_c_binding
    use iso_fortran_env, only: real64

    use input_params

    implicit none
  
    include 'shtns.f03'
  
    integer, parameter :: dp=real64
    real(dp), parameter :: pi=acos(-1.0_dp)
  
    type(shtns_info), pointer :: shtns
    real(dp), pointer :: cosTheta(:), sinTheta(:)
    type(c_ptr) :: shtns_c
  
    !integer :: lmax, mmax, mres, nthreads
    integer :: mmax, mres, nthreads
    integer :: nlat, nphi, n_p, lm, l, m, lmStart, lmStop
    integer :: nlm, norm, nout, layout, i
    real(dp) :: eps_polar, err
    real(dp), allocatable :: Sh(:,:), weights(:)
    complex(dp), allocatable :: Slm(:), ShT(:)
  
    ! my routine
    call runtime_init()
  
    ! set size of transforms, lmax is set in the input file
    mmax = lmax              ! max azimuthal wave number
    mres = 1                 ! azimuthal periodicity is 2*pi/mres
    nlat = int(3*(lmax+1)/2) ! number of points in theta
    nphi = 2*nlat            ! number of points in phi
  
    ! set polar optimisation threshold
    eps_polar = 1.0e-10_dp
  
    !-- Normalization and layout for SHTs
    norm = sht_orthonormal
    layout = sht_gauss + sht_phi_contiguous ! use Legendre grid/quadrature
    !layout = sht_gauss_fly + sht_phi_contiguous ! use on-the-fly Legendre grid/quadrature
  
    ! set output level: 0=none, 1=some, 2=debug
    call shtns_verbose(2)

    ! return number of OpenMP threads that will be used, try to use passed value
    !   send_threads > 0  : set max number of threads to use
    !   send_threads <= 0 : max number of threads is set to number of processors
    !   send_threads = 1  : disable OpenMP
    nthreads = shtns_use_threads(0)
    write(*,*) 'nthreads=', nthreads
  
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
  
    write(*,*) 'cos(theta)'
    do i=1,5
       write(*,*) '  ', costheta(i)
    enddo
  
    ! access various quantities stored in the structure
    write(*,*) 'Total number of modes NLM=', shtns%nlm
    write(*,*) 'l_max =',shtns%lmax
    write(*,*) 'Nth =',shtns%nlat
    write(*,*) 'Nphi =',shtns%nphi
    write(*,*) 'Nspat =',shtns%nspat
    write(*,*) 'Nth*Nphi =',shtns%nlat*shtns%nphi
  
    ! allocate space for some stuff to be used later
    !   sh  : Physical representation, size(nphi,nth)
    !   slm : Spectral representation, size(nlm)
    allocate( sh(1:shtns%nphi,1:shtns%nlat), slm(1:shtns%nlm), sht(1:shtns%nlat) )
    allocate( weights(1:shtns%nlat) )
  
    ! helper function to convert (l,m) pair into mode index lm
    l = 1 
    m = 0
    lm = shtns_lmidx(shtns_c,l,m)
    write(*,*) 'global index to (l,m)=',l,m,' is lm=',lm

    ! set spectral (1,0) coefficient
    lm = shtns_lmidx(shtns_c,1,0)
    slm(:) = 0.0_dp
    slm(lm) = (1.0_dp,0.0_dp)
  
    ! helper function to convert mode index lm into l and m
    lm = 7
    l = shtns_lm2l(shtns_c,lm)
    m = shtns_lm2m(shtns_c,lm)
    write(*,*) 'global index lm=',lm,'leads to (l,m)=',l,m

    ! get the quadrature weights
    call shtns_gauss_wts(shtns_c, weights)
    write(*,*) 'weights(theta)'
    do i=1,5
       write(*,*) '  ', weights(i)
    enddo
  
    ! full SHT from Spectral to Physical
    !   slm : size(nphi,nth)
    !   sh  : size(nlm)
    call sh_to_spat(shtns_c, slm, sh)
  
    err = maxval(abs(sh(1,:) - sqrt(0.75_dp/pi)*costheta(:)))
    write(*,*) 'Spectral-->Physical max error', err
  
    ! full SHT from Physical to Spectral
    !   sh  : size(nlm)
    !   slm : size(nphi,nth)
    call spat_to_sh(shtns_c, sh, slm)
  
    lm = shtns_lmidx(shtns_c,1,0)
    err = abs(slm(lm) - 1.0_dp)
    write(*,*) 'Physical-->Spectral max error', err
  
    ! Legendre transform, i.e., no FFT at a given m value: Spectral-->Physical
    !   get starting l index, i.e., lm index of l start, usually will be 0
    !   get ending l index, i.e., lm index of l stop, usually will be lmax
    !   sht : physical space of size(nth)
    !   0 : the m value
    !   lmax : the max spherical harmonic to use (could do a truncated transform)
    lmstart = shtns_lmidx(shtns_c,0,0)   ! get lm index of (l,m)=(0,0)
    lmstop = shtns_lmidx(shtns_c,lmax,0) ! get lm index of (l,m)=(lmax,0)
    call sh_to_spat_ml(shtns_c, 0, slm(lmstart:lmstop), sht, lmax)
  
    err = maxval(abs(sht(:) - sqrt(0.75_dp/pi)*costheta(:)))
    write(*,*) 'LegendreTransform: Spectral-->Physical max error', err

    ! Legrendre Transform, i.e., no FFT at a given m value: Physical-->Spectral
    call spat_to_sh_ml(shtns_c, 0, sht, slm(lmstart:lmstop), lmax)

    lm = shtns_lmidx(shtns_c,1,0)
    err = abs(slm(lm) - 1.0_dp)
    write(*,*) 'LegendreTransform: Physical-->Spectral max error', err

    ! clean up stuff we made
    deallocate(sh, slm, sht, weights)

    ! clean up SHTns stuff
    call shtns_unset_grid(shtns_c)
    call shtns_destroy(shtns_c)
  
end program sht_example
