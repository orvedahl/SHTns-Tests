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
    !integer :: nlat, nphi, n_p, lm, l, m, lmStart, lmStop
    integer :: nthreads, n_p, lm, l, m, lmStart, lmStop
    integer :: nlm, norm, nout, layout
    real(dp) :: eps_polar
    real(dp), allocatable :: Sh(:,:)
    complex(dp), allocatable :: Slm(:), ShT(:)
  
    ! my routine
    call runtime_init()
  
    !--  Set the size of the transforms
    !lmax = 5
    !mmax = 3
    !mres = 1
    !nphi = 10
    !nlat = 64
  
    !-- Polar optimisation threshold
    eps_polar = 1.0e-10_dp
  
    !-- Norm and layout for SHTs
    norm = sht_orthonormal + sht_no_cs_phase
    layout = sht_gauss + sht_phi_contiguous
  
    call shtns_verbose(2)
    nthreads = shtns_use_threads(0)
    print*, 'nthreads=', nthreads
  
    shtns_c = shtns_create(lmax, mmax, mres, norm)
    call shtns_set_grid(shtns_c, layout, eps_polar, nlat, nphi)
  
    !-- C/Fortran pointer mapping
    call c_f_pointer(cptr=shtns_c, fptr=shtns)
    call c_f_pointer(cptr=shtns%ct, fptr=costheta, shape=[shtns%nlat])
    call c_f_pointer(cptr=shtns%st, fptr=sintheta, shape=[shtns%nlat])
  
    print*, 'cosTheta', costheta
    print*, 'sinTheta', sintheta
  
    print*, 'NLM=', shtns%nlm
  
    allocate( sh(1:shtns%nphi,1:shtns%nlat), slm(1:shtns%nlm), sht(1:shtns%nlat) )
    slm(:)=0.0_dp
  
    !-- Get lm index from the (l,m) pair
    l = 1 
    m = 0
    lm = shtns_lmidx(shtns_c,l,m)
    !-- Set the l=1, m=0 mode to 1.0
    slm(lm)=(1.0_dp,0.0_dp)
  
    !-- Get (l,m) from lm
    lm = 17
    l = shtns_lm2l(shtns_c,lm)
    m = shtns_lm2m(shtns_c,lm)
    print*, '(l,m)=', l,m
  
    !-- Spec -> Spat
    call sh_to_spat(shtns_c, slm, sh)
  
    print*, 'Y(1,0)', sh(1,:) ! It should be sqrt(3/4/pi) * cos(theta)
  
    !-- Spat -> Spec
    call spat_to_sh(shtns_c, sh, slm)
  
    !-- print S(1,0) to check it 1.0 is recovered.
    print*, 'S(1,0)', slm(lm)
  
    !-- Legendre only for m=0
    lmstart = shtns_lmidx(shtns_c,0,0)
    lmstop = shtns_lmidx(shtns_c,lmax,0)
    call sh_to_spat_ml(shtns_c, 0, slm(lmstart:lmstop), sht, lmax)
    print*, 'Legendre only', sht
  
    deallocate(sh, slm, sht)
    call shtns_unset_grid(shtns_c)
    call shtns_destroy(shtns_c)
  
end program sht_example
