module Legendre_Transforms

   use ProblemSize
   !use Legendre_Polynomials
   use iso_c_binding

   implicit none

   include 'shtns.f03'

   type(shtns_info), pointer :: shtns
   type(c_ptr) :: shtns_c

   integer :: m_res, n_lm
   integer :: norm, layout, num_threads

   contains

   subroutine initialize_legendre_transforms(eps_thresh, n_threads_in, on_the_fly, verbose)
      real(kind=dpt), intent(in) :: eps_thresh
      integer, intent(in) :: n_threads_in, verbose
      logical, intent(in) :: on_the_fly

      m_res = 1

      ! choose grid type/layout and normalizations
      if (on_the_fly) then
         layout = sht_gauss + sht_phi_contiguous
      else
         layout = sht_gauss_fly + sht_phi_contiguous
      endif
      norm = sht_orthonormal

      call shtns_verbose(verbose)

      num_threads = shtns_use_threads(n_threads_in)
      if (verbose .gt. 0) then
         write(*,*) 'nthreads=', num_threads
      endif

      ! initialize/allocate transforms and build useful arrays
      shtns_c = shtns_create(l_max, m_max, m_res, norm)

      ! attach a grid to the created SHT
      call shtns_set_grid(shtns_c, layout, eps_thresh, n_theta, n_phi)

      ! map Fortran pointers to the generated/existing C pointers
      call c_f_pointer(cptr=shtns_c, fptr=shtns)
      !call c_f_pointer(cptr=shtns%ct, fptr=costheta, shape=[shtns%nlat])
      !call c_f_pointer(cptr=shtns%st, fptr=sintheta, shape=[shtns%nlat])

      ! get total number of (l,m) modes
      n_lm = shtns%nlm
      write(*,*) 'lmax',shtns%lmax
      write(*,*) 'mmax',shtns%mmax
      write(*,*) 'nth', shtns%nlat
      write(*,*) 'nphi', shtns%nphi
      write(*,*) 'nlm', shtns%nlm

   end subroutine initialize_legendre_transforms

   ! convert (l,m) value to mode index
   function lm_index(l, m) result(lm_idx)
      integer, intent(in) :: l, m
      integer :: lm_idx
      lm_idx = shtns_lmidx(shtns_c,l,m)
   end function lm_index

   ! convert mode index to l value of (l,m) pair
   function l_value(lm_idx) result(l)
      integer, intent(in) :: lm_idx
      integer :: l
      l = shtns_lm2l(shtns_c,lm_idx)
   end function l_value

   ! convert mode index to m value of (l,m) pair
   function m_value(lm_idx) result(m)
      integer, intent(in) :: lm_idx
      integer :: m
      m = shtns_lm2m(shtns_c,lm_idx)
   end function m_value

   ! get quadrature weights
   subroutine gauss_weights(weights)
      real(kind=dpt), intent(inout) :: weights(1:n_theta)
      call shtns_gauss_wts(shtns_c, weights)
   end subroutine gauss_weights

   subroutine finalize_legendre_transforms()
      ! clean up SHTns stuff
      call shtns_unset_grid(shtns_c)
      call shtns_destroy(shtns_c)
   end subroutine finalize_legendre_transforms

   ! Legrendre Transform, i.e., no FFT at a given m value: Physical-->Spectral
   subroutine LT_ToSpectral_single(data_in, data_out, m_val)
      integer, intent(in) :: m_val
      complex(kind=dpt), intent(inout) :: data_in(1:n_theta)
      complex(kind=dpt), intent(inout) :: data_out(1:n_lm)

      integer :: lmstart, lmstop

      lmstart = lm_index(0, m_val)
      lmstop  = lm_index(l_max, m_val)
      call spat_to_sh_ml(shtns_c, m_val, data_in, data_out(lmstart:lmstop), l_max)

   end subroutine LT_ToSpectral_single

    ! Legendre transform, i.e., no FFT at a given m value: Spectral-->Physical
   subroutine LT_ToPhysical_single(data_in, data_out, m_val)
      integer, intent(in) :: m_val
      complex(kind=dpt), intent(inout) :: data_in(1:n_lm)
      complex(kind=dpt), intent(inout) :: data_out(1:n_theta)

      integer :: lmstart, lmstop

      lmstart = lm_index(0, m_val)
      lmstop  = lm_index(l_max, m_val)
      call sh_to_spat_ml(shtns_c, m_val, data_in(lmstart:lmstop), data_out, l_max)

   end subroutine LT_ToPhysical_single

end module Legendre_Transforms
