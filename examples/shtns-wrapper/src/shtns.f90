module Legendre_Transforms

   use iso_c_binding

   use input_params ! access to namelist variables

   implicit none

   include 'shtns.f03'

   type(shtns_info), pointer :: shtns
   type(c_ptr) :: shtns_c

   integer :: n_phi, m_max, m_res, n_l, n_m, n_lm
   integer :: norm, layout

   contains

   subroutine initialize_legendre_transforms
      if ( (l_max .le. 0) .and. (n_l .le. 0) ) then

         if (dealias) then
            l_max = (2*n_theta-1)/3
         else
            l_max = n_theta-1
         endif

      else
         if (l_max .le. 0) l_max = n_l-1

         if (dealias) then
            n_theta = (3*(l_max+1))/2
         else
            n_theta = l_max+1
         endif
      endif

      n_phi = 2*n_theta
      m_max = l_max
      m_res = 1
      n_l = l_max+1
      n_m = m_max+1

      ! choose grid type/layout and normalizations
      if (on_the_fly) then
         layout = sht_gauss + sht_phi_contiguous
      else
         layout = sht_gauss_fly + sht_phi_contiguous
      norm = sht_orthonormal

      call shtns_verbose(verbose)

      n_threads = shtns_use_threads(n_threads)
      if (verbose .gt. 0) then
         write(*,*) 'nthreads=', nthreads
      endif

      ! initialize/allocate transforms and build useful arrays
      shtns_c = shtns_create(l_max, m_max, m_res, norm)

      ! attach a grid to the created SHT
      call shtns_set_grid(shtns_c, layout, eps_polar, n_theta, n_phi)

      ! map Fortran pointers to the generated/existing C pointers
      call c_f_pointer(cptr=shtns_c, fptr=shtns)
      !call c_f_pointer(cptr=shtns%ct, fptr=costheta, shape=[shtns%nlat])
      !call c_f_pointer(cptr=shtns%st, fptr=sintheta, shape=[shtns%nlat])

      ! get total number of (l,m) modes
      n_lm = shtns%nlm

   end subroutine initialize_legendre_transforms

   ! convert (l,m) value to mode index
   subroutine lm_index(l, m, lm_index)
      integer, intent(in) :: l, m
      integer, intent(out) :: lm_index
      lm_index = shtns_lmidx(shtns_c,l,m)
   end subroutine lm_index

   ! convert mode index to l value of (l,m) pair
   subroutine l_value(lm_index, l)
      integer, intent(in) :: lm_index
      integer, intent(out) :: l
      l = shtns_lm2l(shtns_c,lm_index)
   end subroutine l_value

   ! convert mode index to m value of (l,m) pair
   subroutine m_value(lm_index, m)
      integer, intent(in) :: lm_index
      integer, intent(out) :: m
      m = shtns_lm2m(shtns_c,lm_index)
   end subroutine m_value

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
      real(kind=dpt), intent(in) :: data_in(1:n_theta)
      real(kind=dpt), intent(out) :: data_out(1:n_lm)

      integer :: lmstart, lmstop

      call lm_index(0, m_val, lmstart)
      call lm_index(l_max, m_val, lmstop)
      call spat_to_sh_ml(shtns_c, m_val, data_in, data_out(lmstart:lmstop), l_max)

   end subroutine LT_ToSpectral_single

    ! Legendre transform, i.e., no FFT at a given m value: Spectral-->Physical
   subroutine LT_ToSpectral_single(data_in, data_out, m_val)
      integer, intent(in) :: m_val
      real(kind=dpt), intent(in) :: data_in(1:n_lm)
      real(kind=dpt), intent(out) :: data_out(1:n_theta)

      integer :: lmstart, lmstop

      call lm_index(0, m_val, lmstart)
      call lm_index(l_max, m_val, lmstop)
      call sh_to_spat_ml(shtns_c, m_val, data_in(lmstart:lmstop), data_out, l_max)

   end subroutine LT_ToSpectral_single

end module Legendre_Transforms
