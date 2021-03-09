module Legendre_Polynomials ! mimic the Rayleigh module

   implicit none

   integer :: n_theta
   integer :: l_max, n_m
   integer :: n_l      ! this should really be a 1d array...
   integer :: m_values ! this should really be a 1d array...
   integer :: m_max    ! does this belong here?

   ! this module also holds:
   !   x colocation grid, x=cos(th)
   !   weights, which are not renormalized
   !   P_l^m(cos(th))  - evaluate Plm on theta grid
   !   iP_l^m(cos(th)) - evaluate Plm on theta grid, including integration weights

   contains

   subroutine initialize_legendre(n_th, lmax, mval) ! (basically copy ProblemSize values)
      integer, intent(in) :: n_th, lmax
      integer, intent(in) :: mval
      n_theta = n_th
      l_max = lmax
      m_values = mval ! should be array...

      ! derived quantities
      m_max = l_max
      n_l = l_max+1 ! should be array...
      n_m = m_max+1
   end subroutine initialize_legendre

   subroutine finalize_legendre() ! handle deallocations if any
      return
   end subroutine finalize_legendre

end module Legendre_Polynomials
