module ProblemSize ! mimic the Rayleigh module

   use data_types
   use input_params, only: l_max, runtime_init
   use Legendre_Polynomials, only: initialize_legendre
   !use Legendre_Transforms, only: initialize_legendre_transforms

   implicit none

   integer :: n_theta
   integer :: n_phi
   !integer :: l_max   ! l_max is declared in the input_params module
   integer :: n_l      ! this should really be a 1d array...
   integer :: m_max, n_m
   integer :: m_values ! this should really be a 1d array...

   ! this module also holds:
   !   parallel layout information
   !   cos(th),sin(th),cos^2(th),phi,cos(phi),etc.

   contains

   subroutine init_problemsize(eps_thresh, n_threads_in, on_the_fly, verbose)
      real(kind=dpt), intent(in) :: eps_thresh
      integer, intent(in) :: n_threads_in, verbose
      logical, intent(in) :: on_the_fly

      call runtime_init() ! read namelist values

      m_values = l_max ! this is really set by parallel layout and should be an array
      m_max = l_max

      n_theta = 3*(l_max+1)/2

      call initialize_legendre(n_theta, l_max, m_values)
      n_phi = 2*n_theta

      !call initialize_legendre_transforms(eps_thresh, n_threads_in, on_the_fly, verbose)

      ! initialize angular derivative coefficients

      ! build various grid quantities: cos(th),sin(th),cos^2(th),phi,cos(phi),etc.

   end subroutine init_problemsize

end module ProblemSize
