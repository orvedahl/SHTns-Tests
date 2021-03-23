!
! module to return grids used in Rayleigh
!
module cheb_grids

   implicit none

   integer :: Nr
   real*8, allocatable :: x_cheb(:), radius(:)

   contains

   subroutine initialize_grids(Nradius)
      integer, intent(in) :: Nradius

      real*8 :: a,b

      ! assume 0.35 aspect, depth 1
      b = 1.0d0/(1.0d0-0.35d0)
      a = 0.35d0*b

      Nr = Nradius
      allocate(x_cheb(1:Nr),radius(1:Nr))
      call init_chebyshev_grid(a, b)
      write(*,*) 'Chebyshev grid initialized with Nr = ',Nr

   end subroutine initialize_grids

   ! zeros based chebyshev grid
   subroutine init_chebyshev_grid(a, b)
      real*8,intent(in) :: a,b
      integer :: i
      real*8 :: arg, dtheta, xlo, xhi, pi
      pi = acos(-1.0d0)
      dtheta = pi/Nr
      arg = 0.5d0*dtheta
      do i=1,Nr
         x_cheb(Nr-i+1) = cos(arg) ! fancy indexing to reverse grid
         arg = arg + dtheta
      enddo

      xlo = minval(x_cheb)
      xhi = maxval(x_cheb)
      do i=1,Nr
         radius(i) = (x_cheb(i)-xlo)*(b-a)/(xhi-xlo) + a
      enddo
   end subroutine init_chebyshev_grid

   subroutine cleanup_grids()
      if (allocated(x_cheb)) then
         deallocate(x_cheb)
      endif
      if (allocated(radius)) then
         deallocate(radius)
      endif
   end subroutine cleanup_grids

end module cheb_grids
