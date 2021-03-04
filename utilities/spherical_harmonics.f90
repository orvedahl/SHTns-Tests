!
! module to handle simple spherical harmonic tasks
!
module spherical_harmonics

   use data_types
   use variables
   use grids

   implicit none

   contains

   ! find index of theta coordinate
   subroutine theta_lookup(th, ind)
      real(kind=dpt), intent(in) :: th
      integer, intent(out) :: ind

      real(kind=dpt) :: eps
      integer :: i

      eps = 1e-9

      ind = -1
      do i=1,Nth
         if (abs(th-theta(i)) .le. eps) then
             ind = i
             exit
         endif
      enddo

   end subroutine theta_lookup

   ! evaluate spherical harmonic at given theta/phi coordinate
   subroutine Ylm(ell, m, theta_coord, phi_coord, ylm_value_r, ylm_value_i)
      integer, intent(in) :: ell, m
      real(kind=dpt), intent(in) :: theta_coord, phi_coord
      real(kind=dpt), intent(inout) :: ylm_value_r, ylm_value_i

      integer :: ind

      if (m .lt. 0) then
         write(*,*) "Ylm: negative values are not supported, convert them"
         stop
      endif

      call theta_lookup(theta_coord, ind)
      if (ind .lt. 0) then
         write(*,*) "Ylm: theta value not found"
         stop
      endif

      ylm_value_r = Plm(ind,ell,m)*cos(m*phi_coord)
      ylm_value_i = Plm(ind,ell,m)*sin(m*phi_coord)

   end subroutine Ylm

end module spherical_harmonics
