!
! module to return grids used in Rayleigh
!
module grids

   use data_types
   use variables

   implicit none

   ! chebyshev grid
   integer :: Nr
   real(kind=dpt), allocatable :: x_cheb(:), theta_cheb(:)

   ! angular grid
   integer :: l_max, Nth, Nphi
   real(kind=dpt), allocatable :: x_leg(:), theta(:), weights_leg(:), phi(:)
   real(kind=dpt), allocatable :: Plm(:,:,:)

   contains

   subroutine initialize_grids(Nradius, Ntheta)
      integer, intent(in) :: Nradius, Ntheta

      Nr = Nradius
      allocate(x_cheb(1:Nr),theta_cheb(1:Nr))
      call init_chebyshev_grid()
      write(*,*) 'Chebyshev grid initialized with Nr = ',Nr

      Nth = Ntheta
      l_max = int(2*Nth/3.0 - 1)
      allocate(x_leg(1:Nth), theta(1:Nth), weights_leg(1:Nth), Plm(1:Nth,0:l_max,0:l_max))
      call init_legendre_grid()
      write(*,*) 'Legendre grid initialized with Nth = ',Nth

      Nphi = 2*Nth
      allocate(phi(1:Nphi))
      call init_fourier_grid()
      write(*,*) 'Fourier grid initialized with Nphi = ',Nphi

   end subroutine initialize_grids

   ! zeros based chebyshev grid
   subroutine init_chebyshev_grid()
      integer :: i
      real(kind=dpt) :: arg, dtheta
      dtheta = pi/Nr
      arg = 0.5_dpt*dtheta
      do i=1,Nr
         theta_cheb(i) = arg
         x_cheb(i) = cos(arg)
         arg = arg + dtheta
      enddo
   end subroutine init_chebyshev_grid

   subroutine init_legendre_grid()
      call legendre_roots(Nth, x_leg, weights_leg)
      call compute_plm(x_leg, Nth, l_max, Plm)
      theta(:) = acos(x_leg(:))
   end subroutine init_legendre_grid

   subroutine init_fourier_grid()
      real(kind=dpt) :: dphi, twopi
      integer :: i
      twopi = 2.0_dpt*pi
      dphi = twopi/Nphi
      do i=1,Nphi
         phi(i) = (i-1)*dphi
      enddo
   end subroutine init_fourier_grid

   subroutine cleanup_grids()
      if (allocated(x_cheb)) then
         deallocate(x_cheb)
      endif
      if (allocated(theta_cheb)) then
         deallocate(theta_cheb)
      endif
      if (allocated(x_leg)) then
         deallocate(x_leg)
      endif
      if (allocated(theta)) then
         deallocate(theta)
      endif
      if (allocated(weights_leg)) then
         deallocate(weights_leg)
      endif
      if (allocated(Plm)) then
         deallocate(Plm)
      endif
      if (allocated(phi)) then
         deallocate(phi)
      endif
   end subroutine cleanup_grids

   !-----------------------------------------------------------------
   ! compute Alm*Pl(x, l, m)
   !    Alm**2 = (2*l+1)/(4*pi) * (l-m)! / (l+m)!
   ! assumes nm = nl = lmax+1
   ! returns array of shape (x,l,m)
   !-----------------------------------------------------------------
   subroutine compute_plm(x, n, lmax, Pl)
      integer, intent(in) :: n, lmax
      double precision, intent(in) :: x(1:n)

      ! assume nl = nm = lmax+1
      real(kind=dpt), intent(out) :: Pl(1:n,0:lmax,0:lmax) ! size(n, lmax+1, nm)

      integer :: i, m, l, mv
      real(kind=dpt) :: ratio, amp, tmp, xi

      Pl(:,:,:) = 0.0_dpt

      do m=1,lmax+1 ! 1-based array indexing

         mv = m - 1 ! convert to azimuthal integer

         ! first do the l=m & l=m+1 pieces
         ratio = 1.0_dpt
         call compute_factorial_ratio(mv, ratio)
         amp = ((mv+0.5_dpt)/(2.0_dpt*pi))**0.5_dpt
         amp = amp*ratio
         do i=1,n
            xi = x(i)
            tmp = 1.0_dpt - xi*xi
            if (mod(mv,2) .eq. 1) then
               Pl(i,mv,mv) = -amp*tmp**(mv/2+0.5_dpt) ! odd m
            else
               Pl(i,mv,mv) =  amp*tmp**(mv/2)       ! even m
            endif

            ! l=m+1 part
            if (mv .lt. lmax) then
               Pl(i,mv+1,mv) = Pl(i,mv,mv)*xi*(2.0_dpt*mv+3)**0.5_dpt
            endif
         enddo

         ! now do the l>m+1
         do l=mv+2,lmax
            do i=1,n
               xi = x(i)
               amp = (l-1)**2-mv*mv
               amp = amp / (4.0_dpt*(l-1)**2-1.0_dpt)
               amp = amp**0.5_dpt
               tmp = Pl(i,l-1,mv)*xi - amp*Pl(i,l-2,mv)
               amp = (4.0_dpt*l*l-1.0_dpt)/(l*l-mv*mv)
               Pl(i,l,mv) = tmp*amp**0.5_dpt
            enddo
         enddo
      enddo

   end subroutine compute_plm

   !-----------------------------------------------------------------
   ! compute factorial ratio = sqrt( (2m)! / 4**m / m! / m! )
   !-----------------------------------------------------------------
   subroutine compute_factorial_ratio(m, ratio)
      integer, intent(in) :: m
      real(kind=dpt), intent(out) :: ratio
      integer :: i
      ratio = 1.0_dpt
      do i=1,m
         ratio = ratio*((i-0.5_dpt)/i)**0.5_dpt
      enddo
   end subroutine compute_factorial_ratio

   !-----------------------------------------------------------------
   ! Legendre grid points ordered as x(i) < x(i+1) and x in (-1,1)
   !-----------------------------------------------------------------
   subroutine legendre_roots(n, x, weights)
      integer, intent(in) :: n
      real(kind=dpt), intent(out) :: x(1:n), weights(1:n)

      integer :: i, n_roots
      logical :: converged
      real(kind=dpt) :: midpoint, scaling, eps, ith_root(1:1)
      real(kind=dpt) :: Np_half, new_guess, delta, Pl(1:1), dPldx(1:1)

      midpoint = 0.0_dpt
      scaling  = 1.0_dpt
      n_roots  = (n+1)/2
      Np_half = n + 0.5_dpt
      eps = 3.0d-15

      do i=1,n_roots
         ith_root(1) = cos(pi*(i-0.25_dpt)/Np_half)
         converged = .false.
         do
            ! evaluate Pl at single root
            call eval_nth_legendre_quad(ith_root, 1, n, Pl, dPldx)

            new_guess = ith_root(1) - Pl(1)/dPldx(1)
            delta = abs(ith_root(1) - new_guess)
            ith_root(1) = new_guess
            if (delta .le. eps) then
               converged = .true.
            endif

            x(i)     = midpoint - scaling*ith_root(1)
            x(n+1-i) = midpoint + scaling*ith_root(1)

            weights(i)     = 2.0_dpt*scaling/((1.-ith_root(1)*ith_root(1))*dPldx(1)*dPldx(1))
            weights(n+1-i) = weights(i)

            if (converged) exit
         enddo
      enddo

   end subroutine legendre_roots

   !-----------------------------------------------------------------
   ! evaluate P_l and d(P_l)/dx at the given grid points
   !-----------------------------------------------------------------
   subroutine eval_nth_legendre_quad(x, m, l, Pl, dPldx)
      integer, intent(in) :: m, l
      real(kind=dpt), intent(in) :: x(1:m)
      real(kind=dpt), intent(out) :: Pl(1:m), dPldx(1:m)

      integer :: i
      real(kind=dpt) :: Pl_minus1(1:m), Pl_minus2(1:m)

      Pl_minus1(:) = 0.0_dpt
      Pl(:) = 1.0_dpt

      ! use recursion relation
      do i=1,l
         Pl_minus2(:) = Pl_minus1(:)
         Pl_minus1(:) = Pl(:)
         Pl(:) = ((2.0_dpt*i-1.0_dpt)*x(:)*Pl_minus1(:) - (i-1.0_dpt)*Pl_minus2(:))/i
      enddo

      ! get derivative
      dPldx(:) = l*(x(:)*Pl(:) - Pl_minus1(:))/(x(:)*x(:) - 1._dpt)

   end subroutine eval_nth_legendre_quad

end module grids
