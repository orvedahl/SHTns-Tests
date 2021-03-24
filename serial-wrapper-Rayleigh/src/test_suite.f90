
Module test_suite

   use cheb_grids ! holds radial grids

   Use Structures
   Use Legendre_Polynomials
   Use Legendre_Transforms, only: Legendre_Transform

   implicit none

   real*8 :: pi = acos(-1.0d0)

   real*8, allocatable :: costheta(:) ! ideally, this comes in from ProblemSize

   contains

   ! Y_l^m for l=2, does not include the e^(i*m*phi) part
   subroutine Y_2m(costh, Yout, m)
      integer, intent(in) :: m
      real*8, intent(in) :: costh(:)
      real*8, intent(inout) :: Yout(:)

      real*8 :: A

      if (m .eq. 0) then
         A = sqrt(5.0d0 / (4.0d0*pi))
         Yout(:) = A*0.5d0*(3.0d0*costh(:)*costh(:)-1.0d0)

      elseif (m .eq. 1) Then
         A = sqrt(5.0d0 / (4.0d0*pi*6.0d0))
         Yout(:) = -3.0d0*A*costh(:)*sqrt(1.0d0-costh(:)*costh(:))

      elseif (m .eq. 2) Then
         A = sqrt(5.0d0 / (4.0d0*pi*24.0d0))
         Yout(:) = 3.0d0*A*(1.0d0-costh(:)*costh(:))

      else
         Yout(:) = -1.0d0
      endif
   end subroutine Y_2m

   !-------------------------
   ! Rayleigh interfaces
   !-------------------------
   subroutine test_LT(ntest)
      integer, intent(in) :: ntest

      type(rmcontainer4d), allocatable :: spectral(:)
      real*8, allocatable :: physical(:,:,:), true_phys(:)
      character(len=128) :: msg
      integer :: mp, l, r, m, nfields, stat, nq, mind, i
      real*8 :: diff1, diff2, diff3, diff, mxdiff

      allocate(costheta(1:n_theta)) ! only necessary b/c coloc=real*16
      costheta(:) = coloc(:)

      write(*,*)
      write(*,*) '============================'
      write(*,*) 'Using the Rayleigh interface'
      write(*,*) '============================'
      write(*,*)

      nfields = 3

      allocate(spectral(1:n_m), stat=stat, errmsg=msg)
      if (stat .ne. 0) then
         write(*,*) 'allocate spectral(1:nm) failed'
         write(*,*) msg
      endif
      do mp=1,n_m
         m = m_values(mp)
         allocate(spectral(mp)%data(m:l_max,1:nr,1:2,nfields), stat=stat, errmsg=msg)
         if (stat .ne. 0) then
            write(*,*) 'allocate spectral(mp) failed, mp=',mp,'m_value',m
            write(*,*) msg
         endif
      enddo

      nq = 2*nr*nfields
      allocate(physical(1:n_theta,1:nq,1:n_m), stat=stat, errmsg=msg)
      if (stat .ne. 0) then
         write(*,*) 'allocate physical failed'
         write(*,*) msg
      endif
      physical(:,:,:) = 0.0d0

      !--------------------
      ! setup/run the test
      !--------------------
      if (ntest .eq. 1) then

         mind = -1
         do mp = 1, n_m ! modes
            m = m_values(mp)
            if (m .eq. 1) mind = mp
            do l = m, l_max ! \ell values

               spectral(mp)%data(l,:,:,:) = 0.0d0
               if ((l .eq. 2) .and. (m .eq. 1)) then
                  do r = 1, nr ! radius
                     spectral(mp)%data(l,r,1,1) = 1.0d0*radius(r)
                     spectral(mp)%data(l,r,2,2) = 2.0d0*radius(r)**2
                     spectral(mp)%data(l,r,2,3) = 5.0d0*radius(r)**3
                  enddo
               endif
            enddo
         enddo

         allocate(true_phys(1:n_theta)) ! build expected answer
         call Y_2m(costheta, true_phys, 1) ! "1" is m value as used in above initialization

         ! move to physical space
         write(*,*) 'S-->P'
         call Legendre_Transform(spectral, physical)

         ! check error
         write(*,*) 'expected         Rayleigh'
         do i=1,n_theta
            write(*,*) true_phys(i), physical(i,:,mind)
         enddo

         ! move back to spectral space
         write(*,*) 'P-->S'
         call Legendre_Transform(physical, spectral)

         ! zero out l_max mode
         do mp=1,n_m
            spectral(mp)%data(l_max,:,:,:) = 0.0d0
         enddo

         ! compute error
         diff = -1.0d0
         mxdiff = -1.0d0
         do mp = 1, n_m
            m = m_values(mp)
            do l = m, l_max
               if ((l .eq. 2) .and. (m .eq. 1)) then
                   diff1 = maxval(abs(1.0d0*radius(:) - spectral(mp)%data(l,:,1,1)))
                   diff2 = maxval(abs(2.0d0*radius(:)**2 - spectral(mp)%data(l,:,2,2)))
                   diff3 = maxval(abs(5.0d0*radius(:)**3 - spectral(mp)%data(l,:,2,3)))
                   diff = max(diff1,diff2,diff3)
               else
                   diff = maxval(abs(spectral(mp)%data(l,:,:,:)))
               endif
               if (diff .gt. mxdiff) mxdiff = diff
            enddo
         enddo
         write(*,*) 'Spec-->Phys-->Spec max error',mxdiff

         deallocate(true_phys)
      else
         write(*,*) 'Test number not written yet, ntest=',ntest
      endif

      ! cleanup
      do mp=1,n_m
         deallocate(spectral(mp)%data)
      enddo
      deallocate(spectral, physical)

      deallocate(costheta)

   end subroutine test_LT

End Module test_suite

