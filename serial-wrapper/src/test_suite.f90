
Module test_suite

   use cheb_grids ! holds radial grids

   Use Structures
   Use Legendre_Polynomials
#ifdef USE_SHTns
   Use Legendre_Transforms_SHTns
#else
   Use Legendre_Transforms, only: Legendre_Transform
#endif

   implicit none

   real*8 :: pi = acos(-1.0d0)

   real*8, allocatable :: costheta(:) ! ideally, this comes in from ProblemSize

   contains

   ! Y_l^m for l=1, does not include the e^(i*m*phi) part
   subroutine Y_1m(costh, Yout, m)
      integer, intent(in) :: m
      real*8, intent(in) :: costh(:)
      real*8, intent(inout) :: Yout(:)

      real*8 :: A

      if (m .eq. 0) then
         A = sqrt(3.0d0 / (4.0d0*pi))
         Yout(:) = A*costh(:)

      elseif (m .eq. 1) then
         A = sqrt(3.0d0 / (4.0d0*pi*2.0d0))
         Yout(:) = -A*sqrt(1.0d0-costh(:)*costh(:))

      else
         Yout(:) = -1.0d0
      endif
   end subroutine Y_1m

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
      integer :: mp, l, r, m, nfields, stat, nq, mind, i, f, imi, rlo, rhi
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
         imi = 2
         f = 2
         rlo = 1 + (imi-1)*nr + (f-1)*nr*2  ! index = (r-rlo+1) + (imi-1)*Nr + (f-1)*Nr*2
         rhi = nr + (imi-1)*nr + (f-1)*nr*2
         mxdiff = -1.0d0
         do i=1,n_theta
            diff = maxval(abs(physical(i,rlo:rhi,mind)-true_phys(i)*2*radius(:)**2))
            if (diff .gt. mxdiff) mxdiff = diff
         enddo
         write(*,*) 'Spec-->Phys max error',mxdiff

         ! move back to spectral space
         write(*,*)
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

#ifdef USE_SHTns
   !-------------------------
   ! SHTns interfaces
   !-------------------------
   subroutine test_SHTns(ntest)
      integer, intent(in) :: ntest

      type(rmcontainer4d), allocatable :: Rspectral(:)
      real*8, allocatable :: Rphysical(:,:,:)
      complex*16, allocatable :: spectral(:), physical(:)
      integer :: l, m, lm, i, lval, mval, f, imi, mp, nfields, r, rhi, rlo, mind, nq
      real*8 :: diff, diff1, diff2, diff3, mxdiff
      real*8, allocatable :: true_phys(:)

      allocate(costheta(1:n_theta)) ! only necessary b/c coloc=real*16
      costheta(:) = coloc(:)

      write(*,*)
      write(*,*) '============================'
      write(*,*) 'Using the SHTns interface'
      write(*,*) '============================'
      write(*,*)

      allocate(spectral(1:SHTns%nlm), physical(1:SHTns%nlat))
      spectral(:) = 0.0d0
      physical(:) = 0.0d0

      if (ntest .eq. 0) then
         l = 1; m = 1
         lm = SHTns_lm_index(l,m)

         allocate(true_phys(1:SHTns%nlat)) ! build expected answer
         if (l .eq. 2) then
            call Y_2m(costheta, true_phys, m)
         else
            call Y_1m(costheta, true_phys, m)
         endif

         physical(:) = true_phys(:) ! need to convert from real to complex

         ! move to spectral
         call LT_ToSpectral_single(physical, spectral, m)

         ! zero out lmax mode
         do i=1,SHTns%nlm
            lval = SHTns_l_value(i)
            if (lval .eq. l_max) then
               mval = SHTns_m_value(i)
               spectral(i) = 0.0d0
            endif
         enddo

         ! check error
         diff = -1.0d0; mxdiff = -1.0d0
         do i=1,SHTns%nlm
            lval = SHTns_l_value(i)
            mval = SHTns_m_value(i)
            if ((lval .eq. l) .and. (mval .eq. m)) then
               diff = max(diff, abs(spectral(i)-1.0d0))
            else
               mxdiff = max(mxdiff, abs(spectral(i)))
            endif
            if (abs(spectral(i)) .gt. 1e-10) then
               write(*,*) i, lval, mval, spectral(i)
            endif
         enddo
         write(*,*) 'Phys-->Spec max error of expected mode  ',diff
         write(*,*) 'Phys-->Spec max error of all other modes',mxdiff

         deallocate(true_phys)

      elseif (ntest .eq. 1) then
         l = 2; m = 1
         lm = SHTns_lm_index(l,m)
         spectral(lm) = (1.0d0, 0.0d0)

         allocate(true_phys(1:SHTns%nlat)) ! build expected answer
         call Y_2m(costheta, true_phys, m)

         ! move to physical space
         call LT_ToPhysical_single(spectral, physical, m)

         !write(*,*) 'expected   SHTns'
         do i=1,n_theta
            write(*,*) true_phys(i), physical(i)
         enddo

         ! check error
         write(*,*) 'Spec-->Phys max error',maxval(abs(real(physical) - true_phys))

         ! move back to spectral space
         spectral(:) = 0.0d0 ! reset
         call LT_ToSpectral_single(physical, spectral, m)

         ! check error
         diff = abs(spectral(lm) - 1.0d0)
         diff1 = maxval(abs(spectral(1:lm-1))) ! blindly assume lm != n_lm
         diff2 = maxval(abs(spectral(lm+1:)))
         mxdiff = max(diff, diff1, diff2)
         write(*,*) 'Spec-->Phys-->Spec max error',mxdiff

         deallocate(true_phys)

      elseif (ntest .eq. 2) then

         nfields = 3

         allocate(Rspectral(1:n_m))
         do mp=1,n_m
            m = m_values(mp)
            allocate(Rspectral(mp)%data(m:l_max,1:nr,1:2,nfields))
         enddo

         nq = 2*nr*nfields
         allocate(Rphysical(1:n_theta,1:nq,1:n_m))

         mind = -1
         do mp = 1, n_m ! modes
            m = m_values(mp)
            if (m .eq. 1) mind = mp
            do l = m, l_max ! \ell values

               Rspectral(mp)%data(l,:,:,:) = 0.0d0
               if ((l .eq. 2) .and. (m .eq. 1)) then
                  do r = 1, nr ! radius
                     Rspectral(mp)%data(l,r,1,1) = 1.0d0*radius(r)
                     Rspectral(mp)%data(l,r,2,2) = 2.0d0*radius(r)**2
                     Rspectral(mp)%data(l,r,2,3) = 5.0d0*radius(r)**3
                  enddo
               endif
            enddo
         enddo

         allocate(true_phys(1:n_theta)) ! build expected answer
         call Y_2m(costheta, true_phys, 1) ! "1" is m value as used in above initialization

         ! move to physical space
         write(*,*) 'S-->P'
         call Legendre_Transform(Rspectral, Rphysical)

         ! check error
         imi = 2
         f = 2
         rlo = 1 + (imi-1)*nr + (f-1)*nr*2  ! index = (r-rlo+1) + (imi-1)*Nr + (f-1)*Nr*2
         rhi = nr + (imi-1)*nr + (f-1)*nr*2
         mxdiff = -1.0d0
         do i=1,n_theta
            diff = maxval(abs(Rphysical(i,rlo:rhi,mind)-true_phys(i)*2*radius(:)**2))
            if (diff .gt. mxdiff) mxdiff = diff
         enddo
         write(*,*) 'Spec-->Phys max error',mxdiff

         ! move back to spectral space
         write(*,*)
         write(*,*) 'P-->S'
         call Legendre_Transform(Rphysical, Rspectral)

         ! zero out l_max mode
         do mp=1,n_m
            Rspectral(mp)%data(l_max,:,:,:) = 0.0d0
         enddo

         ! compute error
         diff = -1.0d0
         mxdiff = -1.0d0
         do mp = 1, n_m
            m = m_values(mp)
            do l = m, l_max
               if ((l .eq. 2) .and. (m .eq. 1)) then
                   diff1 = maxval(abs(1.0d0*radius(:) - Rspectral(mp)%data(l,:,1,1)))
                   diff2 = maxval(abs(2.0d0*radius(:)**2 - Rspectral(mp)%data(l,:,2,2)))
                   diff3 = maxval(abs(5.0d0*radius(:)**3 - Rspectral(mp)%data(l,:,2,3)))
                   diff = max(diff1,diff2,diff3)
               else
                   diff = maxval(abs(Rspectral(mp)%data(l,:,:,:)))
               endif
               if (diff .gt. mxdiff) mxdiff = diff
            enddo
         enddo
         write(*,*) 'Spec-->Phys-->Spec max error',mxdiff

         ! cleanup
         deallocate(true_phys)
         do mp=1,n_m
            deallocate(Rspectral(mp)%data)
         enddo
         deallocate(Rspectral, Rphysical)

      else
         write(*,*) 'Test number not written yet, ntest=',ntest
      endif

      ! cleanup
      deallocate(spectral, physical)

      deallocate(costheta)

   end subroutine test_SHTns
#endif

End Module test_suite

