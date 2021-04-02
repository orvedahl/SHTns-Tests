!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!

Module Legendre_Transforms_SHTns
    Use Legendre_Polynomials
    Use Structures
    Use iso_c_binding

    Implicit None

    ! include the SHTns interfaces
    Include 'shtns.f03'

    Type(c_ptr) :: SHTns_c             ! main SHTns structure built on initialization
    Type(SHTns_info), Pointer :: SHTns ! the Fortran friendly main SHTns structure

    Integer :: SHTns_nthreads          ! keep track of how many threads SHTns is using

    Real*8, Allocatable :: PTS_normalization(:), STP_normalization(:)

    Interface Legendre_Transform
        Module Procedure SHTns_ToSpectral, SHTns_ToPhysical
    End Interface

Contains

   ! convert (l,m) value to mode index
   function SHTns_lm_index(l, m) result(lm_idx)
      integer, intent(in) :: l, m
      integer :: lm_idx
      lm_idx = SHTns_lmidx(SHTns_c,l,m)
   end function SHTns_lm_index

   ! convert mode index to l value of (l,m) pair
   function SHTns_l_value(lm_idx) result(l)
      integer, intent(in) :: lm_idx
      integer :: l
      l = SHTns_lm2l(SHTns_c,lm_idx)
   end function SHTns_l_value

   ! convert mode index to m value of (l,m) pair
   function SHTns_m_value(lm_idx) result(m)
      integer, intent(in) :: lm_idx
      integer :: m
      m = SHTns_lm2m(SHTns_c,lm_idx)
   end function SHTns_m_value

   ! get quadrature weights
   subroutine SHTns_gauss_weights(weights)
      real*8, intent(inout) :: weights(1:n_theta)
      call SHTns_gauss_wts(SHTns_c, weights)
   end subroutine SHTns_gauss_weights

   Subroutine SHTns_Initialize(n_threads, SHTns_on_the_fly, &
                            SHTns_information, SHTns_polar_threshold, &
                            n_phi, m_max)
       Integer, intent(in) :: n_threads, SHTns_information, n_phi, m_max
       Logical, intent(in) :: SHTns_on_the_fly
       Real*8, intent(in) :: SHTns_polar_threshold

       Integer :: norm, layout, m_res

       m_res = 1 ! 2*pi/m_res is the azimuthal periodicity, m_max*m_res is max azimuthal order

       ! choose grid, data layout, and Y_l^m normalization --- be consistent with Rayleigh
       If (SHTns_on_the_fly) Then
           layout = SHT_gauss_fly
       Else
           layout = SHT_gauss
       Endif
       ! theta_contiguous does not work?
       ! we only need the scalar transforms
       ! Rayleigh has x in (-1,1) & theta in (pi,0) ---> so south pole is first
       layout = layout + SHT_phi_contiguous + SHT_scalar_only + SHT_south_pole_first
       !layout = layout + SHT_theta_contiguous + SHT_scalar_only + SHT_south_pole_first

       ! Rayleigh uses the very sane choice of orthonormal Y_l^m
       norm = SHT_orthonormal

       Call SHTns_verbose(SHTns_information) ! set how much information SHTns will display

       SHTns_nthreads = SHTns_use_threads(n_threads) ! set OpenMP threads
       if (SHTns_nthreads .gt. 1) then
          write(*,*) 'SHTns is using OpenMP threads:',SHTns_nthreads
       endif

       ! initialize/allocate transforms and build useful arrays
       SHTns_c = SHTns_create(l_max, m_max, m_res, norm)

       ! attach a grid to the SHT object & determine optimal algorithm
       Call SHTns_set_grid(SHTns_c, layout, SHTns_polar_threshold, n_theta, n_phi)

       ! map the C SHTns structure to the Fortran one
       Call C_F_pointer(cptr=SHTns_C, fptr=SHTns)
       !Call C_F_pointer(cptr=SHTns%ct, fptr=costheta, shape=[SHTns%nlat])
       !Call C_F_pointer(cptr=SHTns%st, fptr=sintheta, shape=[SHTns%nlat])

       ! apply some m-dependent FFT normalizations during the Legendre Transforms
       Allocate(PTS_normalization(0:l_max), STP_normalization(0:l_max))

       !PTS_normalization(0) = 1.0d0/n_phi    ! m=0
       !STP_normalization(0) = 1.0d0
       !PTS_normalization(1:) = 1.0d0/n_theta ! m/=0
       !STP_normalization(1:) = 0.5d0
       PTS_normalization(:) = 1.0d0
       STP_normalization(:) = 1.0d0

   End Subroutine SHTns_Initialize

   Subroutine SHTns_Finalize()
       Call SHTns_unset_grid(SHTns_c)
       Call SHTns_destroy(SHTns_c)
       DeAllocate(PTS_normalization, STP_normalization)
   End Subroutine SHTns_Finalize

   Subroutine SHTns_ToSpectral(data_in, data_out)
       Real*8, Intent(In) :: data_in(:,:,:)
       Type(rmcontainer4d), Intent(InOut) :: data_out(1:)
       ! ingoing data has shape:
       !   data_out(th,nrhs,lm)
       !       th = theta (in-processor)
       !     nrhs = number of RHS elements, stacked over: radius/real/imag/nfields
       !       lm = mode index (distributed)
       ! nrhs axis is 1-based indexing, except for radius, ordered: radius-real-imag-field
       !     index = (r-rlo+1) + (imi-1)*Nr + (f-1)*Nr*2
       !
       ! outgoing data has shape:
       !   data_out(lm)%data(l,r,imi,nf)
       !       lm = mode index (distributed)
       !        l = spherical harmonic (in-processor)
       !        r = radius (distributed)
       !      imi = real/imag parts
       !       nf = number of fields
       Complex*16, Allocatable :: temp_spec(:), temp_phys(:)
       Complex*16 :: ai, ar
       Integer :: m, f, r, mode
       Integer :: ddims(3), oddims(4)
       Integer :: n_m, nrhs, nfield, rmn, rmx, my_Nr, indr, indi, lind
       Real*8 :: norm

       ai = (0.0d0, 1.0d0) ! imaginary unit
       ar = (1.0d0, 0.0d0) ! regular one

       ! find lower/upper bounds of incoming/outgoing data
       oddims = shape(data_out(1)%data)
       nfield = oddims(4)
       rmn = LBOUND(data_out(1)%data,2)
       rmx = UBOUND(data_out(1)%data,2)

       my_Nr = rmx - rmn + 1

       ddims = shape(data_in)
       n_m = ddims(3)
       nrhs = ddims(2)

       Allocate(temp_spec(1:l_max+1), temp_phys(1:n_theta)) ! SHTns expects complex arrays

       Do mode = 1, n_m ! loop over lm modes

           m = m_values(mode) ! extract actual m value

           norm = PTS_normalization(m) ! extract FFT normalizations, based on m

           lind = l_max - m + 1 ! number of modes for this m-value

           Do f = 1, nfield ! number of fields
               Do r = rmn, rmx ! radius

                   ! package incoming data into complex array at this radius/field
                   indr = (r-rmn+1) + (f-1)*my_Nr*2
                   indi = (r-rmn+1) + my_Nr + (f-1)*my_Nr*2
                   temp_phys(:) = ar*data_in(:,indr,mode) + ai*data_in(:,indi,mode)

                   temp_spec(:) = 0.0d0 ! spectral output is of size (lmax-m+1)
                   Call spat_to_sh_ml(SHTns_c, m, temp_phys, temp_spec(1:lind), l_max)

                   ! extract results and store in output array
                   data_out(mode)%data(m:l_max,r,1,f) = norm*real(temp_spec(1:lind))
                   data_out(mode)%data(m:l_max,r,2,f) = norm*aimag(temp_spec(1:lind))

               Enddo
           Enddo
       Enddo

       DeAllocate(temp_spec, temp_phys)

   End Subroutine SHTns_ToSpectral

   Subroutine SHTns_ToPhysical(data_in, data_out)
       Type(rmcontainer4D), Intent(In) :: data_in(1:)
       Real*8, Intent(InOut) :: data_out(:,:,:)
       ! ingoing data has shape:
       !   data_out(lm)%data(l,r,imi,nf)
       !       lm = mode index (distributed)
       !        l = spherical harmonic (in-processor)
       !        r = radius (distributed)
       !      imi = real/imag parts
       !       nf = number of fields
       !
       ! outgoing data has shape:
       !   data_out(th,nrhs,lm)
       !       th = theta (in-processor)
       !     nrhs = number of RHS elements, stacked over: radius/real/imag/nfields
       !       lm = mode index (distributed)
       Complex*16, Allocatable :: temp_phys(:), temp_spec(:)
       Complex*16 :: ai, ar
       Integer :: odims(3), nm, nrhs, my_Nr
       Integer :: idims(4), nfield, rmn, rmx, mode, m, f, r, ind
       Real*8 :: norm

       ai = (0.0d0, 1.0d0) ! imaginary unit
       ar = (1.0d0, 0.0d0) ! regular one

       odims = shape(data_out)
       nm = odims(3)
       nrhs = odims(2)

       idims = shape(data_in(1)%data)
       nfield = idims(4)
       rmn = LBOUND(data_in(1)%data,2)
       rmx = UBOUND(data_in(1)%data,2)

       my_Nr = rmx - rmn + 1

       ! build temporary storage spaces
       allocate(temp_phys(1:n_theta), temp_spec(0:l_max))

       Do mode = 1, nm ! loop over lm modes

           m = m_values(mode) ! extract actual m value

           norm = STP_normalization(m) ! extract FFT normalizations, based on m

           Do f = 1, nfield ! number of fields
               Do r = rmn, rmx ! radius

                  ! package the input data into a complex array of size lmax-m+1
                  temp_spec(:) = (0.0d0, 0.0d0)
                  temp_spec(m:l_max) = ar*data_in(mode)%data(m:l_max,r,1,f) &
                                     + ai*data_in(mode)%data(m:l_max,r,2,f)

                  temp_phys(:) = (0.0d0, 0.0d0)
                  Call sh_to_spat_ml(SHTns_c, m, temp_spec(m:l_max), temp_phys, l_max)

                  ! 1-based indexing, except for radius, ordered: radius-real-imag-field
                  !     index = (r-rlo+1) + (imi-1)*Nr + (f-1)*Nr*2
                  ind = (r-rmn+1) + (f-1)*my_Nr*2         ! real part for this r/field
                  data_out(1:n_theta,ind,mode) = norm*real(temp_phys(1:n_theta))

                  ind = (r-rmn+1) + my_Nr + (f-1)*my_Nr*2 ! imag part for this r/field
                  data_out(1:n_theta,ind,mode) = norm*aimag(temp_phys(1:n_theta))

               Enddo
           Enddo
       Enddo
       deallocate(temp_phys, temp_spec)

   End Subroutine SHTns_ToPhysical

   ! Legrendre Transform, i.e., no FFT at a given m value: Physical-->Spectral
   subroutine LT_ToSpectral_single(data_in, data_out, m_val)
      integer, intent(in) :: m_val
      complex*16, intent(inout) :: data_in(1:n_theta)
      complex*16, intent(inout) :: data_out(1:SHTns%nlm)

      integer :: lmstart, lmstop

      lmstart = SHTns_lmidx(SHTns_c, m_val, m_val)
      lmstop  = SHTns_lmidx(SHTns_c, l_max, m_val)
      call spat_to_sh_ml(shtns_c, m_val, data_in, data_out(lmstart:lmstop), l_max)

   end subroutine LT_ToSpectral_single

    ! Legendre transform, i.e., no FFT at a given m value: Spectral-->Physical
   subroutine LT_ToPhysical_single(data_in, data_out, m_val)
      integer, intent(in) :: m_val
      complex*16, intent(inout) :: data_in(1:SHTns%nlm)
      complex*16, intent(inout) :: data_out(1:n_theta)

      integer :: lmstart, lmstop

      lmstart = SHTns_lmidx(SHTns_c, m_val, m_val)
      lmstop  = SHTns_lmidx(SHTns_c, l_max, m_val)
      call sh_to_spat_ml(shtns_c, m_val, data_in(lmstart:lmstop), data_out, l_max)

   end subroutine LT_ToPhysical_single

End Module Legendre_Transforms_SHTns

