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
       norm = SHT_orthonormal
       If (SHTns_on_the_fly) Then
           layout = SHT_gauss_fly
       Else
           layout = SHT_gauss
       Endif
       ! theta_contiguous does not work?
       ! we only need the scalar transforms
       ! Rayleigh has x in (-1,1) & theta in (pi,0) ---> so south pole is first
       layout = layout + SHT_phi_contiguous + SHT_scalar_only + SHT_south_pole_first

       Call SHTns_verbose(SHTns_information) ! set how much information SHTns will display

       SHTns_nthreads = SHTns_use_threads(n_threads) ! set OpenMP threads

       ! initialize/allocate transforms and build useful arrays
       SHTns_c = SHTns_create(l_max, m_max, m_res, norm)

       ! attach a grid to the SHT object
       Call SHTns_set_grid(SHTns_c, layout, SHTns_polar_threshold, n_theta, n_phi)

       ! map the C SHTns structure to the Fortran one
       Call C_F_pointer(cptr=SHTns_C, fptr=SHTns)
       !Call C_F_pointer(cptr=SHTns%ct, fptr=costheta, shape=[SHTns%nlat])
       !Call C_F_pointer(cptr=SHTns%st, fptr=sintheta, shape=[SHTns%nlat])

       ! apply some m-dependent FFT normalizations during the Legendre Transforms
       Allocate(PTS_normalization(0:l_max), STP_normalization(0:l_max))

       PTS_normalization(0) = 1.0d0/n_phi  ! m=0
       STP_normalization(0) = 1.0d0
       PTS_normalization(1:) = 1.0d0/n_theta ! m/=0
       STP_normalization(1:) = 0.5d0

   End Subroutine SHTns_Initialize

   Subroutine SHTns_Finalize()
       Call SHTns_unset_grid(SHTns_c)
       Call SHTns_destroy(SHTns_c)
       DeAllocate(PTS_normalization, STP_normalization)
   End Subroutine SHTns_Finalize

   Subroutine SHTns_ToSpectral(data_in, data_out)
       Real*8, Intent(In) :: data_in(:,:,:)
       Type(rmcontainer4d), Intent(InOut) :: data_out(1:)

       Complex*16, Allocatable :: temp_in(:,:,:), temp_out(:,:,:,:)
       Integer :: m, f, imi, r, mode
       Integer :: ddims(3), oddims(4)
       Integer :: n_m, nrhs, nfield, rmn, rmx, lme, lms
       Real*8 :: norm

       ! incoming data has shape:
       !   data_out(lm)%data(th,r,r/i,nfields)
       !       lm = mode index
       !       th = theta coordinate (in-processor)
       !        r = radius (distributed)
       !      r/i = real/imag parts
       !     nfields = number of fields
       !
       ! outgoing data has shape:
       !   data_out(:,nrhs,lm)
       !       lm = mode index
       !     nrhs = number of RHS elements
       !     1st index = the theta & real/imag data is stacked

       ! find lower/upper bounds of incoming/outgoing data
       oddims = shape(data_out(1)%data)
       nfield = oddims(4)
       rmn = LBOUND(data_out(1)%data,2)
       rmx = UBOUND(data_out(1)%data,2)

       ddims = shape(data_in)
       n_m = ddims(3)
       nrhs = ddims(2)

       ! repackage the in/output data as complex arrays for SHTns
       ! SHTns input is 1d data of size (nth), output is size(nlm)
       !
       Allocate(temp_in(2,3,4), temp_out(2,3,4,5))

       ! old way:
       ! do m=1,n_m
       !    nl = l_max - m_values(m) + 1
       !    Call DGEMM('T', 'N', nl, nrhs, nth, 1.0, iP_lm(m)%data, nth, data_in(:,:,m),
       !               nth, 0.0, data_out(m)%data, nl)
       ! enddo

       Do mode = 1, n_m ! loop over lm modes

           m = m_values(mode) ! extract actual m value

           norm = PTS_normalization(m) ! extract FFT normalizations, based on m

           lms = SHTns_lmidx(SHTns_c, 0, m)     ! mode index for (l,m)=(0,m)
           lme = SHTns_lmidx(SHTns_c, l_max, m) ! mode index for (l,m)=(lmax,m)

           Do f = 1, nfield ! number of fields
               Do imi = 1, 2 ! real/imag part
                   Do r = rmn, rmx ! radius
                      !Call spat_to_sh_ml(SHTns_c, m, temp_in, temp_out(lms:lme), l_max)
                   Enddo
               Enddo
           Enddo
       Enddo

       ! un-complex the data into the output format

       DeAllocate(temp_in, temp_out)

   End Subroutine SHTns_ToSpectral

   Subroutine SHTns_ToPhysical(data_in, data_out)
       Type(rmcontainer4D), Intent(In) :: data_in(1:)
       Real*8, Intent(InOut) :: data_out(:,:,:)

       Complex*16, Allocatable :: temp_phys(:), temp_spec(:)
       Complex*16 :: ai, ar
       Real*8, Allocatable :: rpart(:), ipart(:)
       Integer :: ddims(3), n_m, nrhs
       Integer :: oddims(4), nfield, rmn, rmx
       Integer :: mode, m, f, r, ind
       Real*8 :: norm

       ai = (0.0d0, 1.0d0) ! imaginary unit
       ar = (1.0d0, 0.0d0) ! regular one

       ddims = shape(data_out)
       n_m = ddims(3)
       nrhs = ddims(2)

       oddims = shape(data_in(1)%data)
       nfield = oddims(4)
       rmn = LBOUND(data_in(1)%data,2)
       rmx = UBOUND(data_in(1)%data,2)

       ! build temporary storage spaces
       allocate(temp_phys(1:n_theta), temp_spec(0:l_max))
       allocate(rpart(1:n_theta), ipart(1:n_theta))
       temp_phys(:) = (0.0d0, 0.0d0)

       Do mode = 1, n_m ! loop over lm modes

           m = m_values(mode) ! extract actual m value

           norm = STP_normalization(m) ! extract FFT normalizations, based on m

           Do f = 1, nfield ! number of fields
               Do r = rmn, rmx ! radius

                  ! package the input data into a larger array of size(0:l_max)
                  temp_spec(:) = (0.0d0, 0.0d0)
                  temp_spec(m:l_max) = ar*data_in(mode)%data(m:l_max,r,1,f) &
                                     + ai*data_in(mode)%data(m:l_max,r,2,f)

                  Call sh_to_spat_ml(SHTns_c, m, temp_spec, temp_phys, l_max)

                  rpart(:) = norm*temp_phys(1:n_theta)
                  ipart(:) = norm*temp_phys(1:n_theta)

                  ! store the output somewhere
              !    ind = something that stripes the data?
              !    data_out(1:n_theta,ind,mode) = rpart or ipart

               Enddo
           Enddo
       Enddo
       deallocate(temp_phys, temp_spec, rpart, ipart)

   End Subroutine SHTns_ToPhysical

   ! Legrendre Transform, i.e., no FFT at a given m value: Physical-->Spectral
   subroutine LT_ToSpectral(data_in, data_out, m_val)
      integer, intent(in) :: m_val
      complex*16, intent(inout) :: data_in(1:n_theta)
      complex*16, intent(inout) :: data_out(1:SHTns%nlm)

      integer :: lmstart, lmstop

      lmstart = SHTns_lmidx(SHTns_c, 0, m_val)
      lmstop  = SHTns_lmidx(SHTns_c, l_max, m_val)
      call spat_to_sh_ml(shtns_c, m_val, data_in, data_out(lmstart:lmstop), l_max)

   end subroutine LT_ToSpectral

    ! Legendre transform, i.e., no FFT at a given m value: Spectral-->Physical
   subroutine LT_ToPhysical(data_in, data_out, m_val)
      integer, intent(in) :: m_val
      complex*16, intent(inout) :: data_in(1:SHTns%nlm)
      complex*16, intent(inout) :: data_out(1:n_theta)

      integer :: lmstart, lmstop

      lmstart = SHTns_lmidx(SHTns_c, 0, m_val)
      lmstop  = SHTns_lmidx(SHTns_c, l_max, m_val)
      call sh_to_spat_ml(shtns_c, m_val, data_in(lmstart:lmstop), data_out, l_max)

   end subroutine LT_ToPhysical

End Module Legendre_Transforms_SHTns

