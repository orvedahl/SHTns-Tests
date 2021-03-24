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

Module Legendre_Transforms
    Use Legendre_Polynomials
    Use Structures
    Implicit None
    Interface Legendre_Transform
        Module Procedure PtS_4d_dgpv3, StP_4d_dgp2
    End Interface
Contains

!/////////////////////////////////////////////////////////////////////////////
! These 4-D routines use a higher dimension of temp array
Subroutine StP_4d_dgp2(data_in, data_out)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer4d), Intent(In) :: data_in(:)
    Real*8, Intent(InOut) :: data_out(:,:,:)
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),temp2(:,:), temp3(:,:,:,:)
    Integer :: m,nl,nt1,i,j,l, nt2, ddims(3), nrhs, jstart, jend,r,lstart
    Integer :: nfield, rmn, rmx, oddims(4),imi,f,iend,istart

    ddims = shape(data_out)
    n_m = ddims(3)
    nrhs = ddims(2)

    oddims = shape(data_in(1)%data)
    nfield = oddims(4)

    rmn = LBOUND(data_in(1)%data,2)
    rmx = UBOUND(data_in(1)%data,2)



    alpha = 1.0d0
    beta = 0.0d0
    if (.not. parity) then
        Do m = 1, n_m
            nl = l_max-m_values(m)+1
            CALL DGEMM('T','N',n_theta,nrhs,nl, alpha, p_lm(m)%data,  &
                nl,data_in(m)%data , nl, beta,data_out,n_theta)
        Enddo

    else
    !////////////////////////////////////
    ! In progress
    nt1 = n_theta+1
    nt2 = n_theta/2
    data_out(:,:,:) = 0.0d0
    Allocate(temp(1:nt2,1:nrhs))
    ! Solve for odd and even functions
    Do m = 1, n_m

        If (n_l_even(m) .gt. 0) then

            jend = n_l_even(m)
            !write(6,*)jend,rmn,rmx,nfield
            Allocate(temp3(1:jend,rmn:rmx,1:2,1:nfield))
            Do f = 1, nfield
            Do imi = 1, 2
            Do r = rmn, rmx
            Do j = 1,jend
                l =  lvalsi(m)%even(j)
                temp3(j,r,imi,f) = data_in(m)%data(l,r,imi,f)
            Enddo
            Enddo
            Enddo
            Enddo

            CALL DGEMM('T','N',nt2,nrhs,n_l_even(m), alpha, p_lm_even(m)%data, n_l_even(m),temp3 , n_l_even(m), beta,temp,nt2)

            data_out(1:nt2,1:nrhs,m) = temp(1:nt2,1:nrhs)    ! store symmetric part in data_out

            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(nt1-j,i,m) = temp(j,i)    ! reflect even modes about equator
                Enddo
            Enddo
            DeAllocate(temp3)
        Endif

        If (n_l_odd(m) .gt. 0) then

            jend = n_l_odd(m)

            Allocate(temp3(1:jend,rmn:rmx,1:2,1:nfield))
            Do f = 1, nfield
            Do imi = 1, 2
            Do r = rmn, rmx
            Do j = 1,jend
                l =  lvalsi(m)%odd(j)
                temp3(j,r,imi,f) = data_in(m)%data(l,r,imi,f)
            Enddo
            Enddo
            Enddo
            Enddo



            CALL DGEMM('T','N',nt2,nrhs,n_l_odd(m), alpha, p_lm_odd(m)%data, n_l_odd(m),temp3 , n_l_odd(m), beta,temp,nt2)
            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(j,i,m) = data_out(j,i,m)+temp(j,i)
                    data_out(nt1-j,i,m) = data_out(nt1-j,i,m)-temp(j,i)    ! antisymmetric about equator
                Enddo
            Enddo
            DeAllocate(temp3)
        Endif
    Enddo


    ! Note - not sure if it's faster to make a variable named nt2j1 = nt2-j+1 or just let it compute on the fly
    DeAllocate(temp)
    Endif

End Subroutine StP_4d_dgp2


Subroutine PtS_4d_dgpv3(data_in, data_out)
    Implicit None
    Type(rmcontainer4D), Intent(InOut) :: data_out(1:)
    Real*8, Intent(In) :: data_in(:,:,:)
    Integer  :: nrhs
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:,:,:),fodd(:,:,:), feven(:,:,:)
    Integer :: m,nl,nt1,i,j,l, nt2,ddims(3),k
    Integer :: oddims(4), nfield
    Integer :: rmn, rmx, f, imi, istart, iend, jend,r

    oddims = shape(data_out(1)%data)
    nfield = oddims(4)
    rmn = LBOUND(data_out(1)%data,2)
    rmx = UBOUND(data_out(1)%data,2)


    ddims = shape(data_in)
    n_m = ddims(3)
    nrhs = ddims(2)
    alpha = 1.0d0
    beta = 0.0d0
    if (.not. parity) then
        Do m =1, n_m
            nl = l_max-m_values(m)+1

            CALL DGEMM('T','N',nl,nrhs,n_theta, alpha, ip_lm(m)%data, &
                n_theta,data_in(:,:,m) , n_theta, beta,data_out(m)%data,nl)
        Enddo
    else
    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs,1:n_m))
    Allocate(fodd(1:n_theta/2,1:nrhs,1:n_m))
    feven(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    fodd(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    nt1 = n_theta+1
    nt2 = n_theta/2
    Do k = 1, n_m
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i,k) = feven(j,i,k)+data_in(nt1-j,i,k)
             fodd(j,i,k) =  fodd(j,i,k)-data_in(nt1-j,i,k)
        Enddo
    enddo
    Enddo



    Do m = 1, n_m

            If (n_l_even(m) .gt. 0) then
                Allocate(temp(1:n_l_even(m),rmn:rmx,1:2,1:nfield))
                CALL DGEMM('T','N',n_l_even(m),nrhs,nt2, alpha, ip_lm_even(m)%data, nt2,feven(:,:,m) , nt2, beta,temp,n_l_even(m))


                jend = n_l_even(m)
                Do f = 1, nfield
                Do imi = 1, 2
                Do r = rmn, rmx
                Do j = 1,jend
                    l =  lvalsi(m)%even(j)
                    data_out(m)%data(l,r,imi,f)=temp(j,r,imi,f)
                Enddo
                Enddo
                Enddo
                Enddo

                DeAllocate(temp)
            Endif

            If (n_l_odd(m) .gt. 0) then
                Allocate(temp(1:n_l_odd(m),rmn:rmx,1:2,1:nfield))
                CALL DGEMM('T','N',n_l_odd(m),nrhs,nt2, alpha, ip_lm_odd(m)%data, nt2,fodd(:,:,m) , nt2, beta,temp,n_l_odd(m))


                jend = n_l_odd(m)
                Do f = 1, nfield
                Do imi = 1, 2
                Do r = rmn, rmx
                Do j = 1,jend
                    l =  lvalsi(m)%odd(j)
                    data_out(m)%data(l,r,imi,f)=temp(j,r,imi,f)
                Enddo
                Enddo
                Enddo
                Enddo


                DeAllocate(temp)
            Endif


    Enddo

    DeAllocate(feven)
    DeAllocate(fodd)
    endif
End Subroutine PtS_4d_dgpv3

End Module Legendre_Transforms
