! docns.f90:       A Polarisation consistent charge-fitting tool 
!                  A Molecolab Tool www.molecolab.dcci.unipi.it/tools
!
! Copyright (C) 2014, 2015, 2016, 2017
!   S. Caprasecca, C. Curutchet, B. Mennucci
!
! This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
! A copy of the GNU General Public License can be found in LICENSE or at
!   <http://www.gnu.org/licenses/>.
!
subroutine DoCns

  use constants
  use mmpoldata
  use strings
  use constraints
  use time

  implicit real(a-h,o-z)

  call starttime

! Allocate X and B vectors

  NDim = NChg+NCons
  allocate(X(NDim,NDim), B(NDim))

  X  = zero
  B  = zero
  II = 0

! Read again and fill vectors

! 1. CHG

  II = II + 1 
  X(NChg+II,1:NChg) = one
  X(1:NChg,NChg+II) = one
  B(NChg+II) = RCChg

! 2. FRG
 
  do i = 1, NCFrg
    II = II + 1
    do j = 1, ICFrg(i)
      X(NChg+II,VCFrg(i,j)) = one   
      X(VCFrg(i,j),NChg+II) = one
      B(NChg+II) = RCFrg(i)
    enddo
  enddo

! 3. EQV

  do i = 1, NCEqv
    do j = 1, ICEqv(i)-1
      II = II + 1
      X(NChg+II,VCEqv(i,j))   = -one
      X(NChg+II,VCEqv(i,j+1)) =  one
      X(VCEqv(i,j),NChg+II)   = -one
      X(VCEqv(i,j+1),NChg+II) =  one
    enddo
  enddo

! 4. RES does not add constraints

  call gettime('')
  if (iprt.ge.1) call prttime('building constraint matrices')

  return
 
end subroutine
