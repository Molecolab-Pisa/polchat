! potchg.f90:      A Polarisation consistent charge-fitting tool 
!                  A Molecolab Tool www.molecolab.dcci.unipi.it/tools
!
! Copyright (C) 2014, 2015, 2016, 2017
!   S. Caprasecca, C. Curutchet, S. Jurinovich, B. Mennucci
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
subroutine potchg(ng,nq,q,R,V)

  use constants

  implicit real*8(a-h,o-z)

  real*8 :: q(nq), R(ng,nq), V(ng)

  do i = 1, ng
    V(i) = zero
    do j = 1, nq
      V(i) = V(i) + q(j)/R(i,j)
    enddo
  enddo

  if (iprt.ge.2) call PrtMat(iout,ng,1,V,'ESP potential',.false.)

  return

end subroutine

