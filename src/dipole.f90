! dipole.f90       A Polarisation consistent charge-fitting tool 
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
subroutine dipole(nq,q,cq,dip)

  use constants

  implicit real*8(a-h,o-z)
  
  real*8 :: q(nq), cq(3,nq), dip(3)

  dip = zero
  do i = 1, nq
    dip(1) = dip(1) + q(i)*cq(1,i)
    dip(2) = dip(2) + q(i)*cq(2,i)
    dip(3) = dip(3) + q(i)*cq(3,i)
  enddo

  return

end subroutine
