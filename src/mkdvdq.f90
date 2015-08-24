! MkDvDq.f90: A Polarisation consistent charge-fitting tool 
!             A Molecolab Tool www.dcci.unipi.it/molecolab/tools
!
! Copyright (C) 2014, 2015 
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
Subroutine MkDvDq(IOut,IPrint,nch,Rij,DvDq)
  Implicit Real*8 (A-H,O-Z)
  Real*8 :: Rij(4,nch,nch), DvDq(nch,nch)
  1000 format(' dV/dq:')
  1010 format(2(1x,i6),1x,f12.6)

  DvDq = 0.0d0
  do i = 1, nch
    do j = i+1, nch
      r = Rij(4,i,j)
      DvDq(i,j) = 1.0d0/r
      DvDq(j,i) = 1.0d0/r
    enddo
  enddo

  if (IPrint.ge.2) then
    write(IOut,1000)
    do i = 1, nch
      do j = i+1, nch
        write(IOut,1010) i,j,DvDq(i,j)
      enddo
    enddo
  endif

  Return
End Subroutine
