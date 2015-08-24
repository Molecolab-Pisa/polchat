! MkRij.f90: A Polarisation consistent charge-fitting tool 
!            A Molecolab Tool www.dcci.unipi.it/molecolab/tools
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
Subroutine MkRij(IOut,IPrint,nch,CChg,Rij,Rij3)
  Implicit Real*8 (A-H,O-Z)
  Real*8 :: CChg(3,nch), Rij(4,nch,nch), Rij3(nch,nch)
  1000 format(' Distance between charges:')
  1010 format(2(1x,i6),1x,f12.6)

  Rij = 0.0d0
  do i = 1, nch
    do j = i+1, nch
      dx  = CChg(1,i) - CChg(1,j)
      dy  = CChg(2,i) - CChg(2,j)
      dz  = CChg(3,i) - CChg(3,j)
      dd  = sqrt(dx*dx + dy*dy + dz*dz)
      dd3 = dd**3
      Rij(1,i,j) = dx
      Rij(1,j,i) = -dx
      Rij(2,i,j) = dy
      Rij(2,j,i) = -dy
      Rij(3,i,j) = dz
      Rij(3,j,i) = -dz
      Rij(4,i,j) = dd
      Rij(4,j,i) = dd
      Rij3(i,j)  = dd3
      Rij3(j,i)  = dd3
    enddo
  enddo
  
  if (IPrint.ge.2) then
    write(IOut,1000)
    do i = 1, nch
      do j = i+1, nch
        write(IOut,1010) i,j,Rij(4,i,j)
      enddo
    enddo
  endif
  
  Return
End Subroutine
