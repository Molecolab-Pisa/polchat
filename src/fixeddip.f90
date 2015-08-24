! FixedDip.f90: A Polarisation consistent charge-fitting tool 
!               A Molecolab Tool www.dcci.unipi.it/molecolab/tools
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
!
! -------------------------------------------------------------------
!
Subroutine FixedDip(IOut,IPrint,nch,CChg,q,Dipole)
  Implicit Real*8 (A-H,O-Z)
  Dimension CChg(3,nch),q(nch),Dipole(3)
  
  Dipole = 0.0d0
  do i = 1, nch
    Dipole(1) = Dipole(1) + q(i)*CChg(1,i)
    Dipole(2) = Dipole(2) + q(i)*CChg(2,i)
    Dipole(3) = Dipole(3) + q(i)*CChg(3,i)
  enddo
  Return
End Subroutine
