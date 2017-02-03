! rms.f90: A Polarisation consistent charge-fitting tool 
!          A Molecolab Tool www.dcci.unipi.it/molecolab/tools
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
!     Compute the RMS between two quantities 
!
REAL*8 Function rms(n,QRef,Qvar)
implicit real*8 (A-H,O-Z)
dimension Qref(n),Qvar(n)
rms = 0.0d0
Do I=1,n
  rms = rms + (Qref(I)-Qvar(I))**2
Enddo
return 
End Function
