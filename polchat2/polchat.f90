! polchat.f90: A Polarisation consistent charge-fitting tool 
!              A Molecolab Tool www.molecolab.dcci.unipi.it/tools
! version 2.0
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
program polchat

  use constants
  use mmpoldata

! Print header
  call PrtHdr

! Parse input
  call RdOpts

! Read gaussian ESP file, produced using keywords POP=ESP IOP(6/50=1)
  call RdGESP

! Read connectivity from mol2 file
  call RdConn

! Read polarisability file
  call RdPol

! Compute screenings
  call MkScr

! Compute and invert MMPol matrix
  call MMPMat

! Read constraints
  call RdCns

! Print data
  call PrtDat

! Form constraint matrices
  call DoCns

! Compute ESP charges
  call ESP

! Form constraint matrices again
  call DoCns

! Compute pol-ESP charges
  call PESP

! Printout
  call Printout

end program
