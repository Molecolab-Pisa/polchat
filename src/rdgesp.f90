! rdgesp.f90:      A Polarisation consistent charge-fitting tool 
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
subroutine RdGESP

  use constants
  use mmpoldata
  use gespinfo
  use operative
  use time

  implicit real*8 (a-h,o-z)

 1001 format(t46,i9)
 1002 Format(t9,4(d16.8))
 1003 Format(4(d16.8))
 1004 Format(t50,i8)
 1005 Format(4(d16.8))
 1006 Format(3x,d16.8,3x,d16.8,3x,d16.8)
 1007 Format(2x,(A),4(d16.8))
 2000 format(' Reading gesp file.')
 2010 format(' GESP: Number of atoms:      ',i8)
 2020 format(' GESP: Number of gridpoints: ',i8)

  if (iprt .gt. 0) write(iout,2000)
  call starttime
  open(unit=10,file=filenam,status='unknown')

! Skip first 2 lines and read number of atoms
  read(10,*)
  read(10,*)
  read(10,1001) NChg
  allocate (atmnam(NChg))
  allocate (CChg(3,NChg), gesp(NChg))
  allocate (IAnMMP(NChg,LAnMMP),pol(NChg))
  allocate (scrcc(NChg,NChg),scrcp(NChg,NChg))
  allocate (neigh(NChg,NChg))
  allocate (D(3*NChg,3*NChg))

! Read atom coordinates and GESP charges
  if (ldbs.or.lgau) then
    do i = 1, NChg
      read(10,1007) atmnam(i),(CChg(j,i), j=1,3), gesp(i)
    enddo
  else
    do i = 1, NChg
      read(10,1002) (CChg(j,i), j=1,3), gesp(i)
    enddo
  endif

  read(10,*)
  read(10,1006) (DipQM(i),i=1,3)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,1004) NGrd

! Read QM potential and grid coordinates

  allocate (CGrd(3,NGrd),VQM(NGrd))
  do i = 1, NGrd
    read(10,1005) VQM(I), (CGrd(j,i), j=1,3)
  enddo

  close(10)

  call gettime("Reading gesp file")

! Compute distance between charges and gridpoints

  allocate (RChGr(4,NGrd,NChg))
  call dist(NChg,NGrd,.false.,CChg,CGrd,RChGr,RJunk)

  call gettime('Computing chg-grid distance')

! Compute distance between charges

  allocate (Rij(4,NChg,NChg),Rij3(NChg,NChg))
  call dist(NChg,NChg,.true.,CChg,CChg,Rij,Rij3)

  call gettime('Computing chg-chg distance')

! Printout

  if (iprt.ge.1) then
    write(iout,2010) NChg
    write(iout,2020) NGrd
    if (iprt.ge.2) then
      call PrtMat(iout,1,3,DipQM,'Dipole moment (au) from GESP file',.false.)
      call PrtMat(iout,NChg,1,gesp,'GESP charges from GESP file',.false.)
      call PrtMat(iout,3,NChg,CChg,'Atom coordinates (au) from GESP file',.true.)
      call PrtMat(iout,3,NGrd,CGrd,'Gridpoint coordinates (au) from GESP file',.true.)
      call PrtMat(iout,NGrd,1,VQM,'ES potential at gridpoints from GESP file',.false.)
      call PrtMat(iout,NGrd,NChg,RChGr(4,:,:),'Charge-gridpoint distance matrix (au)',.false.)
      call PrtMat(iout,NChg,NChg,Rij(4,:,:),'Charge-charge distance matrix (au)',.false.)
    endif
    call gettime('Debug printout')
  endif

  if (iprt.ge.1) call prttime('reading GESP')

  return
end subroutine


