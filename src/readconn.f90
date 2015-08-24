! readconn.f90: A Polarisation consistent charge-fitting tool 
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
Subroutine readconn(IPrint,IOut,fileconn,nch,IAnMMP,LAnMMP)
IMPLICIT REAL*8 (A-H,O-Z)
INTEGER :: nbond,nch,natoms
INTEGER,ALLOCATABLE :: ib1(:),ib2(:)
INTEGER :: IAnMMP(nch,LAnMMP)
CHARACTER(LEN=50) :: fileconn,what
1101 Format(I5,I6)
1102 Format('ERROR: Number of atom in mol2 is ',I5,' expected number is ',I6)
1103 Format(6X,2(I5))
1104 Format(I5,'|',9(I5))
OPEN(unit=11,file=fileconn,status='unknown')
!
! Skip the first two lines and read the number of atoms and number of 
! bonds. Check if the number of atoms in the mol2 file is the same of
! the number of charge to fit.
!
READ(11,*)
READ(11,*)
READ(11,*) natoms,nbond
If (natoms.ne.nch ) Then
  WRITE(IOut,1102) natoms,nch
  STOP
EndIf
ALLOCATE(ib1(nbond),ib2(nbond))
IAnMMP = 0
!
! Skip the ATOM Section and goes to the connectivity information.
! Read the connectivity and store the information in ib1 and ib2.
!
DO I=1,nch+10
  READ(11,*) what
  If (what.eq.'@<TRIPOS>BOND') GOTO 10
ENDDO
10 continue ! WRITE(IOut,*) 'Found',what
DO I=1,nbond
  READ(11,*) IXXX, ib1(I), ib2(I)
ENDDO
CLOSE(11)
!
! Build the connectivity in the standard MMPol Format and add it
! to the IAnMMP matrix.
!
DO N=1,nch
  K=1
  DO I=1,nbond
    If (N.eq.ib1(I))  Then 
      IAnMMP(N,K) = ib2(I)
!     IAnMMP(K,N) = ib2(I)
      K=K+1
    EndIf
    If (N.eq.ib2(I))  Then
      IAnMMP(N,K) = ib1(I)
!     IAnMMP(K,N) = ib1(I)
      K=K+1
    EndIf
    If (K.gt.LAnMMP) Then
      WRITE(IOut,*) 'ERROR: Exceed in IAnMMP dimension!'
      STOP
    EndIf
  ENDDO
ENDDO
DEALLOCATE(ib1,ib2)

if (IPrint.ge.2) then
 DO I=1,nch
   WRITE(IOut,1104) I,(IAnMMP(I,J),J=1,LAnMMP)
 ENDDO
endif

RETURN
End Subroutine
