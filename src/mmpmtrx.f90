! MMPMTRX.f90: A Polarisation consistent charge-fitting tool 
!              A Molecolab Tool www.dcci.unipi.it/molecolab/tools
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
! This subroutine compute, invert and store the MMPol matrix (A^-1) 
!
Subroutine MMPMTRX(IPrint,IOut,npol,pol,CPol,D,IAnMMP,IMMPCn,LAnMMP)
!
! IMMPCn ..... Treatment of the polarization
! IAnMMP ..... Matrix containing the connectivity information
! LAnMMP ..... Dimension of IAnMMP
! pol ........ isotropic polarizabilities (a.u.)
! D .......... MMPol Matrix (a.u.)
!

IMPLICIT REAL*8(A-H,O-Z)
PARAMETER(LWORK=30000)
LOGICAL   :: DoInter,DoThole,DoCalc
DIMENSION :: pol(npol),CPol(3,npol),WORK(LWORK),IPIV(npol*3)
DIMENSION :: D(3*npol,3*npol),IAnMMP(npol,LAnMMP)
LDA=npol*3

!
! Clear the whole MMPol matrix
!
D = 0.0d0

!
! Fill in diagonal elements of MMPMTRX
!
DO I=1,npol
  D(I,I)               = 1.0d0/pol(I)     ! x-block
  D(I+npol,I+npol)     = 1.0d0/pol(I)     ! y-block
  D(I+2*npol,I+2*npol) = 1.0d0/pol(I)     ! z-block
ENDDO

!
! Fill in non-diagonal elements of MMPTRX
!
DO I=1,npol
  DO J=I+1,npol
    DoCalc = .False.

    ! Thole smeared dipole interaction tensor between all dipoles (Thole option)
    If (IMMPCn.eq.3) Then
      DoThole = .True.
      DoCalc  = .True.
      IScreen = 1

    ! Exclude Amber 1-2 and 1-3 interactions (Amber option) or interactions
    ! inside sites from same group (Groups option).
    ElseIf (IMMPCn.eq.0 .or. IMMPCn.eq.1 .or. IMMPCn.eq.2) then
      DoThole = .False.
      IScreen = 1
      if (DoInter(I,J,IAnMMP,IMMPCn,npol)) DoCalc=.True.

    ! Exclude Amber 1-2 and 1-3 interactions (Amber option) and use Thole
    ! smeared dipole interaction tensor for the others with Wang screening
    ! parameter
    ElseIf (IMMPCn.eq.4) Then
      DoThole = .True.
      IScreen = 2
      If (DoInter(I,J,IAnMMP,IMMPCn,npol)) DoCalc=.True.

    ! DL Amber model (compute all interactions with Thole linear screening)
    ElseIf (IMMPCn.eq.5) Then
      DoThole = .True.
      DoCalc=.True.
      IScreen = 3

    Else
      Write(6,*) 'Confused in polarization treatment.'
      Stop
    EndIf

    If (DoCalc) then
      D(I,J) = Tensormm(1,1,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J,I) = D(I,J)
      D(I,J+npol) = Tensormm(1,2,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J,I+npol) = D(I,J+npol)
      D(I+npol,J) = D(I,J+npol)
      D(J+npol,I) = D(I,J+npol)
      D(I,J+2*npol) = Tensormm(1,3,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J,I+2*npol) = D(I,J+2*npol)
      D(I+2*npol,J) = D(I,J+2*npol)
      D(J+2*npol,I) = D(I,J+2*npol)
      D(I+npol,J+npol) = Tensormm(2,2,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J+npol,I+npol) = D(I+npol,J+npol)
      D(I+npol,J+npol*2) = Tensormm(2,3,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J+npol,I+npol*2) = D(I+npol,J+npol*2)
      D(I+npol*2,J+npol) = D(I+npol,J+npol*2)
      D(J+npol*2,I+npol) = D(I+npol,J+npol*2)
      D(I+npol*2,J+npol*2) = Tensormm(3,3,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J+npol*2,I+npol*2) = D(I+npol*2,J+npol*2)
    EndIf

  ENDDO
ENDDO
 
!
!     Invert the MMPol matrix using lapac library
!
999   Format(6(1x,D15.5))      
If (IPrint.ge.2) then
  write(6,*) 'Matrice non invertita'
  Do i=1,npol*3
    Write(IOut,999) (D(i,j),j=1,npol*3)
  Enddo
EndIf

M = npol*3
N = npol*3

call DGETRF (M,N,D,LDA,IPIV,INFO)
call DGETRI (N,D,LDA,IPIV,WORK,LWORK,INFO)

!
!     Save the inverted MMPol Matrix
!
If (IPrint.ge.2) then
  write(6,*) 'Matrice invertita'
  Do i=1,npol*3
    Write(IOut,999) (D(i,j),j=1,npol*3)
  Enddo
EndIf

End Subroutine
