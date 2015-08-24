! DoInter.f90: A Polarisation consistent charge-fitting tool 
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
!--------------------------------------------------------------------
! Function to decide if two dipole should interact or not
!
Logical Function DoInter(I,J,IAnMMP,IMMPCn,npol)
Implicit Real*8(A-H,O-Z)
COMMON /OPTIONS/ LAnMMP
DIMENSION :: IAnMMP(LAnMMP,*)
!
! QM/MMPol connectivity check:
! Return whether an interaction between two sites has to be computed.
!
DoInter = .True.
!   
! Check if sites belong to same group or to a corresponding excluded group: 
! Groups option.
! 
If(IMMPCn.eq.2.or.IMMPCn.eq.3) Then
  IinJ = 0
  JinI = 0
  Do 50 K=1,LAnMMP
  ! Check if I-group is in J list
    If (IAnMMP(1,I).eq.IAnMMP(K,J)) IinJ = 1
  ! Check if J-group is in I list
    If (IAnMMP(K,I).eq.IAnMMP(1,J)) JinI = 1
50 Continue       
If (IinJ.ne.JinI) Write(6,*) 'Inconsistency error in MultGroups information.'
If (IinJ.eq.1) DoInter = .False.
!
!     Check if sites are 1-2 or 1-3 neighbours: Amber and Wang option
!
Else if(IMMPCn.eq.1.or.IMMPCn.eq.4) Then
  Do 20 K=1,LAnMMP
  Neigh12 = IAnMMP(K,I)
  If (Neigh12.eq.0) Goto 10
  If (Neigh12.eq.J) Then
    DoInter = .False.
    Goto 10
  Endif
  Do 30 L=1,LAnMMP
    Neigh13 = IAnMMP(L,Neigh12)
    If (Neigh13.eq.0) Goto 20
    If (Neigh13.eq.J) Then
      DoInter = .False.
      Goto 10
    Endif

30  Continue
20  Continue
10  Endif

Return
End Function
