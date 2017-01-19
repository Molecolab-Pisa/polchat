! MkScr.f90: A Polarisation consistent charge-fitting tool 
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
!
! -------------------------------------------------------------------
!
Subroutine MkScr(IOut,IPrint,NPol,CPol,Pol,ScrChCh,ScrChPl,IAnMMP,IMMPCn,LAnMMP, &
  LScrChPl)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION :: ScrChCh(npol,npol), ScrChPl(npol,npol), IAnMMP(npol,LAnMMP)
  DIMENSION :: CPol(3,*), R(3), Pol(NPol)
  INTEGER, ALLOCATABLE :: Neigh(:,:)
  LOGICAL :: DoThole, LScrChPl, LN1213
  PARAMETER(a=1.7278d0,b=2.5874d0,c=2.0580d0)
!
! > DoThole .... True if you want to use Thole
! > IScreen .... 1 = use Amber Screening Factor (a)
!           .... 2 = use AL (Wang) Screening Factor (b)
!           .... 3 = use DL (Wang) Screening Factor (c) [NYI]
! Detect intramolecular polarization treatment
    IF (LScrChPl) THEN
      IF (IMMPCn.EQ.3) THEN
        DoThole = .TRUE.
        IScreen = 1
      ELSEIF (IMMPCn.EQ.1) THEN
        DoThole = .FALSE.
        IScreen = 1
      ELSEIF (IMMPCn.EQ.4) THEN
        DoThole = .TRUE.
        IScreen = 2
      ELSEIF (IMMPCn.EQ.5) THEN
        DoThole = .TRUE.
        IScreen = 3
      ELSE
        WRITE(6,*) 'Confused in polarization treatment.'
        STOP
      ENDIF
    ENDIF

  Allocate (Neigh(npol,npol))
  Neigh = 0

! Neigh integer matrix:
!  N(i,i) = 1
!  N(i,j) = 2 if i and j are 1-2 neighbours
!  N(i,j) = 3 if i and j are 1-3 neighbours
!  N(i,j) = 4 if i and j are 1-4 neighbours
!  N(i,j) = 0 if i and j are not within 4 bonds

  do i = 1, npol
    Neigh(i,i) = 1
    do j = 1, LAnMMP
      if (IAnMMP(i,j).eq.0) then
        goto 700
      else
        Neigh(IAnMMP(i,j),i) = 2
        Neigh(i,IAnMMP(i,j)) = 2
      endif
    enddo
 700 continue
  enddo

  do i = 1, npol
    do j = 1, npol
      if (Neigh(i,j).eq.2) then
        do k = 1, npol
          if (Neigh(j,k).eq.2 .and. Neigh(i,k).ne.1 .and. Neigh(i,k).ne.2) then
            Neigh(i,k) = 3
            Neigh(k,i) = 3
          endif
        enddo
      endif
    enddo
  enddo

  do i = 1, npol
    do j = 1, npol
      if (Neigh(i,j).eq.3) then
        do k = 1, npol
          if (Neigh(j,k).eq.2 .and. Neigh(i,k).ne.1 .and. Neigh(i,k).ne.2 .and. Neigh(i,k).ne.3) then
            Neigh(i,k) = 4
            Neigh(k,i) = 4
          endif
        enddo
      endif
    enddo
  enddo

  if (IPrint.ge.2) then
    write(IOut,1000)
    do i = 1, npol
      write(IOut,2000) (Neigh(i,j),j=1,npol)
    enddo
  endif
 1000 format(' Neighbours matrix:')
 2000 format(20(1x,i4))

! Matrix Neigh is complete: make Screening matrices 
!   This is only implemented for Amber12 for the time being
 ScrChCh = 1.0d0
 ScrChPl = 1.0d0

  do i = 1, npol
    ScrChCh(i,i) = 0.0d0
    ScrChPl(i,i) = 0.0d0
    do j = i+1, npol
      LN1213 = (Neigh(i,j) .eq. 2 .or. Neigh(i,j) .eq. 3)
      if (IMMPCn.eq.1.or.IMMPCn.eq.4) then
        if (LN1213) then
          ScrChCh(i,j) = 0.0d0
          ScrChCh(j,i) = 0.0d0
          ScrChPl(i,j) = 0.0d0
          ScrChPl(j,i) = 0.0d0
        endif
      endif
      if (LScrChPl.AND..NOT.LN1213.AND.(IMMPCn.EQ.3.OR.IMMPCn.EQ.4.OR.IMMPCn.EQ.5)) THEN
!
! Compute the distance between two polarizable sites
!
        R(1) = CPol(1,J)-CPol(1,I)
        R(2) = CPol(2,J)-CPol(2,I)
        R(3) = CPol(3,J)-CPol(3,I)
        Rij = Sqrt(R(1)*R(1) + R(2)*R(2) + R(3)*R(3))
        If (IScreen.eq.1) s = a*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
        If (IScreen.eq.2) s = b*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
        If (IScreen.eq.3) s = c*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
!
! Thole linear screening
!
        If (Rij.le.s) Then
          v = Rij/s 
          Scale3 = 4.0d0*(v**3)-3.0d0*(v**4)
          ScrChPl(i,j) = Scale3
          ScrChPl(j,i) = Scale3
        endif
      elseif (IMMPCn.eq.2.or.IMMPCn.eq.0) then
        Write(*,*) 'Code not ready to deal with groups option'
        STOP
      endif
    enddo
  enddo
  
  Deallocate (Neigh)
  Return
End Subroutine
