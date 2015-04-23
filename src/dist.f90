!--------------------------------------------------------------------
!     Subroutine to compute the distances between the gridpoint and
!     the atoms.
!     If symmetric, also compute R^3
!
Subroutine dist(nch,ngr,LSymm,CChg,CGrd,R,R3)
  IMPLICIT REAL*8 (A-H,O-Z)
  LOGICAL :: LSymm
  DIMENSION CChg(3,nch),CGrd(3,ngr),R(4,ngr,nch),R3(nch,nch)
!
  R = 0.0d0
  if (LSymm) then 
    do i = 1, nch
      do j = i+1, nch
        dx  = CChg(1,i) - CChg(1,j)
        dy  = CChg(2,i) - CChg(2,j)
        dz  = CChg(3,i) - CChg(3,j)
        dd  = sqrt(dx*dx + dy*dy + dz*dz)
        dd3 = dd**3
        R(1,i,j) = dx
        R(1,j,i) = -dx
        R(2,i,j) = dy
        R(2,j,i) = -dy
        R(3,i,j) = dz
        R(3,j,i) = -dz
        R(4,i,j) = dd
        R(4,j,i) = dd
        R3(i,j)  = dd3
        R3(j,i)  = dd3
      enddo
    enddo
  else
    do i = 1, ngr
      do j = 1, nch
        dx = CGrd(1,i) - CChg(1,j)
        dy = CGrd(2,i) - CChg(2,j)
        dz = CGrd(3,i) - CChg(3,j)
        dd = sqrt(dx*dx + dy*dy + dz*dz)
        R(1,i,j) = dx
        R(2,i,j) = dy
        R(3,i,j) = dz
        R(4,i,j) = dd
      enddo
    enddo
  endif
        
  RETURN
End Subroutine
