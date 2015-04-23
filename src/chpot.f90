!
! -------------------------------------------------------------------
!     Subroutine to compute the potential generated by a set of
!     charges over the gridpoints.
!
Subroutine chpot(nch,ngr,Chg,R,Vch)
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION Chg(nch),R(4,ngr,nch),Vch(ngr)
DO I=1,ngr
  Vch(I) = 0.0d0
  DO J=1,nch
    Vch(I) = Vch(I) + Chg(J)/R(4,I,J)
  ENDDO
ENDDO
RETURN
End Subroutine
