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
