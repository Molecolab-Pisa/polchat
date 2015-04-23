!--------------------------------------------------------------------
! Function to compute the tensor
!
REAL*8 Function TensorMM(M,N,I,J,CPol,Pol,npol,DoThole,IScreen)
IMPLICIT REAL*8 (A-H,O-Z)
LOGICAL :: DoThole
DIMENSION  :: CPol(3,*),Pol(*),R(3)
PARAMETER(a=1.7278d0,b=2.5874d0,c=2.0580d0)
!
! > DoThole .... True if you want to use Thole  
! > IScreen .... 1 = use Amber Screening Factor (a)
!           .... 2 = use AL (Wang) Screening Factor (b)
!           .... 3 = use DL (Wang) Screening Factor (c)
! > M .......... x,y,z component of site I
! > N .......... x,y,z component of site J
! > I .......... index of I polarizable site
! > J .......... index of J polarizable site
!
!
! Compute the distance between two polarizable sites.
!
R(1) = CPol(1,J)-CPol(1,I)
R(2) = CPol(2,J)-CPol(2,I)
R(3) = CPol(3,J)-CPol(3,I)
Rij = Sqrt(R(1)*R(1) + R(2)*R(2) + R(3)*R(3))
!
If (DoThole) then
  If (IScreen.eq.1) s = a*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
  If (IScreen.eq.2) s = b*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
  If (IScreen.eq.3) s = c*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
EndIf
! 
! Thole smeared dipole interaction tensor 
! (Van Duijnen parameters from JPCA, 102, 2399, 1998)
!
If (DoThole.and.Rij.le.s) Then
  v = Rij/s
  Scale3 = 4.0d0*(v**3)-3.0d0*(v**4)
  Scale5 = v**4
  If (M.eq.N) Then
    TensorMM = (1/(Rij**3))*(Scale3-(Scale5*(3/(Rij**2))*R(M)*R(N)))
  Else
    TensorMM = - Scale5*(3/(Rij**5))*R(M)*R(N)
  Endif
!
! Standard dipole interaction tensor
!
Else
  If (M.eq.N) Then
    TensorMM = (1/(Rij**3))*(1 - ((3/(Rij**2))*R(M)*R(N)))
  Else
    TensorMM = - (3/(Rij**5))*R(M)*R(N)
  Endif
Endif
Return
End Function
