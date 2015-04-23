!
! -------------------------------------------------------------------
!     Compute the RMS between two quantities 
!
REAL*8 Function rms(n,QRef,Qvar)
implicit real*8 (A-H,O-Z)
dimension Qref(n),Qvar(n)
rms = 0.0d0
Do I=1,n
  rms = rms + (Qref(I)-Qvar(I))**2
Enddo
return 
End Function
