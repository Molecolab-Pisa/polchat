!
! -------------------------------------------------------------------
! Compute the fitting error
!
Function FitErr(IOut,IPrint,LPol,nch,ngr,Q,RChGr,RChCh,R3,DInv,Scr,Vqm,DipInd)
  Implicit Real*8(A-H,O-Z)
  Logical LPol
  Real*8 FitErr,Q(nch),RChGr(4,ngr,nch),RChCh(4,nch,nch),R3(nch,nch),Vch(ngr)
  Real*8 DInv(3*nch,3*nch),Scr(nch,nch), Vqm(ngr), DipInd(3)
  if (LPol) then
    call chpotpol(nch,ngr,Q,DInv,Scr,RChGr,RChCh,R3,Vch,DipInd)
  else
    call chpot(nch,ngr,Q,RChGr,Vch)
  endif
  FitErr = rms(nch,Vqm,Vch)
  return
End Function
