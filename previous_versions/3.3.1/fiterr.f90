! FitErr.f90: A Polarisation consistent charge-fitting tool 
!             A Molecolab Tool www.dcci.unipi.it/molecolab/tools
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
    if (IPrint.ge.2) call PrtMat(IOut,nch,1,Q,'ESP chg')
    if (IPrint.ge.2) call PrtMat(IOut,ngr,nch,RChGr(4,:,:),'dist')
    if (IPrint.ge.2) call PrtMat(IOut,ngr,1,Vch,'Chg potential')
  endif
  FitErr = rms(nch,Vqm,Vch)
  return
End Function
