! PESP.f90: A Polarisation consistent charge-fitting tool 
!           A Molecolab Tool www.dcci.unipi.it/molecolab/tools
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
Subroutine PESP(IOut,IPrint,npol,ngr,mcon,RChGr,Rij,Rij3,Vqm,X,B,ScrChPl,ScrChCh,DInv,Qpesp,QSum,restr,irestr,nrestr)
  Implicit Real*8 (A-H,O-Z)
  Dimension RChGr(4,ngr,npol), Vqm(ngr), X(npol+mcon,npol+mcon), B(npol+mcon)
  Dimension Rij(4,npol,npol), Rij3(npol,npol)
  Dimension Qpesp(npol+mcon), ScrChPl(npol,npol), ScrChCh(npol,npol)
  Dimension dVdQ(ngr,npol),dEdQ(3*npol,npol),dMudQ(3*npol,npol),DInv(3*npol,3*npol)
  Dimension RTemp(ngr,3*npol),dVdQChg(ngr,npol),dVdQPol(ngr,npol), E(npol,npol), irestr(nrestr)
  Integer, Allocatable :: IPIV(:)
  Real*8, Allocatable  :: WORK(:)
  Real*8, Allocatable  :: Ax(:,:),BxPol(:,:),BxGrd(:,:),Cx(:,:),Ex(:,:),Fx(:,:),Gx(:,:),Hx(:)
  Parameter (One = 1.0d0, Small = 1.0d-10)
 1000 Format(' polESP matrix: Element ',i6,' has an illegal value',/,' Halting.')
 1010 Format(' polESP matrix: Diagonal element ',i6,' is zero.',/, &
             ' The matrix is singular and cannot be inverted',/,   &
             ' If you have a symmetric molecule and you are imposing a dipole',/,&
             '   constraint,  then  you  may be asking for linearly dependent',/,&
             '   conditions.    Please remove constraint on dipole and see if',/,&
             '   it works. Otherwise contact the author.---------------------',/, &
             ' Halting now.')
 1020 Format(' polESP matrix: Inversion successful.')
 1100 Format(' polESP charges:')
 1110 Format(1x,i6,1x,f12.6)
 1120 Format(' --------------------',/,' Total  ',f12.6)
!
! For definitions and notation, please refer to document by S.Caprasecca
!
  Allocate(Ax(ngr,npol), BxPol(3*npol,npol), BxGrd(ngr,3*npol), Cx(3*npol,npol))
  Allocate(Ex(ngr,npol), Fx(ngr,npol), Gx(npol,npol), Hx(npol))
! --- Matrix Ax
  Ax = One/RChGr(4,:,:)
! --- Sanitise Rij3
  do i = 1, npol
    do j = 1, npol
      if (Rij3(i,j).lt.Small) Rij3(i,j)=Small
    enddo
  enddo
! --- Matrix BxPol
  BxPol(1:npol,:)          = Rij(1,:,:)*(ScrChPl/(Rij3))
  BxPol(npol+1:2*npol,:)   = Rij(2,:,:)*(ScrChPl/(Rij3))
  BxPol(2*npol+1:3*npol,:) = Rij(3,:,:)*(ScrChPl/(Rij3))
! --- Matrix BxGrd
  BxGrd(:,1:npol)          = RChGr(1,:,:)/(RChGr(4,:,:)**3)
  BxGrd(:,npol+1:2*npol)   = RChGr(2,:,:)/(RChGr(4,:,:)**3)
  BxGrd(:,2*npol+1:3*npol) = RChGr(3,:,:)/(RChGr(4,:,:)**3)
! --- Matrix Cx
  Cx = Matmul(DInv,BxPol)
! --- Matrix Ex
  Ex = Matmul(BxGrd,Cx)
! --- Matrix Fx
  Fx = Ax + Ex
! --- Matrix Gx
  Gx = Matmul(Transpose(Fx),Fx)
! --- Vector Hx
  Hx = Matmul(Transpose(Fx),Vqm)
! --- Matrix X
  X(1:npol,1:npol) = Gx
! --- Restraints
  if (restr.gt.1.0d-08) then
    do k = 1, nrestr
      X(irestr(k),irestr(k)) = X(irestr(k),irestr(k)) + restr
    enddo
  endif
! --- Vector B
  B(1:npol) = Hx
!---   
!---   
!---     Rij3 = Rij3 + 1.0d-10
!---   ! Form matrix dVdQ (Chg)
!---     dVdQChg = 1.0d0/RChGr(4,:,:)
!---   ! Form matrix dEdQ
!---     do i = 1, npol
!---       do j = 1, npol
!---         dEdQ(i,j) = ScrChPl(i,j)*Rij(1,i,j)/(Rij3(i,j))
!---         dEdQ(npol+i,j) = ScrChPl(i,j)*Rij(2,i,j)/(Rij3(i,j))
!---         dEdQ(2*npol+i,j) = ScrChPl(i,j)*Rij(3,i,j)/(Rij3(i,j))
!---       enddo
!---     enddo
!---   ! Form matrix dMudQ
!---     dMudQ = matmul(DInv,dEdQ)
!---   ! Form matrix dVdQ (Pol)
!---     RTemp(:,1:npol) = RChGr(1,:,:)
!---     RTemp(:,npol+1:2*npol) = RChGr(2,:,:)
!---     RTemp(:,2*npol+1:3*npol) = RChGr(3,:,:)
!---     dVdQPol = matmul(RTemp,dMudQ)
!---   ! Form matrix dVdQ (Tot)
!---     dVdQ = dVdQChg-dVdQPol
!---   ! Form matrix X
!---     X(1:npol,1:npol) = matmul(transpose(dVdQ),dVdQ)
!---   ! Form vector B
!---     B(1:npol) = matmul(transpose(dVdQ),Vqm)
! Invert matrix X
! Solve system to find ESP charges: 
!  (-) initialise elements
  LWORK = (npol+mcon)*(npol+mcon)
  allocate (IPIV(npol+mcon),WORK(LWORK))
  
!  (-) invert matrix X

  M = npol+mcon
  INFO = 0
  Call DGETRF(M,M,X,M,IPIV,INFO)
  Call DGETRI(M,X,M,IPIV,WORK,LWORK,INFO)

  if (INFO .lt. 0) then
    write(IOut,1000) abs(INFO)
    stop
  elseif (INFO .gt. 0) then
    write(IOut,1010) INFO
  elseif (IPrint.ge.2) then
    write(IOut,1020)
  endif

!   (-) Compute charges

  qpESP = matmul(X,B)

  if (IPrint.ge.2) write(IOut,1100)
  qSum = 0.d0
  do i = 1, npol
    qSum = qSum + qpESP(i)
    if (IPrint.ge.2) write(IOut,1110) i,qpESP(i)
  enddo
  if (IPrint.ge.2) write(IOut,1120) qSum

 Deallocate(IPIV,WORK)

  Return
End Subroutine
