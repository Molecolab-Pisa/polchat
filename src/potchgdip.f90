! potchgdip.f90:   A Polarisation consistent charge-fitting tool 
!                  A Molecolab Tool www.molecolab.dcci.unipi.it/tools
!
! Copyright (C) 2014, 2015, 2016, 2017
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
subroutine potchgdip(ng,nq,R,R3,Rqg,D,scr,q,V,DipInd)

  use constants

  implicit real*8(a-h,o-z)

  real*8 :: R(4,nq,nq), R3(nq,nq), Rqg(4,ng,nq), q(nq), D(3*nq,3*nq), scr(nq,nq)
  real*8 :: Ex(3*nq,nq), Dip(3*nq), BxGrd(ng,3*nq), Vdip(ng), Vq(ng), V(ng)
  real*8 :: E(3*nq), DipInd(3)

! Sanitise Rij3

  do i = 1, nq
    do j = 1, nq
      if (R3(i,j).lt.small) R3(i,j)=small
    enddo
  enddo

! Compute electric field due to charges
  
  Ex(1:nq,:)        = R(1,:,:)*(scr/R3)
  Ex(nq+1:2*nq,:)   = R(2,:,:)*(scr/R3)
  Ex(2*nq+1:3*nq,:) = R(3,:,:)*(scr/R3)
  E = matmul(Ex,q) 

! Compute dipoles induced by electric field

  Dip = matmul(D,E)

! Compute potential on the grid due to dipoles

  BxGrd(:,1:nq)        = Rqg(1,:,:)/(Rqg(4,:,:)**3)
  BxGrd(:,nq+1:2*nq)   = Rqg(2,:,:)/(Rqg(4,:,:)**3)
  BxGrd(:,2*nq+1:3*nq) = Rqg(3,:,:)/(Rqg(4,:,:)**3)
  Vdip = matmul(BxGrd,Dip)

! Compute potential on the grid due to charges

  call potchg(ng,nq,q,Rqg(4,:,:),Vq)

! Compute total potential

  V = Vdip + Vq
  if (iprt.ge.2) then
    call PrtMat(iout,ng,1,Vq,'PESP charge potential',.false.)
    call PrtMat(iout,ng,1,Vdip,'PESP dipole potential',.false.)
    call PrtMat(iout,ng,1,Vdip,'PESP total potential',.false.)
  endif

! Compute total induced dipole

  do i = 1, 3
    DipInd(i) = sum(Dip((i-1)*nq+1:i*nq))
  enddo

  return

end subroutine
