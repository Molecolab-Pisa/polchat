! pesp.f90:        A Polarisation consistent charge-fitting tool 
!                  A Molecolab Tool www.molecolab.dcci.unipi.it/tools
!
! Copyright (C) 2014, 2015, 2016, 2017
!   S. Caprasecca, C. Curutchet, B. Mennucci
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
subroutine pesp

  use constants
  use mmpoldata
  use gespinfo
  use constraints
  use espinfo
  use time

  implicit real*8(a-h,o-z)
  integer, allocatable :: IPIV(:)
  real*8,  allocatable :: WORK(:)
  real*8, allocatable  :: Ax(:,:), BxPol(:,:), BxGrd(:,:), Cx(:,:), Ex(:,:)
  real*8, allocatable  :: Fx(:,:), Gx(:,:), Hx(:)

 2000 format(' PESP matrix: Inversion successful.')
 2010 format(' PESP charges computed.')
 2020 format(' Sum of PESP charges: ',f12.6)
 2030 format(' Fit error of PESP charges: ',f12.6)
 2100 format(' Fitting pol-ESP charges.')
 9000 format(' ERROR',/,&
             ' PESP matrix: Element ',i6,' has an illegal value.')
 9010 format(' ERROR',/,&
             ' PESP matrix: Diagonal element ',i6,' is zero.',/,&
             '              The matrix is singular and cannot be inverted.')

! For definition and notation, please refer to documentation

  if (iprt.ge.0) write(iout,2100)
  call starttime

! Allocation

  allocate(Ax(NGrd,NChg), BxPol(3*NChg,NChg), BxGrd(NGrd,3*NChg), Cx(3*NChg,NChg))
  allocate(Ex(NGrd,NChg), Fx(NGrd,NChg), Gx(NChg,NChg), Hx(NChg))

  call gettime('allocating')

! Matrix Ax

  Ax = one/RChGr(4,:,:)

! Sanitise Rij3

  do i = 1, NChg
    do j = i, NChg
      if (Rij3(i,j).lt.small) then
        Rij3(i,j) = small
        Rij3(j,i) = small
      endif
    enddo
  enddo

! Matrix BxPol

  BxPol(1:NChg,:)          = Rij(1,:,:)*(scrcp/(Rij3))
  BxPol(NChg+1:2*NChg,:)   = Rij(2,:,:)*(scrcp/(Rij3))
  BxPol(2*NChg+1:3*NChg,:) = Rij(3,:,:)*(scrcp/(Rij3))

! Matrix BxGrd
 
  BxGrd(:,1:NChg)          = RChGr(1,:,:)/(RChGr(4,:,:)**3)
  BxGrd(:,NChg+1:2*NChg)   = RChGr(2,:,:)/(RChGr(4,:,:)**3)
  BxGrd(:,2*NChg+1:3*NChg) = RChGr(3,:,:)/(RChGr(4,:,:)**3)

! Matrix Cx

  Cx = matmul(D,BxPol)

! Matrix Ex

  Ex = matmul(BxGrd,Cx)

! Matrix Fx

  Fx = Ax + Ex

! Matrix Gx

  Gx = matmul(transpose(Fx),Fx)

! Vector Hx

  Hx = matmul(transpose(Fx),VQM)

! Matrix X
 
  X(1:NChg,1:NChg) = Gx

! Add restraints

  if (NCRes.ne.0 .and. RCRes.gt.small) then
    do k = 1, MCRes
      X(VCRes(k),VCRes(k)) = X(VCRes(k),VCRes(k)) + RCRes
    enddo
  endif

! Vector B

  B(1:NChg) = Hx

  call gettime('forming matrices')

! Initialise elements

  NDim = NChg+NCons
  LWORK = NDim**2
  allocate (IPIV(NDim), WORK(LWORK))

! Invert X

  INFO = 0
  call DGETRF(NDim,NDim,X,NDim,IPIV,INFO)
  call DGETRI(NDim,X,NDim,IPIV,WORK,LWORK,INFO)

  if ( INFO .lt. 0 ) then
    write(iout,9000) abs(INFO)
    stop
  elseif ( INFO .gt. 0 ) then
    write(iout,9010) INFO
    stop
  elseif (iprt.ge.1) then
    write(iout,2000)
  endif

  call gettime('matrix inversion')
    
! Compute PESP charges and fit error

  allocate (qpesp(NChg))
  qpesp = matmul(X,B)
  spesp = sum(qpesp)
  epesp = error(.true.,qpesp)
  call dipole(NChg,qpesp,CChg,dqpesp)
  dtpesp = dqpesp + ddpesp

  if (iprt.ge.2) then
    write(iout,2010)
    call PrtMat(iout,NChg,1,qpesp,'PESP charges',.false.)
    write(iout,2020) spesp
    write(iout,2030) epesp
  endif

  deallocate(IPIV,WORK)
  deallocate(X,B)

  call gettime('solving for charges and computing fit errors')
  if (iprt.ge.1) call prttime('computing pol-ESP charges')

  return

end subroutine
