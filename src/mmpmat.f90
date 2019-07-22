! mmpmat.f90:      A Polarisation consistent charge-fitting tool 
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
subroutine MMPMat

  use constants
  use mmpoldata
  use time

  implicit real*8(a-h,o-z)

  logical            :: DoCalc, DoInter
  integer, parameter :: lwork=30000
  integer            :: IPIV(3*NChg)
  real*8             :: work(lwork)

 200 format(' Building MMPol matrix.')
 205 format(' MM sites ',i5,1x,i5,' neighb. ',i1,' inteaction ',l1)

  if (iprt.ge.1) write(iout,200)
  call starttime

! Clear matrix
 
  D = 0

! Fill diagonal elements

  do i = 1, NChg
    D(i,       i)        = one/pol(i)
    D(i+NChg,  i+NChg)   = one/pol(i)
    D(i+2*NChg,i+2*NChg) = one/pol(i)
  enddo

! Fill non-diagonal elements 

  do i = 1, NChg
    do j = i+1, NChg
      DoCalc = .false.
      if (IMMPCn.eq.3 .or. IMMPCn.eq.5) then
        DoCalc = .true.
      elseif (IMMPCn.eq.0 .or. IMMPCn.eq.1 .or. IMMPCn.eq.2 .or. IMMPCn.eq.4) then
        DoCalc = DoInter(i,j)
      endif
   
 999 format(3(1x,i5),2(1x,l))
      if (iprt.ge.2) write(iout,205) i,j,neigh(i,j),DoCalc
      if (DoCalc) then
        D(i,j) = Tensor(1,1,i,j)
        D(j,i) = D(i,j)
        D(i,j+NChg) = Tensor(1,2,i,j)
        D(j,i+NChg) = D(i,j+NChg)
        D(i+NChg,j) = D(i,j+NChg)
        D(j+NChg,i) = D(i,j+NChg)
        D(i,j+2*NChg) = Tensor(1,3,i,j)
        D(j,i+2*NChg) = D(i,j+2*NChg)
        D(i+2*NChg,j) = D(i,j+2*NChg)
        D(j+2*NChg,i) = D(i,j+2*NChg)
        D(i+NChg,j+NChg) = Tensor(2,2,i,j)
        D(j+NChg,i+NChg) = D(i+NChg,j+NChg)
        D(i+NChg,j+NChg*2) = Tensor(2,3,i,j)
        D(j+NChg,i+NChg*2) = D(i+NChg,j+NChg*2)
        D(i+NChg*2,j+NChg) = D(i+NChg,j+NChg*2)
        D(j+NChg*2,i+NChg) = D(i+NChg,j+NChg*2)
        D(i+NChg*2,j+NChg*2) = Tensor(3,3,i,j)
        D(j+NChg*2,i+NChg*2) = D(i+NChg*2,j+NChg*2)
      endif
    enddo
  enddo

  call gettime('Forming matrix')

  if (iprt.ge.2) then
    call PrtMat(iout,3*NChg,3*NChg,D,'MMPol matrix before inversion',.false.)
    call gettime('Debug printing 1')
  endif

! Invert the MMPol matrix using lapac library

  LDA = 3*NChg
  call DGETRF(3*NChg,3*NChg,D,LDA,IPIV,info)
  call DGETRI(3*NChg,D,LDA,IPIV,work,lwork,info)

  call gettime('Inverting matrix')
  
  if (iprt.ge.2) then
    call PrtMat(iout,3*NChg,3*NChg,D,'MMPol matrix after inversion',.false.)
    call gettime('Debug printing 2')
  endif

  if (iprt.ge.1) call prttime('forming and inverting MMPol matrix')

  return

end subroutine
