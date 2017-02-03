! dointer.f90:     A Polarisation consistent charge-fitting tool 
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
logical function DoInter(i,j)

  use constants
  use mmpoldata

  implicit real*8(a-h,o-z)
  logical ll

 9000 format(' ERROR',/,&
             ' Possible mismatch in connectivity matrix.')

  DoInter = .true.

! Groups option
  if (IMMPCn.eq.2.or.IMMPCn.eq.3) then
    IinJ = 0
    JinI = 0
    do k = 1, LAnMMP
      if (IAnMMP(1,i).eq.IAnMMP(k,j)) IinJ = 1
      if (IAnMMP(k,i).eq.IAnMMP(1,j)) JinI = 1
    enddo
    if (IinJ.ne.JinI) then
      write(iout,9000)
      stop
    elseif (IinJ.eq.1) then
      DoInter = .false.
    endif

! 1-2, 1-3 neighbour
  elseif (IMMPCn.eq.1.or.IMMPCn.eq.4) then
    if (neigh(i,j) .le. 3 .and. neigh(i,j) .ge. 1) DoInter = .false.
  endif

  return

end function

