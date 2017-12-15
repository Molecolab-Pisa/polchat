! mmpoldata.f90:   A Polarisation consistent charge-fitting tool 
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
module mmpoldata

  use constants 

  implicit none

  logical                            :: lscr, DoThole
  integer                            :: NChg, IScreen
  integer, allocatable               :: IAnMMP(:,:), neigh(:,:)
  real*8, allocatable                :: pol(:), scrcc(:,:), scrcp(:,:)
  real*8, allocatable                :: CChg(:,:), Rij(:,:,:), Rij3(:,:)
  real*8, allocatable                :: D(:,:)
  character(len=typmax), allocatable :: atmnam(:), atmtyp(:), moltyp(:)

end module
