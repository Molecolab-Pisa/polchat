! gespinfo.f90:    A Polarisation consistent charge-fitting tool 
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
module gespinfo

  implicit none

  integer             :: NGrd
  real*8              :: DipQM(3), eGESP
  real*8, allocatable :: gesp(:), CGrd(:,:), VQM(:), RChGr(:,:,:)
  
end module
