! constants.f90    A Polarisation consistent charge-fitting tool 
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
module constants

  implicit none

  character(50)     :: version = "4.1.2"
  integer           :: IMMPCn = 4
  integer,parameter :: LAnMMP = 9
  integer           :: iout = 6
  integer           :: iprt
  integer,parameter :: outmsgmax=40, strmax=500, inpargmax=100, typmax=6
  integer,parameter :: timemax=10
  real*8, parameter :: zero=0.0d0, one=1.0d0, sixth=1.0d0/6.0d0, small=1.0d-8
  real*8, parameter :: three=3.0d0, four=4.0d0
  real*8, parameter :: ang2au=1.889725989d0, deb2au=0.393430201407683d0
  real*8, parameter :: scra=1.7278d0, scrb=2.5874d0, scrc=2.0580d0

end module
