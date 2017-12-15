! tensor.f90:      A Polarisation consistent charge-fitting tool 
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
real*8 function Tensor(m,n,i,j)

  use constants
  use mmpoldata

  if (DoThole) then
    select case (IScreen)
      case (1)
        s = scra*((pol(i)*pol(j))**sixth)
      case (2)
        s = scrb*((pol(i)*pol(j))**sixth)
      case (3)
        s = scrc*((pol(i)*pol(j))**sixth)
    end select
  endif

  dist  = Rij(4,i,j)

! Thole smeared dipole interaction tensor
! van Duijnen, JPCA, 102, 2399 (1998)

  if (DoThole .and. dist.le.s) then
    v = dist/s
    Scale3 = four*(v**3)-three*(v**4)
    Scale5 = v**4

! Standard dipole interaction tensor

  else
    Scale3 = one
    Scale5 = one
  endif

  if (m .eq. n) then
    Tensor = (1/Rij3(i,j))*(Scale3-(Scale5*(3/(dist**2))*Rij(m,i,j)*Rij(n,i,j)))
  else
    Tensor = -Scale5*(3/(dist**5))*Rij(m,i,j)*Rij(n,i,j)
  endif

  return

end function
