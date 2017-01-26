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
