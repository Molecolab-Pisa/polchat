subroutine MkScr

  use constants
  use mmpoldata
  use gespinfo
 
  implicit real*8 (a-h,o-z)

  logical :: LN1213

 9000 format(' ERROR',/,&
             ' Mismatch in type of polarisation scheme.')
 9010 format(' ERROR',/,&
             ' Code does not deal with groups option yet.')

  if (lscr) then
    select case (IMMPCn)
      case (1)
        DoThole = .false.
        IScreen = 1
      case (3)
        DoThole = .true.
        IScreen = 1
      case (4)
        DoThole = .true.
        IScreen = 2
      case (5)
        DoThole = .true.
        IScreen = 3
      case default
        write(iout,9000)
        stop
    end select
  endif

! Build screening matrices (Amber12 only for now)
  scrcc = one
  scrcp = one

  do i = 1, NChg
    scrcc(i,i) = zero
    scrcp(i,i) = zero
    do j = i+1, NChg
      LN1213 = (neigh(i,j).eq.2 .or. neigh(i,j).eq.3)
      if (IMMPCn.eq.1 .or. IMMPCn.eq.4) then
        if (LN1213) then
          scrcc(i,j) = zero
          scrcc(j,i) = zero
          scrcp(i,j) = zero
          scrcp(j,i) = zero
        endif
      endif
      if (lscr .and. .not.LN1213 .and. (IMMPCn.eq.3.or.IMMPCn.eq.4.or.IMMPCn.eq.5)) then
        select case (IScreen)
          case (1)
            s = scra*((pol(i)*pol(j))**sixth)
          case (2)
            s = scrb*((pol(i)*pol(j))**sixth)
          case (3)
            s = scrc*((pol(i)*pol(j))**sixth)
        end select

        dist = Rij(4,i,j)
        if (dist .le. s) then
          v = dist/s
          Scale3 = four*(v**3)-three*(v**4)
          scrcp(i,j) = Scale3
          scrcp(j,i) = Scale3
        endif
      elseif (IMMPCn.eq.2 .or. IMMPCn.eq.0) then
        write(iout,9010)
        stop
      endif
    enddo
  enddo

! Printout

  if (iprt .gt. 1) then
    call PrtMat(iout,NChg,NChg,scrcc,' Charge-charge screening',.false.)
    call PrtMat(iout,NChg,NChg,scrcp,' Charge-dipole screening',.false.)
  endif

  return

end subroutine
