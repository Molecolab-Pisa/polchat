! mkscr.f90:       A Polarisation consistent charge-fitting tool 
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
subroutine MkScr

  use constants
  use mmpoldata
  use gespinfo
  use time
 
  implicit real*8 (a-h,o-z)

  logical :: LN1213

 2000 format(' SCR: Sites ',i5,',',i5,' are 1-2/1-3 neighb.; scr c-c: ',f6.3,', scr c-p: ',f6.3)
 2010 format(' SCR: Sites ',i5,',',i5,' are not close;       scr c-c: ',f6.3,', scr c-p: ',f6.3,'; dist: ',f8.2,' au, s=',f5.2)
 2020 format(' SCR: Sites ',i5,',',i5,' are not close;       scr c-c: ',f6.3,', scr c-p: ',f6.3,'; (screening not active)') 
 2100 format(' Computing chg-chg and chg-dip screening.')
 9000 format(' ERROR',/,&
             ' Mismatch in type of polarisation scheme.')
 9010 format(' ERROR',/,&
             ' Code does not deal with groups option yet.')

  call starttime

  if (iprt.gt.0) write(iout,2100)

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

      if (iprt.ge.2) then
        if (LN1213) then
          write(iout,2000) i,j,scrcc(i,j),scrcp(i,j)
        else
          if (lscr) then 
            write(iout,2010) i,j,scrcc(i,j),scrcp(i,j),dist,s
          else
            write(iout,2020) i,j,scrcc(i,j),scrcp(i,j)
          endif
        endif
      endif

    enddo
  enddo

  call gettime('computing screenings')

! Printout

  if (iprt .gt. 2) then
    call PrtMat(iout,NChg,NChg,scrcc,' Charge-charge screening',.false.)
    call PrtMat(iout,NChg,NChg,scrcp,' Charge-dipole screening',.false.)
    call gettime('debug printout')
  endif

  if (iprt.ge.1) call prttime('computing screenings')

  return

end subroutine
