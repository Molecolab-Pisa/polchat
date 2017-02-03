! time.f90:        A Polarisation consistent charge-fitting tool 
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
module time

  use constants

  integer                                      :: timercnt, timerate, timestart
  integer, dimension(timemax)                  :: timer
  character(len=outmsgmax), dimension(timemax) :: tmess

contains

! inittime

  subroutine inittime
    call system_clock(count_rate=timerate)
    return
  end subroutine

! starttime

  subroutine starttime
    timercnt = 0
    call system_clock(count=timestart)
    return
  end subroutine 

! gettime

  subroutine gettime(message)

    character(*) :: message
 9000 format(' WARNING',/,&
             ' Max number of time records exceeded. Not recording anymore.')

    if (timercnt .gt. timemax) then
      write(iout,9000)
    else
      timercnt = timercnt + 1
      call system_clock(count=timer(timercnt))
      tmess(timercnt) = message
    endif
    return
  end subroutine

! prttime
 
  subroutine prttime(message)

    character(*) :: message
    integer :: i, ttot
 1000 format(1x,55('-'))
 1010 format(1x,' Times for ',(A))
 1020 format(3x,(A),1x,f9.3,' s')
 1030 format(44x,11('-'),/,44x,f9.3,' s')
 9000 format('WARNING: No time to display.')

    if (timercnt.le.0) write(iout,9000)
    write(iout,1000)
    write(iout,1010) message
    write(iout,1020) tmess(1),dble(timer(1)-timestart)/dble(timerate)
    do i = 2, timercnt
      write(iout,1020) tmess(i),dble(timer(i)-timer(i-1))/dble(timerate)
    enddo
    if (timercnt.gt.1) write(iout,1030) dble(timer(timercnt)-timestart)/dble(timerate)
    write(iout,1000)
    return
  end subroutine

end module
