! rdopts.f90:      A Polarisation consistent charge-fitting tool 
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
subroutine RdOpts
!
! This subroutine reads input filenames and printout options from line
!
  use constants
  use mmpoldata
  use operative
  use time

  integer                                             :: iarg,narg,iget
  character(len=inpargmax), dimension(:), allocatable :: args

 1000 format(' Error in input stream: nothing followed option ',(A))
 1001 format(' Error in input stream: option ',(A),' unknown')
 1200 format(' Use the following options to run the program:',/,              &
            '   -g  --gesp     (required) ... Followed by gesp file name',/,              &
            '   -m  --mol2     (required) ... Followed by mol2 file name',/,              &
            '   -p  --pol      (required) ... Followed by polarisability file name',/,    &
            '   -c  --constr   (required) ... Followed by constraints file name',/,       &
            '   -x  --screen   (optional) ... Activate Wang Chg-Pol screening',/,&
            '   -db --database (optional) ... Print database (followed by database file name)',/, &
            '   -h  --help     (optional) ... Get this help message',/,                &
            '   -v  --verbose  (optional) ... Run in debug mode (extra printout)',/,   &
            '   -d  --debug    (optional) ... Run in extra debug mode (lots of printout)',/,&
            '   -s  --silent   (optional) ... Run in silent mode (minimum printout)')

  call starttime

  narg = command_argument_count()
  allocate(args(narg))

  iget = 0
  iprt = 0
  lscr = .false.
  ldbs = .false.
  filenam = ''
  filecon = ''
  filepol = ''
  filecns = ''
  filedbs = ''

  if (narg.eq.0) then
    write(IOut,1200)
    stop
  endif

  do iarg = 1, narg
    call get_command_argument(iarg,args(iarg))
  enddo

  do iarg = 1, narg
    if (iget.ne.0) then
      select case (iget)
        case (1) 
          filenam = args(iarg)
        case (2)
          filecon = args(iarg)
        case (3)
          filepol = args(iarg)
        case (4)
          filecns = args(iarg)
        case (5) 
          filedbs = args(iarg)
      end select
      iget = 0
    else
      select case (args(iarg))
        case ('-h')
          write(iout,1200)
          stop
        case ('--help')
          write(iout,1200)
          stop
        case ('-s')
          iprt = -1
        case ('--silent')
          iprt = -1
        case ('-v')
          iprt = 1
        case ('--verbose')
          iprt = 1
        case ('-d')
          iprt = 2
        case ('--debug')
          iprt = 2
        case ('-x')
          lscr = .true.
        case ('--screen')
          lscr = .true.
        case ('-g')
          iget = 1
        case ('--gesp')
          iget = 1
        case ('-m')
          iget = 2
        case ('--mol2')
          iget = 2
        case ('-p')
          iget = 3
        case ('--pol')
          iget = 3
        case ('-c')
          iget = 4
        case ('--constr')
          iget = 4
        case ('-db')
          ldbs = .true.
          iget = 5
        case ('--database')
          ldbs = .true.
          iget = 5
        case default
          write(iout,1001) trim(args(iarg))
          stop
      end select
      if (iget.ne.0 .and. iarg.eq.narg) then
        write(iout,1000) trim(args(iarg))
        stop
      endif
    endif
  enddo

  call gettime('')
  if (iprt.ge.1) call prttime('reading arguments')

  return

end subroutine
