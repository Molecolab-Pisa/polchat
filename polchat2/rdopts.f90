subroutine RdOpts
!
! This subroutine reads input filenames and printout options from line
!
! Use the following options:
!  -h ... Get help message
!  -d ... Debug mode (extra printout)
!  -s ... Silent mode (minimum printout)
!  -g ... before gesp file name
!  -m ... before mol2 file name
!  -p ... before specifying polarizability file name
!  -c ... before specifying constraint file name
!  -x ... exclude chg-pol screening as in Wang
!
! Example: chpol d g xxx.gesp m xxx.mol2 p pol.in c cnstr.in
!
  ! DECLARATION

  use constants
  use mmpoldata
  use operative

  integer                                       :: iarg,narg,iget
  character(len=100), dimension(:), allocatable :: args

 1000 format(' Error in input stream: nothing followed option ',(A))
 1001 format(' Error in input stream: option ',(A),' unknown')
 1200 format(' Use the following options to run the program:',/,              &
            '   -g  (required) ... Followed by gesp file name',/,              &
            '   -m  (required) ... Followed by mol2 file name',/,              &
            '   -p  (required) ... Followed by polarisability file name',/,    &
            '   -c  (required) ... Followed by constraints file name',/,       &
            '   -x  (optional) ... Do not include Wang Chg-Pol screening --',/,&
            '                      Only use if the MMPol code does not include such screening. If in doubt, do not use.',/, &
            '   -db (optional) ... Print database (followed by database file name)',/, &
            '   -h  (optional) ... Get this help message',/,                &
            '   -d  (optional) ... Run in debug mode (extra printout)',/,   &
            '   -dd (optional) ... Run in extra debug mode (lots of printout)',/,&
            '   -s  (optional) ... Run in silent mode (minimum printout)')

  narg = command_argument_count()
  allocate(args(narg))

  iget = 0
  iprt = 1
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
    elseif (args(iarg) .eq. '-h') then
      write(iout,1200)
      stop
    elseif (args(iarg) .eq. '-d') then
      iprt = 1
    elseif (args(iarg) .eq. '-dd') then
      iprt = 2
    elseif (args(iarg) .eq. '-s') then
      iprt = -1
    elseif (args(iarg) .eq. '-x') then
      lscr = .true.
    elseif (args(iarg) .eq. '-g') then
      iget = 1
      if (iarg .eq. narg) goto 800
    elseif (args(iarg) .eq. '-m') then
      iget = 2
      if (iarg .eq. narg) goto 800
    elseif (args(iarg) .eq. '-p') then
      iget = 3
      if (iarg .eq. narg) goto 800
    elseif (args(iarg) .eq. '-c') then
      iget = 4
      if (iarg .eq. narg) goto 800
    elseif (args(iarg) .eq. '-db') then
      ldbs = .true.
      iget = 5
      if (iarg .eq. narg) goto 800
    else
      goto 801
    endif
  enddo

  goto 802 

 800 write(iout,1000) args(iarg-1)
  stop

 801 write(iout,1001) args(iarg)
  stop

 802 return

end subroutine
