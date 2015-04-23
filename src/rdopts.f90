Subroutine RdOpts(IOut,IPrint,filename,fileconn,filepol,filecnst,LScrChPl)
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
!  -x ... include chg-pol screening as in Wang (experts only)
!
! Example: chpol d g xxx.gesp m xxx.mol2 p pol.in c cnstr.in
!
  integer IOut,IPrint,IArg,num_args
  logical LScrChPl
  CHARACTER(LEN=50) filename,filepol,fileconn,filecnst
  character(len=50), dimension(:), allocatable :: args

1000 format(' Error in input stream: nothing followed option ',(A))
1001 format(' Error in input stream: option ',(A),' unknown')
1002 format(' Error in input stream: not all files have been defined.',/, &
            '   Use g m p c to define all files.')
1200 format(' Use the following options to run the program:',/,            &
            '   -g (required) ... Followed by gesp file name',/,           &
            '   -m (required) ... Followed by mol2 file name',/,           &
            '   -p (required) ... Followed by polarisability file name',/, &
            '   -c (required) ... Followed by constraints file name',/,    &
            '   -x (optional) ... Also include Wang Chg-Pol screening -- WARNING: Expert use only.',/,             &
            '                     Only use if the MMPol code includes such screening. If in doubt, do not use.',/, &
            '   -h (optional) ... Get this help message',/,                &
            '   -d (optional) ... Run in debug mode (extra printout)',/,   &
            '   -s (optional) ... Run in silent mode (minimum printout)')

  num_args = command_argument_count()
  allocate(args(num_args))  

  IPrint   = 1
  IGet     = 0
  LScrChPl = .FALSE.
  filename = ''
  fileconn = ''
  filepol  = ''
  filecnst = ''

  if (num_args.eq.0) then
    write(IOut,1200)
    stop
  endif
  do IArg = 1, num_args
    call get_command_argument(IArg,args(IArg))
  end do
  do IArg = 1, num_args
    if (IGet.ne.0) then
      if (IGet.eq.1) then
        filename = args(IArg)
      elseif(IGet.eq.2) then
        fileconn = args(IArg)
      elseif(IGet.eq.3) then
        filepol = args(IArg)
      elseif(IGet.eq.4) then
        filecnst = args(IArg)
      endif
      IGet = 0
    elseif (args(IArg).eq.'-h') then
      write(IOut,1200)
      stop
    elseif (args(IArg).eq.'-d') then
      IPrint = 2
    elseif (args(IArg).eq.'-s') then
      IPrint = 0
    elseif (args(IArg).eq.'-x') then
      LScrChPl = .TRUE.
    elseif (args(IArg).eq.'-g') then
      IGet = 1
      if (IArg.eq.num_args) goto 800
    elseif (args(IArg).eq.'-m') then
      IGet = 2
      if (IArg.eq.num_args) goto 800
    elseif (args(IArg).eq.'-p') then
      IGet = 3
      if (IArg.eq.num_args) goto 800
    elseif (args(IArg).eq.'-c') then
      IGet = 4
      if (IArg.eq.num_args) goto 800
    else
      goto 801
    endif
  enddo

!  if (filename.eq.''.or.fileconn.eq.''.or.filepol.eq.''.or.filecnst.eq.'') goto 802 
  goto 900

800 write(IOut,1000) args(IArg-1)
  stop

801 write(IOut,1001) args(IArg)
  write(*,*) IArg
  stop
 
802 write(IOut,1002)
  stop
 
900 return

end subroutine
  
subroutine PrtHdr(IOut,VERSION)
  integer IOut
  character*50 :: Version
!
! Print header
!
1000 FORMAT(/,                                                                  &
     '                                                                      ',/,&
     '                                       mm                             ',/,&
     '                                    mMMm                              ',/,&
     '                                  mMMMMm         m                    ',/,&
     '                                 mMMMMm          mMm                  ',/,&
     '                                 mMMMMm          mMm                  ',/,&
     '                                 mMMMMMm        mMMm                  ',/,&
     '                                 MMMMMMMMMMMMMMMMMMm                  ',/,&
     '                                mMMMMMMMMMMMMMMMMMm                   ',/,&
     '       __  ___      __    ____________      __MMMm     __             ',/,&
     '      /  |/  /___  / /   / ____/ ____/___  / /  ____ _/ /_            ',/,&
     '     / /|_/ / __ \/ /   / __/ / /   / __ \/ /  / __ `/ __ \           ',/,&
     '    / /  / / /_/ / /___/ /___/ /___/ /_/ / /__/ /_/ / /_/ /           ',/,&
     '   /_/  /_/\__________/_____/\____/_____/_____|__,_/_.___/            ',/,&
     '           /_  __/ __ \/ __ \/ /  / ___/                              ',/,&
     '            / / / / / / / / / /   \__ \                               ',/,&
     '           / / / /_/ / /_/ / /___ __/ /                               ',/,&
     '          /_/  \____/\____/_____/____/                                ',/,&
     '            mMMMMMMMMMMMMMMMm                                         ',/,&
     '          mMMMMMMMMMMMMMMMm                                           ',/,&
     '        mMMMMMMMMMMMMMMMMM   + ------------------------------------ + ',/,&
     '       mMMMMMMMMMMMMMMMMm    |            P O L C H A T             | ',/,&
     '      mMMMMMMMMMMMMMMMMMm    + ------------------------------------ + ',/,&
     '      mMMMMm       mMMMMMm   | Stefano Caprasecca                   | ',/,&
     '      mMMMm       mMMMMMMm   | Sandro Jurinovich                    | ',/,&
     '       mMm       mMMMMMMm    | Carles Curutchet           ver ',A5,' | ',/,&
     '        m       mMMMMMMm     |          www.dcci.unipi.it/molecolab | ',/,&
     '               mMMMMMm       + ------------------------------------ + ',/)

  write(IOut,1000),VERSION
  return
end subroutine
