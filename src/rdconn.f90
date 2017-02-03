! rdconn.f90:      A Polarisation consistent charge-fitting tool 
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
subroutine RdConn

  use constants
  use mmpoldata
  use gespinfo
  use operative
  use strings
  use time

  implicit real(a-h,o-z)

  logical                              :: donemol, doneatm, donecon
  integer, allocatable                 :: ib1(:), ib2(:)
  character*1                          :: space
  character(len=strmax)                :: what
  character(len=strmax), dimension(20) :: args

 1010 format(i5,' | ',9(i5,1x))
 2000 format(' Reading mol2 file.')
 9000 format(' ERROR',/,&
             ' Mismatch in data read from files:',/,&
             '   Number of atoms read from gesp file: ',i6,/,&
             '   Number of atoms read from mol2 file: ',i6)
 9020 format(' ERROR',/,&
             ' Maximum number of neighbours exceeded.')
 9030 format(' ERROR',/,&
             ' In mol2 file found bond section before molecule section.')
 9040 format(' ERROR',/,&
             ' In mol2 file found atom section before molecule section.')
 9050 format(' ERROR',/,&
             ' Could not find molecule section in mol2 file.')
 9060 format(' ERROR',/,&
             ' Could not find bond section in mol2 file.')
 9070 format(' ERROR',/,&
             ' Could not find atom section in mol2 file.')

  if (iprt .gt. 1) write(iout,2000)
  call starttime
  open(unit=11,file=filecon,status='unknown')

  donemol = .false.
  doneatm = .false.
  donecon = .false.
  IAnMMP  = 0
  space   = ' '

! Find molecule information
  do while (.true.)
    read(11,*,END=800) what

! Read header and get number of atoms and bonds

    if (what .eq. '@<TRIPOS>MOLECULE') then
      read(11,*,END=800) what
      read(11,*,END=800) NAtoms, NBonds
      if (NAtoms .ne. NChg) then
        write(iout,9000) NAtoms, NChg
        stop
      endif
      allocate(ib1(NBonds),ib2(NBonds))
      donemol = .true.

! If rbase requested, read atom and molecule types 

    elseif (what .eq. '@<TRIPOS>ATOM' .and. ldbs) then
      if (.not.donemol) then
        write(iout,9040)
        stop
      endif
      allocate (atmtyp(NChg), moltyp(NChg))
      do i = 1, NAtoms
        read(11,'(A)',END=800) what
        call parse(trim(what),space,args,nargs)
        atmtyp(i) = trim(args(2))
        moltyp(i) = trim(args(8))
      enddo
      doneatm = .true.

! Read connectivity and store information on ib1 and ib2

    elseif (what .eq. '@<TRIPOS>BOND') then
      if (.not.donemol) then
        write(iout,9030)
        stop
      endif
      do i = 1, NBonds
        read(11,*,END=800) ixxx, ib1(i), ib2(i)
      enddo
      donecon = .true.
    endif
  enddo

 800 close(11)

  call gettime('Reading mol2 file')

! Deal with errors

  if (.not.donemol) then
    write(iout,9050)
    stop
  elseif (.not.donecon) then
    write(iout,9060)
    stop
  elseif (ldbs .and. .not.doneatm) then
    write(iout,9070)
    stop
  endif

! Build connectivity in standard MMPol format and store it on IAnMMP

  do n = 1, NChg
    k = 1
    do i = 1, NBonds
      if (n .eq. ib1(i)) then
        IAnMMP(n,k) = ib2(i)
        k = k+1
      endif
      if (n .eq. ib2(i)) then
        IAnMMP(n,k) = ib1(i)
        k = k+1
      endif
      if (k .gt. LAnMMP) then
        write(iout,9020)
        stop
      endif
    enddo
  enddo
  deallocate(ib1,ib2)

  call gettime('Building MMPol connectivity')

! Build neighbourhood matrix and store it

  do i = 1, NChg
    neigh(i,i) = 1
    do j = 1, LAnMMP
      if (IAnMMP(i,j) .eq. 0) then
        goto 700
      else
        neigh(IAnMMP(i,j),i) = 2
        neigh(i,IAnMMP(i,j)) = 2
      endif
    enddo

 700 continue
  enddo

  do i = 1, NChg
    do j = 1, NChg
      if (neigh(i,j) .eq. 2) then
        do k = 1, NChg
          njk = neigh(j,k)
          nik = neigh(i,k)
          if (njk.eq.2 .and. nik.ne.1 .and. nik.ne.2) then
            neigh(i,k) = 3
            neigh(k,i) = 3
          endif
        enddo
      endif
    enddo
  enddo

  do i = 1, NChg
    do j = 1, NChg
      if (neigh(i,j) .eq. 3) then
        do k = 1, NChg
          njk = neigh(j,k)
          nik = neigh(i,k)
          if (njk.eq.2 .and. nik.ne.1 .and. nik.ne.2 .and. nik.ne.3) then
            neigh(i,k) = 4
            neigh(k,i) = 4
          endif
        enddo
      endif
    enddo
  enddo
  call gettime('Building neighbourhood matrix')

  if (iprt .ge. 1) then
    call PrtIMat(iout,NChg,NChg,neigh,'Neighbourhood matrix',.false.)
    call gettime('Debug printout')
  endif

  if (iprt.ge.1) call prttime('reading MOL2')

  return
end subroutine
