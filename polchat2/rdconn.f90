subroutine RdConn

  use constants
  use mmpoldata
  use gespinfo
  use operative

  implicit real(a-h,o-z)

  integer, allocatable          :: ib1(:), ib2(:)
  character(len=strmax)         :: what

 1000 format(' Gaussian-style connectivity following:')
 1010 format(i5,' | ',9(i5,1x))
 2000 format(' Reading mol2 file.')
 9000 format(' ERROR',/,&
             ' Mismatch in data read from files:',/,&
             '   Number of atoms read from gesp file: ',i6,/,&
             '   Number of atoms read from mol2 file: ',i6)
 9010 format(' ERROR',/,&
             ' Could not find connectivity information in mol2 file.')
 9020 format(' ERROR',/,&
             ' Maximum number of neighbours exceeded.')
            

  if (iprt .gt. 1) write(iout,2000)
  open(unit=11,file=filecon,status='unknown')

! Skip first 2 lines
  read(11,*)
  read(11,*)
  read(11,*) NAtoms, NBonds

  if (NAtoms .ne. NChg) then
    write(iout,9000) NAtoms, NChg
    stop
  endif
  allocate(ib1(NBonds),ib2(NBonds))

  IAnMMP = 0

! Skip atom section and go to connectivity

  do i = 1, NChg+10
    read(11,*) what
    if (what .eq. '@<TRIPOS>BOND') goto 10
  enddo
  write(iout,9010)
  stop

! Read connectivity and store information on ib1 and ib2

 10 continue
  do i = 1, NBonds
    read(11,*) ixxx, ib1(i), ib2(i)
  enddo
  
  close(11)

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

! Build neighbourhood matrix

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

! Printout
  if (iprt .ge. 0) then
    write(iout,1000)
    do i = 1, NChg
      write(iout,1010) i, (IAnMMP(i,j), j=1,LAnMMP)
    enddo
  endif

  if (iprt .ge. 2) call PrtIMat(iout,NChg,NChg,neigh,'Neighbourhood matrix',.false.)

  return
end subroutine
