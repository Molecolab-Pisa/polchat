subroutine RdPol

  use constants
  use mmpoldata
  use operative

  implicit real(a-h,o-z)

 2000 format(' Reading polarisability file.')

  if (iprt .gt. 1) write(iout,2000)
  open(unit=12,file=filepol,status='unknown')

! Read polarisabilities
  do i = 1, NChg
    read(12,*) pol(i)
  enddo
  close(12)

! Printout
  if (iprt .gt. 0) call PrtMat(iout,NChg,1,pol,' Polarisabilities',.false.)

  return
end subroutine
