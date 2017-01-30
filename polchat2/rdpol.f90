subroutine RdPol

  use constants
  use mmpoldata
  use operative
  use time

  implicit real(a-h,o-z)

 2000 format(' Reading polarisability file.')

  call starttime 

  if (iprt .gt. 1) write(iout,2000)
  open(unit=12,file=filepol,status='unknown')

! Read polarisabilities
  do i = 1, NChg
    read(12,*) pol(i)
  enddo
  close(12)

! Printout
  if (iprt .ge. 2) call PrtMat(iout,NChg,1,pol,'Polarisabilities',.false.)
  
  call gettime('')
  if (iprt.ge.1) call prttime('reading pol file')

  return
end subroutine
