subroutine RdGESP

  use constants
  use mmpoldata
  use gespinfo
  use operative

  implicit real*8 (a-h,o-z)

 1001 format(t46,i9)
 1002 Format(t9,4(d16.8))
 1003 Format(4(d16.8))
 1004 Format(t50,i8)
 1005 Format(4(d16.8))
 1006 Format(3x,d16.8,3x,d16.8,3x,d16.8)
 2000 format(' Reading gesp file.')
 2010 format(' GESP: Number of atoms:      ',i8)
 2020 format(' GESP: Number of gridpoints: ',i8)

  if (iprt .gt. 1) write(iout,2000)
  open(unit=10,file=filenam,status='unknown')

! Skip first 2 lines and read number of atoms
  read(10,*)
  read(10,*)
  read(10,1001) NChg
  allocate (CChg(3,NChg), gesp(NChg))
  allocate (IAnMMP(NChg,LAnMMP),pol(NChg))
  allocate (scrcc(NChg,NChg),scrcp(NChg,NChg))
  allocate (neigh(NChg,NChg))
  allocate (D(3*NChg,3*NChg))

! Read atom coordinates and GESP charges
  do i = 1, NChg
    read(10,1002) (CChg(j,i), j=1,3), gesp(i)
  enddo

  read(10,*)
  read(10,1006) (DipQM(i),i=1,3)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,1004) NGrd

! Read QM potential and grid coordinates

  allocate (CGrd(3,NGrd),VQM(NGrd))
  do i = 1, NGrd
    read(10,1005) VQM(I), (CGrd(j,i), j=1,3)
  enddo

  close(10)

! Compute dustance between charges and gridpoints

  allocate (RChGr(4,NGrd,NChg))
  call dist(NChg,NGrd,.false.,CChg,CGrd,RChGr,RJunk)

! Compute distance between charges

  allocate (Rij(4,NChg,NChg),Rij3(NChg,NChg))
  call dist(NChg,NChg,.true.,CChg,CChg,Rij,Rij3)

! Printout

  if (iprt.gt.1) then
    write(iout,2010) NChg
    call PrtMat(iout,3,NChg,CChg,' Atom coordinates from GESP file',.true.)
    call PrtMat(iout,NChg,1,gesp,' GESP charges from GESP file',.false.)
    call PrtMat(iout,1,3,DipQM,'Dipole moment from GESP file',.false.)
    write(iout,2020) NGrd
    call PrtMat(iout,3,NGrd,CGrd,' Gridpoint coordinates from GESP file',.true.)
    call PrtMat(iout,NGrd,1,VQM,' ES potential at gridpoints from GESP file',.false.)
    call PrtMat(iout,NGrd,NChg,RChGr(4,:,:),' Charge-gridpoint distance matrix',.false.)
    call PrtMat(iout,NChg,NChg,Rij(4,:,:),' Charge-charge distance matrix',.false.)
  endif

  return
end subroutine


