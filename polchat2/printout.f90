subroutine printout

  use constants
  use mmpoldata
  use gespinfo
  use espinfo
  use operative

  implicit real*8 (a-h,o-z)
  
  integer :: atnum

 1000 format(/,14x,'GESP         ESP        pol-ESP')
 1010 format(1x,i6,3(1x,f12.6))
 1020 format(1x,46('-'))
 1030 format(' Sum    ',3(f12.6,1x))
 1040 format(' Fit Error ',3(es10.3,3x))
 1050 format(/,' Dipole Analysis:',/,1x,69('-'),/,19x,'GESP',9x,'ESP',17x,'pol-ESP',/,&
             6x,'QM calc.',4x,'charges',5x,'charges',5x,'charges',3x,'dipoles',4x,'total',/,&
             5x,3(9('-'),3x),3(9('-'),1x))
 1051 format(1x,'x',4(3x,f9.4),2(1x,f9.4))
 1052 format(1x,'y',4(3x,f9.4),2(1x,f9.4))
 1053 format(1x,'z',4(3x,f9.4),2(1x,f9.4))
 1054 format(1x,69('-'))
 2000 format(' Writing database.')
 2010 format(' Database written.')
 3000 format(1x,(A),1x,(A),2(f12.6,1x),f8.4,1x,i3)

! Print GESP, ESP and pol-ESP to standard output

  write(iout,1000)
  write(iout,1020)
  do i = 1, NChg
    write(iout,1010) i, gesp(i), qesp(i), qpesp(i)
  enddo
  write(iout,1020)
  write(iout,1030) sini, sesp, spesp
  write(iout,1020)
  write(iout,1040) egesp, eesp, epesp
  write(iout,1020)
  write(iout,1050)
  write(iout,1051) dipqm(1), dini(1), desp(1), dqpesp(1), ddpesp(1), dtpesp(1)
  write(iout,1052) dipqm(2), dini(2), desp(2), dqpesp(2), ddpesp(2), dtpesp(2)
  write(iout,1053) dipqm(3), dini(3), desp(3), dqpesp(3), ddpesp(3), dtpesp(3)
  write(iout,1054) 
  
! If required, print database

  if (ldbs) then
    if (iprt.ge.0) write(iout,2000) 
    open(unit=16,file=filedbs,status='unknown')
    do i = 1, NChg
      write(6,*) moltyp(i),atmtyp(i),qesp(i),qpesp(i),pol(i),atnum(atmnam(i))
      write(16,3000) moltyp(i),atmtyp(i),qesp(i),qpesp(i),pol(i),atnum(atmnam(i))
    enddo
    close(16)
    if (iprt.ge.1) write(iout,2010) 
  endif
  
  return

end subroutine
