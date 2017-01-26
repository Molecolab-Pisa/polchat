subroutine printout

  use constants
  use mmpoldata
  use gespinfo
  use espinfo

  implicit real*8 (a-h,o-z)

 1000 format(/,14x,'GESP         ESP        pol-ESP')
 1010 format(1x,i6,3(1x,f12.6))
 1020 format(1x,46('-'))
 1030 format(' Sum    ',3(f12.6,1x))
 1040 format(' Fit Error ',13x,2(es10.3,3x))
 1050 format(/,' Dipole Analysis:',/,22x,'GESP',25x,'ESP',24x,'pol-ESP',/,10x,3(26('-'),3x))
 1060 format('  fixed   ',4(3(f8.4,1x),2x))
 1070 format('  induced ',58x,2(3(f8.4,1x),2x))
 1080 format(10x,3(26('-'),3x),/,'  total   ',3(3(f8.4,1x),2x))

  write(iout,1000)
  write(iout,1020)
  do i = 1, NChg
    write(iout,1010) i, gesp(i), qesp(i), qpesp(i)
  enddo
  write(iout,1020)
  write(iout,1030) sini, sesp, spesp
  write(iout,1020)
  write(iout,1040) eesp, epesp
  write(iout,1020)
  write(iout,1050)
  write(iout,1060) dini(1:3), desp(1:3), dqpesp(1:3)
  write(iout,1070) ddpesp(1:3)
  write(iout,1080) dini(1:3), desp(1:3), dtpesp(1:3)
  
  return

end subroutine
