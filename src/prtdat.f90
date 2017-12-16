! prtdat.f90:      A Polarisation consistent charge-fitting tool 
!                  A Molecolab Tool www.molecolab.dcci.unipi.it/tools
!
! Copyright (C) 2014, 2015, 2016, 2017
!   S. Caprasecca, C. Curutchet, B. Mennucci
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
subroutine PrtDat

  use constants
  use mmpoldata
  use gespinfo
  use constraints
  use operative

 1000 format(/,' PARAMETERS ------------------------------------------------')
 1010 format(  '   > number of charges             > ',i6,/,&
               '   > number of gridpoints          > ',i6)
 1100 format(/,' FILES -----------------------------------------------------')
 1110 format(  '   > GESP file                     > '(A),/,&
               '   > mol2 file                     > '(A),/,&
               '   > polarisability file           > '(A),/,&
               '   > constraint file               > '(A))
 1120 format(  '   > database file (output)        > '(A))
 1200 format(/,' CHARGE-DIPOLE SCREENING -----------------------------------')
 1210 format(  '                                   > ON')
 1220 format(  '                                   > OFF')
 1300 format(/,' PRINTING --------------------------------------------------')
 1310 format(  '   > printout level                > -1 (silent)')
 1320 format(  '   > printout level                >  0 (default)')
 1330 format(  '   > printout level                >  1 (debug)')
 1340 format(  '   > printout level                >  2 (super-debug)')
 1350 format(  '   > Gaussian-style input          > print in log file')
 1355 format(  '   > Gaussian-style input          > do not print')
 1360 format(  '   > connectivity shift            > ',i6)
 1370 format(  '   > residue number                > ',i6)
 1400 format(/,' CONSTRAINTS -----------------------------------------------')
 1410 format(  '   > whole molecule                > charge ',f9.4)
 1420 format(  '   > fragment                      > charge ',f9.4,/,&
               '                                   > atoms ',10(i6,1x),(44x,10(i6,1x)))
 1430 format(  '   > equivalence                   > atoms ',10(i6,1x),(44x,10(i6,1x)))
 1440 format(  '   > restraint                     > type ',i1,' > alpha ',f8.4,/,&
               '                                   > atoms ',10(i6,1x),(44x,10(i6,1x)))
 1500 format(/,' Gaussian-style connectivity following:')
 1510 format(1x,i5,1x,'|',9(1x,i5))

  if (iprt.lt.0) return

! Print working parameters

  write(iout,1000) 
  write(iout,1010) NChg,NGrd

  write(iout,1100)
  write(iout,1110) trim(filenam),trim(filecon),trim(filepol),trim(filecns)
  if (ldbs) write(iout,1120) trim(filedbs)

  write(iout,1200)
  if (lscr) then
    write(iout,1210)
  else
    write(iout,1220)
  endif

  write(iout,1300)
  select case (iprt)
    case (-1)
      write(iout,1310) 
    case (0)
      write(iout,1320) 
    case (1)
      write(iout,1330) 
    case (2)
      write(iout,1340) 
  end select
  if (lgau) then
    write(iout,1350) 
    write(iout,1360) shiftc
    write(iout,1370) resnum
  else 
    write(iout,1355) 
  endif

  write(iout,1400)
  write(iout,1410) RCChg
  do i = 1, NCFrg
    write(iout,1420) RCFrg(i), (VCFrg(i,j), j=1,ICFrg(i))
  enddo
  do i = 1, NCEqv
    write(iout,1430) (VCEqv(i,j), j=1,ICEqv(i))
  enddo
  if (NCRes.ne.0) then
    write(iout,1440) ICRes, RCRes, (VCRes(j), j=1,MCRes)
  endif

! Print connectivity

  if (iprt.ge.0) then
    write(iout,1500)
    do i = 1, NChg
      write(iout,1510) i, (IAnMMP(i,j), j=1,LAnMMP)
    enddo
  endif

  return

end subroutine
