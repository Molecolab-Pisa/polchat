! prtimat.f90:     A Polarisation consistent charge-fitting tool 
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
subroutine PrtIMat(iout,nr,nc,M,Mes,LTrn)
!
! Print to IOut a (nr,nc) integer matrix, eventually transposed.
!

  logical                   :: LTrn
  integer                   :: i, j, k, nr, nc, ncyc, ext1, ext2
  integer, dimension(nr,nc) :: M
  character(*)              :: Mes

 1000 format(1x,A,' following:')
 1010 format(1x,i6,2x,10(i6,1x))
 1020 format(7x,10(i6,1x))

  write(IOut,1000) Mes

! Print transpose
  if (LTrn) then
    if (nr.le.10) then
      write(IOut,1020) (j, j=1,nr)
      do i = 1, nc
        write(IOut,1010) i,(M(j,i), j=1,nr)
      enddo
    else
      ncyc = nr/10
      if (mod(nr,10).ne.0) ncyc = ncyc+1
      do k = 1, ncyc
        ext1 = (k-1)*10 + 1
        ext2 = min(k*10,nr)
        write(IOut,1020) (j, j=ext1,ext2)
        do i = 1, nc
          write(IOut,1010) i, (M(j,i), j=ext1,ext2)
        enddo
      enddo
    endif

! Print original
  else
    if (nc.le.10) then
      write(IOut,1020) (j, j=1,nc)
      do i = 1, nr
        write(IOut,1010) i,(M(i,j), j=1,nc)
      enddo
    else
      ncyc = nc/10
      if (mod(nc,10).ne.0) ncyc = ncyc+1
      do k = 1, ncyc
        ext1 = (k-1)*10 + 1
        ext2 = min(k*10,nc)
        write(IOut,1020) (j, j=ext1,ext2)
        do i = 1, nr
          write(IOut,1010) i, (M(i,j), j=ext1,ext2)
        enddo
      enddo
    endif
  endif

end subroutine PrtIMat
