subroutine PrtMat(iout,nr,nc,M,Mes,LTrn)
!
! Print to IOut a (nr,nc) real matrix, eventually transposed.
!

  logical                  :: LTrn
  integer                  :: i, j, k, nr, nc, ncyc, ext1, ext2
  real*8, dimension(nr,nc) :: M
  character(*)             :: Mes

 1000 format(1x,A,' following:')
 1010 format(1x,i6,2x,5(d13.6,1x))
 1011 format(1x,i6,2x,5(d16.9,1x))
 1020 format(12x,5(i6,8x))
 1021 format(12x,5(i6,11x))

  write(IOut,1000) Mes

! Print transpose
  if (LTrn) then
    if (nr.le.5) then
      if (iprt.lt.2) then
        write(IOut,1020) (j, j=1,nr)
        do i = 1, nc
          write(IOut,1010) i,(M(j,i), j=1,nr)
        enddo
      else
        write(IOut,1021) (j, j=1,nr)
        do i = 1, nc
          write(IOut,1011) i,(M(j,i), j=1,nr)
        enddo
      endif
    else
      ncyc = nr/5
      if (mod(nr,5).ne.0) ncyc = ncyc+1
      do k = 1, ncyc
        ext1 = (k-1)*5 + 1
        ext2 = min(k*5,nr)
        if (iprt.lt.2) then
          write(IOut,1020) (j, j=ext1,ext2)
          do i = 1, nc
            write(IOut,1010) i, (M(j,i), j=ext1,ext2)
          enddo
        else
          write(IOut,1021) (j, j=ext1,ext2)
          do i = 1, nc
            write(IOut,1011) i, (M(j,i), j=ext1,ext2)
          enddo
        endif
      enddo
    endif

! Print original
  else
    if (nc.le.5) then
      if (iprt.lt.2) then
        write(IOut,1020) (j, j=1,nc)
        do i = 1, nr
          write(IOut,1010) i,(M(i,j), j=1,nc)
        enddo
      else
        write(IOut,1021) (j, j=1,nc)
        do i = 1, nr
          write(IOut,1011) i,(M(i,j), j=1,nc)
        enddo
      endif
    else
      ncyc = nc/5
      if (mod(nc,5).ne.0) ncyc = ncyc+1
      do k = 1, ncyc
        ext1 = (k-1)*5 + 1
        ext2 = min(k*5,nc)
        if (iprt.lt.2) then
          write(IOut,1020) (j, j=ext1,ext2)
          do i = 1, nr
            write(IOut,1010) i, (M(i,j), j=ext1,ext2)
          enddo
        else
          write(IOut,1021) (j, j=ext1,ext2)
          do i = 1, nr
            write(IOut,1011) i, (M(i,j), j=ext1,ext2)
          enddo
        endif
      enddo
    endif
  endif

end subroutine PrtMat
