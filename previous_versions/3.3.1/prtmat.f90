  subroutine PrtMat(IOut,nr,nc,M,Mes)

    integer                  :: i, j, k, nr, nc, ncyc, ext1, ext2
    real*8, dimension(nr,nc) :: M
    character(*) Mes

 1000 format(1x,A,' following:')
 1010 format(1x,i6,2x,5(d16.9,1x))
 1020 format(12x,5(i6,11x))

    write(IOut,1000) Mes
    if (nc.le.5) then
      write(IOut,1020) (j, j=1,nc)
      do i = 1, nr
        write(IOut,1010) i,(M(i,j), j=1,nc)
      enddo
    else
      ncyc = nc/5
      if (mod(nc,5).ne.0) ncyc = ncyc+1
      do k = 1, ncyc
        ext1 = (k-1)*5 + 1
        ext2 = min(k*5,nc)
        write(IOut,1020) (j, j=ext1,ext2)
        do i = 1, nr
          write(IOut,1010) i, (M(i,j), j=ext1,ext2)
        enddo
      enddo
    endif

  end subroutine PrtMat
