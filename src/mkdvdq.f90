Subroutine MkDvDq(IOut,IPrint,nch,Rij,DvDq)
  Implicit Real*8 (A-H,O-Z)
  Real*8 :: Rij(4,nch,nch), DvDq(nch,nch)
  1000 format(' dV/dq:')
  1010 format(2(1x,i6),1x,f12.6)

  DvDq = 0.0d0
  do i = 1, nch
    do j = i+1, nch
      r = Rij(4,i,j)
      DvDq(i,j) = 1.0d0/r
      DvDq(j,i) = 1.0d0/r
    enddo
  enddo

  if (IPrint.ge.2) then
    write(IOut,1000)
    do i = 1, nch
      do j = i+1, nch
        write(IOut,1010) i,j,DvDq(i,j)
      enddo
    enddo
  endif

  Return
End Subroutine
