Subroutine MkRij(IOut,IPrint,nch,CChg,Rij,Rij3)
  Implicit Real*8 (A-H,O-Z)
  Real*8 :: CChg(3,nch), Rij(4,nch,nch), Rij3(nch,nch)
  1000 format(' Distance between charges:')
  1010 format(2(1x,i6),1x,f12.6)

  Rij = 0.0d0
  do i = 1, nch
    do j = i+1, nch
      dx  = CChg(1,i) - CChg(1,j)
      dy  = CChg(2,i) - CChg(2,j)
      dz  = CChg(3,i) - CChg(3,j)
      dd  = sqrt(dx*dx + dy*dy + dz*dz)
      dd3 = dd**3
      Rij(1,i,j) = dx
      Rij(1,j,i) = -dx
      Rij(2,i,j) = dy
      Rij(2,j,i) = -dy
      Rij(3,i,j) = dz
      Rij(3,j,i) = -dz
      Rij(4,i,j) = dd
      Rij(4,j,i) = dd
      Rij3(i,j)  = dd3
      Rij3(j,i)  = dd3
    enddo
  enddo
  
  if (IPrint.ge.2) then
    write(IOut,1000)
    do i = 1, nch
      do j = i+1, nch
        write(IOut,1010) i,j,Rij(4,i,j)
      enddo
    enddo
  endif
  
  Return
End Subroutine
