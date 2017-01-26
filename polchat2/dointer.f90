logical function DoInter(i,j)

  use constants
  use mmpoldata

  implicit real*8(a-h,o-z)
  logical ll

 9000 format(' ERROR',/,&
             ' Possible mismatch in connectivity matrix.')

  DoInter = .true.

! Groups option
  if (IMMPCn.eq.2.or.IMMPCn.eq.3) then
    IinJ = 0
    JinI = 0
    do k = 1, LAnMMP
      if (IAnMMP(1,i).eq.IAnMMP(k,j)) IinJ = 1
      if (IAnMMP(k,i).eq.IAnMMP(1,j)) JinI = 1
    enddo
    if (IinJ.ne.JinI) then
      write(iout,9000)
      stop
    elseif (IinJ.eq.1) then
      DoInter = .false.
    endif

! 1-2, 1-3 neighbour
  elseif (IMMPCn.eq.1.or.IMMPCn.eq.4) then
    if (neigh(i,j) .le. 3 .and. neigh(i,j) .ge. 1) DoInter = .false.
  endif

  return

end function

