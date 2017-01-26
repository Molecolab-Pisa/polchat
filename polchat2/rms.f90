real*8 function rms(N,QRef,Q)

  use constants

  implicit real*8(a-h,o-z)

  dimension QRef(N), Q(N)

  rms = zero

  do i = 1, N
    rms = rms + (QRef(i)-Q(i))**2
  enddo

  return

end function
