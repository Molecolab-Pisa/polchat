subroutine potchg(ng,nq,q,R,V)

  use constants

  implicit real*8(a-h,o-z)

  real*8 :: q(nq), R(ng,nq), V(ng)

  do i = 1, ng
    V(i) = zero
    do j = 1, nq
      V(i) = V(i) + q(j)/R(i,j)
    enddo
  enddo

  if (iprt.ge.2) call PrtMat(iout,ng,1,V,'ESP potential',.false.)

  return

end subroutine

