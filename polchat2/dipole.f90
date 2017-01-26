subroutine dipole(nq,q,cq,dip)

  use constants

  implicit real*8(a-h,o-z)
  
  real*8 :: q(nq), cq(3,nq), dip(3)

  dip = zero
  do i = 1, nq
    dip(1) = dip(1) + q(i)*cq(1,i)
    dip(2) = dip(2) + q(i)*cq(2,i)
    dip(3) = dip(3) + q(i)*cq(3,i)
  enddo

  return

end subroutine
