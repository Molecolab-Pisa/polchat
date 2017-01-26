module espinfo

  implicit none

  real*8, allocatable :: qESP(:), qpesp(:)
  real*8              :: eESP, ePESP
  real*8              :: sini, sESP, sPESP
  real*8              :: dini(3), dESP(3), dqPESP(3), ddPESP(3), dtPESP(3)

end module
