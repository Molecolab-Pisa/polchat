module gespinfo

  implicit none

  integer             :: NGrd
  real*8              :: DipQM(3), eGESP
  real*8, allocatable :: gesp(:), CGrd(:,:), VQM(:), RChGr(:,:,:)
  
end module
