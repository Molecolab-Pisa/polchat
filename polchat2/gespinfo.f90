module gespinfo

  implicit none

  integer             :: NGrd
  real*8              :: DipQM(3)
  real*8, allocatable :: gesp(:), CGrd(:,:), VQM(:), RChGr(:,:,:)
  
end module
