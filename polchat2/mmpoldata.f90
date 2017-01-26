module mmpoldata

  implicit none

  logical :: lscr, DoThole
  integer :: NChg, IScreen
  integer, allocatable :: IAnMMP(:,:), neigh(:,:)
  real*8, allocatable  :: pol(:), scrcc(:,:), scrcp(:,:)
  real*8, allocatable  :: CChg(:,:), Rij(:,:,:), Rij3(:,:)
  real*8, allocatable  :: D(:,:)

end module