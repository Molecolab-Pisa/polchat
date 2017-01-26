module constraints

  implicit none

  integer              :: NCons, NCChg, NCFrg, NCEqv, NCRes
  integer              :: ICRes, MCRes
  integer, allocatable :: ICFrg(:), ICEqv(:), VCFrg(:,:), VCEqv(:,:), VCRes(:)
  real*8               :: RCChg, RCRes
  real*8, allocatable  :: RCFrg(:)

  real*8, allocatable  :: X(:,:), B(:)

end module
