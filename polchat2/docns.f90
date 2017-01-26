subroutine DoCns

  use constants
  use mmpoldata
  use strings
  use constraints

  implicit real(a-h,o-z)

! Allocate X and B vectors

  NDim = NChg+NCons
  allocate(X(NDim,NDim), B(NDim))

  X  = zero
  B  = zero
  II = 0

! Read again and fill vectors

! 1. CHG

  II = II + 1 
  X(NChg+II,1:NChg) = one
  X(1:NChg,NChg+II) = one
  B(NChg+II) = RCChg

! 2. FRG
 
  do i = 1, NCFrg
    II = II + 1
    do j = 1, ICFrg(i)
      X(NChg+II,VCFrg(i,j)) = one   
      X(VCFrg(i,j),NChg+II) = one
      B(NChg+II) = RCFrg(i)
    enddo
  enddo

! 3. EQV

  do i = 1, NCEqv
    do j = 1, ICEqv(i)-1
      II = II + 1
      X(NChg+II,VCEqv(i,j))   = -one
      X(NChg+II,VCEqv(i,j+1)) =  one
      X(VCEqv(i,j),NChg+II)   = -one
      X(VCEqv(i,j+1),NChg+II) =  one
    enddo
  enddo

! 4. RES does not add constraints

  return
 
end subroutine
