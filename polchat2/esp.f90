subroutine esp

  use constants
  use mmpoldata
  use gespinfo
  use constraints
  use espinfo

  implicit real*8(a-h,o-z)
  integer, allocatable :: IPIV(:)
  real*8,  allocatable :: WORK(:)

 2000 format(' ESP matrix: Inversion successful.')
 2010 format(' ESP charges computed.')
 2020 format(' Sum of ESP charges: ',f12.6)
 2030 format(' Fit error of ESP charges: ',f12.6)
 9000 format(' ERROR',/,&
             ' ESP matrix: Element ',i6,' has an illegal value.')
 9010 format(' ERROR',/,&
             ' ESP matrix: Diagonal element ',i6,' is zero.',/,&
             '             The matrix is singular and cannot be inverted.')

! Complete matrices X and B

  do k = 1, NChg
    bb = zero
    do i = 1, NGrd
      bb = bb + VQM(i)/RChGr(4,i,k)
    enddo
    B(k) = bb
    do j = 1, NChg
      xx = zero
      do i = 1, NGrd
        xx = xx + one/(RChGr(4,i,j)*RChGr(4,i,k))
      enddo
      X(j,k) = xx
    enddo
  enddo

! Add restraints

  if (NCRes.ne.0 .and. RCRes.gt.small) then
    do k = 1, MCRes
      X(VCRes(k),VCRes(k)) = X(VCRes(k),VCRes(k)) + RCRes
    enddo
  endif

! Initialise elements

  NDim = NChg+NCons
  LWORK = NDim**2
  allocate (IPIV(NDim), WORK(LWORK))
    
! Invert X

  INFO = 0
  call DGETRF(NDim,NDim,X,NDim,IPIV,INFO)
  call DGETRI(NDim,X,NDim,IPIV,WORK,LWORK,INFO)

  if ( INFO .lt. 0 ) then
    write(iout,9000) abs(INFO)
    stop
  elseif ( INFO .gt. 0 ) then
    write(iout,9010) INFO
    stop
  elseif (iprt.ge.2) then
    write(iout,2000)
  endif
    
! Compute ESP charges and fit error

  allocate (qesp(NChg))
  qesp = matmul(X,B)
  sesp = sum(qesp)
  sini = sum(gesp)
  eesp = error(.false.)
  call dipole(NChg,qesp,CChg,desp)
  call dipole(NChg,gesp,CChg,dini)
  if (iprt.ge.2) write(iout,2010)
  if (iprt.ge.2) call PrtMat(iout,NChg,1,qesp,'ESP charges',.false.)
  if (iprt.ge.2) write(iout,2020) sesp
  if (iprt.ge.2) write(iout,2030) eesp

  deallocate(IPIV,WORK)
  deallocate(X,B)

  return

end subroutine
