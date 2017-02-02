subroutine esp

  use constants
  use mmpoldata
  use gespinfo
  use constraints
  use espinfo
  use time

  implicit real*8(a-h,o-z)
  integer, allocatable :: IPIV(:)
  real*8,  allocatable :: WORK(:)

 2000 format(' ESP matrix: Inversion successful.')
 2010 format(' ESP charges computed.')
 2020 format(' Sum of ESP charges: ',f12.6)
 2030 format(' Fit error of ESP charges: ',f12.6)
 2100 format(' Fitting ESP charges.')
 9000 format(' ERROR',/,&
             ' ESP matrix: Element ',i6,' has an illegal value.')
 9010 format(' ERROR',/,&
             ' ESP matrix: Diagonal element ',i6,' is zero.',/,&
             '             The matrix is singular and cannot be inverted.')

  if (iprt.ge.1) write(iout,2100)
  call starttime

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

  call gettime('forming matrices')

  NDim = NChg+NCons

  if (iprt.ge.2) then
    call prtmat(iout,NDim,NDim,X,'X after constraints',.false.)
    call prtmat(iout,NDim,1,B,'B after constraints',.false.)
    call gettime('debug printout')
  endif

! Initialise

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
  elseif (iprt.ge.1) then
    write(iout,2000)
  endif

  call gettime('matrix inversion')
    
! Compute ESP charges and fit error

  allocate (qesp(NChg))
  qesp = matmul(X,B)

  sesp = sum(qesp)
  sini = sum(gesp)
  eesp = error(.false.,qesp)
  egesp = error(.false.,gesp)
  call dipole(NChg,qesp,CChg,desp)
  call dipole(NChg,gesp,CChg,dini)
  if (iprt.ge.2) write(iout,2010)
  if (iprt.ge.2) call PrtMat(iout,NChg,1,qesp,'ESP charges',.false.)
  if (iprt.ge.2) write(iout,2020) sesp
  if (iprt.ge.2) write(iout,2030) eesp

  deallocate(IPIV,WORK)
  deallocate(X,B)

  call gettime('solving for charges and computing fit errors')
  if (iprt.ge.1) call prttime('computing ESP charges')

  return

end subroutine
