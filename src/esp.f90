!
! -------------------------------------------------------------------
!
Subroutine ESP(IOut,IPrint,npol,ngr,mcon,RChGr,Vqm,X,B,Qesp,QSum)
  Implicit Real*8 (A-H,O-Z)
  Dimension RChGr(4,ngr,npol), Vqm(ngr), X(npol+mcon,npol+mcon), B(npol+mcon)
  Dimension Qesp(npol+mcon)
  Integer, Allocatable :: IPIV(:)
  Real*8, Allocatable  :: WORK(:)
 1000 Format(' ESP matrix: Element ',i6,' has an illegal value',/,' Halting.')
 1010 Format(' ESP matrix: Diagonal element ',i6,' is zero.',/, &
             ' The matrix is singular and cannot be inverted',/,' Halting.')
 1020 Format(' ESP matrix: Inversion successful.')
 1100 Format(' ESP charges:')
 1110 Format(1x,i6,1x,f12.6)
 1120 Format(' --------------------',/,' Total  ',f12.6)

! Complete matrices X and B (constraints already done)

  do k = 1, npol
    bb = 0.0d0
    do i = 1, ngr
      bb = bb + Vqm(i)/RChGr(4,i,k)
    enddo
    B(k) = bb
    do j = 1, npol
      xx = 0.0d0
      do i = 1, ngr
        xx = xx + 1.0d0/(RChGr(4,i,j)*RChGr(4,i,k))
      enddo
      X(j,k) = xx
    enddo
  enddo

! Solve system to find ESP charges: 
!  (-) initialise elements

  LWORK = (npol+mcon)*(npol+mcon)
  allocate (IPIV(npol+mcon),WORK(LWORK))
  
!  (-) invert matrix X

  M = npol+mcon
  Call DGETRF(M,M,X,M,IPIV,INFO)
  Call DGETRI(M,X,M,IPIV,WORK,LWORK,INFO)

  if (INFO .lt. 0) then
    write(IOut,1000) abs(INFO)
    stop
  elseif (INFO .gt. 0) then
    write(IOut,1010) INFO
    stop
  elseif (IPrint.ge.2) then
    write(IOut,1020)
  endif

!   (-) Compute charges

  qESP = matmul(X,B)

  if (IPrint.ge.2) write(IOut,1100)
  qSum = 0.d0
  do i = 1, npol
    qSum = qSum + qESP(i)
    if (IPrint.ge.2) write(IOut,1110) i,qESP(i)
  enddo
  if (IPrint.ge.2) write(IOut,1120) qSum

 Deallocate(IPIV,WORK)

  Return
End Subroutine
