subroutine RdCns

  use constants
  use mmpoldata
  use strings
  use constraints
  use operative

  implicit real(a-h,o-z)

  integer                                :: VCns(NChg)
  character*1                            :: del
  character(len=strmax)                  :: line,select
  character(len=strmax), dimension(NChg) :: args
 
 2000 format(' Reading constraint file.')
 2010 format(' Constraint on total charge   : ',f12.4)
 2020 format(' Constraint on fragment charge: ',f12.4,' on atoms:'/,&
             (15x,10(i6,1x)))
 2030 format(' Constraint on equivalence on atoms:'/,&
             (15x,10(i6,1x)))
 2040 format(' Restraint with strength ',f12.4,' of type ',i1,' on atoms:'/,&
             (15x,10(i6,1x)))
 9000 format(' ERROR',/,&
             ' Cannot understand this string from constraint file: ',(A))
 9010 format(' ERROR',/,&
             ' In constraint file expected only one constraint on total charge.')
 9020 format(' ERROR',/,&
             '  In constraint file found more than one restraint.')
 9030 format(' ERROR',/,&
             ' In constraint file found negative restraint strength: ',f12.6)
 9040 format(' ERROR',/,&
             ' In constraint file found unknown restraint type (only 1 and 2 allowed):',i6)
 9100 format(' Found constraint keyword ',A3,' - Adding ',i10,' constraint(s) - Total ',i10)
 9110 format(' NChg=',i6,' NCns=',i6,' Dimension d=',i6,' - X(d,d), B(d)')
 9120 format(' # cns on CHG=',i3,' # cns on FRG=',i3,' # cns on EQV=',i3,' # cns on RES=',i3)
 9200 format(' Parsing line ',(A))

  NCChg = 0
  NCFrg = 0
  NCEqv = 0
  NCRes = 0
  NCons = 0
  del = ' '

  if (iprt .gt. 1) write(iout,2000)
  open(unit=12,file=filecns,status='unknown')
  rewind(12)

! First time just count number of constraints

  do while(.true.)
    read(12,'(A)',END=8000) line
    call parse(trim(line),del,args,nargs)
    select = uppercase(trim(args(1)))
    if (trim(select) .eq. 'CHG') then
      NCChg = NCChg + 1
      NCons = NCons+ 1
      if (iprt.ge.2) write(iout,9100) 'CHG',1,NCons
    elseif (trim(select) .eq. 'FRG') then
      NCFrg = NCFrg + 1
      NCons = NCons+ 1
      if (iprt.ge.2) write(iout,9100) 'FRG',1,NCons
    elseif (trim(select) .eq. 'EQV') then
      call rdrange(args(2),VCns,N)
      NCEqv = NCEqv + 1
      NCons = NCons + N - 1
      if (iprt.ge.2) write(iout,9100) 'EQV',N-1,NCons
    elseif (trim(select) .eq. 'RES') then
      NCRes = NCRes + 1
    else
      write(iout,9000) trim(select)
      stop
    endif
  enddo

  if (iprt.ge.2) then
    write(iout,9110) NChg,NCons,NChg+NCons
    write(iout,9120) NCChg,NCFrg,NCEqv,NCRes
  endif
  
! Consistency check

 8000 continue
  if (NCChg.ne.1) then
    write(iout,9010) 
    stop
  endif
  if (NCRes.gt.1) then
    write(iout,9020) 
    stop
  endif

! Allocate constraints

  allocate ( RCFrg(NCFrg), ICFrg(NCFrg), VCFrg(NCFrg,NChg) )
  allocate ( ICEqv(NCEqv), VCEqv(NCEqv,NChg) )
  allocate ( VCRes(NChg) )

! Read again and record constraint variables and vectors

  rewind(12)
  NCFrg = 0
  NCEqv = 0
  
  do while(.true.)
    read(12,'(A)',END=900) line
    if (iprt.ge.2) write(iout,9200) trim(line)
    call parse(trim(line),del,args,nargs)
    select = uppercase(trim(args(1)))

    if (trim(select) .eq. 'CHG') then
      read(args(2),*) RCChg
      if (iprt.ge.1) write(iout,2010) RCChg

    elseif (trim(select) .eq. 'FRG') then
      NCFrg = NCFrg + 1
      read(args(2),*) RCns
      call rdrange(args(3),VCns,N)
      RCFrg(NCFrg) = RCns
      ICFrg(NCFrg) = N
      VCFrg(NCFrg,1:N) = VCns(1:N)
!     do i = 1, N
!       VCFrg(NCFrg,i) = VCns(i)
!     enddo
      if (iprt.ge.1) write(iout,2020) RCns,(VCns(j),j=1,N)

    elseif (trim(select) .eq. 'EQV') then
      NCEqv = NCEqv + 1
      call rdrange(args(2),VCns,N)
      ICEqv(NCEqv) = N
      VCEqv(NCEqv,1:N) = VCns(1:N)
!     do i = 1, N
!       VCEqv(NCEqv,i) = VCns(i)
!     enddo
      if (iprt.ge.1) write(iout,2030) (VCns(j),j=1,N)

    elseif (trim(select) .eq. 'RES') then
      read(args(2),*) RCns
      read(args(3),*) ICns
      call rdrange(args(4),VCns,N)
      if ( RCns.lt.zero ) then
        write(iout,9030) RCns
        stop
      elseif ( ICns.ne.1 .and. ICns.ne.2 ) then
        write(iout,9040) ICns
        stop
      endif
      if (iprt.ge.1) write(iout,2040) RCns,ICns,(VCns(j),j=1,N)
      ICRes = ICns
      RCRes = RCns
      MCRes = N
      VCRes = VCns 
    endif
  enddo

 900 continue

  close(12)

end subroutine
