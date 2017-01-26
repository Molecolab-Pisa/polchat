subroutine rdrange(string,V,icount)

  use constants, only : iout,strmax
  use mmpoldata, only : NChg
  use strings

  implicit real*8(a-h,o-z)

  integer                                :: V(NChg)
  character(len=1)                       :: comma=',', dash='-'
  character(len=strmax)                  :: string
  character(len=strmax), dimension(NChg) :: rng
  character(len=strmax), dimension(3)    :: ext

 9000 format(' ERROR',/,&
             ' Range values are messed up in constraint file.')
 9010 format(' ERROR',/,&
             ' Range values are messed up in constraint file.',/,&
             ' Found value: ',(A))
 9020 format(' ERROR',/,&
             ' Range values are messed up in constraint file.',/,&
             ' Found interval: ',(A),'-',(A))
 9030 format(' WARNING',/,&
             ' In constraint file found interval: ',i6,'-',i6,/,&
             ' Will nterpret it as interval:      ',i6,'-',i6)
 9040 format(' WARNING',/,&
             ' In constraint file found interval: ',i6,'-',i6,/,&
             ' Will nterpret it as single value:  ',i6)

  V = 0
  icount = 0

  call parse(trim(string),comma,rng,nrng)
  
  do irng = 1, nrng
    call parse(trim(rng(irng)),dash,ext,next)
    select case (next)
      case(1)
        call value(trim(ext(1)),intext1,ios1)
        if (ios1.ne.0) then
          write(iout,9010) trim(ext(1))
          stop
        endif
        icount = icount + 1
        V(icount) = intext1
      case(2)
        call value(trim(ext(1)),intext1,ios1)
        call value(trim(ext(2)),intext2,ios2)
        if (ios1.ne.0 .or. ios2.ne.0) then
          write(iout,9020) trim(ext(1)),trim(ext(2))
          stop
        elseif (intext1.gt.intext2) then
          write(iout,9030) intext1,intext2,intext2,intext1
          itemp   = intext1
          intext1 = intext2
          intext2 = itemp
        elseif (intext1.eq.intext2) then
          write(iout,9040) intext1,intext1,intext1
        endif
        do ii = intext1,intext2
          icount = icount + 1
          V(icount) = ii
        enddo
      case default
        write(iout,9000)
        stop
    end select
  enddo

  return

end subroutine
