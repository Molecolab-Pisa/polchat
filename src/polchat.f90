Program polchat
!
! This program computes the ESP charges consistent with an Induced-dipole
! treatment of polarisability.
!
!
IMPLICIT REAL*8 (A-H,O-Z)
DATA Ang2au,Deb2Au/1.889725989d0,0.393430201407683d0/
CHARACTER(50) :: VERSION = "3.1.1 "
REAL*8, ALLOCATABLE :: CChg(:,:),CGrid(:,:),GESP(:),ChgESP(:),pol(:)
REAL*8, ALLOCATABLE :: Vqm(:),X(:,:),B(:),RChg(:,:,:),D(:,:),R(:,:,:)
REAL*8, ALLOCATABLE :: Rij(:,:,:),Rij3(:,:),RChGr(:,:,:), Qesp(:), Qpesp(:)
REAL*8, ALLOCATABLE :: ScrChCh(:,:),ScrChPl(:,:)
REAL*8 :: DipFixIni(3), DipFixESP(3), DipFixPESP(3), DipIndPESP(3), DipPESP(3)
REAL*8 :: DipQM(3) !DipFixOct(3), DipIndOct(3), DipOct(3)
!-debug
!real*8,allocatable :: QOct(:)
INTEGER, ALLOCATABLE :: IAnMMP(:,:)
Parameter LAnMMP = 9
! CChg ....... coordinates of the atoms (a.u.)
! CGrid ...... coordinates of gridpoints (a.u.)
! Gesp ....... Gaussian ESP charge result (for check)
! ChgESP ..... fitted ESP charges 
! pol ........ isotropic polarizabilities (a.u.)
CHARACTER(LEN=50) filename,filepol,fileconn,filecnst !,fileoct
LOGICAL :: LScrChPl
! filename ... name of the Gaussian ESP file
!
! Declaration of vectors for minimization routine:
!
1001 Format(I6,3(F10.5))
1002 Format(' Vqm / V(Gesp) RMS  :',F10.5)
1003 Format(I6,3(F10.5))
1004 Format('Minimization results',/,'IGO   = ',I1,/,'fnorm =',F10.5,/,'SumCh =',F10.5)
1100 Format(' Reading file       : ',A50)
3000 Format(/,' ---- Input files ----------------------------',/,          &
            ' > gesp file              > ',a50,/, &
            ' > mol2 file              > ',a50,/, &
            ' > polarisability file    > ',a50,/, &
            ' > constraints file       > ',a50)
3010 format(/,' ---- Printout -------------------------------',/, &
              ' > standard               > normal printout')
3020 format(/,' ---- Printout -------------------------------',/, &
              ' > debug                  > extra printout')
3030 format(/,' > Wang Chg-Pol screening > active (expert use only)')
3035 format(/,' > Wang Chg-Pol screening > not active')
5110 format(/,' Computing ESP charges...')
5120 format(' Computing pol-ESP charges...',/)
!
! Read input files and printout options from input stream
!
  IOut = 6
  Call PrtHdr(IOut,VERSION)
  Call RdOpts(IOut,IPrint,filename,fileconn,filepol,filecnst,LScrChPl)

!
! Read the Gaussian ESP file produced using the keywords:
! pop=esp Iop(6/50=1)
!
 Call readgesp(IOut,IPrint,filename,CChg,CGrid,Gesp,Vqm,nch,ngr,Rij,Rij3,RChGr,DipQM)
!
! Read the connectivity from a .mol2 file and store it in IAnMMP
! IAnMMP ..... connectivity matrix.
! Then read isotropic polarizabilities from an external file
!
 ALLOCATE(IAnMMP(nch,LAnMMP),pol(nch))
 if (IPrint.ne.0) then
   write(IOut,3000) filename,fileconn,filepol,filecnst
   if (LScrChPl) write(IOut,3030)
   if (.NOT.LScrChPl) write(IOut,3035)
   if (IPrint.eq.1) write(IOut,3010)
   if (IPrint.eq.2) write(IOut,3020)
 endif


 if (IPrint.ge.2) write(IOut,1100) fileconn
 Call readconn(IPrint,IOut,fileconn,nch,IAnMMP,LAnMMP)
 if (IPrint.ge.2) write(IOut,1100) filepol
 OPEN(unit=12,file=filepol,status='unknown')
 DO I=1,nch
   READ(12,*) pol(I)
 ENDDO
 CLOSE(12)
! Wang AL option IMMPCn = 4
  IMMPCn = 4
! Wang DL option IMMPCn = 5 (NYI)
! IMMPCn = 5
!
! *******************************************************************
! Compute charge-charge and charge-dipole screenings
!
  ALLOCATE(ScrChCh(nch,nch), ScrChPl(nch,nch))
  Call MkScr(IOut,IPrint,nch,CChg,Pol,ScrChCh,ScrChPl,IAnMMP,IMMPCn,LAnMMP, &
    LScrChPl)
!
! *******************************************************************
! Compute, invert and store the MMPol matrix
!
 ALLOCATE(D(3*nch,3*nch))
 Call MMPMTRX(IPrint,IOut,nch,pol,CChg,D,IAnMMP,IMMPCn,LAnMMP)
!
! *******************************************************************
! Compute constraint conditions
!
 if (IPrint.ge.2) write(IOut,1100) filecnst
 call FixedDip(IOut,IPrint,nch,CChg,Gesp,DipFixIni)
 call RdConstr(IOut,IPrint,filecnst,nch,mcon,X,B,DipFixIni,CChg,DipQM)
 Allocate (Qesp(nch+mcon))
!
! Compute X and B matrices for ESP, and solve for ESP charges
!
 if (IPrint.ge.2) write(IOut,5110) 
 CALL CPU_TIME(T1)
 call ESP(IOut,IPrint,nch,ngr,mcon,RChGr,Vqm,X,B,Qesp,QSumE)
 deallocate(X,B)
 ErrESP  = FitErr(IOut,IPrint,.false.,nch,ngr,Qesp,RChGr,Rij,Rij3,RJunk,RJunk,Vqm,RJunk)

  call RdConstrPol(IOut,IPrint,filecnst,nch,mcon,X,B,DipFixIni,Rij,Rij3,ScrChPl,D,CChg,DipQM)
  Allocate (Qpesp(nch+mcon))
  if (IPrint.ge.2) write(IOut,5120) 
  call PESP(IOut,IPrint,nch,ngr,mcon,RChGr,Rij,Rij3,Vqm,X,B,ScrChPl,ScrChCh,D,Qpesp,QSumP)
  ErrPESP = FitErr(IOut,IPrint,.true.,nch,ngr,Qpesp,RChGr,Rij,Rij3,D,ScrChPl,Vqm,DipIndPESP)
!
! *******************************************************************
! Compute initial and final dipoles
!
  call FixedDip(IOut,IPrint,nch,CChg,Gesp,DipFixIni)
  call FixedDip(IOut,IPrint,nch,CChg,Qesp,DipFixESP)
  call FixedDip(IOut,IPrint,nch,CChg,Qpesp,DipFixPESP)
  DipPESP = DipFixPESP + DipIndPESP

  write(IOut,6000)
  write(IOut,6020)
  do i = 1, nch
    write(IOut,6010) i,Qesp(i),Qpesp(i) !,QOct(i)
  enddo
  write(IOut,6020)
  write(IOut,6030) QSumE,QSumP !,QSumO
  write(IOut,6020) 
  write(IOut,6040) ErrESP,ErrPESP !,ErrOct
  write(IOut,6020) 
  write(IOut,6050)
  write(IOut,6060) DipFixIni(1:3),DipFixESP(1:3),DipFixPESP(1:3) !,DipFixOct(1:3)
  write(IOut,6070) DipIndPESP(1:3) !,DipIndOct(1:3)
  write(IOut,6080) DipFixIni(1:3),DipFixESP(1:3),DipPESP(1:3) !,DipOct(1:3)
    
 6000 Format(' ESP charges:',/,14x,'ESP        pol-ESP') !    Comparison')
 6010 Format(1x,i6,3(1x,f12.6))
 6020 Format(' --------------------------------------')
 6030 Format(' Total  ',3(f12.6,1x))
 6040 Format(' Fit Error ',3(es10.3,3x))
 6050 Format(/,' Dipole Analysis:',/,20x,'Initial',24x,'ESP',24x,'pol-ESP',/,10x,3(26('-'),3x))
 6060 Format('  fixed   ',4(3(f8.4,1x),2x))
 6070 Format('  induced ',58x,2(3(f8.4,1x),2x))
 6080 Format(10x,3(26('-'),3x),/,'  total   ',3(3(f8.4,1x),2x))
 
!
! *******************************************************************
! Compute the distance between atom and gridpoints and check if
! the Gaussian ESP fitted charges are correct
!
!ALLOCATE (Vch(ngr))
!  Call chpot(nch,ngr,Gesp,R,Vch)
!  WRITE(6,1002) rms(ngr,Vqm,Vch)
!
! *******************************************************************
!
! Compute r_ij and r_ij^3 matrices !! Done already
!
!  Allocate (Rij(4,nch,nch),Rij3(nch,nch))
!  Call MkRij(IOut,IPrint,nch,CChg,Rij,Rij3)
!
! Compute dV/dq
!
!Allocate (DvDq(nch,nch))
!Call MkDvDq(IOut,IPrint,nch,Rij,DvDq)
!  !
!  ! Compute dE/dq
!  !
!  Allocate (DeDq(3,nch))
!  Call MkDeDq(IOut,IPrint,nch,Chg,Rij,Rij3,DeDq)
!  
!  !mcon = 1
!  ldfjac = mcon+ngr
!  ldfj   = mcon+ngr
!  !ALLOCATE (ind(nch+mcon),bl(nch+mcon),bu(nch+mcon),Chg(nch),ChgESP(nch),fjac(ldfjac,nch+1))
!  ALLOCATE (Chg(nch),ChgESP(nch),fjac(ldfjac,nch+1))
!  ! Fuck
!  NT = 5             ! because I'm not using option 15
!  NA = mcon+2*nch+NT ! because I'm not using option 14
!  MX = max(ngr,nch)
!  LIWORK = 2*(3*mcon+9*nch+4*NT+NA+10)
!  LWORK  = 2*(NA*(NA+4)+nch*(NT+33)+(ngr+MX+14)*NT+9*mcon+26)
!  ALLOCATE (iwa(2*LIWORK),wa(2*LWORK),iopt(17),ropt(1))
!  !ind = 4
!  !bl = 0.0d0
!  !bu = 0.0d0
!  !ind(nch+1) = 3
!  !bl(nch+1) = 0.00
!  !bu(nch+1) = 0.00
!  Chg = 0.0d0
!  iopt(01)=99
!  iwa(1) = LWORK
!  iwa(2) = LIWORK
!  Call dqed (dqedev_esp,ngr,nch,mcon,ind,bl,bu,Chg,fjac,ldfjac,fnorm,igo,iopt,ropt,iwa,wa)
!  write(6,*) '******************** BL **********************'
!  write(6,*) bl
!  write(6,*) '******************** BU **********************'
!  write(6,*) bu
!  write(6,*) '******************** IND **********************'
!  write(6,*) ind
!  write(6,*) '******************** FJAC **********************'
!  write(6,*) fjac
!  IGOESP = IGO
!  fnrESP = fnorm
!  ChgESP = Chg
!  
!  write(6,5120) 
!  ALLOCATE(ECh(3,nch),dip(3,nch),Vdip(ngr))
!  Call EChCalc(1,6,nch,ngr,IAnMMP,IMMPCn,RChg,Gesp,ECh)
!  Call DipCalc(1,6,D,nch,ECh,dip)
!  Call dipot(nch,ngr,R,dip,Vdip)
!  ALLOCATE(dVdip(ngr,nch))
!  Call ddipot(nch,ngr,IAnMMP,IMMPCn,RChg,R,D,dVdip)
!  !WRITE(6,*) dVdip
!  
!  Call dqed (dqedev_pesp,ngr,nch,mcon,ind,bl,bu,Chg,fjac,ldfjac,fnorm,igo,iopt,ropt,iwa,wa)
!  !WRITE(IOut,1004) IGO,fnorm,sum(chg)
!  
!  Call chpot(nch,ngr,Chg,R,Vch)
!  Call EChCalc(1,6,nch,ngr,IAnMMP,IMMPCn,RChg,Chg,ECh)
!  Call DipCalc(1,6,D,nch,ECh,dip)
!  Call dipot(nch,ngr,R,dip,Vdip)
!  
!  !
!  ! Printout
!  !
!  write(6,5000) 
!  write(6,5001) IGOESP,IGO,fnrESP,fnorm
!  DO I=1,nch
!    WRITE(6,5002) I,Gesp(I),ChgESP(I),Chg(I),ChgESP(I)-Chg(I)
!  ENDDO
!  write(6,5004) sum(Gesp),sum(ChgESP),sum(Chg)
!  
!  5000 format(/,1x,'       Gaussian      ESP      pol-ESP      diff',/,1x,50('-'))
!  5001 format(' IGO:  ',3x,2(1x,i10),/,' fnorm:',9x,2(1x,f10.5),/,1x,50('-'))
!  5002 format(1x,i4,4(1x,f10.5))
!  5004 format(50('-'),/,' Sum:',4(1x,f10.5))
!  5120 format(' Computing pol-ESP charges...')
Contains
! -------------------------------------------------------------------
!     Subroutine to read Gaussian ESP file
!
Subroutine readgesp(IOut,IPrint,filename,CChg,CGrid,Gesp,Vqm,nch,ngr,   &
  Rij,R3,RChGr,DipQM)
  IMPLICIT REAL*8 (A-H,O-Z)
  REAL*8 , ALLOCATABLE :: CChg(:,:),CGrid(:,:),Gesp(:),Vqm(:),Rij(:,:,:)
  REAL*8 , ALLOCATABLE :: R3(:,:), RChGr(:,:,:)
  REAL*8 :: DipQM(3)
  CHARACTER(LEN=50) filename
  1001 Format(T46,I9)
  1002 Format(T9,4(D16.8))
  1003 Format(4(D16.8))
  1004 Format(T50,I8)
  1005 Format(4(D16.8))
  1006 Format(3x,D16.8,3x,D16.8,3x,D16.8)
  1010 Format(' ---- Parameters -----------------------------',/, &
              ' > number of charges      > ',I8)
  1011 Format(' > grid points            > ',I8)
  OPEN(unit=10,file=filename,status='unknown')
! 
! Skip the first two black lines and read the number of atoms
!
  READ(10,*) 
  READ(10,*) 
  READ(10,1001) nch
! 
! Read the atom coordinates (C) and the calculated 
! Gaussian esp charges (Gesp).
! 
  ALLOCATE (CChg(3,nch),Gesp(nch))
  DO I=1,nch
    READ(10,1002) (CChg(J,I),J=1,3),Gesp(I)
  ENDDO
!
! Skip the first two black lines and read the number of atoms
! 
  READ(10,*) 
  READ(10,1006) (DipQM(i),i=1,3)
  READ(10,*) 
  READ(10,*) 
  READ(10,*) 
  READ(10,1004) ngr
! 
! Read the atom Quantum-mechanical potential (Vqm) and 
! the gridpoint coordinates (CGrid).
! 
  ALLOCATE (CGrid(3,ngr),Vqm(ngr)) 
  DO I=1,ngr
    READ(10,1005) Vqm(I),(CGrid(J,I),J=1,3)
  ENDDO
  CLOSE(10)
!
! Compute the distance between charges and gridpoints
!
  ALLOCATE(RChGr(4,ngr,nch))
  call dist(nch,ngr,.false.,CChg,CGrid,RChGr,RJunk)
!
! Compute the distance between charges
!
  ALLOCATE(Rij(4,nch,nch),R3(nch,nch))
  call dist(nch,nch,.true.,CChg,CChg,Rij,R3)
!
  WRITE(IOut,1010) nch
  WRITE(IOut,1011) ngr

  RETURN
End Subroutine
! -------------------------------------------------------------------
Subroutine RdConstr(IOut,IPrint,filecnst,nch,mcon,X,B,ESPDip,CChg,QMDip)
!
! Read the constraints
! 
implicit real*8 (A-H,O-Z)
character(LEN=50) info,filecnst
Real*8, Allocatable :: X(:,:),B(:)
Real*8 :: RdDip(3), ESPDip(3), QMDip(3), CChg(3,nch)
integer, allocatable :: temp(:)
Logical :: LDoDip
!
1000 format(/,' ---- Constraints ----------------------------')
1010 format(' > total charge           > ',f6.3)
1020 format(' > fragment               > charge ',f6.3,'    > atoms ',(15(1x,i4)))
1030 format(' > equivalence            >                  > atoms ',(15(1x,i4)))
1040 format(' > total dipole           > read input       > ',3(f7.4,1x))
1050 format(' > total dipole           > get from ESP     > ',3(f7.4,1x))
1055 format(' > total dipole           > use QM dipole    > ',3(f7.4,1x))
1060 format(' > total dipole           > no constraint')
1070 format(/)
1100 format(' ERROR IN FILE')
2000 format(' X matrix (constraints only):')
2010 format(' B vector (constraints only):')
2020 format(1x,i6,15(1x,f10.4))

open(unit=15,file=filecnst,status='unknown')
if (IPrint.ge.1) write(IOut,1000)
rewind(15)
! Just loop and count number of constraints
read(15,*,END=8000) info
mcon = 1   ! There will be 1 constraint on the total charge 
           ! only---
LDoDip = .false.
do while (1.eq.1)
  read(15,*,END=9000) info
  if (info.eq.'fragm') then    
    read(15,*,END=8000) info
    read(15,*,END=8000) info
    mcon = mcon + 1
  elseif (info.eq.'equiv') then
    read(15,*,END=8000) it
    read(15,*,END=8000) info
    mcon = mcon + it - 1
  elseif (info.eq.'dipole') then
    LDoDip = .true.
    mcon = mcon + 3
    read(15,*,END=8000) info
    if (info.ne.'esp'.and.info.ne.'qm') read(15,*,END=8000) info !read extra line
  else
    goto 8000
  endif
enddo
! Errors go here
8000 write(IOut,1100)
stop
! Continue here: allocate constraints
9000 continue
allocate(temp(nch),X(nch+mcon,nch+mcon),B(nch+mcon))
X = 0.0d0
B = 0.0d0
! Read again: fill in vectors
rewind(15)
icon = 1
!    Total charge
read(15,*) ch
X(nch+icon,1:nch) = 1.0d0
X(1:nch,nch+icon) = 1.0d0
B(nch+icon) = ch
if (IPrint.ge.1) write(IOut,1010) ch
!
icon = icon + 1
do while (1.eq.1)
  read(15,*,END=9100) info
  if (info.eq.'fragm') then
    read(15,*) it,ch
    read(15,*) (temp(i),i=1,it)
    do i = 1, it
      X(nch+icon,temp(i)) = 1.0d0
      X(temp(i),nch+icon) = 1.0d0
      B(nch+icon) = ch
    enddo
    icon = icon+1
    if (IPrint.ge.1) write(IOut,1020) ch,(temp(i),i=1,it)
  elseif (info.eq.'equiv') then
    read(15,*) it
    read(15,*) (temp(i),i=1,it)
    if (IPrint.ge.1) write(IOut,1030) (temp(i),i=1,it)
    do i = 1, it-1
      X(nch+icon,temp(i))  = -1.0d0
      X(nch+icon,temp(i+1)) = 1.0d0
      X(temp(i),nch+icon)  = -1.0d0
      X(temp(i+1),nch+icon) = 1.0d0
      icon = icon+1
    enddo
  elseif (info.eq.'dipole') then 
    read(15,*) info
    if (info.eq.'esp') then
      RdDip = ESPDip
      if (IPrint.ge.1) write(IOut,1050) (RdDip(i), i=1,3)
    elseif (info.eq.'qm') then
      RdDip = QMDip
      if (IPrint.ge.1) write(IOut,1055) (RdDip(i), i=1,3)
    elseif (info.eq.'read') then
      read(15,*) (RdDip(i), i=1,3)
      if (IPrint.ge.1) write(IOut,1040) (RdDip(i), i=1,3)
    else
      goto 8000
    endif
  endif
enddo
! Continue here
9100 continue
if (.not.LDoDip .and. IPrint.ge.1) write(IOut,1060)
if (LDoDip) then
  nStr = nch+mcon-2
  nEnd = nch+mcon
  X(nStr:nEnd,1:nch) = CChg
  X(1:nch,nStr:nEnd) = Transpose(CChg)
  B(nStr:nEnd) = RdDip
endif

if (IPrint.ge.1) write(IOut,1070)
!  ---
if (IPrint.ge.2) then
  write(IOut,2000)
  do i = 1, nch+mcon
    write(IOut,2020) i,X(i,:)
  enddo
  write(IOut,2010)
  do i = 1, nch+mcon
    write(IOut,2020) i,B(i)
  enddo
endif
! End
close(15)
return
end subroutine
! -------------------------------------------------------------------
Subroutine RdConstrPol(IOut,IPrint,filecnst,nch,mcon,X,B,ESPDip,Rij,Rij3,Scr,DInv,CChg,QMDip)
!
! Read the constraints
! 
implicit real*8 (A-H,O-Z)
character(LEN=50) info,filecnst
Real*8, Allocatable :: X(:,:),B(:)
Real*8 :: RdDip(3), ESPDip(3), Rij(4,nch,nch), Rij3(nch,nch), QMDip(3)
Real*8 :: Scr(nch,nch), DInv(3*nch,3*nch), CChg(3,nch)
Real*8 :: Bx(3*nch,nch), Cx(3*nch,nch), Fx(nch,3), Gx(nch,3)
integer, allocatable :: temp(:)
Logical :: LDoDip
Parameter Small=1.0d-10
!
1000 format(/,'    Summary of constraints')
1010 format(' Total charge [',f6.3,']')
1020 format(' Fragment     [ch ',f6.3,'] :',(15(1x,i4)))
1030 format(' Equivalence:            atoms:',(15(1x,i4)))
1040 format(' Total dipole: read input    [',3(f7.4,1x),']')
1050 format(' Total dipole: get from ESP  [',3(f7.4,1x),']')
1055 format(' Total dipole: use QM dipole [',3(f7.4,1x),']')
1060 format(' Total dipole: no constraint')
1100 format(' ERROR IN FILE')
2000 format(' X matrix (constraints only):')
2010 format(' B vector (constraints only):')
2020 format(1x,i6,15(1x,f10.4))

open(unit=15,file=filecnst,status='unknown')
if (IPrint.ge.2) write(IOut,1000)
rewind(15)
! Just loop and count number of constraints
read(15,*,END=8000) info
mcon = 1   ! There will be 1 constraint on the total charge 
           ! only--- Constraint on total dipole will be im-
           ! posed only if dipole keyword is read.
           ! Subkeywords: esp :: use total dipole from ESP
           !              read x y z :: read dipole constr.
LDoDip = .false.
do while (1.eq.1)
  read(15,*,END=9000) info
  if (info.eq.'fragm') then    
    read(15,*,END=8000) info
    read(15,*,END=8000) info
    mcon = mcon + 1
  elseif (info.eq.'equiv') then
    read(15,*,END=8000) it
    read(15,*,END=8000) info
    mcon = mcon + it - 1
  elseif (info.eq.'dipole') then
    LDoDip = .true.
    mcon = mcon + 3
    read(15,*,END=8000) info
    if (info.ne.'esp'.and.info.ne.'qm') read(15,*,END=8000) info !read extra line
  else
    goto 8000
  endif
enddo
! Errors go here
8000 write(IOut,1100)
stop
! Continue here: allocate constraints
9000 continue
allocate(temp(nch),X(nch+mcon,nch+mcon),B(nch+mcon))
X = 0.0d0
B = 0.0d0
! Read again: fill in vectors
rewind(15)
icon = 1
!    Total charge
read(15,*) ch
X(nch+icon,1:nch) = 1.0d0
X(1:nch,nch+icon) = 1.0d0
B(nch+icon) = ch
if (IPrint.ge.2) write(IOut,1010) ch
!
icon = icon + 1
do while (1.eq.1)
  read(15,*,END=9100) info
  if (info.eq.'fragm') then
    read(15,*) it,ch
    read(15,*) (temp(i),i=1,it)
    do i = 1, it
      X(nch+icon,temp(i)) = 1.0d0
      X(temp(i),nch+icon) = 1.0d0
      B(nch+icon) = ch
    enddo
    icon = icon+1
    if (IPrint.ge.2) write(IOut,1020) ch,(temp(i),i=1,it)
  elseif (info.eq.'equiv') then
    read(15,*) it
    read(15,*) (temp(i),i=1,it)
    if (IPrint.ge.2) write(IOut,1030) (temp(i),i=1,it)
    do i = 1, it-1
      X(nch+icon,temp(i))  = -1.0d0
      X(nch+icon,temp(i+1)) = 1.0d0
      X(temp(i),nch+icon)  = -1.0d0
      X(temp(i+1),nch+icon) = 1.0d0
      icon = icon+1
    enddo
  elseif (info.eq.'dipole') then
    read(15,*) info
    if (info.eq.'esp') then
      RdDip = ESPDip
      if (IPrint.ge.2) write(IOut,1050) (RdDip(i), i=1,3)
    elseif (info.eq.'qm') then
      RdDip = QMDip
      if (IPrint.ge.2) write(IOut,1055) (RdDip(i), i=1,3)
    elseif (info.eq.'read') then
      read(15,*) (RdDip(i), i=1,3)
      if (IPrint.ge.2) write(IOut,1040) (RdDip(i), i=1,3)
    else
      goto 8000
    endif
  endif
enddo
! Continue here
9100 continue

! Take care of dipole constraint
if (.not.LDoDip) then
  if (IPrint.ge.2) write(IOut,1060)
else
! Calculate F elements
!   --sanitise Rij3
  do i = 1, nch
    do j = 1, nch
      if (Rij3(i,j).lt.Small) Rij3(i,j)=Small
    enddo
  enddo
!   --compute Bx
  Bx(1:nch,:)         = Rij(1,:,:)*(Scr/(Rij3))
  Bx(nch+1:2*nch,:)   = Rij(2,:,:)*(Scr/(Rij3))
  Bx(2*nch+1:3*nch,:) = Rij(3,:,:)*(Scr/(Rij3))
!   --compute Cx
  Cx = Matmul(DInv,Bx)
!   --compute Fx
  do i = 1, nch
    do k = 1, 3
      nStr = (k-1)*nch + 1
      nEnd = k*nch
      Fx(i,k) = sum(Cx(nStr:nEnd,i))
!err      nRow = (k-1)*nch + i
!err      Fx(i,k) = sum(Cx(nRow,:))
    enddo
  enddo
!   --compute Gx
  Gx = Fx + Transpose(CChg)
!   --fill X and B
  nStr = nch+mcon-2
  nEnd = nch+mcon
  X(nStr:nEnd,1:nch) = Transpose(Gx)
  X(1:nch,nStr:nEnd) = Gx
  B(nStr:nEnd) = RdDip
!  ---
endif
if (IPrint.ge.2) then
  write(IOut,2000)
  do i = 1, nch+mcon
    write(IOut,2020) i,X(i,:)
  enddo
  write(IOut,2010)
  do i = 1, nch+mcon
    write(IOut,2020) i,B(i)
  enddo
endif
! End
close(15)
return
end subroutine
END PROGRAM polchat
