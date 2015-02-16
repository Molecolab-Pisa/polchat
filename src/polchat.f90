Program polchat
!
! This program computes the ESP charges consistent with an Induced-dipole
! treatment of polarisability.
!
!
IMPLICIT REAL*8 (A-H,O-Z)
DATA Ang2au,Deb2Au/1.889725989d0,0.393430201407683d0/
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
  Call PrtHdr(IOut)
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
 6040 Format(' Fit Error  ',3(es9.3,4x))
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
!------------------------------------------------------------------------------------   
!------------------------------------------------------------------------------------   
!------------------------------------------------------------------------------------   
Subroutine MkRij(IOut,IPrint,nch,CChg,Rij,Rij3)
  Implicit Real*8 (A-H,O-Z)
  Real*8 :: CChg(3,nch), Rij(4,nch,nch), Rij3(nch,nch)
  1000 format(' Distance between charges:')
  1010 format(2(1x,i6),1x,f12.6)

  Rij = 0.0d0
  do i = 1, nch
    do j = i+1, nch
      dx  = CChg(1,i) - CChg(1,j)
      dy  = CChg(2,i) - CChg(2,j)
      dz  = CChg(3,i) - CChg(3,j)
      dd  = sqrt(dx*dx + dy*dy + dz*dz)
      dd3 = dd**3
      Rij(1,i,j) = dx
      Rij(1,j,i) = -dx
      Rij(2,i,j) = dy
      Rij(2,j,i) = -dy
      Rij(3,i,j) = dz
      Rij(3,j,i) = -dz
      Rij(4,i,j) = dd
      Rij(4,j,i) = dd
      Rij3(i,j)  = dd3
      Rij3(j,i)  = dd3
    enddo
  enddo
  
  if (IPrint.ge.2) then
    write(IOut,1000)
    do i = 1, nch
      do j = i+1, nch
        write(IOut,1010) i,j,Rij(4,i,j)
      enddo
    enddo
  endif
  
  Return
End Subroutine
!------------------------------------------------------------------------------------   
Subroutine MkDvDq(IOut,IPrint,nch,Rij,DvDq)
  Implicit Real*8 (A-H,O-Z)
  Real*8 :: Rij(4,nch,nch), DvDq(nch,nch)
  1000 format(' dV/dq:')
  1010 format(2(1x,i6),1x,f12.6)

  DvDq = 0.0d0
  do i = 1, nch
    do j = i+1, nch
      r = Rij(4,i,j)
      DvDq(i,j) = 1.0d0/r
      DvDq(j,i) = 1.0d0/r
    enddo
  enddo

  if (IPrint.ge.2) then
    write(IOut,1000)
    do i = 1, nch
      do j = i+1, nch
        write(IOut,1010) i,j,DvDq(i,j)
      enddo
    enddo
  endif

  Return
End Subroutine
!  !------------------------------------------------------------------------------------   
!  Subroutine MkDeDq(IOut,IPrint,nch,Chg,Rij,Rij3,DeDq)
!    Implicit Real*8 (A-H,O-Z)
!    Real*8 :: Chg(nch), Rij(4,nch,nch), Rij3(nch,nch), DeDq(3,nch)
!  
!    DeDq = 0.0d0
!    do i = 1, nch
!      do j = 1, nch
!        
!  
!  
!  
!  subroutine dqedev_esp ( x, fj, ldfj, igo, iopt, ropt )
!  USE POTENTIAL
!  USE CONSTRAINTS
!    implicit none
!    integer ( kind = 4 ) ldfj
!    real ( kind = 8 ) fj(ldfj,*)
!    integer ( kind = 4 ) igo,i,j
!    integer ( kind = 4 ) iopt(*)
!    real ( kind = 8 ) ropt(*)
!    real ( kind = 8 ) x(nch)
!    real*8 chtot
!  !
!  write(6,*) 'LDFJ = ',ldfj
!  write(6,*) 'MCON = ',mcon
!  write(6,*) 'NCH  = ',nch
!  write(6,*) 'IGO  = ',igo
!  ! Compute F
!    DO I=mcon+1,ldfj
!      fj(I,nch+1) = Vqm(I-mcon)
!      DO J=1,nch
!        fj(I,nch+1) = fj(I,nch+1) - x(J)/R(4,I-mcon,J)
!      ENDDO
!    ENDDO
!  ! Compute G
!    fj(1:mcon,nch+1) = matmul(dG,x)
!    write(6,*) 'dG'
!    write(6,1111) dG
!    write(6,*) 'x'
!    write(6,1111) x
!    write(6,*) 'fj'
!    write(6,1111) fj(1:mcon,nch+1)
!  1111 format(10(1x,f12.6))
!  !DEBUG
!  ! fj(1,nch+1) = fj(1,nch+1)+2.00
!  !ENDDEBUG
!    write(6,*) fj(1:mcon,nch+1)
!    
!    
!  ! chtot = 0.0d0
!  ! do i = 1, nch
!  !   chtot = chtot + x(i)
!  ! enddo
!  ! fj(1,nch+1) = chtot
!  
!    if (igo.ne.0) then
!  
!  ! Compute dF
!      DO I=mcon+1,ldfj
!        DO J=1,nch
!          fj(I,J) = -1.0d0/R(4,I-mcon,J)
!        ENDDO
!      ENDDO
!  
!  ! Compute dG
!      fj(1:mcon,1:nch) = dG
!  !   DO J=1,nch
!  !     fj(1,J) = 1.0d0
!  !   ENDDO
!    endif
!  
!  !write(6,*) 'ldfj =',ldfj,nch
!  !do i = 1, mcon
!  !  write(6,*) (FJ(i,j),j=1,NCh+1)
!  !enddo
!  
!  write(6,*) 'cazz'
!  do i = 1, ldfj
!    write(6,7777) i,(fj(I,J),j=1,nch+1)
!  enddo
!  7777 format(i4,(17(f10.5)))
!  RETURN
!  end subroutine
!  !
!  ! ---------------------------------------------------------------
!  !
!  subroutine dqedev_pesp ( x, fj, ldfj, igo, iopt, ropt )
!  USE POTENTIAL
!  USE MMPol
!  USE CONSTRAINTS
!  implicit none
!  integer ( kind = 4 ) ldfj
!  real ( kind = 8 ) fj(ldfj,*)
!  integer ( kind = 4 ) igo,i,j
!  integer ( kind = 4 ) iopt(*)
!  real ( kind = 8 ) ropt(*)
!  real ( kind = 8 ) x(nch)
!  real*8 chtot,rmx,rmy,rmz,d3
!  real*8  VCh(ngr),Vdip(ngr),dip(3,nch),dVch(ngr,nch),dVdip(ngr,nch)
!  real*8  ECh(3,nch)
!  
!  !
!  ! Compute potential due to charges (Vch)
!  !
!  call chpot(nch,ngr,x,R,Vch)
!  !
!  ! Compute the potential due to the dipoles (Vdip)
!  !
!  Call EChCalc(1,6,nch,ngr,IAnMMP,IMMPCn,RChg,x,ECh)
!  Call DipCalc(1,6,D,nch,ECh,dip)
!  Call dipot(nch,ngr,R,Dip,Vdip)
!  !
!  ! Compute fj
!  !
!  DO I=mcon+1,ldfj
!  !  fj(I,nch+1) = Vqm(I-mcon) - VCh(I-mcon) 
!     fj(I,nch+1) = Vqm(I-mcon) - VCh(I-mcon) - Vdip(I-mcon)
!  ENDDO
!  !
!  ! Compute partial derivatives of the potential due to charges (dVch)
!  !
!  Call dchpot(nch,ngr,R,dVch)
!  ! >>> dVch corretto
!  !
!  ! Compute partial derivatives of the potential due to induced dipoles (dVdip)
!  !
!  Call ddipot(nch,ngr,IAnMMP,IMMPCn,RChg,R,D,dVdip)
!  
!  ! Compute G
!    fj(1:mcon,nch+1) = matmul(dG,x)
!  ! chtot = 0.0d0
!  ! do i = 1, nch
!  !   chtot = chtot + x(i)
!  ! enddo
!  ! fj(1,nch+1) = chtot
!  
!    if (igo.ne.0) then
!  
!  ! Compute dF
!      DO I=mcon+1,ldfj
!        DO J=1,nch
!  !        fj(I,J) = -1.0d0/R(4,I-mcon,J)
!           fj(I,J) = -dVch(I-mcon,J) - dVdip(I-mcon,J)
!        ENDDO
!      ENDDO
!  write(6,*) 'cazz 2'
!  do i = 1+mcon, ldfj
!    write(6,7777) i,(fj(I,J),j=1,nch)
!  enddo
!  7777 format(i4,(16(f10.5)))
!  ! Compute dG
!      fj(1:mcon,1:nch) = dG
!  !   DO J=1,nch
!  !     fj(1,J) = 1.0d0
!  !   ENDDO
!    endif
!  
!  ! write(6,*) 'ldfj =',ldfj
!  ! do i = mcon, mcon+1
!  !   write(6,*) (FJ(i,j),j=1,NCh+1)
!  ! enddo
!  
!  RETURN
!  end subroutine
!  
!  !
!  ! *******************************************************************
!  ! Compute the potential due to induced dipoles
!  !
!  Subroutine dchpot(nch,ngr,R,dVch)
!  implicit real*8 (A-H,O-Z)
!  real*8 dVch(ngr,nch),R(4,ngr,nch)
!  DO I=1,ngr
!    DO J=1,nch
!      dVCh(I,J) = 1.0d0/R(4,I,J)
!    ENDDO
!  ENDDO
!  RETURN
!  End Subroutine
!  
!  !
!  ! *******************************************************************
!  ! Compute the partial derivative of the induced dipoles
!  !
!  Subroutine ddipot(nch,ngr,IAnMMP,IMMPCn,RChg,R,D,dVdip)
!  implicit real*8 (A-H,O-Z)
!  COMMON /OPTIONS/ LAnMMP
!  real*8 RChg(4,nch,nch),dE(3*nch,nch),dVdip(ngr,nch),ddip(3*nch,nch),D(3*nch,3*nch),R(4,ngr,nch)
!  INTEGER IAnMMP(LAnMMP,*)
!  LOGICAL DoInter
!  
!  !
!  ! Compute the partial derivative of electric field
!  !
!  dE = 0.0d0
!  DO I = 1,nch
!    DO J = 1,nch
!    If (I.ne.J) then
!      If (DoInter(I,J,IAnMMP,IMMPCn,nch)) then
!        dx    = RChg(1,I,J)
!        dy    = RChg(2,I,J)
!        dz    = RChg(3,I,J)
!        dist3 = RChg(4,I,J)**3
!        dE(I,J)       = dx/dist3
!        dE(I+nch,J)   = dy/dist3
!        dE(I+2*nch,J) = dz/dist3
!       Endif
!      Endif
!    EndDo
!  EndDo
!  !
!  ! Multiply D*dE
!  !
!  ddip = MatMul(D,dE)
!  Call dVdipmake(nch,ngr,ddip,R,dVdip)
!  !write(6,*) dVdip
!  
!  RETURN
!  End Subroutine
!  
!  
!  
!  Subroutine dVdipmake(nch,ngr,ddip,R,dVdip)
!  implicit real*8 (A-H,O-Z)
!  real*8 ddip(3*nch,nch),R(4,ngr,nch),dVdip(ngr,nch),R1(ngr,3*nch)
!  
!  R1 = 0.0d0
!  
!  DO I=1,ngr
!    DO J=1,nch
!      R1(I,J)       = R(1,I,J)/(R(4,I,J)**3)
!      R1(I,J+nch)   = R(2,I,J)/(R(4,I,J)**3)
!      R1(I,J+nch*2) = R(3,I,J)/(R(4,I,J)**3)
!    ENDDO
!  ENDDO
!  
!  !DO I=1,ngr
!  !  DO J=1,3*nch
!  !    WRITE(6,*) I,J,R1(I,J)
!  !  ENDDO
!  !ENDDO
!  
!  dVdip = MatMul(R1,ddip)
!  
!  RETURN
!  End Subroutine
!  
!  
!  
!  !
!  ! *******************************************************************
!  ! Compute the potential due to induced dipoles
!  !
!  Subroutine dipot(nch,ngr,R,dip,Vdip)
!  implicit real*8 (A-H,O-Z)
!  dimension Vdip(ngr),R(4,ngr,nch),dip(3,nch)
!  Do i=1,ngr
!    Vdip(i) = 0
!    Do j=1,nch
!      rmx = R(1,I,J)*dip(1,j)
!      rmy = R(2,I,J)*dip(2,j)
!      rmz = R(3,I,J)*dip(3,j)
!      d3 = R(4,I,J)**3
!      Vdip(i) = Vdip(i) + (1.0d0/d3)*(rmx+rmy+rmz)
!    Enddo
!  EndDo
!  Return
!  End Subroutine
!  
!  !
!  ! *******************************************************************
!  ! Compute electric field due to charges on the charge site by conidering
!  ! the connectivity
!  !
!  Subroutine EChCalc(IPrint,IOut,NumPol,NumGrid,IAnMMP,IMMPCn,RPol,ChInit,Ei0)
!  implicit real*8 (A-H,O-Z)
!  COMMON /OPTIONS/ LAnMMP
!  logical :: DoInter
!  DIMENSION :: IAnMMP(LAnMMP,*),RPol(4,NumPol,NumPol),ChInit(NumPol),Ei0(3,*)
!  
!  Do i = 1,NumPol
!    Ex = 0.0d0
!    Ey = 0.0d0
!    Ez = 0.0d0
!    Do j = 1,NumPol
!    if (i.ne.j) then
!      if (DoInter(I,J,IAnMMP,IMMPCn,NumPol)) then
!        dx = RPol(1,I,J)
!        dy = RPol(2,I,J)
!        dz = RPol(3,I,J)
!        dist = RPol(4,I,J)
!        dist3 = dist**3
!        !WRITE(6,*) I,J,dist3
!        Ex = Ex + ChInit(j)*(dx/dist3)
!        Ey = Ey + ChInit(j)*(dy/dist3)
!        Ez = Ez + ChInit(j)*(dz/dist3)
!       endif
!      endif
!    EndDo
!    Ei0(1,i) = Ex
!    Ei0(2,i) = Ey
!    Ei0(3,i) = Ez
!  EndDo
!  RETURN
!  End Subroutine
!  !
!  !--------------------------------------------------------------------
!  ! Compute the induced dipole (Dip = D * Ei0) 
!  !
!  subroutine DipCalc (IPrint,IOut,D,NumPol,Ei0,Dip)
!  implicit real*8 (A-H,O-Z)
!  dimension :: Ei0(3,NumPol),D(3*NumPol,3*NumPol),Dip(3,NumPol)
!  Do I=1,NumPol
!    DipX = 0.d0
!    DipY = 0.d0
!    DipZ = 0.d0
!    Do J=1,NumPol
!      DipX = DipX + D(I,J)*Ei0(1,J) + D(I,J+NumPol)*Ei0(2,J) + D(I,J+NumPol*2)*Ei0(3,J)
!      DipY = DipY + D(I+NumPol,J)*Ei0(1,J) + D(I+NumPol,J+NumPol)*Ei0(2,J) + D(I+NumPol,J+NumPol*2)*Ei0(3,J)
!      DipZ = DipZ + D(I+NumPol*2,J)*Ei0(1,J) + D(I+NumPol*2,J+NumPol)*Ei0(2,J) + D(I+NumPol*2,J+NumPol*2)*Ei0(3,J)
!    EndDo
!    Dip(1,I) = DipX
!    Dip(2,I) = DipY
!    Dip(3,I) = DipZ
!  EndDo
!  RETURN
!  End Subroutine
!  
!
!--------------------------------------------------------------------
! Read the connectivity information from a .mol2 file
!
Subroutine readconn(IPrint,IOut,fileconn,nch,IAnMMP,LAnMMP)
IMPLICIT REAL*8 (A-H,O-Z)
INTEGER :: nbond,nch,natoms
INTEGER,ALLOCATABLE :: ib1(:),ib2(:)
INTEGER :: IAnMMP(nch,LAnMMP)
CHARACTER(LEN=50) :: fileconn,what
1101 Format(I5,I6)
1102 Format('ERROR: Number of atom in mol2 is ',I5,' expected number is ',I6)
1103 Format(6X,2(I5))
1104 Format(I5,'|',9(I5))
OPEN(unit=11,file=fileconn,status='unknown')
!
! Skip the first two lines and read the number of atoms and number of 
! bonds. Check if the number of atoms in the mol2 file is the same of
! the number of charge to fit.
!
READ(11,*)
READ(11,*)
READ(11,*) natoms,nbond
If (natoms.ne.nch ) Then
  WRITE(IOut,1102) natoms,nch
  STOP
EndIf
ALLOCATE(ib1(nbond),ib2(nbond))
IAnMMP = 0
!
! Skip the ATOM Section and goes to the connectivity information.
! Read the connectivity and store the information in ib1 and ib2.
!
DO I=1,nch+10
  READ(11,*) what
  If (what.eq.'@<TRIPOS>BOND') GOTO 10
ENDDO
10 continue ! WRITE(IOut,*) 'Found',what
DO I=1,nbond
  READ(11,*) IXXX, ib1(I), ib2(I)
ENDDO
CLOSE(11)
!
! Build the connectivity in the standard MMPol Format and add it
! to the IAnMMP matrix.
!
DO N=1,nch
  K=1
  DO I=1,nbond
    If (N.eq.ib1(I))  Then 
      IAnMMP(N,K) = ib2(I)
!     IAnMMP(K,N) = ib2(I)
      K=K+1
    EndIf
    If (N.eq.ib2(I))  Then
      IAnMMP(N,K) = ib1(I)
!     IAnMMP(K,N) = ib1(I)
      K=K+1
    EndIf
    If (K.gt.LAnMMP) Then
      WRITE(IOut,*) 'ERROR: Exceed in IAnMMP dimension!'
      STOP
    EndIf
  ENDDO
ENDDO
DEALLOCATE(ib1,ib2)

if (IPrint.ge.2) then
 DO I=1,nch
   WRITE(IOut,1104) I,(IAnMMP(I,J),J=1,LAnMMP)
 ENDDO
endif

RETURN
End Subroutine
!--------------------------------------------------------------------
!     Subroutine to compute the distances between the gridpoint and
!     the atoms.
!     If symmetric, also compute R^3
!
Subroutine dist(nch,ngr,LSymm,CChg,CGrd,R,R3)
  IMPLICIT REAL*8 (A-H,O-Z)
  LOGICAL :: LSymm
  DIMENSION CChg(3,nch),CGrd(3,ngr),R(4,ngr,nch),R3(nch,nch)
!
  R = 0.0d0
  if (LSymm) then 
    do i = 1, nch
      do j = i+1, nch
        dx  = CChg(1,i) - CChg(1,j)
        dy  = CChg(2,i) - CChg(2,j)
        dz  = CChg(3,i) - CChg(3,j)
        dd  = sqrt(dx*dx + dy*dy + dz*dz)
        dd3 = dd**3
        R(1,i,j) = dx
        R(1,j,i) = -dx
        R(2,i,j) = dy
        R(2,j,i) = -dy
        R(3,i,j) = dz
        R(3,j,i) = -dz
        R(4,i,j) = dd
        R(4,j,i) = dd
        R3(i,j)  = dd3
        R3(j,i)  = dd3
      enddo
    enddo
  else
    do i = 1, ngr
      do j = 1, nch
        dx = CGrd(1,i) - CChg(1,j)
        dy = CGrd(2,i) - CChg(2,j)
        dz = CGrd(3,i) - CChg(3,j)
        dd = sqrt(dx*dx + dy*dy + dz*dz)
        R(1,i,j) = dx
        R(2,i,j) = dy
        R(3,i,j) = dz
        R(4,i,j) = dd
      enddo
    enddo
  endif
        
  RETURN
End Subroutine
!
! -------------------------------------------------------------------
!     Subroutine to compute the potential generated by a set of
!     charges over the gridpoints.
!
Subroutine chpot(nch,ngr,Chg,R,Vch)
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION Chg(nch),R(4,ngr,nch),Vch(ngr)
DO I=1,ngr
  Vch(I) = 0.0d0
  DO J=1,nch
    Vch(I) = Vch(I) + Chg(J)/R(4,I,J)
  ENDDO
ENDDO
RETURN
End Subroutine
!
! -------------------------------------------------------------------
!     Subroutine to compute the potential generated by a set of
!     charges and induced dipoles over the gridpoints.
!
Subroutine chpotpol(nch,ngr,Chg,DInv,Scr,RChGr,RChCh,R3,Vtot,DipInd)
  Implicit Real*8 (A-H,O-Z)
  DIMENSION Chg(nch),RChGr(4,ngr,nch),RChCh(4,nch,nch),Vtot(ngr),DInv(3*nch,3*nch)
  Dimension Scr(nch,nch), R3(nch,nch), DipInd(3)
  Dimension Ex(3*nch,nch),E(3*nch),Dip(3*nch),BxGrd(ngr,3*nch),Vdip(ngr),Vch(ngr)
  Parameter (Zero = 0.0d0, Small = 1.0d-10)
! Sanitise Rij3
  do i = 1, nch
    do j = 1, nch
      if (R3(i,j).lt.Small) R3(i,j)=Small
    enddo
  enddo
! Compute electric field due to charges
  Ex(1:nch,:)         = RChCh(1,:,:)*(Scr/(R3))
  Ex(nch+1:2*nch,:)   = RChCh(2,:,:)*(Scr/(R3))
  Ex(2*nch+1:3*nch,:) = RChCh(3,:,:)*(Scr/(R3))
  E = Matmul(Ex,Chg)
! Compute dipoles generated by the electric field
  Dip = Matmul(DInv,E)
! Compute potential on the grid generated by the dipoles
  BxGrd(:,1:nch)         = RChGr(1,:,:)/(RChGr(4,:,:)**3)
  BxGrd(:,nch+1:2*nch)   = RChGr(2,:,:)/(RChGr(4,:,:)**3)
  BxGrd(:,2*nch+1:3*nch) = RChGr(3,:,:)/(RChGr(4,:,:)**3)
  Vdip = Matmul(BxGrd,Dip)
! Compute potential on the grid generated by the charges
  Call chpot(nch,ngr,Chg,RChGr,Vch)
  Vtot = Vch + Vdip
! Compute total induced dipole
  do i = 1, 3
    DipInd(i) = sum(Dip((i-1)*nch+1:i*nch))
  enddo
!
  Return
End Subroutine
!
! -------------------------------------------------------------------
!     Compute the RMS between two quantities 
!
REAL*8 Function rms(n,QRef,Qvar)
implicit real*8 (A-H,O-Z)
dimension Qref(n),Qvar(n)
rms = 0.0d0
Do I=1,n
  rms = rms + (Qref(I)-Qvar(I))**2
Enddo
return 
End Function
!
! -------------------------------------------------------------------
! Compute the fitting error
!
Function FitErr(IOut,IPrint,LPol,nch,ngr,Q,RChGr,RChCh,R3,DInv,Scr,Vqm,DipInd)
  Implicit Real*8(A-H,O-Z)
  Logical LPol
  Real*8 FitErr,Q(nch),RChGr(4,ngr,nch),RChCh(4,nch,nch),R3(nch,nch),Vch(ngr)
  Real*8 DInv(3*nch,3*nch),Scr(nch,nch), Vqm(ngr), DipInd(3)
  if (LPol) then
    call chpotpol(nch,ngr,Q,DInv,Scr,RChGr,RChCh,R3,Vch,DipInd)
  else
    call chpot(nch,ngr,Q,RChGr,Vch)
  endif
  FitErr = rms(nch,Vqm,Vch)
  return
End Function
!
! -------------------------------------------------------------------
!
Subroutine MkScr(IOut,IPrint,NPol,CPol,Pol,ScrChCh,ScrChPl,IAnMMP,IMMPCn,LAnMMP, &
  LScrChPl)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION :: ScrChCh(npol,npol), ScrChPl(npol,npol), IAnMMP(npol,LAnMMP)
  DIMENSION :: CPol(3,*), R(3), Pol(NPol)
  INTEGER, ALLOCATABLE :: Neigh(:,:)
  LOGICAL :: DoThole, LScrChPl
  PARAMETER(a=1.7278d0,b=2.5874d0,c=2.0580d0)
!
! > DoThole .... True if you want to use Thole
! > IScreen .... 1 = use Amber Screening Factor (a)
!           .... 2 = use AL (Wang) Screening Factor (b)
!           .... 3 = use DL (Wang) Screening Factor (c) [NYI]
! Detect intramolecular polarization treatment
    IF (LScrChPl) THEN
      IF (IMMPCn.EQ.3) THEN
        DoThole = .TRUE.
        IScreen = 1
      ELSEIF (IMMPCn.EQ.1) THEN
        DoThole = .FALSE.
        IScreen = 1
      ELSEIF (IMMPCn.EQ.4) THEN
        DoThole = .TRUE.
        IScreen = 2
      ELSEIF (IMMPCn.EQ.5) THEN
        DoThole = .TRUE.
        IScreen = 3
      ELSE
        WRITE(6,*) 'Confused in polarization treatment.'
        STOP
      ENDIF
    ENDIF

  Allocate (Neigh(npol,npol))
  Neigh = 0

! Neigh integer matrix:
!  N(i,i) = 1
!  N(i,j) = 2 if i and j are 1-2 neighbours
!  N(i,j) = 3 if i and j are 1-3 neighbours
!  N(i,j) = 4 if i and j are 1-4 neighbours
!  N(i,j) = 0 if i and j are not within 4 bonds

  do i = 1, npol
    Neigh(i,i) = 1
    do j = 1, LAnMMP
      if (IAnMMP(i,j).eq.0) then
        goto 700
      else
        Neigh(IAnMMP(i,j),i) = 2
        Neigh(i,IAnMMP(i,j)) = 2
      endif
    enddo
 700 continue
  enddo

  do i = 1, npol
    do j = 1, npol
      if (Neigh(i,j).eq.2) then
        do k = 1, npol
          if (Neigh(j,k).eq.2 .and. Neigh(i,k).ne.1 .and. Neigh(i,k).ne.2) then
            Neigh(i,k) = 3
            Neigh(k,i) = 3
          endif
        enddo
      endif
    enddo
  enddo

  do i = 1, npol
    do j = 1, npol
      if (Neigh(i,j).eq.3) then
        do k = 1, npol
          if (Neigh(j,k).eq.2 .and. Neigh(i,k).ne.1 .and. Neigh(i,k).ne.2 .and. Neigh(i,k).ne.3) then
            Neigh(i,k) = 4
            Neigh(k,i) = 4
          endif
        enddo
      endif
    enddo
  enddo

  if (IPrint.ge.2) then
    write(IOut,1000)
    do i = 1, npol
      write(IOut,2000) (Neigh(i,j),j=1,npol)
    enddo
  endif
 1000 format(' Neighbours matrix:')
 2000 format(20(1x,i4))

! Matrix Neigh is complete: make Screening matrices 
!   This is only implemented for Amber12 for the time being
 ScrChCh = 1.0d0
 ScrChPl = 1.0d0

  do i = 1, npol
    ScrChCh(i,i) = 0.0d0
    ScrChPl(i,i) = 0.0d0
    do j = i+1, npol
      if (IMMPCn.eq.1.or.IMMPCn.eq.4) then
        if (Neigh(i,j) .eq. 2 .or. Neigh(i,j) .eq. 3) then
          ScrChCh(i,j) = 0.0d0
          ScrChCh(j,i) = 0.0d0
          ScrChPl(i,j) = 0.0d0
          ScrChPl(j,i) = 0.0d0
        endif
      endif
      if (LScrChPl.AND.(IMMPCn.EQ.3.OR.IMMPCn.EQ.4.OR.IMMPCn.EQ.5)) THEN
!
! Compute the distance between two polarizable sites
!
        R(1) = CPol(1,J)-CPol(1,I)
        R(2) = CPol(2,J)-CPol(2,I)
        R(3) = CPol(3,J)-CPol(3,I)
        Rij = Sqrt(R(1)*R(1) + R(2)*R(2) + R(3)*R(3))
        If (IScreen.eq.1) s = a*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
        If (IScreen.eq.2) s = b*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
        If (IScreen.eq.3) s = c*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
!
! Thole linear screening
!
        If (Rij.le.s) Then
          v = Rij/s 
          Scale3 = 4.0d0*(v**3)-3.0d0*(v**4)
          ScrChPl(i,j) = Scale3
          ScrChPl(j,i) = Scale3
        endif
      elseif (IMMPCn.eq.2.or.IMMPCn.eq.0) then
        Write(*,*) 'Code not ready to deal with groups option'
        STOP
      endif
    enddo
  enddo
  
  Deallocate (Neigh)
  Return
End Subroutine
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
!
! -------------------------------------------------------------------
!
Subroutine PESP(IOut,IPrint,npol,ngr,mcon,RChGr,Rij,Rij3,Vqm,X,B,ScrChPl,ScrChCh,DInv,Qpesp,QSum)
  Implicit Real*8 (A-H,O-Z)
  Dimension RChGr(4,ngr,npol), Vqm(ngr), X(npol+mcon,npol+mcon), B(npol+mcon)
  Dimension Rij(4,npol,npol), Rij3(npol,npol)
  Dimension Qpesp(npol+mcon), ScrChPl(npol,npol), ScrChCh(npol,npol)
  Dimension dVdQ(ngr,npol),dEdQ(3*npol,npol),dMudQ(3*npol,npol),DInv(3*npol,3*npol)
  Dimension RTemp(ngr,3*npol),dVdQChg(ngr,npol),dVdQPol(ngr,npol), E(npol,npol)
  Integer, Allocatable :: IPIV(:)
  Real*8, Allocatable  :: WORK(:)
  Real*8, Allocatable  :: Ax(:,:),BxPol(:,:),BxGrd(:,:),Cx(:,:),Ex(:,:),Fx(:,:),Gx(:,:),Hx(:)
  Parameter (One = 1.0d0, Small = 1.0d-10)
 1000 Format(' polESP matrix: Element ',i6,' has an illegal value',/,' Halting.')
 1010 Format(' polESP matrix: Diagonal element ',i6,' is zero.',/, &
             ' The matrix is singular and cannot be inverted',/,   &
             ' If you have a symmetric molecule and you are imposing a dipole',/,&
             '   constraint,  then  you  may be asking for linearly dependent',/,&
             '   conditions.    Please remove constraint on dipole and see if',/,&
             '   it works. Otherwise contact the author.---------------------',/, &
             ' Halting now.')
 1020 Format(' polESP matrix: Inversion successful.')
 1100 Format(' polESP charges:')
 1110 Format(1x,i6,1x,f12.6)
 1120 Format(' --------------------',/,' Total  ',f12.6)
!
! For definitions and notation, please refer to document by S.Caprasecca
!
  Allocate(Ax(ngr,npol), BxPol(3*npol,npol), BxGrd(ngr,3*npol), Cx(3*npol,npol))
  Allocate(Ex(ngr,npol), Fx(ngr,npol), Gx(npol,npol), Hx(npol))
! --- Matrix Ax
  Ax = One/RChGr(4,:,:)
! --- Sanitise Rij3
  do i = 1, npol
    do j = 1, npol
      if (Rij3(i,j).lt.Small) Rij3(i,j)=Small
    enddo
  enddo
! --- Matrix BxPol
  BxPol(1:npol,:)          = Rij(1,:,:)*(ScrChPl/(Rij3))
  BxPol(npol+1:2*npol,:)   = Rij(2,:,:)*(ScrChPl/(Rij3))
  BxPol(2*npol+1:3*npol,:) = Rij(3,:,:)*(ScrChPl/(Rij3))
! --- Matrix BxGrd
  BxGrd(:,1:npol)          = RChGr(1,:,:)/(RChGr(4,:,:)**3)
  BxGrd(:,npol+1:2*npol)   = RChGr(2,:,:)/(RChGr(4,:,:)**3)
  BxGrd(:,2*npol+1:3*npol) = RChGr(3,:,:)/(RChGr(4,:,:)**3)
! --- Matrix Cx
  Cx = Matmul(DInv,BxPol)
! --- Matrix Ex
  Ex = Matmul(BxGrd,Cx)
! --- Matrix Fx
  Fx = Ax + Ex
! --- Matrix Gx
  Gx = Matmul(Transpose(Fx),Fx)
! --- Vector Hx
  Hx = Matmul(Transpose(Fx),Vqm)
! --- Matrix X
  X(1:npol,1:npol) = Gx
! --- Vector B
  B(1:npol) = Hx
!---   
!---   
!---     Rij3 = Rij3 + 1.0d-10
!---   ! Form matrix dVdQ (Chg)
!---     dVdQChg = 1.0d0/RChGr(4,:,:)
!---   ! Form matrix dEdQ
!---     do i = 1, npol
!---       do j = 1, npol
!---         dEdQ(i,j) = ScrChPl(i,j)*Rij(1,i,j)/(Rij3(i,j))
!---         dEdQ(npol+i,j) = ScrChPl(i,j)*Rij(2,i,j)/(Rij3(i,j))
!---         dEdQ(2*npol+i,j) = ScrChPl(i,j)*Rij(3,i,j)/(Rij3(i,j))
!---       enddo
!---     enddo
!---   ! Form matrix dMudQ
!---     dMudQ = matmul(DInv,dEdQ)
!---   ! Form matrix dVdQ (Pol)
!---     RTemp(:,1:npol) = RChGr(1,:,:)
!---     RTemp(:,npol+1:2*npol) = RChGr(2,:,:)
!---     RTemp(:,2*npol+1:3*npol) = RChGr(3,:,:)
!---     dVdQPol = matmul(RTemp,dMudQ)
!---   ! Form matrix dVdQ (Tot)
!---     dVdQ = dVdQChg-dVdQPol
!---   ! Form matrix X
!---     X(1:npol,1:npol) = matmul(transpose(dVdQ),dVdQ)
!---   ! Form vector B
!---     B(1:npol) = matmul(transpose(dVdQ),Vqm)
! Invert matrix X
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
  elseif (IPrint.ge.2) then
    write(IOut,1020)
  endif

!   (-) Compute charges

  qpESP = matmul(X,B)

  if (IPrint.ge.2) write(IOut,1100)
  qSum = 0.d0
  do i = 1, npol
    qSum = qSum + qpESP(i)
    if (IPrint.ge.2) write(IOut,1110) i,qpESP(i)
  enddo
  if (IPrint.ge.2) write(IOut,1120) qSum

 Deallocate(IPIV,WORK)

  Return
End Subroutine
!
! -------------------------------------------------------------------
!
Subroutine FixedDip(IOut,IPrint,nch,CChg,q,Dipole)
  Implicit Real*8 (A-H,O-Z)
  Dimension CChg(3,nch),q(nch),Dipole(3)
  
  Dipole = 0.0d0
  do i = 1, nch
    Dipole(1) = Dipole(1) + q(i)*CChg(1,i)
    Dipole(2) = Dipole(2) + q(i)*CChg(2,i)
    Dipole(3) = Dipole(3) + q(i)*CChg(3,i)
  enddo
  Return
End Subroutine

!
! This subroutine compute, invert and store the MMPol matrix (A^-1) 
!
Subroutine MMPMTRX(IPrint,IOut,npol,pol,CPol,D,IAnMMP,IMMPCn,LAnMMP)
!
! IMMPCn ..... Treatment of the polarization
! IAnMMP ..... Matrix containing the connectivity information
! LAnMMP ..... Dimension of IAnMMP
! pol ........ isotropic polarizabilities (a.u.)
! D .......... MMPol Matrix (a.u.)
!

IMPLICIT REAL*8(A-H,O-Z)
PARAMETER(LWORK=30000)
LOGICAL   :: DoInter,DoThole,DoCalc
DIMENSION :: pol(npol),CPol(3,npol),WORK(LWORK),IPIV(npol*3)
DIMENSION :: D(3*npol,3*npol),IAnMMP(npol,LAnMMP)
LDA=npol*3

!
! Clear the whole MMPol matrix
!
D = 0.0d0

!
! Fill in diagonal elements of MMPMTRX
!
DO I=1,npol
  D(I,I)               = 1.0d0/pol(I)     ! x-block
  D(I+npol,I+npol)     = 1.0d0/pol(I)     ! y-block
  D(I+2*npol,I+2*npol) = 1.0d0/pol(I)     ! z-block
ENDDO

!
! Fill in non-diagonal elements of MMPTRX
!
DO I=1,npol
  DO J=I+1,npol
    DoCalc = .False.

    ! Thole smeared dipole interaction tensor between all dipoles (Thole option)
    If (IMMPCn.eq.3) Then
      DoThole = .True.
      DoCalc  = .True.
      IScreen = 1

    ! Exclude Amber 1-2 and 1-3 interactions (Amber option) or interactions
    ! inside sites from same group (Groups option).
    ElseIf (IMMPCn.eq.0 .or. IMMPCn.eq.1 .or. IMMPCn.eq.2) then
      DoThole = .False.
      IScreen = 1
      if (DoInter(I,J,IAnMMP,IMMPCn,npol)) DoCalc=.True.

    ! Exclude Amber 1-2 and 1-3 interactions (Amber option) and use Thole
    ! smeared dipole interaction tensor for the others with Wang screening
    ! parameter
    ElseIf (IMMPCn.eq.4) Then
      DoThole = .True.
      IScreen = 2
      If (DoInter(I,J,IAnMMP,IMMPCn,npol)) DoCalc=.True.

    ! DL Amber model (compute all interactions with Thole linear screening)
    ElseIf (IMMPCn.eq.5) Then
      DoThole = .True.
      DoCalc=.True.
      IScreen = 3

    Else
      Write(6,*) 'Confused in polarization treatment.'
      Stop
    EndIf

    If (DoCalc) then
      D(I,J) = Tensormm(1,1,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J,I) = D(I,J)
      D(I,J+npol) = Tensormm(1,2,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J,I+npol) = D(I,J+npol)
      D(I+npol,J) = D(I,J+npol)
      D(J+npol,I) = D(I,J+npol)
      D(I,J+2*npol) = Tensormm(1,3,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J,I+2*npol) = D(I,J+2*npol)
      D(I+2*npol,J) = D(I,J+2*npol)
      D(J+2*npol,I) = D(I,J+2*npol)
      D(I+npol,J+npol) = Tensormm(2,2,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J+npol,I+npol) = D(I+npol,J+npol)
      D(I+npol,J+npol*2) = Tensormm(2,3,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J+npol,I+npol*2) = D(I+npol,J+npol*2)
      D(I+npol*2,J+npol) = D(I+npol,J+npol*2)
      D(J+npol*2,I+npol) = D(I+npol,J+npol*2)
      D(I+npol*2,J+npol*2) = Tensormm(3,3,I,J,CPol,Pol,npol,DoThole,IScreen)
      D(J+npol*2,I+npol*2) = D(I+npol*2,J+npol*2)
    EndIf

  ENDDO
ENDDO
 
!
!     Invert the MMPol matrix using lapac library
!
999   Format(6(1x,D15.5))      
If (IPrint.ge.2) then
  write(6,*) 'Matrice non invertita'
  Do i=1,npol*3
    Write(IOut,999) (D(i,j),j=1,npol*3)
  Enddo
EndIf

M = npol*3
N = npol*3

call DGETRF (M,N,D,LDA,IPIV,INFO)
call DGETRI (N,D,LDA,IPIV,WORK,LWORK,INFO)

!
!     Save the inverted MMPol Matrix
!
If (IPrint.ge.2) then
  write(6,*) 'Matrice invertita'
  Do i=1,npol*3
    Write(IOut,999) (D(i,j),j=1,npol*3)
  Enddo
EndIf

End Subroutine

!--------------------------------------------------------------------
! Function to decide if two dipole should interact or not
!
Logical Function DoInter(I,J,IAnMMP,IMMPCn,npol)
Implicit Real*8(A-H,O-Z)
COMMON /OPTIONS/ LAnMMP
DIMENSION :: IAnMMP(LAnMMP,*)
!
! QM/MMPol connectivity check:
! Return whether an interaction between two sites has to be computed.
!
DoInter = .True.
!   
! Check if sites belong to same group or to a corresponding excluded group: 
! Groups option.
! 
If(IMMPCn.eq.2.or.IMMPCn.eq.3) Then
  IinJ = 0
  JinI = 0
  Do 50 K=1,LAnMMP
  ! Check if I-group is in J list
    If (IAnMMP(1,I).eq.IAnMMP(K,J)) IinJ = 1
  ! Check if J-group is in I list
    If (IAnMMP(K,I).eq.IAnMMP(1,J)) JinI = 1
50 Continue       
If (IinJ.ne.JinI) Write(6,*) 'Inconsistency error in MultGroups information.'
If (IinJ.eq.1) DoInter = .False.
!
!     Check if sites are 1-2 or 1-3 neighbours: Amber and Wang option
!
Else if(IMMPCn.eq.1.or.IMMPCn.eq.4) Then
  Do 20 K=1,LAnMMP
  Neigh12 = IAnMMP(K,I)
  If (Neigh12.eq.0) Goto 10
  If (Neigh12.eq.J) Then
    DoInter = .False.
    Goto 10
  Endif
  Do 30 L=1,LAnMMP
    Neigh13 = IAnMMP(L,Neigh12)
    If (Neigh13.eq.0) Goto 20
    If (Neigh13.eq.J) Then
      DoInter = .False.
      Goto 10
    Endif

30  Continue
20  Continue
10  Endif

Return
End Function

!--------------------------------------------------------------------
! Function to compute the tensor
!
REAL*8 Function TensorMM(M,N,I,J,CPol,Pol,npol,DoThole,IScreen)
IMPLICIT REAL*8 (A-H,O-Z)
LOGICAL :: DoThole
DIMENSION  :: CPol(3,*),Pol(*),R(3)
PARAMETER(a=1.7278d0,b=2.5874d0,c=2.0580d0)
!
! > DoThole .... True if you want to use Thole  
! > IScreen .... 1 = use Amber Screening Factor (a)
!           .... 2 = use AL (Wang) Screening Factor (b)
!           .... 3 = use DL (Wang) Screening Factor (c)
! > M .......... x,y,z component of site I
! > N .......... x,y,z component of site J
! > I .......... index of I polarizable site
! > J .......... index of J polarizable site
!
!
! Compute the distance between two polarizable sites.
!
R(1) = CPol(1,J)-CPol(1,I)
R(2) = CPol(2,J)-CPol(2,I)
R(3) = CPol(3,J)-CPol(3,I)
Rij = Sqrt(R(1)*R(1) + R(2)*R(2) + R(3)*R(3))
!
If (DoThole) then
  If (IScreen.eq.1) s = a*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
  If (IScreen.eq.2) s = b*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
  If (IScreen.eq.3) s = c*((Pol(I)*Pol(J))**(1.0d0/6.0d0))
EndIf
! 
! Thole smeared dipole interaction tensor 
! (Van Duijnen parameters from JPCA, 102, 2399, 1998)
!
If (DoThole.and.Rij.le.s) Then
  v = Rij/s
  Scale3 = 4.0d0*(v**3)-3.0d0*(v**4)
  Scale5 = v**4
  If (M.eq.N) Then
    TensorMM = (1/(Rij**3))*(Scale3-(Scale5*(3/(Rij**2))*R(M)*R(N)))
  Else
    TensorMM = - Scale5*(3/(Rij**5))*R(M)*R(N)
  Endif
!
! Standard dipole interaction tensor
!
Else
  If (M.eq.N) Then
    TensorMM = (1/(Rij**3))*(1 - ((3/(Rij**2))*R(M)*R(N)))
  Else
    TensorMM = - (3/(Rij**5))*R(M)*R(N)
  Endif
Endif
Return
End Function




Subroutine RdOpts(IOut,IPrint,filename,fileconn,filepol,filecnst,LScrChPl)
!
! This subroutine reads input filenames and printout options from line
!
! Use the following options:
!  -h ... Get help message
!  -d ... Debug mode (extra printout)
!  -s ... Silent mode (minimum printout)
!  -g ... before gesp file name
!  -m ... before mol2 file name
!  -p ... before specifying polarizability file name
!  -c ... before specifying constraint file name
!  -x ... include chg-pol screening as in Wang (experts only)
!
! Example: chpol d g xxx.gesp m xxx.mol2 p pol.in c cnstr.in
!
  integer IOut,IPrint,IArg,num_args
  logical LScrChPl
  CHARACTER(LEN=50) filename,filepol,fileconn,filecnst
  character(len=50), dimension(:), allocatable :: args

1000 format(' Error in input stream: nothing followed option ',(A))
1001 format(' Error in input stream: option ',(A),' unknown')
1002 format(' Error in input stream: not all files have been defined.',/, &
            '   Use g m p c to define all files.')
1200 format(' Use the following options to run the program:',/,            &
            '   -g (required) ... Followed by gesp file name',/,           &
            '   -m (required) ... Followed by mol2 file name',/,           &
            '   -p (required) ... Followed by polarisability file name',/, &
            '   -c (required) ... Followed by constraints file name',/,    &
            '   -x (optional) ... Also include Wang Chg-Pol screening -- WARNING: Expert use only.',/,             &
            '                     Only use if the MMPol code includes such screening. If in doubt, do not use.',/, &
            '   -h (optional) ... Get this help message',/,                &
            '   -d (optional) ... Run in debug mode (extra printout)',/,   &
            '   -s (optional) ... Run in silent mode (minimum printout)')

  num_args = command_argument_count()
  allocate(args(num_args))  

  IPrint   = 1
  IGet     = 0
  LScrChPl = .FALSE.
  filename = ''
  fileconn = ''
  filepol  = ''
  filecnst = ''

  if (num_args.eq.0) then
    write(IOut,1200)
    stop
  endif
  do IArg = 1, num_args
    call get_command_argument(IArg,args(IArg))
  end do
  do IArg = 1, num_args
    if (IGet.ne.0) then
      if (IGet.eq.1) then
        filename = args(IArg)
      elseif(IGet.eq.2) then
        fileconn = args(IArg)
      elseif(IGet.eq.3) then
        filepol = args(IArg)
      elseif(IGet.eq.4) then
        filecnst = args(IArg)
      endif
      IGet = 0
    elseif (args(IArg).eq.'-h') then
      write(IOut,1200)
      stop
    elseif (args(IArg).eq.'-d') then
      IPrint = 2
    elseif (args(IArg).eq.'-s') then
      IPrint = 0
    elseif (args(IArg).eq.'-x') then
      LScrChPl = .TRUE.
    elseif (args(IArg).eq.'-g') then
      IGet = 1
      if (IArg.eq.num_args) goto 800
    elseif (args(IArg).eq.'-m') then
      IGet = 2
      if (IArg.eq.num_args) goto 800
    elseif (args(IArg).eq.'-p') then
      IGet = 3
      if (IArg.eq.num_args) goto 800
    elseif (args(IArg).eq.'-c') then
      IGet = 4
      if (IArg.eq.num_args) goto 800
    else
      goto 801
    endif
  enddo

!  if (filename.eq.''.or.fileconn.eq.''.or.filepol.eq.''.or.filecnst.eq.'') goto 802 
  goto 900

800 write(IOut,1000) args(IArg-1)
  stop

801 write(IOut,1001) args(IArg)
  write(*,*) IArg
  stop
 
802 write(IOut,1002)
  stop
 
900 return

end subroutine
  
subroutine PrtHdr(IOut)
  integer IOut
!
! Print header
!
1000 FORMAT(/,                                                                  &
     '                                                                      ',/,&
     '                                       mm                             ',/,&
     '                                    mMMm                              ',/,&
     '                                  mMMMMm         m                    ',/,&
     '                                 mMMMMm          mMm                  ',/,&
     '                                 mMMMMm          mMm                  ',/,&
     '                                 mMMMMMm        mMMm                  ',/,&
     '                                 MMMMMMMMMMMMMMMMMMm                  ',/,&
     '                                mMMMMMMMMMMMMMMMMMm                   ',/,&
     '       __  ___      __    ____________      __MMMm     __             ',/,&
     '      /  |/  /___  / /   / ____/ ____/___  / /  ____ _/ /_            ',/,&
     '     / /|_/ / __ \/ /   / __/ / /   / __ \/ /  / __ `/ __ \           ',/,&
     '    / /  / / /_/ / /___/ /___/ /___/ /_/ / /__/ /_/ / /_/ /           ',/,&
     '   /_/  /_/\__________/_____/\____/_____/_____|__,_/_.___/            ',/,&
     '           /_  __/ __ \/ __ \/ /  / ___/                              ',/,&
     '            / / / / / / / / / /   \__ \                               ',/,&
     '           / / / /_/ / /_/ / /___ __/ /                               ',/,&
     '          /_/  \____/\____/_____/____/                                ',/,&
     '            mMMMMMMMMMMMMMMMm                                         ',/,&
     '          mMMMMMMMMMMMMMMMm                                           ',/,&
     '        mMMMMMMMMMMMMMMMMM   + ------------------------------------ + ',/,&
     '       mMMMMMMMMMMMMMMMMm    |    P O L C H A T                     | ',/,&
     '      mMMMMMMMMMMMMMMMMMm    + ------------------------------------ + ',/,&
     '      mMMMMm       mMMMMMm   | Stefano Caprasecca                   | ',/,&
     '      mMMMm       mMMMMMMm   | Sandro Jurinovich                    | ',/,&
     '       mMm       mMMMMMMm    | Carles Curutchet           ver 3.1.0 | ',/,&
     '        m       mMMMMMMm     |          www.dcci.unipi.it/molecolab | ',/,&
     '               mMMMMMm       + ------------------------------------ + ',/)
  write(IOut,1000)
  return
end subroutine
