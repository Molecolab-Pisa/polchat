module constants

  implicit none

  character(50)     :: version = "4.0.0"
  integer           :: IMMPCn = 4
  integer,parameter :: LAnMMP = 9
  integer           :: iout = 6
  integer           :: iprt
  integer,parameter :: outmsgmax=40, strmax=500, inpargmax=100, typmax=6
  integer,parameter :: timemax=10
  real*8, parameter :: zero=0.0d0, one=1.0d0, sixth=1.0d0/6.0d0, small=1.0d-8
  real*8, parameter :: three=3.0d0, four=4.0d0
  real*8, parameter :: ang2au=1.889725989d0, deb2au=0.393430201407683d0
  real*8, parameter :: scra=1.7278d0, scrb=2.5874d0, scrc=2.0580d0

end module
