#pgf77 dpmpar.f  enorm.f  lmder1.f  lmder.f  lmpar.f  qrfac.f  qrsolv.f chpolfit.f -o chpolfit.exe -llapack -lblas
rm *.o *.exe
pgf90 chpol.f90 -o chpol.exe -llapack -lblas
#ifort chpol2.f90 MMPMTRX.f90 -o chpolfit.exe -llapack -lblas
#ifort  dqed.f90 -o chpolfit.exe -llapack -lblas
