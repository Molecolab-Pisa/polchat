
# Select fortran compiler here, or
# make FC an env variable
# here we select standard compiler
compiler = ifort

##################################################

ROOTDIR  = $(shell pwd)
OPSYS    = $(shell uname -s )

# override Make's defaults that we don't need here
ifeq ($(FC),f77)
 FC = $(compiler)
endif

# Select compiler-dependent flags
# ... Intel
ifeq ($(FC),ifort)
 # Use sequential Math Kernel Library
 # (should work w. all intel compilers) 
 FFLAGS = -i8  -mkl=sequential -O3 -xHost
 LFLAGS = -lpthread -lm
endif
# ... GNU
ifeq ($(FC),gfortran)
  FFLAGS =  -ffree-line-length-none -O3 
  LFLAGS = -lpthread -lm -llapack -lblas
endif
ifeq ($(FC),pgf90)
  FFLAGS = -O3
  LFLAGS = -lpthread -lm -llapack -lblas
endif
# TODO: PGI

ifndef LFLAGS
all:
	@echo Compiler $(FC) not recognized!
	@echo 'Use make FC=(ifort|gfortran|pgf90)'
endif


OBJECTS = constants.o mmpoldata.o constraints.o precmod.o stringmod.o gespinfo.o espinfo.o operative.o time.o dist.o prthdr.o prtmat.o prtimat.o rdopts.o rdgesp.o rdconn.o rdpol.o tensor.o mkscr.o mmpmat.o dointer.o rdcns.o docns.o rdrange.o esp.o pesp.o error.o potchg.o potchgdip.o rms.o dipole.o prtdat.o printout.o atnum.o polchat.o 
.PHONY: clean clobber

polchat.exe: $(OBJECTS)
	$(FC)  $(FFLAGS) $(OBJECTS) -o polchat.exe $(LFLAGS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean: 
	rm -f $(OBJECTS) *.mod 

clobber:
	rm -f $(OBJECTS) *.mod polchat.exe

