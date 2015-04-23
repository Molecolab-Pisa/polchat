FC = ifort
FFLAGS = -i8 -llapack -lblas -O2 -xHost
LFLAGS = -llapack -lblas
OBJECTS = chpot.o chpotpol.o dist.o dointer.o esp.o fiterr.o fixeddip.o mkdvdq.o mkrij.o mkscr.o mmpmtrx.o pesp.o rdopts.o readconn.o rms.o tensormm.o polchat.o 
.PHONY: clean clobber

polchat.exe: $(OBJECTS)
	$(FC)  $(FFLAGS) $(OBJECTS) -o polchat.exe $(LFLAGS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean: 
	rm -f $(OBJECTS)

clobber:
	rm -f $(OBJECTS) polchat.exe

