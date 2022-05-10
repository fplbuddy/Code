#####################################################
#   Makefile
#
#   Author:: Vassilios Dallas
#   			 Laboratory of Statistical Physics
#   			 Department of Physics
#   			 Ecole Normale Superieure
#   			 e-mail: vdallas@lps.ens.fr
#
#   Date:: 6 June 2012
#####################################################

### Compiler
## For gauss 
#IFORT = mpif90 -f90=ifort
#IFORT = mpif90
IFORT = mpiifort

## For jxb
#IFORT = ifort 

### Math Libraries
## For gauss 
# PRECISION: SINGLE or DOUBLE
PRECISION = DOUBLE
# PRECISION: yes or no
HYBRID = no
FFTW3SINGLE = f
FFTW3DOUBLE =
FFTW3HYByes = -lfftw3$(FFTW3$(PRECISION))_threads
FFTW3HYBno  =
FFTW3FLAGS  = $(FFTW3HYB$(HYBRID)) -lfftw3$(FFTW3$(PRECISION)) -lm
FFTWLDIR = /mi/share/scratch/dallas/Coding/FFTW3lib/lib
LIBS = -L $(FFTWLDIR) $(FFTW3FLAGS) 

## For jxb
#FFTWLDIR = /home/vdallas/FFTWlibv2/lib
#LIBS = -limf -lmpi -L $(FFTWLDIR) -lsrfftw -lsfftw

### Optimisation compiler options
#ifdef OPT

## IA-32 & Intel 64
#OPTIFORT = -traceback -xP -O3 -Zp16 -align 
#OPTIFORT = -traceback -ftz
#OPTIFORT = -traceback -fpe:0 -check all -fp-stack-check -openmp
#OPTIFORT = -traceback -openmp
#OPTIFORT = -traceback -fpe:0 -check all -fp-stack-check
OPTIFORT = -traceback

## IA-64 (Itanium)
#OPTIFORT = -O3 -ip -pad -align -i-static

## Opts used in ghost
#OPTIFORT = -O3 -ip -pad -unroll -align -i-static

#endif

### Flug enabling most warning messages
ifdef MORE
  IFORTFLUGS = -w
endif

### Debugger
ifdef DBG
  #DBGFLUGS = -g
  #DBGFLUGS = -O0 -g -traceback -fpe:0 -check all -fp-stack-check
  DBGFLUGS = -traceback -fpe:0 -check all -fp-stack-check
endif

### Filename modification
PERL     = /usr/bin/perl
UNDER    = 0
EXTEN    = 1
APPEND   = _
F90EXT   = f90
PREC     = 1

### Object files
OBJS = fftp.o fprecision_mod.o fftp_mod.o fftp2D.o \
  	    pseudospec2D_mod.o pseudospec2D_inc.o \
		 mainHD2D.o

### Exe file
EXE = run.exe

### Forced remove
RM = rm -f

### Compiles multiple source files and links them together to an executable file named 'run'.
all : $(OBJS)
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -o $(EXE) $(IFORTFLUGS) $(OBJS) $(LIBS)

### Object compilation rules
fftp.o : fftp.f90
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -c fftp.f90
fprecision_mod.o : fprecision_mod.f90
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -c fprecision_mod.f90
fftp_mod.o : fftp_mod.f90
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -c fftp_mod.f90
fftp2D.o : fftp2D.f90
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -c fftp2D.f90
pseudospec2D_mod.o : pseudospec2D_mod.f90
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -c pseudospec2D_mod.f90
pseudospec2D_inc.o : pseudospec2D_inc.f90
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -c pseudospec2D_inc.f90
mainHD2D.o : mainHD2D.f90
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -c mainHD2D.f90

### make clean
clean: 
	$(RM) $(EXE) *.o *.mod

### make allclean
cleanall:
	$(RM) $(EXE) *.o *.mod *~ *.txt *.dat

### make modify
modify:
	if [ ${UNDER} -eq 1 ]; then \
	$(PERL) -i.bak -pe 's/fftw_f77\(/fftw_f77${APPEND}\(/g' *.f90 ; \
	for item in `cat external`; do \
	  `echo '$(PERL) -i -pe 's/$$item/$${item}${APPEND}/g' *.f90'`;\
	done; fi
	if [ ${EXTEN} -eq 1 ]; then \
        for file in *.f90; do \
          if [ ${UNDER} -eq 0 ]; then \
             cp $$file $${file}.bak; \
          fi; \
          mv $$file $${file%.f90}.$(F90EXT); \
        done; fi

### make undo
undo:
	rm -f *.f
	if [ -e fftp.f90.bak ]; then \
        for file in *.f90.bak; do \
          mv $$file $${file%.bak}; \
        done; fi
