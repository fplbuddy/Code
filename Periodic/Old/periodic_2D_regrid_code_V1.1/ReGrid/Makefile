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
IFORT = ifort

### Math Libraries
## For jxb
#FFTWLDIR = /home/vdallas/FFTWlibv2/lib
#LIBS = -limf -lmpi -L $(FFTWLDIR) -lsrfftw -lsfftw

### Optimisation compiler options
#ifdef OPT
#OPTIFORT = -p -ftrace -C hopt

## IA-32 & Intel 64
#OPTIFORT = -traceback -xP -O3 -Zp16 -align 
#OPTIFORT = -traceback -ftz
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
EXTEN    = 0
APPEND   = _
F90EXT   = f90

### Object files
OBJS = regrid_slices.o

### Exe file
EXE = run.exe

### Forced remove
RM = rm -f

### Compiles multiple source files and links them together to an executable file named 'run'.
all : $(OBJS)
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -o $(EXE) $(IFORTFLUGS) $(OBJS) $(LIBS)

### Object compilation rules
regrid_slices.o : regrid_slices.f90
	$(IFORT) $(DBGFLUGS) $(OPTIFORT) -c regrid_slices.f90

### make clean
clean: 
	$(RM) $(EXE) *.o *.mod

### make allclean
cleanall:
	$(RM) $(EXE) *.o *.mod *~ *.float

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
