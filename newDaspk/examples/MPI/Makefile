#
# This makefile compiles and loads the DDASPK example program for SGI
# machine.
# If necessary, change the constants COMP and FFLAGS below for the
# compiler to be used.

.SUFFIXES : .f .c

COMP = f77
FFLAGS = -O 

SOLVR = ../../solver

PRECON = ../../preconds

#-------------------------------------------------------------------
# DASPK using MPI packages (without ADIFOR packages)
#
OBJS = $(SOLVR)/ddaspk.o $(SOLVR)/daux.o $(SOLVR)/dlinpk.o \
	$(SOLVR)/adf_dummy.o

# MPI fortran compile comands
MPIF = ${MPIF77}

#MPI compile options
MPIFLAGS = -lmpi

HEAT = dheatmpi.o $(OBJS) $(PRECON)/dbanpre.o

HEATILU = dheatilumpi.o $(OBJS) \
          $(PRECON)/dilupre.o $(PRECON)/dsparsk.o

heat : $(HEAT)
	$(MPIF) -o heat $(HEAT)
heatilu : $(HEATILU) 
	$(MPIF) -o heatilu $(HEATILU)

clean : 
	@rm -rf core *~

new :
	@rm -rf *.o heat heatilu web webilu; 


# Rule for compiling a Fortran source file:
.f.o: ;	$(MPIF) -c $*.f -o $*.o
