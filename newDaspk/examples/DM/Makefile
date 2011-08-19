#
# This makefile compiles and loads the DDASPK example program for SGI
# machine.
# If necessary, change the constants COMP and FFLAGS below for the
# compiler to be used.

COMP = gfortran #f77
FFLAGS = -O 

SOLVR = ../../solver

PRECON = ../../preconds

#----------------------------------------------------------------
# DASPK without involving the ADIFOR and MPI packages
#
OBJS = $(SOLVR)/ddaspk.o $(SOLVR)/daux.o $(SOLVR)/dlinpk.o \
       $(SOLVR)/mpi_dummy.o $(SOLVR)/adf_dummy.o
 
# derived informations for heat equation and ILU preconditioner
SENSD = $(SOLVR)/dsensd.o reshiludr.o

HEAT = dheat.o $(OBJS) $(PRECON)/dbanpre.o


HEATILU = dheatilu.o $(OBJS) $(SENSD) $(PRECON)/dilupre.o \
          $(PRECON)/dsparsk.o 

WEB = dweb.o $(OBJS) $(PRECON)/drbdpre.o $(PRECON)/drbgpre.o 

WEBILU = dwebilu.o $(OBJS) $(PRECON)/dilupre.o $(PRECON)/dsparsk.o 

heat :  $(HEAT) 
	$(COMP) $(FFLAGS) -o heat $(HEAT)

# test with sensitivity to the derived information
heatilu: $(HEATILU)
	$(COMP) $(FFLAGS) -o heatilu $(HEATILU) 

web :  $(WEB) 
	$(COMP) $(FFLAGS) -o web $(WEB)

webilu :  $(WEBILU) 
	$(COMP) $(FFLAGS) -o webilu $(WEBILU)

clean : 
	@rm -rf core *~

new :
	@rm -rf *.o heat heatilu web webilu; 


# Rule for compiling a Fortran source file:
.f.o: ;	$(COMP) $(FFLAGS) -c $*.f -o $*.o
