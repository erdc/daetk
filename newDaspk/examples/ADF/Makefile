#
# This makefile compiles and loads the DDASPK example program for SGI
# machine.
# If necessary, change the constants COMP and FFLAGS below for the
# compiler to be used.

COMP = f77
FFLAGS = -O 

SOLVR = ../../solver

PRECON = ../../preconds

#--------------------------------------------------------------------
# DASPK using ADIFOR packages (no MPI packages)
#
OBJS = $(SOLVR)/ddaspk.o $(SOLVR)/daux.o $(SOLVR)/dlinpk.o \
       $(SOLVR)/mpi_dummy.o

# only for IRIX -o32 compiling option. For others, remove -o32 from 
# the following lines.
ADLIB    = $(AD_LIB)/lib/ReqADIntrinsics-$(AD_OS).o \
	   $(AD_LIB)/lib/libADIntrinsics-$(AD_OS).a \
	   $(AD_LIB)/lib/libSparsLinC-$(AD_OS).a 

# Adifor preprocessor command
ADF   = Adifor

# derived informations for heat equation and ILU preconditioner
SENSD = $(SOLVR)/dsensd.o reshiludr.o


# Adifor-generated routines for different examples
#
AD_HEATD   = res4sen/h_resh.o res4sen/g_resh.o jacsp/i_resh.o jac/j_resh.o
AD_HEATILU = res4sen/h_reshilu.o res4sen/g_reshilu.o jacsp/i_reshilu.o \
	     res4sen/h_reshiludr.o
AD_WEBD = res4sen/h_resweb.o res4sen/g_resweb.o jacsp/i_resweb.o
AD_WEBILU = res4sen/h_reswebilu.o res4sen/g_reswebilu.o jacsp/i_reswebilu.o
AD_PEND = res4sen/h_respend.o jacsp/i_respend.o jac/j_respend.o
AD_RAHA = res4sen/h_resraha.o jacsp/i_resraha.o

HEAT = dheat.o $(AD_HEATD) $(OBJS) $(PRECON)/dbanpre.o

HEATILU = dheatilu.o $(AD_HEATILU) $(OBJS) $(SENSD) \
          $(PRECON)/dilupre_adf.o $(PRECON)/dsparsk.o

WEB = dweb.o $(AD_WEBD) $(OBJS) $(PRECON)/drbdpre.o $(PRECON)/drbgpre.o 

WEBILU = dwebilu.o $(AD_WEBILU) $(OBJS) $(PRECON)/dilupre_adf.o \
         $(PRECON)/dsparsk.o 

PEND = pend.o $(AD_PEND) $(OBJS)

heat :  $(HEAT) 
	$(COMP) $(FFLAGS) -o heat $(HEAT) $(ADLIB)  -lm

res4sen/h_resh.f : resh.f
	$(ADF) AD_SCRIPT=resjac_heat.adf; \
	$(ADF) AD_SCRIPT=res4sensm_heat.adf; \
	$(ADF) AD_SCRIPT=res4senmv_heat.adf; \
	$(ADF) AD_SCRIPT=resjacsp_heat.adf

heatilu: $(HEATILU)
	$(COMP) $(FFLAGS) -o heatilu $(HEATILU) $(ADLIB)  -lm

res4sen/h_reshilu.f : reshilu.f
	$(ADF) AD_SCRIPT=res4sensm_heat.adf AD_PROG=heatilu.cmp; \
	$(ADF) AD_SCRIPT=res4senmv_heat.adf AD_PROG=heatilu.cmp; \
	$(ADF) AD_SCRIPT=resjacsp_heat.adf  AD_PROG=heatilu.cmp
res4sen/h_reshiludr.f : reshiludr.f
	$(ADF) AD_SCRIPT=res4sensm_heatdr.adf AD_PROG=heatiludr.cmp
web :  $(WEB)
	$(COMP) $(FFLAGS) -o web $(WEB) $(ADLIB)  -lm 

res4sen/h_resweb.f : resweb.f
	$(ADF) AD_SCRIPT=res4sensm_web.adf; \
	$(ADF) AD_SCRIPT=res4senmv_web.adf; \
	$(ADF) AD_SCRIPT=resjacsp_web.adf 

webilu : $(WEBILU)
	$(COMP) $(FFLAGS) -o webilu $(WEBILU) $(ADLIB)  -lm

res4sen/h_reswebilu.f : reswebilu.f
	$(ADF) AD_SCRIPT=res4sensm_web.adf; \
	$(ADF) AD_SCRIPT=resjacsp_web.adf 

pend :  $(PEND)
	$(COMP) $(FFLAGS) -o pend $(PEND) $(ADLIB)  -lm 

res4sen/h_respend.f : respend.f
	$(ADF) AD_SCRIPT=res4sensm_pend.adf; \
	$(ADF) AD_SCRIPT=resjac_pend.adf; \
	$(ADF) AD_SCRIPT=resjacsp_pend.adf 

clean : 
	@rm -rf core *~

new :
	@rm -rf *.o heat heatilu web webilu jac jacsp res4sen AD_cache; \
         cd ../../solver; rm -rf *.o *~ \
	 cd ../../preconds; rm -rf *.o *~



# Rule for compiling a Fortran source file:
.f.o: ;	$(COMP) $(FFLAGS) -c $*.f -o $*.o
