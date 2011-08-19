SRCS  = BrooksCorey2p.cpp PskHistory_Vec.cpp  \
	BrooksCorey2pSpline.cpp ParkerHysteresis.cpp \
	ConstBC.cpp SlCompress.cpp PetscSecondOrderFd.cpp\
	Stencil.cpp GasDensity.cpp PetscStencilMM.cpp TestPostProc.cpp   \
	Poisson2p.cpp VanGenuchten2p.cpp                  \
	MualemVanGenuchten2p.cpp Psk2p.cpp                           \
	MualemVanGenuchten2pSpline.cpp PskHistory.cpp SaturatedPsk2p.cpp \
        TwopCDM.cpp TwopData.cpp GerhardKueper2p.cpp TwopUtils.cpp BiodegradationData.cpp

#ConstBCECDM.cpp 

OBJS =  BrooksCorey2p.o PskHistory_Vec.o  \
	BrooksCorey2pSpline.o ParkerHysteresis.o \
	ConstBC.o SlCompress.o PetscSecondOrderFd.o\
	Stencil.o GasDensity.o PetscStencilMM.o TestPostProc.o   \
	Poisson2p.o VanGenuchten2p.o                  \
	MualemVanGenuchten2p.o Psk2p.o                           \
	MualemVanGenuchten2pSpline.o PskHistory.o SaturatedPsk2p.o \
        TwopCDM.o TwopData.o GerhardKueper2p.o TwopUtils.o BiodegradationData.o

#ConstBCECDM.o 
