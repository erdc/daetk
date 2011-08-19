SRCS  = AnalyticalJacobian.cpp CMRVecIndex.cpp  Integrator.cpp Jacobian.cpp \
	LinearOperator.cpp LinearSolver.cpp mvvind.cpp NonlinearSolver.cpp  \
	NumericalJacobian.cpp Preconditioner.cpp TimingDataFile.cpp       \
	VectorNorm.cpp DaeDefinition.cpp Chronograph.cpp DaspkPsol.cpp      \
	FullDataFile.cpp LaFullCholeskySolver.cpp Vec.cpp BandColMat.cpp    \
	DaspkRes.cpp FullDirectSolver.cpp LaFullDirectSolver.cpp          \
	mvm.cpp BiCGstab.cpp DASPK.cpp DataCollector.cpp     \
	GlobalDaeDef.cpp LiBandedDirectSolver.cpp PreCG.cpp               \
	WeightedRMSNorm.cpp mvv.cpp CMRVecBlasd.cpp  DataFile.cpp           \
	ParameterDatabase.cpp IntVec.cpp ModifiedNewton.cpp Newton.cpp      \
	PTC.cpp PTCDAE.cpp CMRVecBlass.cpp DaspkJac.cpp FLCBDF.cpp            \
	LaBandedDirectSolver.cpp mvblas.cpp daux.f dlinpk.f ddaspk.f    \
	VectorFunction.cpp DaetkPetscVec.cpp DaetkPetscSys.cpp            \
	DaetkPetscVecBlas.cpp DaetkPetscMat.cpp                         \
	DaetkPetscLinearSolver.cpp DaetkPetscNumericalJacobian.cpp      \
	DaetkPetscAnalyticalJacobian.cpp TexDataFile.cpp FLCBDF_lite.cpp \
        FullDataFileForRichExtrap.cpp 

OBJS =  AnalyticalJacobian.o CMRVecIndex.o  Integrator.o Jacobian.o \
	LinearOperator.o LinearSolver.o mvvind.o NonlinearSolver.o  \
	NumericalJacobian.o Preconditioner.o TimingDataFile.o       \
	VectorNorm.o DaeDefinition.o Chronograph.o DaspkPsol.o      \
	FullDataFile.o LaFullCholeskySolver.o Vec.o BandColMat.o    \
	DaspkRes.o FullDirectSolver.o LaFullDirectSolver.o          \
	mvm.o BiCGstab.o DASPK.o DataCollector.o     \
	GlobalDaeDef.o LiBandedDirectSolver.o PreCG.o               \
	WeightedRMSNorm.o mvv.o CMRVecBlasd.o  DataFile.o           \
	ParameterDatabase.o IntVec.o ModifiedNewton.o Newton.o      \
	PTC.o PTCDAE.o CMRVecBlass.o DaspkJac.o FLCBDF.o            \
	LaBandedDirectSolver.o mvblas.o daux.o dlinpk.o ddaspk.o    \
	VectorFunction.o DaetkPetscVec.o DaetkPetscSys.o            \
	DaetkPetscVecBlas.o DaetkPetscMat.o                         \
	DaetkPetscLinearSolver.o DaetkPetscNumericalJacobian.o      \
	DaetkPetscAnalyticalJacobian.o TexDataFile.o FLCBDF_lite.o  \
        FullDataFileForRichExtrap.o
