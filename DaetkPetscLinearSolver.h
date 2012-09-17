#ifndef DAETKPETSCLINEARSOLVER_H
#define DAETKPETSCLINEARSOLVER_H
#include "LinearSolver.h"
#include "Preconditioner.h"
#include "DataCollector.h"
#include "VectorNorm.h"
#include "DaetkPetscSys.h"
#include "DaetkPetscMat.h"
#include "Vec.h"
#include "IntVec.h"
#include <vector>

namespace Daetk 
{
namespace Petsc
{
  namespace cc
    {
      //mwf 090104 PETSc 2.2.0 got rid of SLES completely
      //struct _p_SLES;
  struct _p_PC;
  struct _p_KSP;
  struct _p_PetscViewer;
      struct _p_IS;
    }
class LinearSolver : public Daetk::LinearSolver, public Daetk::Preconditioner
{
public:
  LinearSolver();
  LinearSolver(DataCollector& dataIn);
  LinearSolver(Petsc::Mat& matIn, DataCollector& dataIn, VectorNorm& normIn);
  LinearSolver(Petsc::Mat& matIn, Petsc::Mat& precIn, DataCollector& dataIn,VectorNorm& normIn);
  virtual ~LinearSolver();
  //int KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal atol,PetscReal dtol,int maxits)
  void useInitialGuess();
  void setTolerances(real rtol=1.e-5, real atol=1.e-5, real dtol=1.e5, int maxits=100000);
  void useScaling(bool scale=true);
  void attachMat(Petsc::Mat& matIn, bool clear=false);
  void attachSerialMat(Petsc::Mat& matIn, bool clear=false);
  void attachPreconditioner(Petsc::Mat& precIn);
  void attachData(DataCollector& dataIn);
  bool prepare(); 
  bool solve(const Vec& bIn,Vec& xIn);
  bool apply(const Vec& x, Vec& Mx);
  void printMatrices(const char* filename);
  void calculateCondition(DataCollector& d);
  void setUp(Vec& v);

  void useTimeIntegrationConvergenceTest(bool forceOneIterationIn=true);
  void useCoarseSolverConvergenceTest(bool forceOneIterationIn=true);
  void useFixedIterationConvergence();
  bool forceOneIteration();
  //Krylov subspace methods
  enum Method {RICHARDSON, CHEBYSHEV, CG, GMRES, TCQMR, BCGS, CGS, 
               TFQMR, CR, LSQR, BICG, PREONLY};

  void setMethod(Method m=GMRES);

  void setRichardsonScale(double damping_factor=1.0); 
  void setChebyshevEigenvalues(double emax=0.01,double emin=100.0); 
  void setGmresRestart(int max_steps=30);
  enum Orthogonalization {MODIFIED,UNMODIFIED};
  void setGmresOrthogonalization(Orthogonalization m=MODIFIED);

  //how to apply preconditioner
  enum Preconditioning {LEFT,RIGHT,SYMMETRIC};
  void setPreconditioning(Preconditioning p=LEFT);
  void usePreconditionedResidual();

  //type of preconditioner
  enum Preconditioner { NONE, JACOBI, SOR, LU, SHELL, BJACOBI, MG, EISENSTAT,
                        ILU, ICC, ASM, SLES, COMPOSITE, REDUNDANT, SPAI, MILU, NN, CHOLESKY};
  void setRedundantPreconditioner(Preconditioner p);
  void setPreconditioner(Preconditioner p);

  struct Precon
  {
    Err ierr;
    Petsc::cc::_p_PC* pc;
  };

  struct Ilu : public Precon 
  {
    void setLevels(int levels);
    void reuseOrdering();
    void setDropTolerance(double dt, double dtcol, int dtcount);
    void dropToleranceReuseFill();
    void inPlace();
    void allowDiagonalFill();
  };

  struct Sor: public Precon
  {
    void setOmega(double omega);
    void setIterations(int itsIn);
    
    enum Sweep {FORWARD_SWEEP, BACKWARD_SWEEP, SYMMETRIC_SWEEP, LOCAL_FORWARD_SWEEP ,
              LOCAL_BACKWARD_SWEEP ,LOCAL_SYMMETRIC_SWEEP};
    
    void setSymmetric(Sweep s);
  };

  struct BlockJacobi: public Precon
  {
    void setTotalBlocks(int blocks,int* size);
  };
  
  class AdditiveSchwarz: public Precon
  {
  public:
    AdditiveSchwarz();
    ~AdditiveSchwarz();
    void setTotalSubdomains(int n,int* indices);
    enum IntergridTransfer {BASIC, INTERPOLATE, RESTRICT, NONE};
    void setType(IntergridTransfer i);
    void setOverlap(int overlap);
    typedef std::vector<int> IndexSet;
    void setSubdomainIndexSets(std::vector< IndexSet  >& nodes, int blocksize=1);
    void setSubdomainSolvers();
    Petsc::cc::_p_IS** is;
    //mwf 090104 PETSc 2.2.0 got rid of SLES completely
    //mwf was Petsc::cc::_p_SLES**   sleses;
    Petsc::cc::_p_KSP** sleses;
    int nsdLocal;
  };
  
  struct Multigrid: public Precon
  {
    //setLevels(int levels);
//      enum multigrid {MULTIPLICATIVE, ADDITIVE, FULL, KASKADE}
//      setType(multigrid m);
//      setCycles(int cycles);
//      setNumberSmoothUp(int m);
//      setNumberSmoothDown(int n);
  };
  
  Ilu ilu;
  Sor sor;
  BlockJacobi blockJacobi;
  AdditiveSchwarz additiveSchwarz;
  Multigrid multigrid;
private:
  Err ierr;
  bool PRINT_MATRICES,CALCULATE_CONDITION, FORCE_ONE_ITERATION, SCALE_SYSTEM;
  int its,eigenvalueArrayDim;
  real *eigenvalueRealPart,
    *eigenvalueComplexPart;
  AttacheVec xAttache,bAttache,scaleAttache;
  Petsc::Mat* mat,*prec;
  DataCollector* data;
  VectorNorm* norm;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf was Petsc::cc::_p_SLES* sles;
  Petsc::cc::_p_KSP* sles;
  Petsc::cc::_p_PC* pc;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles
  //Petsc::cc::_p_KSP* ksp;
  Petsc::cc::_p_PetscViewer* viewer;
};

}//Petsc
}//Daetk
#endif
