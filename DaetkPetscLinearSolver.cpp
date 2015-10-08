#include "DaetkPetscLinearSolver.h"
#include <cstdio>
#include "Tracer.h"
//user defined Prec
//     ierr = PCShellSetApply(PC pc,int (*apply)(void *ctx,Vec,Vec),void *ctx); 
//     ierr = PCShellSetSetUp(PC pc,int (*setup)(void *ctx)); 

namespace Daetk 
{
namespace Petsc
{
  namespace cc
  {
    extern "C"
    {
#ifndef DAETK_DEF_CPLUSPLUS_FOR_PETSC_H
#undef __cplusplus
#endif
#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
      //mwf 090104 PETSc 2.2.0 got rid of SLES completely
      //mwf was #include "petscsles.h"
#include "petscksp.h"
#include "petscpc.h" //cek added for development version of PETSC
#ifndef DAETK_DEF_CPLUSPLUS_FOR_PETSC_H
#define __cplusplus
#endif
    }
  }


using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ostream;
using std::istream;

  void LinearSolver::useScaling(bool scale){SCALE_SYSTEM=scale;}

  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles below

LinearSolver::LinearSolver():
  PRINT_MATRICES(0),
  CALCULATE_CONDITION(0),
  FORCE_ONE_ITERATION(false),
  SCALE_SYSTEM(false),
  xAttache(this),
  bAttache(this),
  scaleAttache(this),
  mat(0),
  prec(0),
  data(0),
  norm(0),
  sles(0),
  pc(0)
  //,ksp(0)
{}
  
LinearSolver::LinearSolver(DataCollector& dataIn):
  PRINT_MATRICES(0),
  CALCULATE_CONDITION(0),
  FORCE_ONE_ITERATION(false),
  SCALE_SYSTEM(false),
  xAttache(this),
  bAttache(this),
  scaleAttache(this),
  mat(0),
  prec(0),
  data(&dataIn),
  norm(0),
  sles(0),
  pc(0)
  //,ksp(0)
{}

LinearSolver::LinearSolver(Mat& matIn, DataCollector& dataIn, VectorNorm& normIn):
  PRINT_MATRICES(0),
  CALCULATE_CONDITION(0),
  FORCE_ONE_ITERATION(false),
  SCALE_SYSTEM(false),
  xAttache(this,matIn.referenceVec.dim()),
  bAttache(this,matIn.referenceVec.dim()),
  scaleAttache(this,matIn.referenceVec.dim()),
  mat(&matIn),
  prec(&matIn),
  data(&dataIn),
  norm(&normIn),
  sles(0),
  pc(0)
  //,ksp(0)
{
  using namespace cc;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replaced SLES prefix with KSP
  ierr = KSPCreate(PETSC_COMM_WORLD,&sles);  
  ierr = KSPSetFromOptions(sles); 
  ierr = KSPGetPC(sles,&pc); 
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf redundant now 
  //ierr = SLESGetKSP(sles,&ksp); 
}
 
LinearSolver::LinearSolver(Mat& matIn, Mat& precIn, DataCollector& dataIn,VectorNorm& normIn):
  PRINT_MATRICES(0),
  CALCULATE_CONDITION(0),
  SCALE_SYSTEM(false),
  xAttache(this,matIn.referenceVec.dim()),
  bAttache(this,matIn.referenceVec.dim()),
  scaleAttache(this,matIn.referenceVec.dim()),
  mat(&matIn),
  prec(&precIn),
  data(&dataIn),
  norm(&normIn),
  sles(0),
  pc(0)
  //,ksp(0)
{
  using namespace cc;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replaced SLES prefix with KSP
  ierr = KSPCreate(PETSC_COMM_WORLD,&sles);  
  ierr = KSPSetFromOptions(sles); 
  ierr = KSPGetPC(sles,&pc); 
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf redundant now
  //ierr = SLESGetKSP(sles,&ksp); 
}
 
LinearSolver::~LinearSolver()
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replaced SLES prefix with KSP
  if (sles)
    ierr = cc::KSPDestroy(&sles);
}

void LinearSolver::attachMat(Mat& matIn, bool clear)
{
  using namespace cc;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replaced SLES prefix with KSP
  if (clear ||  !sles)
    {
      if (sles)
        ierr = KSPDestroy(&sles);
      ierr = KSPCreate(PETSC_COMM_WORLD,&sles);  
      ierr = KSPGetPC(sles,&pc); 
      //mwf 090104 PETSc 2.2.0 got rid of SLES completely
      //mwf redundant now
      //ierr = SLESGetKSP(sles,&ksp); 
    }
  mat  = &matIn;
  xAttache.v_.newsize(mat->referenceVec.dim());
  bAttache.v_.newsize(mat->referenceVec.dim());
  scaleAttache.v_.newsize(mat->referenceVec.dim());
}

void LinearSolver::attachSerialMat(Mat& matIn, bool clear)
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replaced SLES prefix with KSP
  using namespace cc;
  if (clear ||  !sles)
    {
      if (sles)
        ierr = KSPDestroy(&sles);
      ierr = KSPCreate(PETSC_COMM_SELF,&sles);  
      ierr = KSPGetPC(sles,&pc); 
      //mwf 090104 PETSc 2.2.0 got rid of SLES completely
      //mwf redundant now
      //ierr = SLESGetKSP(sles,&ksp); 
    }
  mat  = &matIn;
  xAttache.v_.newsizeSerial(mat->referenceVec.dim());
  bAttache.v_.newsizeSerial(mat->referenceVec.dim());
  scaleAttache.v_.newsizeSerial(mat->referenceVec.dim());
}

void LinearSolver::attachPreconditioner(Mat& precIn)
{
  prec  = &precIn;
}

void LinearSolver::attachData(DataCollector& dataIn)
{
  data = &dataIn;
}

bool LinearSolver::prepare()
{ 
  using namespace cc;
  int i;

  if (PRINT_MATRICES)
    MatView(prec->castToPetsc(),viewer);

  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replaced SLES with KSP 

  i = KSPSetOperators(sles,mat->castToPetsc(),prec->castToPetsc());
  ierr = i;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace with KSPSetUp(ksp) following man pages? web page is wrong
  //ierr = SLESSetUp(sles,mat->referenceVec.castToPetsc(),prec->referenceVec.castToPetsc());
  //mwf nonsuch ierr = KSPSetRhs(sles,mat->referenceVec.castToPetsc());
  //mwf nonsuch ierr = KSPSetSolution(sles,prec->referenceVec.castToPetsc());
  ierr = KSPSetUp(sles);

  return i;
}
  
//int KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal atol,PetscReal dtol,int maxits)
void LinearSolver::setTolerances(real rtol, real atol, real dtol,int maxits)
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles
  cc::KSPSetTolerances(sles,rtol,atol,dtol,maxits);
}

bool LinearSolver::solve(const Vec& bIn,Vec& xIn)
{
  using namespace cc;
  bAttache.attachToTarget(bIn);
  xAttache.attachToTarget(xIn);

  //  int neig;

  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles
  // if (CALCULATE_CONDITION)
  //ierr = KSPSetComputeEigenvalues(sles);

  KSPType ksptype;
  const KSPType ksptype2=KSPPREONLY;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles 
  ierr = KSPGetType(sles,&ksptype);
  if (SCALE_SYSTEM && strcmp(ksptype,ksptype2))
    {
      if (norm)
        {
          scaleAttache.attachToTarget(norm->getScaling());
          PCSetDiagonalScale(pc,const_cast<_p_Vec*>(scaleAttache.v_.castToConstPetsc()));
        }
      else
        {
          cerr<<"You must supply a norm with scalings if you wish to scale the system in PetscLinearSolver"<<endl
              <<"Continuing without scaling"<<endl;
          SCALE_SYSTEM = false;
        }
    }
//      PCDiagonalScaleSet(pc,const_cast<_p_Vec*>(norm->getScaling().castToConstPetsc()));

  bAttache.v_.restoreLocal();
  xAttache.v_.restoreLocal();

  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace with KSPSetRhs(ksp,b),etc following web directions
  //ierr = SLESSetUp(sles,const_cast<_p_Vec*>(bAttache.v_.castToConstPetsc()),
  //xAttache.v_.castToPetsc());
  //mwf Petsc web page wrong, try following web pages
  //mwf nonsuc ierr = KSPSetRhs(sles,const_cast<_p_Vec*>(bAttache.v_.castToConstPetsc()));
  //mwf nonsuch ierr = KSPSetSolution(sles,xAttache.v_.castToPetsc());
  ierr = KSPSetUp(sles);

  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace with KSPSolve(ksp,b,x),etc following web directions
  //ierr  = SLESSolve(sles,const_cast<_p_Vec*>(bAttache.v_.castToConstPetsc()),
  //xAttache.v_.castToPetsc(),&its); 
  ierr   = KSPSolve(sles,
		    const_cast<_p_Vec*>(bAttache.v_.castToConstPetsc()),
		    xAttache.v_.castToPetsc());
  ierr   = KSPGetIterationNumber(sles,&its);


  data->setLinearSolverIterations(its);
  xAttache.v_.getLocal();
  bAttache.v_.getLocal();

  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  

  //    norm->deScale(x,x);
  //if (CALCULATE_CONDITION)
  //ierr = KSPComputeEigenvalues(sles, eigenvalueArrayDim,eigenvalueRealPart,eigenvalueComplexPart,&neig);
  //cout<<"out of solve"<<endl<<flush;
//    Vec tmp(x);
//    x.restoreLocal();
//    tmp.restoreLocal();
//    MatMult(mat->castToPetsc(),x.castToPetsc(),tmp.castToPetsc());
//    x.getLocal();
//    tmp.getLocal();
//    tmp-=b;
//    cout<<norm(tmp)/norm(x)<<'\t'<<norm(tmp)/norm(b)<<endl;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  
  KSPConvergedReason reason;
  KSPGetConvergedReason(sles,&reason);

  
  if (int(reason) > 0)
    {
      if (SOLVE_SUB)
        xIn=bIn;
      xAttache.restoreToTarget();
      return 0;
    }
  else
    {
      //mwf 090104 PETSc 2.2.0 got rid of SLES completely
      //mwf replace ksp with sles  
	KSPBuildSolution(sles,xAttache.v_.castToPetsc(),(Daetk::Petsc::cc::Vec*)(NULL));
      if (SOLVE_SUB)
        xIn=bIn;
      xAttache.restoreToTarget();
      //std::cerr<<"Petsc Iterative Method failed with reason="<<reason<<std::endl;
      return 1;
    }
}

bool LinearSolver::apply(const Vec& x, Vec& Mx)
{
  using namespace cc;
  x.restoreLocal();
  Mx.restoreLocal();
  //mwf 090104 PETSc 2.2.1 got rid of PCSetVector
  //ierr =  PCSetVector(pc,Mx.castToPetsc());
  ierr =  PCApply(pc,const_cast<_p_Vec*>(x.castToConstPetsc()),Mx.castToPetsc());
  x.getLocal();
  Mx.getLocal();
  return 0;
}

void LinearSolver::setRichardsonScale(double damping_factor)
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  
  ierr = cc::KSPRichardsonSetScale(sles,damping_factor); 
}

void LinearSolver::setChebyshevEigenvalues(double emax,double emin)
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  
  ierr = cc::KSPChebyshevSetEigenvalues(sles,emax,emin); 
}

void LinearSolver::setGmresRestart(int max_steps)
{
  std::cerr<<"setGmresRestart doesn't work with Petsc 2.1.1"<<std::endl;
  exit(1);
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  
  //  using namespace cc;
  //ierr = KSPGMRESSetRestart(sles,max_steps); 
}

void LinearSolver::setGmresOrthogonalization(Orthogonalization m)
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  
  //mwf changed UnmodifiedGramSchmidt to ClassicalGramSchmidt
  using namespace cc;
  if (m==MODIFIED)
    ierr = KSPGMRESSetOrthogonalization(sles, 
                                        KSPGMRESModifiedGramSchmidtOrthogonalization); 
  else
    ierr = KSPGMRESSetOrthogonalization(sles, 
                                        KSPGMRESClassicalGramSchmidtOrthogonalization); 
}

void LinearSolver::setMethod(Method m)
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  
  using namespace cc;
  switch (m)
    {
    case RICHARDSON: 
      ierr = KSPSetType(sles,KSPRICHARDSON); 
      break; 
    case CHEBYSHEV: 
      ierr = KSPSetType(sles,KSPCHEBYSHEV); 
      break; 
    case CG: 
      ierr = KSPSetType(sles,KSPCG); 
      break; 
    case GMRES: 
      ierr = KSPSetType(sles,KSPGMRES); 
      break; 
    case TCQMR: 
      ierr = KSPSetType(sles,KSPTCQMR); 
      break; 
    case BCGS: 
      ierr = KSPSetType(sles,KSPBCGS); 
      break; 
    case CGS: 
      ierr = KSPSetType(sles,KSPCGS); 
      break; 
    case TFQMR: 
      ierr = KSPSetType(sles,KSPTFQMR); 
      break; 
    case CR: 
      ierr = KSPSetType(sles,KSPCR); 
      break; 
    case LSQR: 
      ierr = KSPSetType(sles,KSPLSQR); 
      break; 
    case BICG: 
      ierr = KSPSetType(sles,KSPBICG); 
      break; 
    case PREONLY:
      ierr = KSPSetType(sles,KSPPREONLY);
      break;
    }
}
  
void LinearSolver::setPreconditioning(Preconditioning p)
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  
  using namespace cc;
  switch(p)
    {
    case LEFT:
      ierr = KSPSetPCSide(sles,PC_LEFT);
      break; 
    case RIGHT:
      ierr = KSPSetPCSide(sles,PC_RIGHT); 
      break; 
    case SYMMETRIC:
      ierr = KSPSetPCSide(sles,PC_SYMMETRIC);
      break; 
    }
}
  
void LinearSolver::setPreconditioner(Preconditioner p)
{
  using namespace cc;
  switch(p)
    {
    case NONE:
      PCSetType(pc,PCNONE); 
      break; 
    case JACOBI: 
      PCSetType(pc,PCJACOBI); 
      break; 
    case SOR: 
      PCSetType(pc,PCSOR); 
      break; 
    case LU: 
      PCSetType(pc,PCLU); 
      break; 
    case SHELL: 
      PCSetType(pc,PCSHELL); 
      break; 
    case BJACOBI: 
      PCSetType(pc,PCBJACOBI);
      blockJacobi.pc=pc;
      break; 
    case MG: 
      PCSetType(pc,PCMG);
      multigrid.pc=pc;
      break; 
    case EISENSTAT: 
      PCSetType(pc,PCEISENSTAT); 
      break; 
    case ILU:
      PCSetType(pc,PCILU); 
      ilu.pc=pc;
      break; 
    case ICC: 
      PCSetType(pc,PCICC); 
      break; 
    case ASM: 
      PCSetType(pc,PCASM); 
      additiveSchwarz.pc = pc;
      break; 
    case SLES: 
      //mwf 090104 PETSc got rid of this sometime
      //PCSetType(pc,PCSLES); 
      std::cerr<<"PCSLES no longer implemented!!"<<std::endl;
      assert(0);
      break; 
    case COMPOSITE: 
      PCSetType(pc,PCCOMPOSITE); 
      break; 
    case REDUNDANT: 
      PCSetType(pc,PCREDUNDANT); 
      break; 
    case SPAI: 
      PCSetType(pc,PCSPAI); 
      break; 
    case MILU: 
      //mwf 090104 PETSc got rid of this sometime
      //PCSetType(pc,PCMILU); 
      std::cerr<<"PCMILU no longer implemented!!"<<std::endl;
      assert(0);
      break; 
    case NN: 
      PCSetType(pc,PCNN); 
      break; 
    case CHOLESKY:
      PCSetType(pc,PCCHOLESKY); 
      break; 
    }
}

//this is sloppy
void LinearSolver::setRedundantPreconditioner(Preconditioner p)
{
  using namespace cc;
  _p_KSP* kspr;
  _p_PC* pcr;
  PCRedundantGetKSP(pc,&kspr);
  KSPGetPC(kspr,&pcr);
  switch(p)
    {
    case NONE:
      PCSetType(pcr,PCNONE); 
      break; 
    case JACOBI: 
      PCSetType(pcr,PCJACOBI); 
      break; 
    case SOR: 
      PCSetType(pcr,PCSOR); 
      break; 
    case LU: 
      PCSetType(pcr,PCLU); 
      break; 
    case SHELL: 
      PCSetType(pcr,PCSHELL); 
      break; 
    case BJACOBI: 
      PCSetType(pcr,PCBJACOBI);
      blockJacobi.pc=pc;
      break; 
    case MG: 
      PCSetType(pcr,PCMG);
      multigrid.pc=pc;
      break; 
    case EISENSTAT: 
      PCSetType(pcr,PCEISENSTAT); 
      break; 
    case ILU:
      PCSetType(pcr,PCILU); 
      ilu.pc=pc;
      break; 
    case ICC: 
      PCSetType(pcr,PCICC); 
      break; 
    case ASM: 
      PCSetType(pcr,PCASM); 
      additiveSchwarz.pc = pc;
      break; 
    case SLES: 
      //mwf 090104 PETSc got rid of this sometime
      //PCSetType(pcr,PCSLES); 
      std::cerr<<"PCSLES no longer implemented!!!"<<std::endl;
      assert(0);
      break; 
    case COMPOSITE: 
      PCSetType(pcr,PCCOMPOSITE); 
      break; 
    case REDUNDANT: 
      PCSetType(pcr,PCREDUNDANT); 
      break; 
    case SPAI: 
      PCSetType(pcr,PCSPAI); 
      break; 
    case MILU: 
      //mwf 090104 PETSc got rid of this sometime
      //PCSetType(pcr,PCMILU); 
      std::cerr<<"PCMILU no longer implemented!!!"<<std::endl;
      assert(0);
      break; 
    case NN: 
      PCSetType(pcr,PCNN); 
      break; 
    case CHOLESKY:
      PCSetType(pcr,PCCHOLESKY); 
      break; 
    }
}

void LinearSolver::usePreconditionedResidual()
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles  
  using namespace cc;
//   ierr = cc::KSPSetNormType(sles,KSP_PRECONDITIONED_NORM);
  ierr = cc::KSPSetNormType(sles,KSP_NORM_PRECONDITIONED);
}

void LinearSolver::calculateCondition(DataCollector& d)
{
  data = &d;
  CALCULATE_CONDITION=true; 
  eigenvalueRealPart = new real[eigenvalueArrayDim];
  eigenvalueComplexPart = new real[eigenvalueArrayDim];
}

void  LinearSolver::Ilu::setLevels(int levels)
{
  ierr = cc::PCFactorSetLevels(pc,levels); 
}

void LinearSolver::Ilu::reuseOrdering()
{
  using namespace cc;
  ierr = PCFactorSetReuseOrdering(pc,PETSC_TRUE);
}
   
void LinearSolver::Ilu::setDropTolerance(double dt, double dtcol, int dtcount)
{
  ierr = cc::PCFactorSetDropTolerance(pc,dt,dtcol,dtcount); 
}
   
void LinearSolver::Ilu::dropToleranceReuseFill()
{
  using namespace cc;
  ierr = PCFactorSetReuseFill(pc,PETSC_TRUE); 
}

void LinearSolver::Ilu::inPlace()
{
  using namespace cc;
  ierr = PCFactorSetUseInPlace(pc,PETSC_TRUE);   
}
 
void LinearSolver::Ilu::allowDiagonalFill()
{
  using namespace cc;
  ierr = PCFactorSetAllowDiagonalFill(pc,PETSC_TRUE); 
}

void LinearSolver::Sor::setOmega(double omega)
{
  ierr = cc::PCSORSetOmega(pc,omega);
}

void LinearSolver::Sor::setIterations(int itsIn)
{
  //  ierr = cc::PCSORSetIterations(pc,itsIn);
}

void LinearSolver::Sor::setSymmetric(Sweep s)
{
  using namespace cc;
  switch(s)
    {
    case FORWARD_SWEEP: 
      ierr = PCSORSetSymmetric(pc,SOR_FORWARD_SWEEP);
      break; 
    case  BACKWARD_SWEEP: 
      ierr = PCSORSetSymmetric(pc, SOR_BACKWARD_SWEEP);
      break; 
    case  SYMMETRIC_SWEEP: 
      ierr = PCSORSetSymmetric(pc,SOR_SYMMETRIC_SWEEP);
      break; 
    case  LOCAL_FORWARD_SWEEP: 
      ierr = PCSORSetSymmetric(pc,SOR_LOCAL_FORWARD_SWEEP);
      break; 
    case LOCAL_BACKWARD_SWEEP:
      ierr = PCSORSetSymmetric(pc,SOR_LOCAL_BACKWARD_SWEEP);
      break; 
    case LOCAL_SYMMETRIC_SWEEP:
      ierr = PCSORSetSymmetric(pc,SOR_LOCAL_SYMMETRIC_SWEEP);
      break; 
    }
}

void LinearSolver::BlockJacobi::setTotalBlocks(int blocks,int* size)
{
  ierr = cc::PCBJacobiSetTotalBlocks(pc,blocks,size);
} 

    LinearSolver::AdditiveSchwarz::AdditiveSchwarz():is(0),nsdLocal(0){}

void LinearSolver::AdditiveSchwarz::setTotalSubdomains(int n,int* indices)
{
  cerr<<"not yet implemented"<<endl;
  //  ierr = PCASMSetTotalSubdomains(pc,n,IS *is); 
}

void LinearSolver::AdditiveSchwarz::setType(IntergridTransfer a)
{
  using namespace cc;
  switch (a)
    {
    case BASIC:
      ierr = PCASMSetType(pc,PC_ASM_BASIC); 
      break;
    case INTERPOLATE:
      ierr = PCASMSetType(pc,PC_ASM_INTERPOLATE); 
      break;
    case RESTRICT:
      ierr = PCASMSetType(pc,PC_ASM_RESTRICT); 
      break;
    case NONE:
      ierr = PCASMSetType(pc,PC_ASM_NONE); 
      break;
    }
}

void LinearSolver::AdditiveSchwarz::setOverlap(int overlap)
{
  ierr = cc::PCASMSetOverlap(pc,overlap); 
}

void LinearSolver::AdditiveSchwarz::setSubdomainIndexSets(std::vector< IndexSet >& nodes, int blocksize)
{
  using namespace cc;
  nsdLocal=nodes.size();
  is = new IS[blocksize*nodes.size()]; 
  for(unsigned int i=0;i<nodes.size();i++)
    {
      int* int_is = new int[nodes[i].size()];
      for (unsigned int j=0;j<nodes[i].size();j++)
        int_is[j]=nodes[i][j]*blocksize;
      
      ISCreateBlock(PETSC_COMM_SELF,blocksize,nodes[i].size(),int_is,PETSC_COPY_VALUES,&(is[i]));
      delete [] int_is;
    }
  ierr = PCASMSetLocalSubdomains(pc,nodes.size(),is,(IS*)(NULL));
}

LinearSolver::AdditiveSchwarz::~AdditiveSchwarz()
{
  for (int i=nsdLocal-1;i>=0;i--)
    cc::ISDestroy(&is[i]);
  delete [] is;
}

void LinearSolver::AdditiveSchwarz::setSubdomainSolvers()
{
  using namespace cc;
  //mwf 090104 PETSc 2.2.0 removed SLES completely
  //mwf replace sles with KSP
  ierr = PCASMGetSubKSP(pc,
                        reinterpret_cast<PetscInt*>(NULL),
                        reinterpret_cast<PetscInt*>(NULL),
                        &sleses);
  PC subpc;
  //mwf 090104 PETSc 2.2.0 removed SLES completely
  //mwf so now redundant
  //KSP subksp;
  for(int i=0;i<nsdLocal;i++)
    {
      ierr = KSPGetPC(sleses[i],&subpc); 
      //mwf 090104 PETSc 2.2.0 removed SLES completely
      //mwf so now redundant
      //ierr = SLESGetKSP(sleses[i],&subksp); 
      ierr = KSPSetType(sleses[i],KSPPREONLY);
      ierr = PCSetType(subpc,PCILU); 
    }
}

void LinearSolver::setUp(Vec& v)
{
  using namespace cc;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace SLESSetUp with KSPSetRhs(sles,b) etc from web directions
  //ierr = SLESSetUp(sles,const_cast<_p_Vec*>(v.castToConstPetsc()),
  //v.castToPetsc());
  //mwf Petsc web page wrong try to follow man pages
  //mwf nonsuch ierr = KSPSetRhs(sles,const_cast<_p_Vec*>(v.castToConstPetsc()));
  //mwf nonsuch ierr = KSPSetSolution(sles,v.castToPetsc());
  ierr = KSPSetUp(sles);
}

//  //multigrid
//     ierr = PCSetType(PC pc,PCMG); 
//     ierr = MGSetLevels(pc,int levels,MPI_Comm *comms); 
//     ierr = MGSetType(PC pc,MGType mode);   MGMULTIPLICATIVE, MGADDITIVE, MGFULL, MGKASKADE
//     ierr = MGSetCycles(PC pc,int cycles); 
//     ierr = MGSetNumberSmoothUp(PC pc,int m); 
//     ierr = MGSetNumberSmoothDown(PC pc,int n); 
//    ierr = MGGetCoarseSolve(PC pc,SLES *sles); 
//    ierr = MGGetSmoother(PC pc,int level,SLES *sles); 
//  int MGSetResidual(PC pc,int l,int (*residual)(Mat,Vec,Vec,Vec),Mat mat) 
//  int MGDefaultResidual(Mat mat,Vec b,Vec x,Vec r)
//    //need all per level except last only r
//     ierr = MGSetRhs(PC pc,int level,Vec b); 
//     ierr = MGSetX(PC pc,int level,Vec x); 
//     ierr = MGSetR(PC pc,int level,Vec r); 

void LinearSolver::printMatrices(const char* filename)
{
  using namespace cc;
  PRINT_MATRICES = true;
  ierr= PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
}

bool LinearSolver::forceOneIteration()
{
  return FORCE_ONE_ITERATION;
}

extern "C" {
int TIConverged(cc::KSP kspIn,int it,double rnorm, cc::KSPConvergedReason *reason,void *ctx)
{
  using namespace cc;
  // PetscFunctionBegin;
//    PetscValidHeaderSpecific( 
//                             kspIn , 
//                             KSP_COOKIE );
  PetscReal rtol;
  PetscReal atol;
  PetscReal dtol;
  int maxits;
  using namespace cc;
  KSPGetTolerances(kspIn,&rtol,&atol,&dtol,&maxits);
  *reason = KSP_CONVERGED_ITERATING;
  LinearSolver* thisLS=static_cast<LinearSolver*>(ctx);
  if (it || !thisLS->forceOneIteration()) 
    {
      if (rnorm < atol) 
        *reason = KSP_CONVERGED_ATOL;
    } 
  PetscFunctionReturn(0);
}
}

extern "C" 
{
int CoarseSolverConverged(cc::KSP kspIn,int it,double rnorm, cc::KSPConvergedReason *reason,void *ctx)
{
  static real ttol;
  using namespace cc;
  // PetscFunctionBegin;
  //    PetscValidHeaderSpecific( 
  //                             kspIn , 
  //                             KSP_COOKIE );
  PetscReal rtol;
  PetscReal atol;
  PetscReal dtol;
  int maxits;
  using namespace cc;
  KSPGetTolerances(kspIn,&rtol,&atol,&dtol,&maxits);
  *reason = KSP_CONVERGED_ITERATING;
  if (it==0)
    ttol   = rtol*rnorm + atol;
  else if (it && rnorm <= ttol ) 
    {
      if (rnorm < atol) 
        *reason = KSP_CONVERGED_ATOL;
      else                   
        *reason = KSP_CONVERGED_RTOL;
    } 
  PetscFunctionReturn(0);
}
}
cc::PetscErrorCode destroy(void*)
{
  return cc::PetscErrorCode(0);
}
void LinearSolver::useFixedIterationConvergence()
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles
  cc::KSPSetConvergenceTest(sles,cc::KSPConvergedSkip,this,destroy);
}

void LinearSolver::useTimeIntegrationConvergenceTest(bool forceOneIterationIn)
{
  using namespace cc;
  FORCE_ONE_ITERATION = forceOneIterationIn;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles
  ierr = KSPSetConvergenceTest(sles, 
                               TIConverged, 
                               this,
			       destroy); 
}

void LinearSolver::useCoarseSolverConvergenceTest(bool forceOneIterationIn)
{
  using namespace cc;
  FORCE_ONE_ITERATION = forceOneIterationIn;
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles
  ierr = KSPSetConvergenceTest(sles, 
                               CoarseSolverConverged, 
                               this,
			       destroy); 
}

void LinearSolver::useInitialGuess()
{
  //mwf 090104 PETSc 2.2.0 got rid of SLES completely
  //mwf replace ksp with sles
  //  cc::KSPSetInitialGuessNonzero(sles);
}
}//Petsc
}//Daetk


