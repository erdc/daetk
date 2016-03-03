#include "DaetkPetscAnalyticalJacobian.h"
#include <cstdio>

namespace Daetk 
{
namespace Petsc
{
  namespace cc
  {
#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
#include "petscmat.h"
  }
  
AnalyticalJacobian::~AnalyticalJacobian(){}

Mat& AnalyticalJacobian::getMatShell(){return matShell;}

AnalyticalJacobian* AnalyticalJacobian::theJacVec=0;

int AnalyticalJacobian::petscMatVec(cc::_p_Mat* A, cc::_p_Vec* x, cc::_p_Vec* Ax)
{
  using namespace cc;
  static Err ierr;
  void *context;
  ierr =  MatShellGetContext(A,&context);
  theJacVec = static_cast<AnalyticalJacobian*>(context);
  theJacVec->xVec.attachToPetscRepMulti(x);
  theJacVec->AxVec.attachToPetscRepMulti(Ax);
  theJacVec->apply(theJacVec->xVec,theJacVec->AxVec);
  theJacVec->xVec.detachFromPetscRepMulti();
  theJacVec->AxVec.detachFromPetscRepMulti();
  return 0; //? need to find the petsc error codes
}

AnalyticalJacobian::AnalyticalJacobian(Mat& matrixIn, VectorFunction& F):
  Daetk::AnalyticalJacobian(F)
{
  using namespace cc;
  Err ierr;
  void *context(static_cast<void*>(this));
  int m,n,M,N;
  ierr =  MatGetLocalSize(matrixIn.castToPetsc(),&m,&n);
  ierr =  MatGetSize(matrixIn.castToPetsc(),&M,&N);
  ierr =  MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,context,matShell.castToPetscLValue());
  ierr =  MatShellSetOperation(matShell.castToPetsc(),MATOP_MULT, reinterpret_cast<void(*)()>(petscMatVec) );
  matShell.beginAssembly();
  matShell.endAssembly();
  matShell.referenceVec.newsize(matrixIn.dimDomain());
}

}//Petsc
}//Daetk
