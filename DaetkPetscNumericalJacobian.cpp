#include "DaetkPetscNumericalJacobian.h"
#include <cstdio>
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
#include "petscmat.h"
#ifndef DAETK_DEF_CPLUSPLUS_FOR_PETSC_H
#define __cplusplus
#endif
     //   int MatShellSetOperation(Mat,MatOperation,JacobianBase::PetscMatVecType);  
   }
 }
  
JacobianBase::~JacobianBase(){}

  Mat& JacobianBase::getMatShell()
  {
    if (USE_ANALYTICAL_JACOBIAN)
      return pajac.getMatShell();
    return matShell;
  }

JacobianBase* JacobianBase::theJacVec=0;

int JacobianBase::petscMatVec(cc::_p_Mat* A, cc::_p_Vec* x, cc::_p_Vec* Ax)
{
  static Err ierr;
  void *context;
  ierr =  cc::MatShellGetContext(A,&context);
  theJacVec = static_cast<JacobianBase*>(context);
  theJacVec->xVec.attachToPetscRepMulti(x);
  theJacVec->AxVec.attachToPetscRepMulti(Ax);
  theJacVec->apply(theJacVec->xVec,theJacVec->AxVec);
  theJacVec->xVec.detachFromPetscRepMulti();
  theJacVec->AxVec.detachFromPetscRepMulti();
  return 0; //? need to find the petsc error codes
}

JacobianBase::JacobianBase(Petsc::Mat& matrixIn, VectorFunction& F):
  Daetk::NumericalJacobian(F),
  pajac(matrixIn,F)
{
  using namespace  cc;
  Err ierr;
  void *context(static_cast<void*>(this));
  int m,n,M,N;
  ierr =  MatGetLocalSize(matrixIn.castToPetsc(),&m,&n);
  ierr =  MatGetSize(matrixIn.castToPetsc(),&M,&N);
  ierr =  MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,context,matShell.castToPetscLValue());

  ierr =  MatShellSetOperation(matShell.castToPetsc(),MATOP_MULT,reinterpret_cast<void (*)()>(petscMatVec));
  matShell.beginAssembly();
  matShell.endAssembly();
  matShell.referenceVec.newsize(matrixIn.dimDomain());
}

  
}//Petsc
}//Daetk
