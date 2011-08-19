#ifndef DAETKPETSCANALYTICALJACOBIAN_H
#define DAETKPETSCANALYTICALJACOBIAN_H

#include <vector>
#include "Definitions.h"
#include "IntVec.h"
#include "Utilities.h"
#include "AnalyticalJacobian.h"
#include "DaetkPetscMat.h"
#include "DaetkPetscSys.h"

namespace Daetk 
{
namespace Petsc
{
  namespace cc
    {
  struct _p_Mat;
  struct _p_Vec;
    }

class AnalyticalJacobian : public Daetk::AnalyticalJacobian
{
public:
  AnalyticalJacobian(Petsc::Mat& M, VectorFunction& F);
  virtual ~AnalyticalJacobian();
  virtual Petsc::Mat& getMatShell();
  static int petscMatVec(Petsc::cc::_p_Mat* A, Petsc::cc::_p_Vec* x, Petsc::cc::_p_Vec* Ax);
  static AnalyticalJacobian* theJacVec;
  Petsc::Mat matShell;
  Vec xVec,AxVec;
};
}//Petscn
}//Daetk
#endif
