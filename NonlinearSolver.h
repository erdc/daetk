#ifndef NONLINEARSOLVER_H
#define NONLINEARSOLVER_H

#include "Definitions.h"
#include "Utilities.h"
#include "VectorFunction.h"

namespace Daetk 
{
class NonlinearSolver
{
 public:
  virtual bool solve(Vec& correction,VectorFunction& f)=0;
  virtual void setConvergenceFactor(const real& cf)=0;
  virtual void recomputeConvergenceRate()=0;
  //mwf added just to try and stick in line search
  virtual void useGlobalSearch();
  NonlinearSolver();
  virtual ~NonlinearSolver();
};
}//Daetk
#endif













