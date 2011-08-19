#ifndef NUMERICALJACOBIAN_H
#define NUMERICALJACOBIAN_H

#include "Definitions.h"
#include "Utilities.h"
#include "Jacobian.h"
#include "AnalyticalJacobian.h"

namespace Daetk 
{
class NumericalJacobian : public Jacobian
{
 public:
  NumericalJacobian(VectorFunction& F);

  virtual void setFunction(VectorFunction& F);

  virtual~NumericalJacobian();
  
  bool apply(const Vec& x, Vec& Ax);

  virtual bool evaluate(const Vec& x,const Vec& F)=0;

  virtual void solveSubSystem(int start,int end,int stride,
                              int dimLS=0);

  virtual void attachToSubSystem(VectorFunction& F,const Vec& Fatx);

  virtual void useAnalyticalJacobian(bool flag);
protected:
  bool SOLVE_SUB,USE_ANALYTICAL_JACOBIAN;
  int str;
  VecIndex index;
  Vec tempDelta,FatxPdelta;
  Vec tempDeltaAttache,deltaAttache,FatxPdeltaAttache,FatxAttache;
  AnalyticalJacobian ajac;
};
}//Daetk
#endif
