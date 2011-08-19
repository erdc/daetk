#include "AnalyticalJacobian.h"

namespace Daetk 
{

  AnalyticalJacobian::AnalyticalJacobian(VectorFunction& F):
    Jacobian(F)
  {}
  
  AnalyticalJacobian::~AnalyticalJacobian(){}
  
  bool AnalyticalJacobian::evaluate(const Vec&,const Vec& )
  {
      return Fp->evaluateAnalyticalJacobian();
  }
  
  bool AnalyticalJacobian::apply(const Vec& x,Vec& Jx)
  {
    return Fp->analyticalJacVec(x,Jx);
  }
  
}//Daetk
