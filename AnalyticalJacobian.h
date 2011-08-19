#ifndef ANALYTICALJACOBIAN_H
#define ANALYTICALJACOBIAN_H

#include "Definitions.h"
#include "Utilities.h"
#include "Jacobian.h"

namespace Daetk 
{

class AnalyticalJacobian : public Jacobian
{
public:
  AnalyticalJacobian(VectorFunction& F);
  virtual ~AnalyticalJacobian();
  virtual bool evaluate(const Vec&,const Vec& );
  bool apply(const Vec& x,Vec& Jx);
};
}//Daetk
#endif








