#ifndef  JACOBIAN_H
#define  JACOBIAN_H

#include "Definitions.h"
#include "Utilities.h"
#include "Vec.h"
#include "VectorFunction.h"
#include "LinearOperator.h"

namespace Daetk 
{
class Jacobian : public LinearOperator
{
public:
  virtual bool evaluate(const Vec& x, const Vec& F)=0;
  Jacobian(VectorFunction& F);
  virtual ~Jacobian();
  virtual void setFunction(VectorFunction& F);
  virtual void solveSubSystem(int start,int end,int stride,
                              int dimLS=0);
protected:
  VectorFunction* Fp;
};
}//Daetk
#endif
