#include "Jacobian.h"

namespace Daetk
{

Jacobian::Jacobian(VectorFunction& F):
  LinearOperator(F.dimDomain(),F.dimRange())
{
  Tracer tr("Jacobian::~Jacobian()");
  Fp=&F;
}
  
Jacobian::~Jacobian()
{
  Tracer tr("Jacobian::~Jacobian()");
}

void Jacobian::setFunction(VectorFunction& F)
{
  Fp=&F;
}

void Jacobian::solveSubSystem(int start,int end,int stride,
                              int dimLS)
{
}

}//Daetk
