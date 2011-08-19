#include "VectorFunction.h"

namespace Daetk 
{

VectorFunction::VectorFunction(unsigned int m,unsigned int n):
  dimDomain_(m),
  dimRange_(n),
  deltaVF(m,12345)
{
  Tracer tr("VectorFunction::VectorFunction(unsigned int m,unsigned int n)");
}

VectorFunction::~VectorFunction()
{
  Tracer tr("VectorFunction::~VectorFunction()");
}

int VectorFunction::dimRange()
{
  return dimRange_;
}

int VectorFunction::dimDomain()
{
  return dimDomain_;
}

}//Daetk
