#include "ErrorNorms.h"

namespace Daetk
{

WeightedMaxNorm::WeightedMaxNorm(const int& Neq):
  WeightedRMSNorm(Neq)
{
}

WeightedMaxNorm::~WeightedMaxNorm()
{
}

real WeightedMaxNorm::operator()(const Vec& y)
{
  int ldim = weight.ldim();
  real errMax(-1.0);

  for (int i=0; i < ldim; i++)
    {
      real ei = std::fabs(y[i])*weight[i];
      errMax  = std::max(errMax,ei);
    }

#ifdef DEBUG_ERROR_CONTROLLERS
  if (errMax <= 0.0)
    {
      std::cerr<<"WARNING WRMaxNorm errMax= "<<errMax
	       <<" y = "<<y<<" weight= "<<weight<<std::endl;
    }
#endif

  return errMax;
}

}; //Daetk
