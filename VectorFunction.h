#ifndef VECTORFUNCTION_H
#define VECTORFUNCTION_H

#include "Definitions.h"
#include "Utilities.h"
#include "Vec.h"

namespace Daetk 
{
class VectorFunction
{
public:
  VectorFunction(unsigned int domain=0, unsigned int range=0);
  virtual void usePicardApproximation(bool yesNo){std::cerr<<"no picard approximation available"<<std::endl;}
  virtual ~VectorFunction();
  virtual const Vec& argument()=0;
  virtual const Vec& value(bool& evalError)=0;
  virtual void correctArgument(Vec& correction)=0;
  virtual void unCorrect()=0;
  virtual bool evaluateAnalyticalJacobian()=0;
  virtual bool numericalJacVec(const Vec& v, Vec& Jv)=0;
  virtual bool analyticalJacVec(const Vec& v, Vec& Jv)=0;
  //mwf added for modified newton ss?
  virtual void computeDeltaForJacobian() {}
  int dimDomain();
  int dimRange();
  unsigned int dimDomain_, dimRange_;
  Vec deltaVF;
};
}//Daetk
#endif
