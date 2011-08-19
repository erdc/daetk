#ifndef FULLNUMERICALJACOBIAN_H
#define FULLNUMERICALJACOBIAN_H

#include "Definitions.h"
#include "Utilities.h"
#include "NumericalJacobian.h"

namespace Daetk 
{
template<class FullMatrix>
class FullNumericalJacobian : public NumericalJacobian
{
 public:
  FullNumericalJacobian(FullMatrix& M,VectorFunction& F);
  virtual ~FullNumericalJacobian();
  bool evaluate(const Vec& x,const Vec& F);
private:
  Vec tempDelta;
  FullMatrix& matrix;
};
  
template<class FullMatrix>
FullNumericalJacobian<FullMatrix>::FullNumericalJacobian(FullMatrix& M,VectorFunction& F):
  NumericalJacobian(F),
  tempDelta(F.dimDomain()),
  matrix(M)
{}

template<class FullMatrix>
FullNumericalJacobian<FullMatrix>::~FullNumericalJacobian(){}

template<class FullMatrix>
bool FullNumericalJacobian<FullMatrix>::evaluate(const Vec& x,const Vec& Fatx)
{
  matrix=0.0;
  int neq=x.dim();
  real tempDelInverse;
  for (int j=0;j<neq;j++)
    {
      tempDelta=0.0;
      tempDelta(j) = Fp->deltaVF(j);
      Fp->correctArgument(tempDelta);
      bool evalError=false;
      FatxPdelta = Fp->value(evalError);
      if (evalError)
        {
          Fp->unCorrect();
          std::cerr<<"Delta in numerical jacobian is causing S or P to be out of range"<<std::endl;
          return true;
        }
      tempDelInverse = -1.0/Fp->deltaVF(j);
      for (int i=0;i<neq;i++)
	{
	  matrix(i,j)=tempDelInverse*(FatxPdelta(i)-Fatx(i));
	}
      Fp->unCorrect();
    }
  return false;
};
}//Daetk
#endif

