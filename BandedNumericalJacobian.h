#ifndef BANDEDNUMERICALJACOBIAN_H
#define BANDEDNUMERICALJACOBIAN_H

#include "Definitions.h"
#include "Utilities.h"
#include "NumericalJacobian.h"

namespace Daetk 
{
template <class BandedMatrix>
class BandedNumericalJacobian : public NumericalJacobian
{
 public:
  BandedNumericalJacobian(BandedMatrix& M,VectorFunction& F);
  virtual ~BandedNumericalJacobian();
  bool evaluate(const Vec& x,const Vec& F);
  //mwf dummy for same interface
  void setStride(int str) {}
private:
  BandedMatrix& matrix;
};

template <class BandedMatrix>
BandedNumericalJacobian<BandedMatrix>::BandedNumericalJacobian(BandedMatrix& M,VectorFunction& F):
  NumericalJacobian(F),
  matrix(M)
{};


template <class BandedMatrix>
BandedNumericalJacobian<BandedMatrix>::~BandedNumericalJacobian(){}

template <class BandedMatrix>
bool BandedNumericalJacobian<BandedMatrix>::evaluate(const Vec& x,const Vec& Fatx)
{
  if (USE_ANALYTICAL_JACOBIAN)
    return ajac.evaluate(x,Fatx);
  matrix.zeroAll();
  int i,j,k,neq=matrix.getNeq(),
    kl=matrix.getLowerBandWidth(),
    ku=matrix.getUpperBandWidth(),
    bandWidth=ku + kl + 1,begin,end;
  
  attachToSubSystem(*Fp,Fatx);

  for (j=0;j<bandWidth;j++)    
    {
      tempDelta=0.0;
      for (k=j;k<neq;k+=bandWidth)
	{
          tempDeltaAttache(k) = deltaAttache(k);
        }
      Fp->correctArgument(tempDelta);
      bool evalError=false;
      FatxPdelta = Fp->value(evalError);
      if (evalError)
        {
          Fp->unCorrect();
          std::cerr<<"Delta in numerical jacobian is causing S or P to be out of range"<<std::endl;
          return true;
        }
      for (k=j;k<neq;k+=bandWidth)
	{
	  real tempDelInverse;
	  tempDelInverse = -1.0/ deltaAttache(k);
          begin=std::max(0,k-ku);
          end=std::min(neq,k+kl+1);
          for (i=begin;i<end;i++)
            matrix(i,k)=tempDelInverse*(FatxPdeltaAttache(i)-FatxAttache(i));
	}
      Fp->unCorrect();
    }
  return false;
}
}//Daetk
#endif










