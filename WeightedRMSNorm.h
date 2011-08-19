#ifndef WEIGHTEDRMSNORM_H
#define WEIGHTEDRMSNORM_H

#include "Definitions.h"
#include "VectorNorm.h"
#include "Utilities.h"
#include "Vec.h"
#include "VecOperators.h"

namespace Daetk 
{
class WeightedRMSNorm : public VectorNorm
{
public:
  WeightedRMSNorm(const int& Neq=0);
  WeightedRMSNorm(WeightedRMSNorm& WN);
  virtual ~WeightedRMSNorm();
  real operator()(const Vec& y);
  void setWeight(const Vec& y);
  const Vec& getWeight() const;
  const Vec&  getScaling();
  void setTolerances(const Vec& atolIn,const Vec& rtolIn);
  void setTolerances(const real& atolR,const real& rtolR);
  real* getWeightBegin();
  real* getWeightEnd();
  void scale(const Vec& x,Vec& y);
  void deScale(const Vec& x, Vec& y);
  void print();
protected: //mwf was private
  int i,neq;
  real rneq,oneOver_sqrt_neq;
  Vec weight,scaling,atol,rtol,tmp;
};
   
inline const Vec& WeightedRMSNorm::getWeight() const
{
  return weight;
}

  class WeightedL2Norm : public WeightedRMSNorm
  {
  public:
    WeightedL2Norm(const int& Neq=0);
    void setTolerances(const real& atolR,const real& rtolR, const Vec& dV);
    real operator()(const Vec& y);
    void setWeight(const Vec& y);
  protected:
    double V;
    Vec dV;
  };
  
     
}//Daetk
#endif



