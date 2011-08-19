#ifndef VECTORNORM_H
#define VECTORNORM_H

#include "Definitions.h"
#include "Utilities.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
class VectorNorm
{
public:
  VectorNorm();
  virtual ~VectorNorm();
  virtual real operator()(const Vec& y)=0;
  virtual void setWeight(const Vec& y)=0;
  virtual const Vec& getWeight() const=0;
  virtual const Vec& getScaling()=0;
  virtual void setTolerances(const Vec& atol, const Vec& rtol)=0;
  virtual void setTolerances(const real& atol,const real& rtol)=0;
  virtual real* getWeightBegin()=0;
  virtual real* getWeightEnd()=0;
  virtual void scale(const Vec& x,Vec& y)=0;
  virtual void deScale(const Vec& x, Vec& y)=0;
};
    
class Norm2 : public VectorNorm
{
public:
  Norm2();
  ~Norm2();
  virtual real operator()(const Vec& y);
  virtual void setWeight(const Vec& y);
  virtual const Vec& getWeight() const;
  virtual const Vec& getScaling();
  virtual void setTolerances(const Vec& atol, const Vec& rtol);
  virtual void setTolerances(const real& atol,const real& rtol);
  virtual real* getWeightBegin();
  virtual real* getWeightEnd();
  virtual void scale(const Vec& x,Vec& y);
  virtual void deScale(const Vec& x, Vec& y);
};
}//Daetk    
#endif
