#ifndef DAEDEFINITION_H
#define DAEDEFINITION_H

#include "Definitions.h"
#include "Utilities.h"
#include "DataCollector.h"
#include "Vec.h"
#include "VectorFunction.h"
#include "VecOperators.h"
#include "VecBlas.h"
//mwf added
#include "VectorNorm.h"

namespace Daetk 
{

class DaeDefinition : public VectorFunction
{ 
public:
  DaeDefinition(int dim=0,DataCollector* dataIn=0);
  virtual bool residual(const real& t,const Vec& y,const Vec& yp,Vec& res)=0;
  virtual bool yPrimeValue(const real& t,const Vec& y,Vec& yp);
  virtual bool evaluateDaeJacobian(const real& t,const Vec& y,const Vec& yp, 
                                   const real& alpha);
  virtual bool jacVec(const Vec& x, Vec& Jx);
  virtual const real& getT0()=0;
  virtual const Vec& getY0()=0;
  virtual const Vec& getY0prime()=0;
  virtual ~DaeDefinition();
  virtual void stepTaken(){}
  const Vec& argument();
  const Vec& value(bool& evalError);
  virtual void correctArgument(Vec& correction);
  virtual void unCorrect();
  void resetFunction();
  virtual void initializeFunction(const Vec& y, const Vec& yp, const Vec& f);
  bool evaluateAnalyticalJacobian();
  virtual bool numericalJacVec(const Vec& v, Vec& Jv);
  bool analyticalJacVec(const Vec& v, Vec& Jv);

  //mwf added to keep track of weighted norm for nonlinear system
  virtual void setWeightedNorm(VectorNorm & wNormIn)
    { weightedNorm = &wNormIn; }
  virtual VectorNorm* getWeightedNorm() 
    { return weightedNorm; }
  virtual Petsc::VecIndex getTimeIndex() const
    { return timeIndex; }
  virtual void setTimeIndex(const Petsc::VecIndex& tIndexIn)
    { timeIndex = tIndexIn; }
   //mwf added so that can compute own Jacobian for other integrators
  virtual void computeDeltaForJacobian();

  //mwf

  bool updateF,updateJac;
  real alphaDaeDef,
    tDaeDef;
  Vec yDaeDef,
    ypDaeDef,
    Fcurrent,
    yLast,
    ypLast,
    Flast,
    del,
    betaDaeDef;
  DataCollector* data;

  //mwf added
  VectorNorm * weightedNorm; 
  Petsc::VecIndex timeIndex;
  //allow DaeDefinition to compute its own delta for jacobian 
  //by default if integrator hasn't done it
  bool computeOwnJacobianDelta;
};

}//Daetk
#endif
