#ifndef BICGSTAB_H
#define BICGSTAB_H

#include <bitset>
#include <fstream>
#include "Definitions.h"
#include "IntVec.h"
#include "LinearSolver.h"
#include "LaFullDirectSolver.h"
#include "LaBandedDirectSolver.h"
#include "DataCollector.h"
#include "Definitions.h" 
#include "VectorNorm.h"
#include "VectorFunction.h"
#include "LinearOperator.h"
#include "Preconditioner.h"
#include "ParameterDatabase.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
class BiCGstab : public LinearSolver
{ 
public:
  BiCGstab();
  BiCGstab(LinearOperator& linOpIn, Preconditioner& precIn, 
           VectorNorm& normIn, 
           DataCollector& dataIn,
           ParameterDatabase& pd);
  BiCGstab(LinearOperator& linOpIn, Preconditioner& precIn, 
           VectorNorm& normIn, 
           DataCollector& dataIn,
           real tolIn=1.65e-4,
           int maxIts=100,int neq=0);

  virtual ~BiCGstab();
  bool prepare(); 
  void useInitialGuess();
  bool solve(const Vec& bIn,Vec& xIn);
  void checkTrueResidual(real resTol);
  void dump(std::ostream& os);
private:
  bool RESTART, CHECK_RESIDUAL,error,USE_GUESS;
  int maxIterations;
  real tol,trueResTol,rhoK,rhoKm1,alpha,omega,k,beta,realTemp,rNorm,xNorm,r0;
  Vec r,rHat,v,p,s,unscaled_p,unscaled_s,t,vecTemp;
  AttacheVec x,b;
  VectorNorm* vectorNorm;
  DataCollector* data;
  LinearOperator* linearOperator;
  Preconditioner* preconditioner;
};
} //Daetk 
#endif








