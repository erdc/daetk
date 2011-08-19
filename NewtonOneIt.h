#ifndef NEWTON_ONE_IT_H
#define NEWTON_ONE_IT_H

#include "Newton.h"

namespace Daetk 
{
class NewtonOneIt : public Newton
{
public:
  NewtonOneIt();

  NewtonOneIt(LinearSolver& linearSolverIn, Jacobian& jacIn,VectorNorm& W,
	      DataCollector& dataIn, int neq, real lTol=0.005*0.33, 
	      real nlTol=0.33, int maxit=4, 
	      LineSearchType lsType = poorMansLS); 

  virtual ~NewtonOneIt();
  bool solve(Vec& correction,VectorFunction& F);

protected:

};
}//Daetk
#endif
