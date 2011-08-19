#ifndef PRECG_H
#define PRECG_H

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
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
class PreCG : public LinearSolver
{ 
public:
  PreCG();
  PreCG(LinearOperator& linOpIn, Preconditioner& precIn, 
	VectorNorm& normIn, 
	DataCollector& dataIn,
	real tolIn=1.65e-4,
	int maxIts=100);

  virtual ~PreCG();
  bool prepare(); 
  bool solve(const Vec& bIn,Vec& xIn);
  void useInitialGuess();
private:
  bool error,USE_GUESS;
  int maxIterations,k;
  real tol,
    rhoK,
    alpha,
    beta,
    tauKm1,
    tauKm2,
    bNorm,
    pDotw,
    zNorm,
    xNorm;
  Vec r,
    p,
    Ax,
    z,
    w;
  AttacheVec x,b;
  VectorNorm* vectorNorm;
  DataCollector* data;
  LinearOperator* linearOperator;
  Preconditioner* preconditioner;
};
}//Daetk  
#endif
