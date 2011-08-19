#ifndef NEWTON_H
#define NEWTON_H

#include <fstream>
#include "Definitions.h"
#include "NonlinearSolver.h"
#include "Utilities.h"
#include "LinearSolver.h"
#include "Jacobian.h"
#include "VectorNorm.h"
#include "DataCollector.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
class Newton : public NonlinearSolver
{
public:
  enum LineSearchType { poorMansLS, cubicLS, armijoLS };
  void setConvergenceFactor(const real&){}
  void recomputeConvergenceRate(){}
  
  Newton();

  Newton(LinearSolver& linearSolverIn, Jacobian& jacIn,VectorNorm& W,
                   DataCollector& dataIn, int neq, real lTol=0.005*0.33, 
 		   real nlTol=0.33, int maxit=4, 
		   LineSearchType lsType = poorMansLS);

  void testResidual(real tol);
  void testConvergenceRate(bool flag=true);
  void linearSolverIsInexact();
  virtual ~Newton();
  bool solve(Vec& correction,VectorFunction& F);
  void computeRate();
  bool converged(VectorFunction& F);
  virtual void useLineSearch(int nls=10);
  virtual void usePicardHybrid(){USE_PICARD_HYBRID=true;}
  virtual void setLineSearchFact(real lsRedIn);
  virtual void setLineSearchMethod(LineSearchType lsType);
  bool logSteps(bool yesNo=true){LOG_STEPS=yesNo; return false;}
  //mwf added for NonIterative time integrators
  virtual void setMaxIterations(const int& maxItIn)
  { maxIterations = maxItIn; }
  virtual int getMaxIterations()
  { return maxIterations; }
protected:
  bool USE_LINE_SEARCH,LOG_STEPS,USE_PICARD_HYBRID;
  LineSearchType lineSearchMethod;
  int  nLineSearches;
  real lsRedFact;
  //
  bool INEXACT_LINEAR_SOLVER,
    TEST_RESIDUAL,
    TEST_RATE,
    evalError;
  int iterations,
    maxIterations,
    lin_it;

  real s,
    rate,
    normOfInitialGuess,
    norm2OfInitialGuess,
    normOfLastCorrection,
    normOfCorrection,
    nonlinearTolerance,
    roundOffTolerance,
    r0,
    resTol,
    linearTolerance,FnNew,FnOld;

  Vec p,pLS,p0,tmpArg,tmpF,argPrev,funPrev,
    residual;

  VectorNorm* weightedNorm;
  DataCollector* data;
  LinearSolver* linearSolver;
  std::ofstream mnout;

  //go through and perform switch to see which linesearch to use
  bool lineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
		  bool& evalFailed,
		  VectorFunction& F);

  bool cubicLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
		       bool& evalFailed,
		       VectorFunction& F);

  bool poorMansLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
			  bool& evalFailed,
			  VectorFunction& F);
  
  bool armijoLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
			bool& evalFailed,
			VectorFunction& F);
  Jacobian* jac;

};
}//Daetk
#endif
