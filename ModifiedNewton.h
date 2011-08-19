#ifndef MODIFIEDNEWTON_H
#define MODIFIEDNEWTON_H

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
class ModifiedNewton : public NonlinearSolver
{
public:
  //mwf add these as part of base?
  enum LineSearchType { poorMansLS, cubicLS, armijoLS,simpleLS };

  ModifiedNewton();
  ModifiedNewton(LinearSolver& linearSolverIn, VectorNorm& W,
                 DataCollector& dataIn, int neq, real lTol=0.005*0.33, real nlTol=0.33, int maxit=4);
  virtual ~ModifiedNewton();
  void doChordIteration();
  void doFullNewton();
  void testResidual(real tol);
  void testConvergenceRate(bool flag=true);
  void linearSolverIsInexact();
  bool solve(Vec& correction,VectorFunction& F);
  void setConvergenceFactor(const real& cf);
  void computeRate();
  void recomputeConvergenceRate();
  bool converged(VectorFunction& F);
  //mwf added just to see if I could turn line search on and off
  virtual void useLineSearch(int nls=10);
  virtual void setLineSearchFact(real lsRedIn);
  virtual void setLineSearchMethod(LineSearchType lsType);
  //mwf also added dimension here
  virtual void solveSubSystem(int start,int end,int stride,
			      int dimLS=0);
 
  //mwf changed to protected so I could use in derived class
protected:
  bool SOLVE_SUB;
  VecIndex index;
  Vec attache;
  int str;
  //mwf added just to see if I could turn line search on and off
  bool USE_LINE_SEARCH;
  LineSearchType lineSearchMethod;
  int  nLineSearches;
  real lsRedFact;
  //
  bool CHORD_ITERATION, 
    INEXACT_LINEAR_SOLVER,
    TEST_RESIDUAL,
    TEST_RATE,
    RECOMPUTE_RATE,
    evalError;
  int iterations,
    maxIterations;

  const real *convergenceFactorp;

  real s,
    rate,
    normOfPredictor,
    normOfFirstCorrection,
    normOfCorrection,
    nonlinearTolerance,
    roundOffTolerance,
    r0,
    resTol;

  Vec p,
    residual;

  VectorNorm* weightedNorm;
  DataCollector* data;
  LinearSolver* linearSolver;
  std::ofstream mnout;
};

class ModifiedNewtonMM : public ModifiedNewton
{
public:
  enum LineSearchType { poorMansLS, cubicLS, armijoLS, simpleLS};
  ModifiedNewtonMM();

  ModifiedNewtonMM(LinearSolver& linearSolverIn, VectorNorm& W,
		   DataCollector& dataIn, int neq, 
		   real lTol=0.005*0.33, 
		   real nlTol=0.33, int maxit=4);

  ModifiedNewtonMM(LinearSolver& linearSolverIn, VectorNorm& W,
		   DataCollector& dataIn, int neq,
		   VectorNorm& Wlin,
		   real lTol=0.005*0.33, 
		   real nlTol=0.33, int maxit=4);
  virtual ~ModifiedNewtonMM();

  bool solve(Vec& correction,VectorFunction& F);
  virtual void solveSubSystem(int start,int end,int stride,
			      int dimLS=0);

bool lineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
				  bool& evalFailed,
				  Vec& corr,
				  VectorFunction& F);
bool cubicLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
				     bool& evalFailed,
				     Vec& corr,
                                       VectorFunction& F);
bool simpleLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
				      bool& evalFailed,
				      Vec& corr,
                                        VectorFunction& F);
bool poorMansLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
					  bool& evalFailed,
					  Vec& corr,
					  VectorFunction& F);
bool armijoLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
					bool& evalFailed,
					Vec& corr,
					VectorFunction& F);

protected:
  Vec attache2,attache3;
  Vec resForLSolve,pForLSolve;
  VectorNorm* weightedNormLinSys;
};

class ModifiedNewtonSS : public ModifiedNewton
{
public:

  ModifiedNewtonSS();

  ModifiedNewtonSS(LinearSolver& linearSolverIn, Jacobian& jacIn,VectorNorm& W,
                   DataCollector& dataIn, int neq, real lTol=0.005*0.33, 
 		   real nlTol=0.33, int maxit=4, 
		   LineSearchType lsType = poorMansLS);

  virtual ~ModifiedNewtonSS();
  bool solve(Vec& correction,VectorFunction& F);
  bool converged(VectorFunction& F);
protected:
  //go through and perform switch to see which linesearch to use
  bool lineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
		  bool& evalFailed,
		  Vec& corr,
		  VectorFunction& F);

  bool cubicLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
		       bool& evalFailed,
		       Vec& corr,
		       VectorFunction& F);

  bool poorMansLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
			  bool& evalFailed,
			  Vec& corr,
			  VectorFunction& F);
  
  bool armijoLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
			bool& evalFailed,
			Vec& corr,
			VectorFunction& F);
  Jacobian* jac;

};
}//Daetk
#endif
