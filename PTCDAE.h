#ifndef PTCDAE_H
#define PTCDAE_H

#include <fstream>
#include <iostream>
#include "Definitions.h"
#include "Integrator.h"
#include "Utilities.h"
#include "DaeDefinition.h"
#include "LinearSolver.h"
#include "ModifiedNewton.h"
#include "Jacobian.h"
#include "VectorNorm.h"
#include "DataCollector.h"
#include "ParameterDatabase.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
class PTCDAE : public Integrator
{
public:
  PTCDAE();
  virtual ~PTCDAE();
  PTCDAE(DaeDefinition& daeIn,LinearSolver& lsIn,Jacobian& jIn,
         VectorNorm& nmIn,DataCollector& dataIn,real lTol=0.001);
  void readParameters(ParameterDatabase& pd);
  bool step(const real& tout,real& tStep,Vec& solutionAtTStep, Vec& sp);
  bool calculateSolution(const real& tout,Vec& solutionAtTout, Vec& sp);
  void linearSolverIsInexact();
  void reset();
  const Vec& Fvalue(bool& errFlag);
protected:
  const real STEADYSTATE;
  int nIts,delform,maxPTCits,nCorrectorIts;
  bool err,evalFailed,useWRMSstop,useL2stop,INEXACT_LINEAR_SOLVER;
  real dt,alpha,dyNorm,dyNormLast,dy2Norm,dy2NormLast,res2Norm,res2NormLast,res0,tau,ssTol,ss2Tol,
    dtMax,dtMin,t,ypNorm,yp2Norm,userDt0,ypNormLast,yp2NormLast,dtLast,dtNew,
    roundOffTolerance,linearTolerance;
  Vec residual,residualLast,solutionLast,solutionBeforeLast,solutionPrimeLast,correction,correctionLast,totalCorrection,dyDiff,xtt,xtt_a;
  std::ofstream iterout;
//  ("iterations.grf");
//    iterout.setf(std::ios::scientific);
//    iterout.precision(10);
  DaeDefinition* dae;
  DataCollector* data;
  VectorNorm* weightedNorm;
  Jacobian* jacobian;
  LinearSolver* linearSolver;
  void chooseDt(const Vec& solution, const Vec& solutionPrime);
  bool evaluateJacobian();
  void computeDeltaForJacobian();
  bool converged(const Vec& solution, const Vec& solutionPrime);
};

}//Daetk
#endif
