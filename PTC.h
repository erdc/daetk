#ifndef PTC_H
#define PTC_H

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
class PTC : public Integrator
{
public:
  PTC();
  virtual ~PTC();
  PTC(DaeDefinition& daeIn,LinearSolver& lsIn,Jacobian& jIn,
         VectorNorm& nmIn,DataCollector& dataIn,real lTol=0.001);
  void readParameters(ParameterDatabase& pd);
  bool step(const real& tout,real& tStep,Vec& solutionAtTStep, Vec& sp);
  bool calculateSolution(const real& tout,Vec& solutionAtTout, Vec& sp);
  void linearSolverIsInexact();
  void reset();
protected:
  const real STEADYSTATE;
  int nIts,delform,maxPTCits;
  bool err,evalFailed,useWRMSstop,useL2stop,INEXACT_LINEAR_SOLVER;
  real dt,alpha,dyNorm,dyNormOld,dy2Norm,dy2NormOld,res2Norm,res2NormOld,res0,tau,ssTol,ss2Tol,
    dtMax,dtMin,t,ypNorm,yp2Norm,userDt0,ypNormOld,yp2NormOld,dtOld,dtNew,
    roundOffTolerance,linearTolerance;
  Vec residual,solutionOld,solutionPrev,correction,solutionPrimeSave,residualSave,dyDiff,xtt,xtt_a;
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
  bool evaluateJacobian(const Vec& solution, const Vec& solutionPrime);
  void computeDeltaForJacobian(const Vec& solution, const Vec& solutionPrime);
  bool converged(const Vec& solution, const Vec& solutionPrime);
};

}//Daetk
#endif
