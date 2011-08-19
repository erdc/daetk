#ifndef FLCBDF_lite_H
#define FLCBDF_lite_H

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
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
class FLCBDF_lite
{
public:
  FLCBDF_lite();

  virtual ~FLCBDF_lite();

  FLCBDF_lite(int dim,
              VectorNorm&,
              DataCollector&);

  bool takeInitialSteps(real& t,real& tStep,Vec& solutionAtTStep, 
			Vec& sp);
  void useFixedStep(const real step);
  void useFixedOrder(const int order);
  void updateJacobians();
  void useInterpolant();
  void reset();
  real chooseDT(const Daetk::real&, const Daetk::real&);
  real setDT(const real& hin);
  real chooseInitialStepSize(const real& t, 
                             const real& tout,
                             const Vec& y,
                             const Vec& yPrime);
  bool initializationPhaseStep();
  bool step();
  real estimateError(const Vec& y);
  bool initializationPhaseSolverFailure(int&);
  bool initializationPhaseCheckError(const Daetk::Petsc::Vec&);
  bool errorForStepTooLarge(const Daetk::Petsc::Vec&);
  bool checkError(const Vec& y);
  bool solverFailure();
  bool calculate_yprime(const Vec& y,
                        const Vec& Dy,
                        Vec& yprime,
                        Vec& Dyprime);
  bool initializeTimeHistory(const Vec& yIn, const Vec& yPrimeIn);
  double retryStep_solverFailure();
  double retryStep_errorFailure();
  double getCurrentAlpha() const;//mwf for building off-diagonal jacobians, not necessarily best idea
  Vec yn,ynprime,yCminusyP,tempvec;
private:
  
  bool firstStep,
    prepareJacobian,                  
    predictorAlreadyUpdated,                      
    firstCallToSolver,
    jacobianIsOld,
    correctorStepFailed,
    kRaisedOnLastStep,
    USE_FIXED_STEP,
    USE_FIXED_ORDER,
    UPDATE_JACOBIANS,
    STEP_EXACT,
    inInitializationPhase,
    recomputeNonlinearConvergenceRate,
    stepFailed,
    useOldChooseStepSize;

  int neq,                               //number of equtions    
    k,                                   //order
    kLast,
    constStepsTaken,
    constStepsTakenSave,
    kLastUpdate,
    constStepsTakenLastUpdate,
    constStepsTakenLastUpdateSave,
    stepsFailed,
    errorFailures,
    MAX_FAILURES,
    shortCutFactor,
    kFixed,
    iStep;
  
  real h,                                 //stepsize
    hLast,
    hLastUpdate,
    tn,                                   //time of previous step
    tnPlusH,
    error,
    errorEstimate,
    norm_yCminusyP,
    convergenceFactor,                    //for nonlinear solver
    hFixed,
    alpha,
    stepIncreaseCeiling,
    stepIncreaseFloor,
    stepDecreaseCeiling,
    stepDecreaseFloor;
							    
  //coefficient data

  real psiNp1[7],
    alphaNp1[6],
    betaNp1[6],
    sigmaNp1[7],
    gammaNp1[6],
    alphaS,
    alphaoNp1,
    alphaOld,
    alphaLast;
  
  Vec phiN[7];
  
  DataCollector* data;
 
  VectorNorm* weightedNorm;

  //private member functions 

  bool analyzeJacobian();
  bool evaluateJacobian();
  void computeDeltaForJacobian(Vec& deltaVF);
  void resetCoefficients();
  void updateCoefficients();
  void chooseOrderForStep();
  void predictor();
  void  interpolant(const real& T,Vec& yAtT,Vec& yPAtT);
  void chooseStepSize();
  //mwf added so that I can try and eliminate need to interpolate to
  //mwf time for final soluiton
  void chooseStepSize(const real& tin, const real& tout);
  bool corrector();
  bool errorForStepTooLarge();
  void checkStepSize();

};
}//Daetk
#endif
